import logging
from collections import defaultdict, Counter

import numpy as np
from ete3 import TreeNode

from mtbd_models import BirthDeathModel, BirthDeathExposedInfectiousModel
from treesimulator import STATE, save_forest, TIME
from treesimulator.generator import check_skyline_times, simulate_epidemic_gillespie, reconstruct_ltt, \
    annotate_tree_with_time, name_forest


def remove_certain_leaves(tr, to_remove=lambda node: False):
    """
    Removes all the branches leading to leaves identified positively by to_remove function.
    :param tr: the tree of interest (ete3 Tree)
    :param to_remove: a method to check is a leaf should be removed.
    :return: void, modifies the initial tree.
    """

    tips = [tip for tip in tr if to_remove(tip)]
    for node in tips:
        if node.is_root():
            return None
        parent = node.up
        parent.remove_child(node)
        # If the parent node has only one child now, merge them.
        if len(parent.children) == 1:
            brother = parent.children[0]
            brother.dist += parent.dist
            if parent.is_root():
                brother.up = None
                tr = brother
            else:
                grandparent = parent.up
                grandparent.remove_child(parent)
                grandparent.add_child(brother)
    return tr

def get_n2avg_time(model, max_time, root_state, sampled=False):
    n2times = defaultdict(list)

    for _ in range(10):
        id2parent_id, id2node_time, sampled_ids, id2state, time, observed_nums = (
            simulate_epidemic_gillespie([model], skyline_times=None,
                                    max_time=max_time,
                                    max_sampled=100,
                                    state_frequencies=None,
                                    max_notified_contacts=0, notify_at_removal=False,
                                    root_state=root_state))
        if sampled:
            for i, tid in enumerate(sorted(sampled_ids, key=lambda _: id2node_time[_]), start=1):
                n2times[i].append(id2node_time[tid])
        else:
            time2num = reconstruct_ltt(id2parent_id, id2node_time)
            num2time = defaultdict(lambda: np.inf)
            for time in sorted(time2num.keys()):
                n = time2num[time]
                if n >= 100:
                    num2time[100] = min(num2time[100], time)
                    break
                num2time[n] = min(num2time[n], time)
                for _ in range(n - 1, 0, -1):
                    if num2time[_] > time:
                        num2time[_] = time
                    else:
                        break
            first_sampling_time = min(id2node_time[tid] for tid in sampled_ids) if sampled_ids else -np.inf
            for n, t in num2time.items():
                n2times[n].append(max(t, first_sampling_time))

    n2times = {n: np.mean(times) for n, times in n2times.items()}
    mtime = max(n2times.values())
    res = defaultdict(lambda: mtime)
    res.update(n2times)
    return res

def generate_hierarchically(models, skyline_times=None, max_time=np.inf, max_sampled=np.inf,
                            state_frequencies=None, root_state=None):
    if max_sampled <= 100:
        return simulate_epidemic_gillespie(models, skyline_times=skyline_times, max_time=max_time,
                                           max_sampled=max_sampled, state_frequencies=state_frequencies,
                                           max_notified_contacts=0,  notify_at_removal=False,
                                           root_state=root_state)

    if max_sampled >= np.inf and max_time >= np.inf:
        raise ValueError(
            'A limit on tree generation must be specified: it can be either the maximum number of sampled tips (max_sampled) '
            'or the maximum length of the sampling period (max_time). '
            'If both are specified, then the generation will stop when the first of them is reached.')

    # Check if we're using skyline models (list with more than one model)
    use_skyline = isinstance(models, list) and len(models) > 1

    # Determine next model change time (first element of skyline_times or infinity if only one model)
    next_model_change = np.inf if not use_skyline else skyline_times[0]

    # Validate skyline times if using skyline
    if use_skyline:
        check_skyline_times(skyline_times, len(models))

    # Initialize variables to track the current model
    current_model = models[0] if isinstance(models, list) else models
    current_model_index = 0

    num_states = len(current_model.states)
    state_indices = np.arange(num_states)
    root_states = state_indices[current_model.states == root_state]
    if len(root_states) > 0:
        root_state = root_states[0]
    else:
        if state_frequencies is None:
            state_frequencies = current_model.state_frequencies
        root_state = np.random.choice(state_indices, size=1, p=state_frequencies)[0]

    n2avg_time = get_n2avg_time(model=current_model, max_time=min(next_model_change, max_time), root_state=root_state)
    time_step = n2avg_time[100]

    root = None
    while root is None:
        cur_time = 0

        while root is None:
            root = generate_tree_over_time(current_model, root_state, time_step)
        cur_time += time_step
        n_sampled = sum(1 for tip in root if getattr(tip, 'sampled', 'no') == 'yes')
        while n_sampled < max_sampled and cur_time < max_time:
            unfinished_tips = [tip for tip in root if getattr(tip, 'sampled', 'no') == 'no']
            if not unfinished_tips:
                root = None
                break

            if cur_time >= next_model_change:
                # Move to the next model in the skyline
                current_model_index += 1
                current_model = models[current_model_index]
                next_model_change = skyline_times[current_model_index] if current_model_index < len(skyline_times) else np.inf
                n2avg_time = get_n2avg_time(model=current_model, max_time=min(next_model_change, max_time),
                                            root_state=Counter(getattr(tip, STATE) for tip in unfinished_tips)\
                                            .most_common(n=1)[0][0])

            n_unfinished_tips = len(unfinished_tips)
            time_step = n2avg_time[min(100, max(2, int((max_sampled - n_sampled) / n_unfinished_tips))) \
                if max_sampled < np.inf else 100]
            time_step = min(time_step, min(next_model_change, max_time) - cur_time)

            cur_tips = sum(1 for tip in root if getattr(tip, 'sampled', 'no') == 'yes')
            logging.info(f'Have {cur_tips} sampled tips at time {cur_time:.2f}, '
                         f'extending by {time_step:.2f} time units each of {len(unfinished_tips)} unsampled tips')
            for tip in unfinished_tips:
                state = getattr(tip, STATE)
                tree = generate_tree_over_time(current_model, state, time_step)
                if tree:
                    if tip.is_root():
                        tree.dist += root.dist
                        root = tree
                    else:
                        parent = tip.up
                        tip.up.remove_child(tip)
                        parent.add_child(tree, dist=tip.dist + tree.dist)
            cur_time += time_step
            n_sampled = sum(1 for tip in root if getattr(tip, 'sampled', 'no') == 'yes')

    logging.info(f'Generated a tree with {len(root)} sampled or unfinished tips.')
    root = remove_certain_leaves(root, lambda n: getattr(n, 'sampled', 'no') == 'no')
    name_forest([root])
    logging.info(f'Pruned unfinished tips, kept {len(root)} sampled tips.')
    if len(root) > max_sampled:
        annotate_tree_with_time(root)
        first_tips = {_.name for _ in sorted([tip for tip in root], key=lambda n: getattr(n, TIME))[:max_sampled]}
        root = remove_certain_leaves(root, lambda n: n.name not in first_tips)
        logging.info(f'Pruned extra tips, kept {len(root)} sampled tips.')
    return root



def generate_tree_over_time(model, root_state, max_time):
    id2parent_id, id2node_time, sampled_ids, id2state, time, observed_nums = (
        simulate_epidemic_gillespie(model, skyline_times=None,
                                    max_time=max_time,
                                    max_sampled=np.inf,
                                    state_frequencies=None,
                                    max_notified_contacts=0, notify_at_removal=False,
                                    root_state=root_state))
    id2children_ids = defaultdict(set)

    # Figure out the last branch id for each individual
    individual_id2last_id = defaultdict(lambda: 0)
    for (individual_id, branch_id) in id2parent_id.keys():
        individual_id2last_id[individual_id] = max(branch_id, individual_id2last_id[individual_id])
    tip_ids = list(individual_id2last_id.items())

    # if there was no transmission, then root is also a tip
    if not tip_ids:
        if (0, 0) in sampled_ids:
            root = TreeNode(name='0_0', dist=id2node_time[(0, 0)])
            root.add_feature(STATE, id2state[0])
            return root
        else:
            return None

    for id, parent_id in id2parent_id.items():
        id2children_ids[parent_id].add(id)

    root_id = (0, 0)
    for tip_id in tip_ids:
        # unsampled but removed tip -- prune it
        if tip_id not in sampled_ids and tip_id in id2node_time:
            if tip_id == root_id:
                return None
            del id2node_time[tip_id]
            del id2state[tip_id[0]]
            parent_id = id2parent_id[tip_id]
            del id2parent_id[tip_id]
            id2children_ids[parent_id].remove(tip_id)
            siblings = id2children_ids[parent_id]
            # if parent has only one child left, branch it to the grandparent directly if there is one
            if len(siblings) == 1:
                only_child = list(siblings)[0]
                del id2children_ids[parent_id]
                if parent_id in id2parent_id:
                    grandparent_id = id2parent_id[parent_id]
                    del id2parent_id[parent_id]
                    id2parent_id[only_child] = grandparent_id
                    id2children_ids[grandparent_id].remove(parent_id)
                    id2children_ids[grandparent_id].add(only_child)
                else:
                    root_id = only_child
                    del id2parent_id[only_child]

    id2node = {}
    sampled_ids = set(sampled_ids)

    def get_node(id):
        if id not in id2node:
            pid = id2parent_id[id] if id in id2parent_id else None
            parent_time = id2node_time[pid] if pid in id2node_time else 0
            node_time = id2node_time[id] if id in id2node_time else max_time
            node_dist = node_time - parent_time
            node = TreeNode(dist=node_dist)
            if id in sampled_ids:
                node.add_feature('sampled', 'yes')
            if id in tip_ids:
                node.add_feature(STATE, id2state[id[0]])
            if pid is not None:
                parent = get_node(pid)
                parent.add_child(node, dist=node_dist)
            id2node[id] = node
        return id2node[id]

    for id in id2parent_id.keys():
        get_node(id)

    root = get_node(root_id)
    # annotate time till now
    for node in root.traverse('postorder'):
        if node.is_leaf():
            node.add_feature('sampled', getattr(node, 'sampled', 'no'))
            node.add_feature('visible', getattr(node, 'sampled'))
        else:
            node.add_feature('visible',
                             'yes' if any(getattr(child, 'visible') == 'yes' for child in node.children) else 'no')
            node.add_feature('sampled', 'yes' if sum(
                getattr(child, 'visible') == 'yes' for child in node.children) > 1 else 'no')
    return root


if '__main__' == __name__:
    logging.getLogger().handlers = []
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s | %(levelname)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
    root = generate_hierarchically([BirthDeathExposedInfectiousModel(p=0.5, la=2, psi=1, mu=0.5)], max_sampled=1000000)
    save_forest([root], 'hierarchical_tree.nwk', format=3)