import logging
import random
from collections import Counter

import numpy as np
import scipy
from ete3 import TreeNode

from treesimulator import STATE, DIST_TO_START, TIME_TILL_NOW
from treesimulator.mtbd_models import CTModel

TRANSITION = 0
TRANSMISSION = 1
REMOVAL = 2
EVENT_TYPES = np.array([TRANSITION, TRANSMISSION, REMOVAL])


def simulate_tree_gillespie(model, max_time=np.inf, min_sampled=0, max_sampled=np.inf,
                            state_feature=STATE, state_frequencies=None, ltt=False, max_notified_contacts=1,
                            root_state=None):
    """
    Simulates the tree evolution from a root over the given time based on the given model.

    :param state_feature: a name of the tree feature which will store node states
    :param state_frequencies: array of equilibrium frequencies of model states
        (to be used to draw the root state). If not given, will be taken from the model (by default all equal).
    :param min_sampled: minimal number of sampled nodes (when reached, the simulation could stop), by default 0
    :param max_sampled: maximal number of sampled nodes (when reached, the simulation stops), by default infinity
    :type max_sampled: int
    :param max_time: time over which we generate a tree, by default infinity
    :type max_time: float
    :param model: MTBD model used to simulate the tree
    :type model: treesimulator.mtbd_models.Model
    :param ltt: whether to return LTT values as well (default False)
    :type ltt: bool
    :param max_notified_contacts: maximum notified contacts for -CT models (by default 1, meaning the most recent contact)
    :type max_notified_contacts: int
    :return: the simulated tree and the time it covers
    :rtype: ete3.Tree
    """
    num_states = len(model.states)
    if state_frequencies is None:
        state_frequencies = model.state_frequencies
    state_indices = np.arange(num_states)
    root_states = state_indices[model.states == root_state]
    root_state = np.random.choice(state_indices, size=1, p=state_frequencies)[0] \
        if len(root_states) == 0 else root_states[0]
    # evolve till the time is up, following Gillespie
    time = 0
    infectious_nums = np.zeros(num_states, dtype=np.int64)
    infectious_nums[root_state] = 1
    sampled_nums = np.zeros(num_states, dtype=np.int64)

    infectious_state2id = [set() for _ in model.states]
    cur_id = 0, 0
    infectious_state2id[root_state].add(cur_id)
    id2time = {}
    id2parent_id = {}
    sampled_id2state = {}
    donor_id2recipient_id = {}
    id2current_id = {0: 0}
    id2state = {0: root_state}

    target_sampled = np.round(np.random.uniform(low=min_sampled, high=max_sampled, size=1)[0], 0) \
        if max_sampled < np.inf else np.inf

    logging.debug('Aiming for {} sampled cases over time {}'
                  .format(target_sampled if target_sampled < np.inf else 'any number of', max_time))

    num_states_squared = np.power(num_states, 2)
    rate_vector = np.concatenate([model.transition_rates.reshape(num_states_squared),
                                  model.transmission_rates.reshape(num_states_squared),
                                  model.removal_rates])
    index_vector = np.arange(rate_vector.shape[0])
    transmission_rates_per_state = model.transmission_rates.sum(axis=1)
    # to avoid dividion by zero
    transmission_rates_per_state[transmission_rates_per_state == 0] = 1
    transmission_probs = model.transmission_rates / transmission_rates_per_state.reshape((num_states, 1))

    while infectious_nums.sum() and sampled_nums.sum() < target_sampled and time < max_time:
        infectious_num_vector = np.concatenate([np.tile(infectious_nums.reshape((num_states, 1)), (2, num_states))
                                               .reshape(num_states_squared * 2), infectious_nums])
        total_rate_vector = rate_vector * infectious_num_vector
        total_rate = total_rate_vector.sum()

        # check if we passed the time limit
        time += np.random.exponential(1 / total_rate, 1)[0]
        if time >= max_time:
            time = max_time
            break

        # now let us see which event will happen
        random_event_index = np.random.choice(index_vector, p=total_rate_vector / total_rate_vector.sum(),
                                              replace=False, size=1)[0]
        # case 1: state transition
        if random_event_index < num_states_squared:
            i, j = random_event_index // num_states, random_event_index % num_states
            infectious_nums[i] -= 1
            infectious_nums[j] += 1
            state_changing_id = random_pop(infectious_state2id[i])
            id2state[state_changing_id[0]] = j
            infectious_state2id[j].add(state_changing_id)
            logging.debug('Time {}:\t{} changed state from {} to {}'
                          .format(time, state_changing_id, model.states[i], model.states[j]))

            continue
        random_event_index -= num_states_squared

        # case 2: transmission
        if random_event_index < num_states_squared:
            i, js = random_event_index // num_states, [random_event_index % num_states]
            if model.n_recipients[i] > 1:
                extra_recipients = scipy.stats.poisson.rvs(model.n_recipients[i] - 1, size=1)[0]
                if extra_recipients:
                    js.extend(np.random.choice(state_indices, p=transmission_probs[i, :], replace=True,
                                               size=extra_recipients))
            parent_id = random_pop(infectious_state2id[i])
            id2time[parent_id] = time
            donor_id = parent_id[0], parent_id[1] + 1
            id2current_id[donor_id[0]] = donor_id[1]
            infectious_state2id[i].add(donor_id)
            id2parent_id[donor_id] = parent_id

            recipient_ids = []
            for j in js:
                infectious_nums[j] += 1
                cur_id = cur_id[0] + 1, 0
                id2current_id[cur_id[0]] = 0
                id2state[cur_id[0]] = j
                infectious_state2id[j].add(cur_id)
                id2parent_id[cur_id] = parent_id
                recipient_ids.append(cur_id)
            donor_id2recipient_id[parent_id] = recipient_ids
            logging.debug('Time {}:\t{} in state {} transmitted to {} in state(s) {}'
                          .format(time, parent_id, model.states[i], recipient_ids, [model.states[j] for j in js]))
            continue
        random_event_index -= num_states_squared

        # case 3: removal
        i = random_event_index
        infectious_nums[i] -= 1
        removed_id = random_pop(infectious_state2id[i])
        id2time[removed_id] = time
        msg = 'Time {}:\t{} in state {} got removed'.format(time, removed_id, model.states[i])

        if np.random.uniform(0, 1, 1)[0] < model.ps[i]:
            sampled_id2state[removed_id] = model.states[i]
            sampled_nums[i] += 1
            msg += ' and sampled'

            # contact tracing
            if isinstance(model, CTModel):
                contact_n = max_notified_contacts
                # if the max number of notified contacts per person allows and it is not the root
                while contact_n > 0 and removed_id in id2parent_id:
                    parent_id = id2parent_id[removed_id]
                    donor_id = (parent_id[0], parent_id[1] + 1)
                    contact_ids = donor_id2recipient_id[parent_id] if removed_id == donor_id else [donor_id]
                    contact_ids = [contact_ids[_]
                                   for _ in np.random.choice(np.arange(len(contact_ids)), replace=False,
                                                             size=min(len(contact_ids), contact_n))]
                    contact_n -= len(contact_ids)
                    for contact_id in contact_ids:
                        contact_id = contact_id[0], id2current_id[contact_id[0]]
                        if removed_id != donor_id:
                            # means this is the last contact this node has
                            contact_n = 0
                        else:
                            removed_id = parent_id
                        # If the notifier agrees and the contact is not yet removed, notify them
                        if (np.random.uniform(0, 1, 1)[0] < model.upsilon
                                and contact_id not in id2time):
                            unnotified_contact_i = id2state[contact_id[0]]
                            # The ids are organised as follows: s1, s2, ..., sm, s1-n, s2-n, ..., sm-n
                            # Hence if we have an id >= m then the contact was already notified by someone else
                            if unnotified_contact_i < num_states // 2:
                                notified_contact_i = num_states // 2 + unnotified_contact_i
                                id2state[contact_id[0]] = notified_contact_i
                                infectious_state2id[unnotified_contact_i].remove(contact_id)
                                infectious_state2id[notified_contact_i].add(contact_id)
                                infectious_nums[unnotified_contact_i] -= 1
                                infectious_nums[notified_contact_i] += 1
                            msg += ' and notified {} in state {}' \
                                .format(contact_id, model.states[unnotified_contact_i])
        logging.debug(msg)

    if max_time == np.inf:
        max_time = time

    root = reconstruct_tree(id2parent_id, id2time, sampled_id2state, max_time, state_feature=state_feature)
    if ltt:
        return root, reconstruct_ltt(id2parent_id, id2time), max_time
    return root, max_time


def reconstruct_tree(id2parent, id2time, sampled_id2state, max_time, state_feature=STATE):
    if not sampled_id2state:
        return None

    root = None
    id2node = {}
    for node_id, state in sampled_id2state.items():
        time = id2time[node_id]
        # if there is PN it could be
        # that we decided to stop at an earlier time when the unsampled contact proportion was lower
        if time > max_time:
            continue
        node = TreeNode(dist=time - (0 if node_id not in id2parent else id2time[id2parent[node_id]]), name=node_id[0])
        node.add_feature(DIST_TO_START, time)
        node.add_feature(state_feature, state)
        id2node[node_id] = node
        while node is not None and node_id in id2parent:
            parent_id = id2parent[node_id]
            if parent_id in id2node:
                id2node[parent_id].add_child(node)
                break
            parent_time = id2time[parent_id]
            parent = TreeNode(dist=parent_time - (0 if parent_id not in id2parent else id2time[id2parent[parent_id]]))
            parent.add_feature(DIST_TO_START, parent_time)
            id2node[parent_id] = parent
            parent.add_child(node)
            node, node_id = parent, parent_id
        if node_id not in id2parent:
            root = node
    # remove internal nodes with just one child
    for node in root.traverse('postorder'):
        if len(node.children) == 1:
            child = node.children[0]
            child.dist += node.dist
            if not node.is_root():
                parent = node.up
                parent.remove_child(node)
                parent.add_child(child)
            else:
                root = child
                child.up = None
    # annotate time till now
    for node in root.traverse('postorder'):
        node.add_feature(TIME_TILL_NOW, max_time - getattr(node, DIST_TO_START))
    return root


def reconstruct_ltt(id2parent_id, id2time):
    time2num = Counter()
    time2num[0] = 1
    parents = set(id2parent_id.values())
    for id, time in id2time.items():
        if id in parents:
            time2num[time] += 1
        else:
            time2num[time] -= 1
    total = 0
    for time in sorted(time2num.keys()):
        total += time2num[time]
        time2num[time] = total
    return time2num


def observed_ltt(forest, T):
    time2num = Counter()
    time2num[0] = len(forest)
    for tree in forest:
        for n in tree.traverse():
            if n.is_leaf():
                time2num[T - getattr(n, TIME_TILL_NOW)] -= 1
            else:
                time2num[T - getattr(n, TIME_TILL_NOW)] += 1
    total = 0
    for time in sorted(time2num.keys()):
        total += time2num[time]
        time2num[time] = total
    return time2num


def random_pop(elements):
    """
    Removes a random element from a list and returns it.
    :param elements: list of elements
    :return: the selected element
    """
    element = random.sample(list(elements), 1)[0]
    elements.remove(element)
    return element


def generate_forest(model, max_time=np.inf, min_tips=1000, keep_nones=False, state_feature=STATE,
                    state_frequencies=None, ltt=False, max_notified_contacts=1):
    total_n_tips = 0
    forest = []
    total_trees = 0
    sampled_trees = 0
    res_ltt = None
    while total_n_tips < min_tips:
        if ltt:
            tree, cur_ltt, _ = simulate_tree_gillespie(model, max_time=max_time, ltt=True,
                                                       state_feature=state_feature, state_frequencies=state_frequencies,
                                                       max_notified_contacts=max_notified_contacts)
            if res_ltt is None:
                res_ltt = cur_ltt
            else:
                total = 0
                prev_res, prev_cur = 0, 0
                for time in sorted(set(cur_ltt.keys()) | set(res_ltt.keys())):
                    if time in cur_ltt:
                        total += cur_ltt[time] - prev_cur
                        prev_cur = cur_ltt[time]
                    if time in res_ltt:
                        total += res_ltt[time] - prev_res
                        prev_res = res_ltt[time]
                    res_ltt[time] = total
        else:
            tree, _ = simulate_tree_gillespie(model, max_time=max_time,
                                              state_feature=state_feature, state_frequencies=state_frequencies,
                                              max_notified_contacts=max_notified_contacts)
        total_trees += 1
        if tree:
            total_n_tips += len(tree)
            sampled_trees += 1
        if tree or keep_nones:
            forest.append(tree)
    if ltt:
        return forest, res_ltt
    else:
        return forest


def generate(model, min_tips, max_tips, T=np.inf, state_frequencies=None, max_notified_contacts=1):
    """
    Simulates a tree (or a forest of trees, if --T is specified) for given MTBD model parameters.

    If a simulation leads to less than --min_tips tips, it is repeated.
    For a tree simulation, if --min_tips and --max_tips are equal, exactly that number of tips will be simulated.
    If --min_tips is less than --max_tips, a value randomly drawn between one and another will be simulated.

    :param model: MTBD model to use for tree/forest generation
    :type model: treesimulator.mtbd_models.Model
    :param min_tips: desired minimal bound on the total number of simulated leaves
    :type min_tips: int
    :param max_tips: desired maximal bound on the total number of simulated leaves
    :type max_tips: int
    :param T: total simulation time.
        If specified, a forest will be simulated instead of one tree.
        The trees in this forest will be simulated during the given time,
        till the --min_tips number is reached.
        If after simulating the last tree, the forest exceeds the --max_tips number, the process will be restarted.
    :type T: float
    :param state_frequencies: array of model state frequencies,
        to be used to draw the root states. If not given, will be taken from the model (by default all equal).
    :type state_frequencies: list(float)
    :param max_notified_contacts: maximum notified contacts for -CT models (by default 1, meaning the most recent contact)
    :type max_notified_contacts: int
    :return: the simulated forest (containing only one tree in case of a tree simulation),
        stats on total number of tips, on the number of hidden trees (0  case of a tree simulation), and on total time T,
        and the LTT numbers as a mapping between times and numbers of infected individuals
    :rtype: tuple(list(ete3.Tree), (int, int, float), dict(float, int))
    """

    if max_tips < min_tips:
        raise ValueError('--max_tips cannot be smaller than --min_tips')

    if T < np.inf:
        while True:
            forest, ltt = generate_forest(model, max_time=T, min_tips=min_tips,
                                          keep_nones=True, state_frequencies=state_frequencies, ltt=True,
                                          max_notified_contacts=max_notified_contacts)
            total_trees = len(forest)
            forest = [tree for tree in forest if tree is not None]
            fl = len(forest)
            u = total_trees - fl
            total_tips = sum(len(list(t.iter_leaves())) for t in forest)
            if total_tips <= max_tips:
                logging.info('Generated a forest of {} visible and {} hidden trees with {} sampled tips over time T={}.'
                             .format(fl, u, total_tips, T))
                return forest, (total_tips, u, T), ltt
            logging.debug('Generated a forest of {} visible and {} hidden trees with {} sampled tips over time T={}.'
                          .format(fl, u, total_tips, T))
    else:
        while True:
            tree, ltt, max_time = simulate_tree_gillespie(model, max_time=np.inf,
                                                          max_sampled=max_tips, min_sampled=min_tips,
                                                          state_frequencies=state_frequencies, ltt=True,
                                                          max_notified_contacts = max_notified_contacts)
            total_tips = len(tree) if tree else 0
            if total_tips >= min_tips:
                logging.info('Generated a tree with {} sampled tips over time T={}.'.format(total_tips, max_time))
                return [tree], (total_tips, 0, max_time), ltt
            logging.debug('Generated a tree with {} sampled tips over time T={}.'.format(total_tips, max_time))
