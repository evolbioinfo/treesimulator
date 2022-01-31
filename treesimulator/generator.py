import logging
import random

import numpy as np
from ete3 import TreeNode

from treesimulator import STATE, DIST_TO_START, TIME_TILL_NOW


def simulate_tree_gillespie(model, max_time=np.inf, max_sampled=np.inf,
                            state_feature=STATE, state_frequencies=None):
    """
    Simulates the tree evolution from a root over the given time based on the given model.

    :param state_feature: a name of the tree feature which will store node states
    :param state_frequencies: array of equilibrium frequencies of model states
        (to be used to draw the root state). If not given, will be taken from the model (by default all equal).
    :param max_sampled: maximal number of sampling node (when reached, the simulation stops), by default infinity
    :type max_sampled: int
    :param max_time: time over which we generate a tree, by default infinity
    :type max_time: float
    :param model: MTBD model used to simulate the tree
    :type model: treesimulator.mtbd_models.Model
    :return: the simulated tree
    :rtype: ete3.Tree
    """
    num_states = len(model.states)
    if state_frequencies is None:
        state_frequencies = model.state_frequencies
    root_state = np.random.choice(np.arange(num_states), size=1, p=state_frequencies)[0]
    # evolve till the time is up, following Gillespie
    time = 0
    infectious_nums = np.zeros(num_states, dtype=np.int)
    infectious_nums[root_state] = 1
    sampled_nums = np.zeros(num_states, dtype=np.int)

    infectious_state2id = [set() for _ in model.states]
    cur_id = 0, 0
    infectious_state2id[root_state].add(cur_id)
    id2time = {}
    id2parent = {}
    sampled_id2state = {}

    while infectious_nums.sum() and sampled_nums.sum() < max_sampled and time < max_time:
        # first we need to calculate rate sum
        transmission_rate_sums = model.transmission_rates.sum(axis=1) * infectious_nums
        transition_rate_sums = model.transition_rates.sum(axis=1) * infectious_nums
        removal_rate_sums = model.removal_rates * infectious_nums
        total_transmission_rate = transmission_rate_sums.sum()
        total_transition_rate = transition_rate_sums.sum()
        total_removal_rate = removal_rate_sums.sum()
        total_rate = total_transmission_rate + total_transition_rate + total_removal_rate

        # now let us see when next event takes place
        time += np.random.exponential(1 / total_rate, 1)[0]

        # Check if the time is up
        if time > max_time:
            break

        # now let us see which event will happen
        random_event = np.random.uniform(0, 1, 1)[0] * total_rate

        # case 1: state transition
        if random_event < total_transition_rate:
            for i in range(num_states):
                if random_event < transition_rate_sums[i]:
                    infectious_nums[i] -= 1
                    random_event /= infectious_nums[i]
                    for j in range(num_states):
                        if random_event < model.transition_rates[i, j]:
                            infectious_nums[j] += 1
                            infectious_state2id[j].add(random_pop(infectious_state2id[i]))
                            break
                        random_event -= model.transition_rates[i, j]
                    break
                random_event -= transition_rate_sums[i]
            continue
        random_event -= total_transition_rate

        # case 2: transmission
        if random_event < total_transmission_rate:
            for i in range(num_states):
                if random_event < transmission_rate_sums[i]:
                    random_event /= infectious_nums[i]
                    for j in range(num_states):
                        if random_event < model.transmission_rates[i, j]:
                            infectious_nums[j] += 1
                            cur_id = cur_id[0] + 1, 0
                            parent_id = random_pop(infectious_state2id[i])
                            donor_id = parent_id[0], parent_id[1] + 1
                            infectious_state2id[i].add(donor_id)
                            infectious_state2id[j].add(cur_id)
                            id2parent[cur_id] = parent_id
                            id2parent[donor_id] = parent_id
                            id2time[parent_id] = time
                            break
                        random_event -= model.transmission_rates[i, j]
                    break
                random_event -= transmission_rate_sums[i]
            continue
        random_event -= total_transmission_rate

        # case 3: sampling
        for i in range(num_states):
            if random_event < removal_rate_sums[i]:
                infectious_nums[i] -= 1

                sampled_id = random_pop(infectious_state2id[i])
                id2time[sampled_id] = time

                if np.random.uniform(0, 1, 1)[0] < model.ps[i]:
                    sampled_id2state[sampled_id] = model.states[i]
                    sampled_nums[i] += 1

                break
            random_event -= total_removal_rate
    if max_time == np.inf:
        max_time = time

    return reconstruct_tree(id2parent, id2time, sampled_id2state, max_time, state_feature=state_feature)


def reconstruct_tree(id2parent, id2time, sampled_id2state, max_time, state_feature=STATE):
    if not sampled_id2state:
        return None

    root = None
    id2node = {}
    for node_id, state in sampled_id2state.items():
        time = id2time[node_id]
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


def random_pop(elements):
    """
    Removes a random element from a list and returns it.
    :param elements: list of elemetns
    :return: the selected element
    """
    element = random.sample(elements, 1)[0]
    elements.remove(element)
    return element


def generate_forest(model, max_time=np.inf, min_tips=1000, max_sampled=np.inf, keep_nones=False, state_feature=STATE,
                    state_frequencies=None):
    total_n_tips = 0
    forest = []
    total_trees = 0
    sampled_trees = 0
    while total_n_tips < min_tips:
        tree = simulate_tree_gillespie(model, max_time=max_time, max_sampled=max_sampled,
                                       state_feature=state_feature, state_frequencies=state_frequencies)
        total_trees += 1
        if tree:
            total_n_tips += len(tree)
            sampled_trees += 1
        if tree or keep_nones:
            forest.append(tree)
    return forest


def generate(model, min_tips, max_tips, T=np.inf, state_frequencies=None):
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
    :type: int
    :return: the simulated forest (containing only one tree in case of a tree simulation)
        and stats on total number of tips, on the number of hidden trees (0  case of a tree simulation),
        and on total time T
    :rtype: tuple(list(ete3.Tree), (int, int, float))
    """

    if max_tips < min_tips:
        raise ValueError('--max_tips cannot be smaller than --min_tips')

    if T < np.inf:
        while True:
            forest = generate_forest(model, max_time=T, min_tips=min_tips,
                                     keep_nones=True, state_frequencies=state_frequencies)
            total_trees = len(forest)
            forest = [tree for tree in forest if tree is not None]
            fl = len(forest)
            u = total_trees - fl
            total_tips = sum(len(list(t.iter_leaves())) for t in forest)
            if total_tips <= max_tips:
                logging.info('Generated a forest of {} visible and {} hidden trees with {} sampled tips over time {}.'
                             .format(fl, u, total_tips, T))
                return forest, (total_tips, u, T)
            logging.debug('Generated a forest of {} visible and {} hidden trees with {} sampled tips over time {}.'
                          .format(fl, u, total_tips, T))
    else:
        while True:
            max_sampled = int(min_tips + np.random.random() * (max_tips + 1 - min_tips))
            tree = simulate_tree_gillespie(model, max_time=np.inf, max_sampled=max_sampled,
                                           state_frequencies=state_frequencies)
            total_tips = len(tree) if tree else 0
            if total_tips >= min_tips:
                logging.info('Generated a tree with {} sampled tips over time {}.'
                             .format(total_tips, getattr(tree, TIME_TILL_NOW) if tree else np.inf))
                return [tree], (total_tips, 0, getattr(tree, TIME_TILL_NOW) + tree.dist)
            logging.debug('Generated a tree with {} sampled tips over time {}.'
                          .format(total_tips, getattr(tree, TIME_TILL_NOW) if tree else np.inf))
