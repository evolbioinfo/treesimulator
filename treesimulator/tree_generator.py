import logging
import random

import numpy as np
from ete3 import TreeNode

from treesimulator import DIST_TO_START, TIME_TILL_NOW, STATE


def simulate_tree_gillespie(states, rates, max_time, max_sampled=np.inf, root_state=None,
                            state_feature=STATE, return_tree=True):
    """
    Simulates the tree evolution from a root over the given time
    based on the given sampling, state transition, transmission, death and partner notification rates.
    :param root_state: root state, if set to None a random state will be chosen according to equilibrium frequencies.
    :param max_sampled: maximal number of sampling node (when reached, the simulation stops)
    :param max_time: float, time over which we generate a tree.
    :param states: possible states
    :param rates: numpy.matrix of rates, where columns correspond to model states,
        and rows correspond to state transitions, transmissions, sampling and equilibrium frequencies.
    :return: the simulated tree (ete3.Tree).
    """
    if root_state is None:
        root_state = np.random.choice(states, size=1, p=rates[-1, :])[0]

    # evolve till the time is up, following Gillespie
    time = 0
    infectious_nums = np.zeros(len(states), dtype=np.int)
    infectious_nums[root_state.index] = 1
    sampled_nums = np.zeros(len(states), dtype=np.int)

    infectious_state2id = [set() for _ in states]
    cur_id = 0, 0
    infectious_state2id[root_state.index].add(cur_id)
    id2time = {}
    id2parent = {}
    sampled_id2state = {}

    while infectious_nums.sum() and sampled_nums.sum() < max_sampled and time < max_time:
        # first we need to calculate rate sum
        rate_sums = rates[:-1, :].dot(infectious_nums)
        total_rate = rate_sums.sum()

        # now let us see when next event takes place
        time += np.random.exponential(1 / total_rate, 1)[0]

        # Check if the time is up
        if time > max_time:
            break

        # now let us see which event will happen
        random_event = np.random.uniform(0, 1, 1)[0] * total_rate

        # case 1: state transition
        if random_event < rate_sums[0]:
            transition_rates = rates[0, :] * infectious_nums
            for i in range(len(states)):
                if random_event < transition_rates[i]:
                    state = states[i]
                    infectious_nums[state.index] -= 1
                    infectious_nums[state.next_state.index] += 1

                    infectious_state2id[state.next_state.index].add(random_pop(infectious_state2id[state.index]))
                    break
                random_event -= transition_rates[i]
            continue
        random_event -= rate_sums[0]

        # case 2: transmission
        if random_event < rate_sums[1]:
            transmission_rates = rates[1, :] * infectious_nums
            for i in range(len(states)):
                if random_event < transmission_rates[i]:
                    state = states[i]
                    infectious_nums[state.recipient.index] += 1

                    cur_id = cur_id[0] + 1, 0
                    parent_id = random_pop(infectious_state2id[state.index])
                    donor_id = parent_id[0], parent_id[1] + 1
                    infectious_state2id[state.index].add(donor_id)
                    infectious_state2id[state.recipient.index].add(cur_id)
                    id2parent[cur_id] = parent_id
                    id2parent[donor_id] = parent_id
                    id2time[parent_id] = time
                    break
                random_event -= transmission_rates[i]
            continue
        random_event -= rate_sums[1]

        # case 3: sampling
        sampling_rates = rates[2, :] * infectious_nums
        for i in range(len(states)):
            if random_event < sampling_rates[i]:
                state = states[i]
                infectious_nums[state.index] -= 1
                sampled_nums[state.index] += 1

                sampled_id = random_pop(infectious_state2id[state.index])
                id2time[sampled_id] = time
                sampled_id2state[sampled_id] = state
                break
            random_event -= sampling_rates[i]

    if return_tree:
        return reconstruct_tree(id2parent, id2time, sampled_id2state, max_time, state_feature=state_feature)
    return [id2time[_] for _ in sampled_id2state.keys()], infectious_nums


def reconstruct_tree(id2parent, id2time, sampled_id2state, max_time, state_feature=STATE):
    if not sampled_id2state:
        return None

    root = None
    id2node = {}
    for id, state in sampled_id2state.items():
        time = id2time[id]
        node = TreeNode(dist=time - (0 if not id in id2parent else id2time[id2parent[id]]), name=id[0])
        node.add_feature(DIST_TO_START, time)
        node.add_feature(state_feature, state)
        id2node[id] = node
        while node is not None and id in id2parent:
            parent_id = id2parent[id]
            if parent_id in id2node:
                id2node[parent_id].add_child(node)
                break
            parent_time = id2time[parent_id]
            parent = TreeNode(dist=parent_time - (0 if not parent_id in id2parent else id2time[id2parent[parent_id]]))
            parent.add_feature(DIST_TO_START, parent_time)
            id2node[parent_id] = parent
            parent.add_child(node)
            node, id = parent, parent_id
        if id not in id2parent:
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
        node.add_feature(TIME_TILL_NOW, max_time - getattr(node, DIST_TO_START) + node.dist / 2)
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


def generate_forest(states, rates, max_time, min_tips=1_000, keep_nones=False, state_feature=STATE):
    total_n_tips = 0
    forest = []
    total_trees = 0
    sampled_trees = 0
    while total_n_tips < min_tips:
        tree = simulate_tree_gillespie(states, rates, max_time=max_time, state_feature=state_feature)
        total_trees += 1
        if tree:
            total_n_tips += len(tree)
            sampled_trees += 1
        if tree or keep_nones:
            forest.append(tree)

    logging.info('# sampled trees in the forest = {} (out of {}), total # of tips = {}, time = {}.'
                 .format(sampled_trees, total_trees, total_n_tips, max_time))
    return forest


def subforest(forest, filter, percentile=.5):
    trees_to_take = []
    for tree in forest:
        for node in tree.traverse('postorder'):
            if node.is_leaf():
                node.add_feature('to_take', filter(node))
            else:
                to_take = all(getattr(_, 'to_take') for _ in node.children) #and filter(getattr(node, STATE).name)
                node.add_feature('to_take', to_take)
                if not to_take:
                    for _ in node.children:
                        # if not (_.is_leaf()) and getattr(_, 'to_take'):
                        if getattr(_, 'to_take'):
                            trees_to_take.extend(_.children)
                elif node.is_root():
                    trees_to_take.extend(node.children)
            for child in node.children:
                delattr(child, 'to_take')
            if node.is_root():
                delattr(node, 'to_take')

    time = np.percentile([getattr(t, TIME_TILL_NOW) + t.dist / 2 for t in trees_to_take], q=percentile)

    subforest = []
    for tree in trees_to_take:
        if getattr(tree, TIME_TILL_NOW) + tree.dist / 2 < time:
            continue

        todo = [tree]
        while todo:
            node = todo.pop()
            time_at_node = getattr(node, TIME_TILL_NOW) - node.dist / 2
            if time_at_node >= time:
                todo.extend(node.children)
            else:
                node = node.copy("deepcopy")
                node.dist = time - time_at_node
                node.add_feature(TIME_TILL_NOW, time_at_node + node.dist / 2)
                subforest.append(node)

    num_trees_at_time = 0
    for tree in forest:
        todo = [tree]
        while todo:
            node = todo.pop()
            time_at_node = getattr(node, TIME_TILL_NOW) - node.dist / 2
            if time_at_node >= time:
                todo.extend(node.children)
            else:
                num_trees_at_time += 1

    return subforest, num_trees_at_time, time
