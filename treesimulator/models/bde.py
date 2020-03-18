import logging

import numpy as np

from likelihood.convolution import get_convolution
from treesimulator.models import Model, SAMPLING, TRANSMISSION, TRANSITION, State

EXPOSED = 'e'
INFECTED = 'i'


class BirthDeathExposedModel(Model):

    def __init__(self, p=0.5, *args, **kwargs):
        Model.__init__(self, ps=[0, p], minus_avg_sigma=False, *args, **kwargs)

    def num_params(self, type=None):
        if SAMPLING == type:
            return 2
        if TRANSMISSION == type:
            return 1
        if TRANSITION == type:
            return 1
        return 4

    def _get_states(self):
        infected_state = State(name=INFECTED, index=1)
        exposed_state = State(name=EXPOSED, index=0, next_state=infected_state)

        # we'll put exposed as a recipient in order not to break the simulator,
        # but as the transmission rate from this state is 0, it will be ok anyway
        exposed_state.recipient = exposed_state
        infected_state.recipient = exposed_state

        return np.array([exposed_state, infected_state])

    def params2rates(self, ps, **kwargs):
        """
        Converts parameters into a rate array.

        :param ps: parameters, in the following order:
            transition E->I, transmission from I to E, sampling of I
        """
        mu, lambda_i, psi_e, psi_i = ps
        self.rates[0, 0] = mu
        self.rates[1, 1] = lambda_i
        self.rates[2, :] = psi_e, psi_i
        self.rates[3, :] = self.get_sd()

    def get_bounds(self, lb, ub, **kwargs):
        return np.array([[lb, ub]] * self.num_params())

    def _get_sd_formulas(self):
        mu = self.rates[0, 0]
        lambda_i = self.rates[1, 1]
        psi_e, psi_i = self.rates[2, :]

        p = [
            lambda_i - psi_i + psi_e,
            psi_i - psi_e + mu,
            -mu
        ]
        roots = np.roots(p)
        roots = roots[np.where(np.isreal(roots) & (~np.isnan(roots)) & (roots >= 0) & (roots <= 1))]
        pi_i = roots[0]

        return np.array([1 - pi_i, pi_i])

    def map_state(self, state):
        state = str(state)
        if state.startswith('e'):
            return self.states[0]
        if state.startswith('i'):
            return self.states[1]
        return self.states

    def get_name(self):
        return 'BDE'

    def get_unsampled_p(self, time, root_state=None, **kwargs):
        """
        The probability for a tree with a root in a specified state to evolve unsampled over specified time,
        given the rates.

        :param time: time
        :type time: float
        :param root_state: (optional) root state, if None, the root state will be averaged using equilibrium frequencies
        :type root_state: treesimulator.models.State

        :return: probability for a tree to evolve unsampled over time t
        :rtype: float
        """

        c = 1 / 2
        for _ in range(2):
             c = self.__get_unsampled_p_const(time / 2, self.states[0], c=c)

        return self.__get_unsampled_p_const(time, root_state, c=c)

    def get_prob(self, time, T, start_state, end_state, **kwargs):
        if EXPOSED == end_state.name:
            return 0

        m = self.rates[0, 0]
        lambda_i = self.rates[1, 1]
        sigma_e, sigma_i = self.sigmas

        avg_u_i = self.get_avg_unsampled_p(time + T, T, self.states[1])
        avg_u_e = self.get_avg_unsampled_p(time + T, T, self.states[0])
        if INFECTED == start_state.name:
            same_patient_p = np.exp(-(sigma_i - lambda_i * avg_u_i) * time)
            res = same_patient_p
            paths = [[1, 0, 1], [1, 1, 0, 1], [1, 0, 1, 1], [1, 1, 0, 1, 1], [1, 0, 1, 1, 1], [1, 0, 1, 0, 1], [1, 0, 1, 0, 1], [1, 1, 1, 0, 1]]
            for path in paths:
                p = get_convolution(time, np.around(self.sigmas[path], 5))
                p *= np.power(m, sum(1 for _ in path if _ == 0))
                p *= np.power(avg_u_i, sum(1 for i in range(len(path) - 1) if path[i] == 1 and path[i + 1] == 0))
                p *= np.power(avg_u_e, sum(1 for i in range(len(path) - 1) if path[i] == 1 and path[i + 1] == 1))
                res += p
        else:
            transition_then_same_patient_p = m * np.exp(-(sigma_i - lambda_i * avg_u_i) * time) \
                                             * (np.exp(-(sigma_e - sigma_i + lambda_i * avg_u_i) * time) - 1) \
                                             / (sigma_i - sigma_e - lambda_i * avg_u_i)
            res = transition_then_same_patient_p
            paths = [[0, 1, 0, 1], [0, 1, 0, 1, 1], [0, 1, 1, 0, 1]]
            for path in paths:
                p = get_convolution(time, np.around(self.sigmas[path], 5))
                p *= np.power(m, sum(1 for _ in path if _ == 0))
                p *= np.power(avg_u_i, sum(1 for i in range(len(path) - 1) if path[i] == 1 and path[i + 1] == 0))
                p *= np.power(avg_u_e, sum(1 for i in range(len(path) - 1) if path[i] == 1 and path[i + 1] == 1))
                res += p
        return res

    def get_tree_log_likelihood(self, tree, root_state=None):
        """
        Computes the likelihoods of possible states for each non-leaf node of the given tree
        with the leaf states specified in their 'state_feature' features, as well as the log likelihood of the whole tree.
        :param root_state: state of the root, if known
        :param tree: the tree of interest (ete3 Tree)
        """
        STATE_PS = 'state_ps'

        lambda_i = self.rates[1, 1]
        psi_e, psi_i = self.rates[1, :]
        p = self.ps[1]

        res = len(tree) * (np.log(psi_i) + np.log(p)) + (len(tree) - 1) * np.log(lambda_i)
        for node in tree.traverse("postorder"):
            if node.is_leaf():
                # already added the sampling rate contribution above
                pass
            elif len(node.children) != 2:
                raise ValueError(
                    'Your tree is not binary (e.g. see node {}) and it causes ML calculation problems!'.format(
                        node.name))
            else:
                left, right = node.get_children()
                # probabilities of evolving the left branch from each state
                left_ps, right_ps = getattr(left, STATE_PS), getattr(right, STATE_PS)
                res += np.log(left_ps.dot(right_ps[[1, 0]]))
            p_values = np.array([self.get_prob(time=node.dist, T=getattr(node, 'T'),
                                               start_state=state, end_state=self.states[1])
                                 for state in self.states])
            # if we got an overflow somewhere, and therefore a negative probability, let's reset it to zero
            if np.any(p_values < 0):
                logging.error('Got an over/underflow: {}'.format(p_values))
            p_values = np.fmax(p_values, 0)

            # If the likelihood of this subtree is 0, the likelihood of the whole tree will be 0
            if np.all(p_values == 0):
                logging.error('P_values: {}'.format(p_values))
                logging.error('An unlikely subtree: {}'
                              .format(node.get_ascii(attributes=['state', 'dist', STATE_PS])))
                if hasattr(tree, STATE_PS):
                    delattr(tree, STATE_PS)
                return None

            node.add_feature(STATE_PS, p_values)

        state_ps = getattr(tree, STATE_PS)
        res += np.log(self.rates[-1, :].dot(state_ps) if root_state is None else state_ps[root_state.index])
        return res

    def __get_unsampled_p_const(self, time, root_state=None, c=0, **kwargs):
        """
        The differential equations are the following (mu == m; lambda_i == l, psi_e == s1, psi_i == s2):
        (1) x'(t) = -(m + s1) x(t) + m y(t) + s1
        (2) y'(t) = -(l + s2) y(t) + s2 (1 - p) + l x(t) y(t)
        with the initial condition: x(0) = y(0) = 1,
        where x(t) = U_e(t), y(t) = U_i(t).

        Replacing x(t) in (2) with a constant 0 <= c <= 1, we can obtain the following solution (Mathematica):

        x(t) = (e^(t (-(-c l + l + s2))) (c^2 l^2 m e^(t (-c l + l + s2) + t (-m - s1)) + c^2 l^2 s1 e^(t (-c l + l + s2)) - 2 c l^2 m e^(t (-c l + l + s2) + t (-m - s1)) + l^2 m e^(t (-c l + l + s2) + t (-m - s1)) - 2 c l^2 s1 e^(t (-c l + l + s2)) + l^2 s1 e^(t (-c l + l + s2)) + m^2 p s2 e^(t (-c l + l + s2)) - m^2 s2 e^(t (-c l + l + s2)) + c l m^2 + m p s2^2 e^(t (-c l + l + s2) + t (-m - s1)) - c l m p s2 e^(t (-c l + l + s2) + t (-m - s1)) + l m p s2 e^(t (-c l + l + s2) + t (-m - s1)) + m p s1 s2 e^(t (-c l + l + s2)) - m p s2^2 e^(t (-c l + l + s2)) + c l m p s2 e^(t (-c l + l + s2)) - l m p s2 e^(t (-c l + l + s2)) + c l m s1 e^(t (-c l + l + s2)) - l m s1 e^(t (-c l + l + s2)) - c l m s2 e^(t (-c l + l + s2) + t (-m - s1)) + l m s2 e^(t (-c l + l + s2) + t (-m - s1)) - 2 m s1 s2 e^(t (-c l + l + s2)) + c l m s1 + m s2^2 e^(t (-c l + l + s2)) - c l m s2 e^(t (-c l + l + s2)) + l m s2 e^(t (-c l + l + s2)) + c l s1^2 e^(t (-c l + l + s2)) - l s1^2 e^(t (-c l + l + s2)) - s1^2 s2 e^(t (-c l + l + s2)) + s1 s2^2 e^(t (-c l + l + s2)) - 2 c l s1 s2 e^(t (-c l + l + s2)) + 2 l s1 s2 e^(t (-c l + l + s2)) - l m^2 - l m s1 - m^2 p s2 - m p s1 s2))/((m + s1) (c l - l - s2) (c l - l + m + s1 - s2))
        y(t) = (e^(t (-(-c l + l + s2))) (p s2 e^(t (-c l + l + s2)) - s2 e^(t (-c l + l + s2)) + c l - l - p s2))/(c l - l - s2)

        :param time:
        :param root_state:
        :param kwargs:
        :return:
        """
        m = self.rates[0, 0]
        l = self.rates[1, 1]
        s1, s2 = self.rates[2, :]
        p = self.ps[1]

        def __U_e(t):
            return (np.exp(t * (-(-c * l + l + s2))) * (np.power(c, 2) * np.power(l, 2) * m * np.exp(t * (-c * l + l + s2) + t * (-m - s1)) + np.power(c, 2) * np.power(l, 2) * s1 * np.exp(t * (-c * l + l + s2)) - 2 * c * np.power(l, 2) * m * np.exp(t * (-c * l + l + s2) + t * (-m - s1)) + np.power(l, 2) * m * np.exp(t * (-c * l + l + s2) + t * (-m - s1)) - 2 * c * np.power(l, 2) * s1 * np.exp(t * (-c * l + l + s2)) + np.power(l, 2) * s1 * np.exp(t * (-c * l + l + s2)) + np.power(m, 2) * p * s2 * np.exp(t * (-c * l + l + s2)) - np.power(m, 2) * s2 * np.exp(t * (-c * l + l + s2)) + c * l * np.power(m, 2) + m * p * np.power(s2, 2) * np.exp(t * (-c * l + l + s2) + t * (-m - s1)) - c * l * m * p * s2 * np.exp(t * (-c * l + l + s2) + t * (-m - s1)) + l * m * p * s2 * np.exp(t * (-c * l + l + s2) + t * (-m - s1)) + m * p * s1 * s2 * np.exp(t * (-c * l + l + s2)) - m * p * np.power(s2, 2) * np.exp(t * (-c * l + l + s2)) + c * l * m * p * s2 * np.exp(t * (-c * l + l + s2)) - l * m * p * s2 * np.exp(t * (-c * l + l + s2)) + c * l * m * s1 * np.exp(t * (-c * l + l + s2)) - l * m * s1 * np.exp(t * (-c * l + l + s2)) - c * l * m * s2 * np.exp(t * (-c * l + l + s2) + t * (-m - s1)) + l * m * s2 * np.exp(t * (-c * l + l + s2) + t * (-m - s1)) - 2 * m * s1 * s2 * np.exp(t * (-c * l + l + s2)) + c * l * m * s1 + m * np.power(s2, 2) * np.exp(t * (-c * l + l + s2)) - c * l * m * s2 * np.exp(t * (-c * l + l + s2)) + l * m * s2 * np.exp(t * (-c * l + l + s2)) + c * l * np.power(s1, 2) * np.exp(t * (-c * l + l + s2)) - l * np.power(s1, 2) * np.exp(t * (-c * l + l + s2)) - np.power(s1, 2) * s2 * np.exp(t * (-c * l + l + s2)) + s1 * np.power(s2, 2) * np.exp(t * (-c * l + l + s2)) - 2 * c * l * s1 * s2 * np.exp(t * (-c * l + l + s2)) + 2 * l * s1 * s2 * np.exp(t * (-c * l + l + s2)) - l * np.power(m, 2) - l * m * s1 - np.power(m, 2) * p * s2 - m * p * s1 * s2))/((m + s1) * (c * l - l - s2) * (c * l - l + m + s1 - s2))

        def __U_i(t):
            return (np.exp(t * (-(-c * l + l + s2))) * (p * s2 * np.exp(t * (-c * l + l + s2)) - s2 * np.exp(t * (-c * l + l + s2)) + c * l - l - p * s2))/(c * l - l - s2)

        if root_state is None:
            return self.rates[-1, :].dot([__U_e(time), __U_i(time)])
        if EXPOSED == root_state.name:
            return __U_e(time)
        if INFECTED == root_state.name:
            return __U_i(time)
        raise ValueError('Unknown root state {}'.format(root_state.name))
