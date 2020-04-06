import logging

import numpy as np

from likelihood.convolution import get_convolution
from treesimulator.models import Model, SAMPLING, TRANSMISSION, TRANSITION, State

M = 5
# 0.1%
FRACTION_OF_SIGNIFICANT_ADDITION = 0.01

EXPOSED = 'e'
INFECTED = 'i'


class BirthDeathExposedModel(Model):

    def __init__(self, p=0.5, *args, **kwargs):
        Model.__init__(self, ps=[0, p], minus_avg_sigma=False, *args, **kwargs)

    def num_params(self, type=None):
        if SAMPLING == type:
            return 1
        if TRANSMISSION == type:
            return 1
        if TRANSITION == type:
            return 1
        return 3

    def clone(self):
        model = BirthDeathExposedModel(p=self.ps[1])
        model.params2rates([self.rates[0, 0], self.rates[1, 1], self.rates[2, 1]])
        return model

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
        mu, lambda_i, psi_i = ps
        self.rates[0, 0] = mu
        self.rates[1, 1] = lambda_i
        self.rates[2, 1] = psi_i
        self.rates[3, :] = self.get_sd()

    def get_bounds(self, lb, ub, **kwargs):
        return np.array([[lb, ub]] * self.num_params())

    def _get_sd_formulas(self):
        mu = self.rates[0, 0]
        lambda_i, psi_i = self.rates[1: -1, 1]

        p = [
            lambda_i - psi_i,
            psi_i + mu,
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

    def get_avg_unsampled_p(self, time_start, time_end, root_state=None, m=M, **kwargs):
        """
        Calculates the probability for a tree to evolve unsampled, averaged over times between time_end and time_start.

        :param time_start: maximum time
        :type time_start: float
        :param time_end: minimum time
        :type time_end: float
        :param root_state: (optional) root state, if None, the root state will be averaged using equilibrium frequencies
        :type root_state: treesimulator.models.State
        :return: probability for a tree to evolve unsampled, averaged over times between time_end and time_start
        :rtype: float
        """
        if time_end == time_start:
            return 1

        time = time_end - time_start
        times = [time_start + (i + 1 / 2) * time / (m + 1) for i in range(0, m + 1)]
        return sum(self.get_unsampled_p(t, root_state, **kwargs) for t in times) / len(times)

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

        def get_iterative_c(t):
            p = max(4, int(np.log2(16 * max(1, t * max([self.rates[0, 0], self.rates[1, 1], self.rates[2, 1]])))))
            c = 1
            tau = t / (2 << (p - 1))
            while tau < t:
                c = self.__get_unsampled_p_const(tau, self.states[0], c=c)
                tau *= 2
            return c
        c = get_iterative_c(time)
        return self.__get_unsampled_p_const(time, root_state, c=c)

    def get_prob(self, time, T, start_state, bother=True, **kwargs):
        n = max(1, int(time * 4 * max([self.rates[0, 0], self.rates[1, 1], self.rates[2, 1]])))
        if bother and n > 1:
            t = time / n
            ttn = T + time - t
            return sum(self._get_prob(t, ttn, start_state, state) * self.get_prob(time - t, T, state)
                       for state in self.states)
        else:
            return self._get_prob(time, T, start_state, self.states[1])

    def _get_prob(self, time, T, start_state, end_state, **kwargs):
        mu = self.rates[0, 0]
        lambda_i = self.rates[1, 1]
        s_i = self.rates[2, 1]

        avg_u_i = self.get_avg_unsampled_p(time + T, T, self.states[1])
        avg_u_e = self.get_avg_unsampled_p(time + T, T, self.states[0])

        # def get_sigma(avg_u_e):
        #     return lambda_i + s_i - lambda_i * avg_u_e
        #
        # mu_la_la = mu * lambda_i * lambda_i
        #
        # def get_avg_us(m):
        #     times = [i * time / m for i in range(1, m + 1)]
        #     return [(self.get_avg_unsampled_p((time - t - time / m) + T, (time - t) + T, self.states[0]),
        #              self.get_avg_unsampled_p((time - t - time / m) + T, (time - t) + T, self.states[1]))
        #             for t in times]

        sigma = lambda_i + s_i - lambda_i * avg_u_e
        coeff = mu * lambda_i * lambda_i * avg_u_i * avg_u_e

        if INFECTED == start_state.name and INFECTED == end_state.name:
            sigmas = [sigma]
            multiplier = 1
            res = np.exp(-sigma * time)
            i = 0
            while True:
                i += 1
                multiplier *= coeff
                sigmas.extend([sigma, mu])
                addition = multiplier * get_convolution(time, sigmas)
                res += addition
                if addition <= 0 or addition < res * FRACTION_OF_SIGNIFICANT_ADDITION:
                    if i > 15:
                        print('I->I: after {} iterations: from {} to {}\n'.format(i, np.exp(-sigma * time), res))
                    break
        elif INFECTED == start_state.name and EXPOSED == end_state.name:
            sigmas = [mu, sigma]
            multiplier = lambda_i * avg_u_i
            res = 0
            i = 0
            while True:
                i += 1
                addition = multiplier * get_convolution(time, sigmas)
                res += addition
                if addition <= 0 or addition < res * FRACTION_OF_SIGNIFICANT_ADDITION:
                    if i > 15:
                        print('I->E: after {} iterations: {}\n'.format(i, res))
                    break
                multiplier *= coeff
                sigmas.extend([sigma, mu])
        elif EXPOSED == start_state.name and INFECTED == end_state.name:
            sigmas = [sigma, mu]
            multiplier = mu
            res = 0
            i = 0
            while True:
                i += 1
                addition = multiplier * get_convolution(time, sigmas)
                res += addition
                if addition <= 0 or addition < res * FRACTION_OF_SIGNIFICANT_ADDITION:
                    if i > 15:
                        print('E->I: after {} iterations: {}\n'.format(i, res))
                    break
                multiplier *= coeff
                sigmas.extend([sigma, mu])
        elif EXPOSED == start_state.name and EXPOSED == end_state.name:
            res = np.exp(-mu * time)
            multiplier = mu * lambda_i * avg_u_i
            sigmas = [mu, mu]
            i = 0
            while True:
                i += 1
                addition = multiplier * get_convolution(time, sigmas)
                res += addition
                if addition <= 0 or addition < res * FRACTION_OF_SIGNIFICANT_ADDITION:
                    if i > 15:
                        print('E->E: after {} iterations: from {} to {}\n'.format(i, np.exp(-mu * time), res))
                    break
                multiplier *= coeff
                sigmas.extend([sigma, mu])
        return res

    def get_tree_log_likelihood(self, tree, root_state=None):
        """
        Computes the likelihoods of possible states for each non-leaf node of the given tree
        with the leaf states specified in their 'state_feature' features, as well as the log likelihood of the whole tree.
        :param root_state: state of the root, if known
        :param tree: the tree of interest (ete3 Tree)
        """
        STATE_PS = 'state_ps_{}'.format(self.__hash__())

        lambda_i, psi_i = self.rates[1:-1, 1]
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
            p_values = np.array([self.get_prob(time=node.dist, T=getattr(node, 'T'), start_state=state)
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
        The differential equations are the following (mu == m; lambda_i == l, psi_i == s):
        (1) x'(t) = -m x(t) + m y(t)
        (2) y'(t) = -(l + s) y(t) + s (1 - p) + l x(t) y(t)
        with the initial condition: x(0) = y(0) = 1,
        where x(t) = U_e(t), y(t) = U_i(t).

        Replacing x(t) in (2) with a constant 0 <= c <= 1, we can obtain the following solution (Mathematica):

        x(t) = (exp(-t ((c - 1) l + m - s) - t (-c l + l + s) - m t) (c^2 l^2 exp(t ((c - 1) l + m - s) + t (-c l + l + s)) - 2 c l^2 exp(t ((c - 1) l + m - s) + t (-c l + l + s)) + l^2 exp(t ((c - 1) l + m - s) + t (-c l + l + s)) + p s^2 exp(t ((c - 1) l + m - s) + t (-c l + l + s)) - p s^2 exp(t ((c - 1) l + m - s) + t (-c l + l + s) + m t) - c l p s exp(t ((c - 1) l + m - s) + t (-c l + l + s)) + l p s exp(t ((c - 1) l + m - s) + t (-c l + l + s)) + c l p s exp(t ((c - 1) l + m - s) + t (-c l + l + s) + m t) - l p s exp(t ((c - 1) l + m - s) + t (-c l + l + s) + m t) + m p s exp(t ((c - 1) l + m - s) + t (-c l + l + s) + m t) - m p s exp(t ((c - 1) l + m - s) + 2 t (-c l + l + s)) + s^2 exp(t ((c - 1) l + m - s) + t (-c l + l + s) + m t) - c l s exp(t ((c - 1) l + m - s) + t (-c l + l + s)) + l s exp(t ((c - 1) l + m - s) + t (-c l + l + s)) - c l s exp(t ((c - 1) l + m - s) + t (-c l + l + s) + m t) + l s exp(t ((c - 1) l + m - s) + t (-c l + l + s) + m t) - m s exp(t ((c - 1) l + m - s) + t (-c l + l + s) + m t) + m s exp(t ((c - 1) l + m - s) + 2 t (-c l + l + s)) - m p s e^(t ((c - 1) l + m - s) + m t) + m p s e^(t (-c l + l + s) + m t) + c l m e^(t ((c - 1) l + m - s) + m t) - l m e^(t ((c - 1) l + m - s) + m t) - m s e^(t (-c l + l + s) + m t)))/((c l - l - s) (c l - l + m - s))
        y(t) = (e^(t (-(-c l + l + s))) (p s e^(t (-c l + l + s)) - s e^(t (-c l + l + s)) + c l - l - p s))/(c l - l - s)

        :param time:
        :param root_state:
        :param kwargs:
        :return:
        """
        m = self.rates[0, 0]
        l, s = self.rates[1: -1, 1]
        p = self.ps[1]

        def __U_e(t):
            res = (np.exp(-t * ((c - 1) * l + m - s) - t * (-c * l + l + s) - m * t)
                    * (np.power(c, 2) * np.power(l, 2) * np.exp(t * ((c - 1) * l + m - s) + t * (-c * l + l + s)) - 2 * c * np.power(l, 2) * np.exp(t * ((c - 1) * l + m - s) + t * (-c * l + l + s)) + np.power(l, 2) * np.exp(t * ((c - 1) * l + m - s) + t * (-c * l + l + s)) + p * np.power(s, 2) * np.exp(t * ((c - 1) * l + m - s) + t * (-c * l + l + s)) - p * np.power(s, 2) * np.exp(t * ((c - 1) * l + m - s) + t * (-c * l + l + s) + m * t) - c * l * p * s * np.exp(t * ((c - 1) * l + m - s) + t * (-c * l + l + s)) + l * p * s * np.exp(t * ((c - 1) * l + m - s) + t * (-c * l + l + s)) + c * l * p * s * np.exp(t * ((c - 1) * l + m - s) + t * (-c * l + l + s) + m * t) - l * p * s * np.exp(t * ((c - 1) * l + m - s) + t * (-c * l + l + s) + m * t) + m * p * s * np.exp(t * ((c - 1) * l + m - s) + t * (-c * l + l + s) + m * t) - m * p * s * np.exp(t * ((c - 1) * l + m - s) + 2 * t * (-c * l + l + s)) + np.power(s, 2) * np.exp(t * ((c - 1) * l + m - s) + t * (-c * l + l + s) + m * t) - c * l * s * np.exp(t * ((c - 1) * l + m - s) + t * (-c * l + l + s)) + l * s * np.exp(t * ((c - 1) * l + m - s) + t * (-c * l + l + s)) - c * l * s * np.exp(t * ((c - 1) * l + m - s) + t * (-c * l + l + s) + m * t) + l * s * np.exp(t * ((c - 1) * l + m - s) + t * (-c * l + l + s) + m * t) - m * s * np.exp(t * ((c - 1) * l + m - s) + t * (-c * l + l + s) + m * t) + m * s * np.exp(t * ((c - 1) * l + m - s) + 2 * t * (-c * l + l + s)) - m * p * s * np.exp(t * ((c - 1) * l + m - s) + m * t) + m * p * s * np.exp(t * (-c * l + l + s) + m * t) + c * l * m * np.exp(t * ((c - 1) * l + m - s) + m * t) - l * m * np.exp(t * ((c - 1) * l + m - s) + m * t) - m * s * np.exp(t * (-c * l + l + s) + m * t)))\
                   / ((c * l - l - s) * (c * l - l + m - s))
            return min(1, max(0, res))

        def __U_i(t):
            res = (np.exp(t * (-(-c * l + l + s)))
                    * (p * s * np.exp(t * (-c * l + l + s)) - s * np.exp(t * (-c * l + l + s)) + c * l - l - p * s))\
                   / (c * l - l - s)
            return min(1, max(0, res))
            
        if root_state is None:
            return self.rates[-1, :].dot([__U_e(time), __U_i(time)])
        if EXPOSED == root_state.name:
            return __U_e(time)
        if INFECTED == root_state.name:
            return __U_i(time)
        raise ValueError('Unknown root state {}'.format(root_state.name))
