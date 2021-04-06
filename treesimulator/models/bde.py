import logging

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special._ufuncs import gammainc, gamma

from likelihood.convolution import get_convolution
from treesimulator.models import Model, SAMPLING, TRANSMISSION, TRANSITION, State

M = 5
# 10%
FRACTION_OF_SIGNIFICANT_ADDITION = 0.1

EXPOSED = 'e'
INFECTED = 'i'


class BirthDeathExposedModel(Model):

    def __init__(self, p=0.5, *args, **kwargs):
        Model.__init__(self, ps=[0, p], minus_avg_sigma=False, *args, **kwargs)
        self.ues, self.uis = np.ones(1000, dtype=np.float128), np.ones(1000, dtype=np.float128)
        self.dt = None

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

    def check_rates(self):
        return np.all(self.rates >= 0) #and np.all(rates[-1, :] <= 1)


    def get_best_params(self, forest, num_trees, lb=1 / 100, ub=2):
        """
        Computes the best (in terms of max likelihood) transition and transmission rates
        for the given tree(s) with the leaf states specified in their 'state_feature' features.
        :param forest: list of trees of interest (ete3 Tree)
        :return: a tuple (transition_rates, transmission_rates), each in a form {from_state_i: (to_state_i, rate_i)},
        with the best estimated parameters, or (None, None) if no combination of parameters was found.
        """

        total_time_till_now = next((getattr(tree, 'T') + tree.dist) for tree in forest)

        logging.info('# trees in the forest = {} (out of {}).'.format(len(forest), num_trees))

        bounds = self.get_bounds(lb, ub)

        def get_v(ps):
            if np.any(pd.isnull(ps)):
                return np.nan
            self.params2rates(ps)
            self.calc_us(total_time_till_now)
            logging.info('rates:\n{rates}'.format(params=ps, rates=self.rates))

            if not self.check_rates():
                res = -np.inf
            else:
                res = self.get_forest_log_likelihood(forest=forest, num_trees=num_trees,
                                                     total_time_till_now=total_time_till_now)

            logging.info('log likelihood: {result}\n'.format(result=res))
            # containing nans means that the rate parameters contain nans
            return np.inf if pd.isnull(res) else -res

        for _ in range(5):
            vs = np.around(np.random.uniform(bounds[:, 0], bounds[:, 1]), 2)
            fres = minimize(get_v, x0=vs, method='L-BFGS-B', bounds=bounds, options={'eps': 1 / 365.})
            if fres.success and not np.any(np.isnan(fres.x)):
                self.params2rates(fres.x)
                logging.info(fres)
                return True
        return False

    def calc_us(self, T):
        mu, lmbda, psi, p = self.rates[0, 0], self.rates[1, 1], self.rates[2, 1], self.ps[1]
        self.dt = T / 1000
        change_prob = mu * self.dt
        transmission_prob = lmbda * self.dt
        removal_prob = psi * self.dt
        unsampled_removal_prob = removal_prob * (1 - p)
        for i in range(1000 - 1):
            self.ues[i + 1] = (self.ues[i] + change_prob * self.uis[i]) / (1 + change_prob)
            self.uis[i + 1] = (self.uis[i] + unsampled_removal_prob) \
                              / (1 + removal_prob + transmission_prob * (1 - self.ues[i + 1]))

    def calc_ps(self, t, T):
        mu, lmbda, psi, p = self.rates[0, 0], self.rates[1, 1], self.rates[2, 1], self.ps[1]
        dt = min(self.dt / 100, t)
        tau = 0
        pei, pie = 0, 0
        pii, pee = 1, 1
        while tau < t:
            pei_next = pei + dt * mu * (pii - pei)
            pee_next = pee + dt * mu * (pie - pee)
            pii_next = pii + dt * (-(lmbda + mu) * pii + lmbda * (pii * self.get_ue(t - tau + T) + pei * self.get_ui(t - tau + T)))
            pie_next = pie + dt * (-(lmbda + mu) * pie + lmbda * (pie * self.get_ue(t - tau + T) + pee * self.get_ui(t - tau + T)))
            tau = min(t, tau + dt)
            pei, pie, pii, pee = pei_next, pie_next, pii_next, pee_next
        return pei, pii

    def get_ue(self, t):
        i = int(t / self.dt)
        if i < 999:
            return self.ues[i] + (self.ues[i + 1] - self.ues[i]) * (t - self.dt * i) / self.dt
        else:
            return self.ues[-1]

    def get_ui(self, t):
        i = int(t / self.dt)
        if i < 999:
            return self.uis[i] + (self.uis[i + 1] - self.uis[i]) * (t - self.dt * i) / self.dt
        else:
            return self.uis[-1]

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
        if root_state is None:
            return self.rates[-1, :].dot([self.get_ue(time), self.get_ui(time)])
        return self.get_ui(time) if INFECTED == root_state.name else self.get_ue(time)


        # def get_iterative_c(t):
        #     p = max(4, int(np.log2(16 * max(1, t * max([self.rates[0, 0], self.rates[1, 1], self.rates[2, 1]])))))
        #     c = 1
        #     tau = t / (2 << (p - 1))
        #     while tau < t:
        #         c = self.__get_unsampled_p_const(tau, self.states[0], c=c)
        #         tau *= 2
        #     return c
        # c = get_iterative_c(time)
        # return self.__get_unsampled_p_const(time, root_state, c=c)

    def get_prob(self, time, T, start_state, frac_sign=0, min_i=64, **kwargs):
        # n = max(1, int(time * 2 * max([self.rates[0, 0], self.rates[1, 1], self.rates[2, 1]])))
        # if bother and n > 1:
        #     t = time / n
        #     ttn = T + time - t
        #     return sum(self._get_prob(t, ttn, start_state, state, frac_sign=frac_sign) * self.get_prob(time - t, T, state)
        #                for state in self.states)
        # else:
        # if gamma:
        #     return self.__get_prob(time, T, start_state, self.states[1], frac_sign=frac_sign, max_i=max_i)
        # else:
        #     return self._get_prob(time, T, start_state, self.states[1], frac_sign=frac_sign, min_i=max_i)

        return self.__get_prob_formulas(time, T, start_state, frac_sign=frac_sign, min_i=min_i)

    def __get_prob_formulas(self, time, T, start_state, frac_sign=FRACTION_OF_SIGNIFICANT_ADDITION, min_i=20, **kwargs):
        mu = self.rates[0, 0]
        lambda_i = self.rates[1, 1]
        s_i = self.rates[2, 1]

        def get_spaced_times(j):
            return ((T + (p - 1 / 2) * time / j) for p in range(1, j + 1))

        def get_avg_u_e_s(times):
            return (self.get_unsampled_p(t, self.states[0]) for t in times)

        def get_avg_u_i_s(times):
            return (self.get_unsampled_p(t, self.states[1]) for t in times)

        def get_i_sigma(k, j):
            return lambda_i + s_i - lambda_i * self.get_ue(T + (k - 1 / 2) * time / (j + 1))

        def get_i_exp(k, j):
            return np.exp(-get_i_sigma(k, j) * time)

        e_exp = np.exp(-mu * time)

        mu_lambda = mu * lambda_i

        if EXPOSED == start_state.name:
            res = (e_exp - get_i_exp(k=1, j=0)) / (get_i_sigma(k=1, j=0) - mu)
            j = 1
            multiplier = 1
            while True:
                multiplier *= mu_lambda
                crazy_sum = np.zeros(j + 1, dtype=np.float128)
                ues = np.array(list(get_avg_u_e_s(get_spaced_times(j + 1))), dtype=np.float128)
                for k in range(1, j + 2):
                    fact_sum = np.ones(j + 1, dtype=np.float128)
                    i_sigma = get_i_sigma(k, j)
                    for p in range(1, j + 1):
                        fact_sum[p] = fact_sum[p - 1] * i_sigma * time / p
                    e_sums = np.ones(j + 1, dtype=np.float128) * ues[k - 1] - ues
                    e_sums[k - 1] = 1
                    fact_sum.sort()

                    crazy_sum[k - 1] = (get_i_exp(k, j) - e_exp * fact_sum.sum()) / np.power(i_sigma, j) / prod(e_sums)
                times_j = np.array(list(get_spaced_times(j)), dtype=np.float128)
                crazy_sum.sort()
                addition = multiplier * crazy_sum.sum() * prod(get_avg_u_e_s(times_j)) * prod(get_avg_u_i_s(times_j))

                if addition < 0 or (res + addition) > 1 / mu or addition < res * frac_sign or min_i < j:
                    # print('E: not adding {} to {} for j={}'.format(addition * mu, res * mu, j))
                    break
                j += 1
                res += addition
            return mu * res

        elif INFECTED == start_state.name:
            res = get_i_exp(k=1, j=0)
            j = 1
            multiplier = 1
            while True:
                multiplier *= mu_lambda
                crazy_sum = np.zeros(j + 1, dtype=np.float128)
                ues = np.array(list(get_avg_u_e_s(get_spaced_times(j + 1))), dtype=np.float128)
                for k in range(1, j + 2):
                    fact_sum = np.ones(j, dtype=np.float128)
                    i_sigma = get_i_sigma(k, j)
                    for p in range(1, j):
                        fact_sum[p] = fact_sum[p - 1] * i_sigma * time / p
                    e_sums = np.ones(j + 1, dtype=np.float128) * ues[k - 1] - ues
                    e_sums[k - 1] = 1
                    fact_sum.sort()

                    crazy_sum[k - 1] = (get_i_exp(k, j) - e_exp * fact_sum.sum()) / np.power(i_sigma, j - 1) / prod(e_sums)
                times_j = np.array(list(get_spaced_times(j)), dtype=np.float128)
                crazy_sum.sort()
                addition = multiplier * crazy_sum.sum() * prod(get_avg_u_e_s(times_j)) * prod(get_avg_u_i_s(times_j))

                if addition < 0 or (res + addition) > 1 or addition < res * frac_sign or min_i < j:
                    # print('I: not adding {} to {} for j={}'.format(addition, res, j))
                    break
                res += addition
                j += 1
            return res

    def _get_prob_gamma(self, time, T, start_state, end_state, frac_sign=FRACTION_OF_SIGNIFICANT_ADDITION, max_i=8, **kwargs):
        mu = self.rates[0, 0]
        lambda_i = self.rates[1, 1]
        s_i = self.rates[2, 1]

        def get_avg_u_e_s(m):
            if m == 0:
                return []
            times = (i * time / m for i in range(1, m + 1))
            return (self.get_unsampled_p((time - t - time / 2 / m) + T, self.states[0]) for t in times)

        def get_avg_u_i_s(m):
            if m == 0:
                return []
            times = (i * time / m for i in range(1, m + 1))
            return (self.get_unsampled_p((time - t - time / 2 / m) + T, self.states[1]) for t in times)

        if INFECTED == start_state.name and INFECTED == end_state.name:
            multiplier = 1
            res = np.exp(-(lambda_i + s_i - lambda_i * next(get_avg_u_e_s(1))) * time)
            i = 1
            while True:
                multiplier *= mu * lambda_i

                u_es = np.array(list(get_avg_u_e_s(i + 1)), dtype=np.float128)
                sigmas = lambda_i + s_i - lambda_i * u_es
                divisors = mu - lambda_i - s_i + lambda_i * u_es
                exps = np.exp(-sigmas * time)
                lambda_u_es = lambda_i * u_es
                prods = np.ones(i + 1)
                for k in range(i + 1):
                    p = lambda_u_es[k] - lambda_u_es
                    p[k] = 1
                    prods[k] = np.prod(p)
                gammas = gammainc(i, divisors * time)
                conv = gamma(i) * np.sum(exps * gammas / prods / np.power(divisors, i))

                addition = max(0, multiplier * prod(get_avg_u_e_s(i)) * prod(get_avg_u_i_s(i)) * conv)
                if addition > 1:
                    break
                res += addition
                i += 1
                multiplier /= i
                if addition < res * frac_sign or max_i < i:
                    break
        elif EXPOSED == start_state.name and INFECTED == end_state.name:
            multiplier = mu
            res = 0
            i = 0
            while True:
                u_es = np.array(list(get_avg_u_e_s(i + 1)), dtype=np.float128)
                sigmas = lambda_i + s_i - lambda_i * u_es
                divisors = mu - lambda_i - s_i + lambda_i * u_es
                exps = np.exp(-sigmas * time)
                lambda_u_es = lambda_i * u_es
                prods = np.ones(i + 1)
                for k in range(i + 1):
                    p = lambda_u_es[k] - lambda_u_es
                    p[k] = 1
                    prods[k] = np.prod(p)
                gammas = gammainc(i + 1, divisors * time)
                conv = gamma(i + 1) * np.sum(exps * gammas / prods / np.power(divisors, i + 1))

                addition = max(0, multiplier * prod(get_avg_u_e_s(i)) * prod(get_avg_u_i_s(i)) * conv)
                if addition > 1:
                    break
                res += addition
                i += 1
                multiplier *= mu * lambda_i / i
                if addition < res * frac_sign or max_i < i:
                    break
        return res

    def _get_prob(self, time, T, start_state, end_state, frac_sign=FRACTION_OF_SIGNIFICANT_ADDITION, min_i=1, **kwargs):
        mu = self.rates[0, 0]
        lambda_i = self.rates[1, 1]
        s_i = self.rates[2, 1]

        def get_sigma(avg_u_e):
            return lambda_i + s_i - lambda_i * avg_u_e

        mu_la_la = mu * lambda_i * lambda_i

        def get_avg_u_e_s(m):
            if m == 0:
                return []
            times = (i * time / m for i in range(1, m + 1))
            return (self.get_unsampled_p((time - t - time / 2 / m) + T, self.states[0]) for t in times)

        def get_avg_u_i_s(m):
            if m == 0:
                return []
            times = (i * time / m for i in range(1, m + 1))
            return (self.get_unsampled_p((time - t - time / 2 / m) + T, self.states[1]) for t in times)

        if INFECTED == start_state.name and INFECTED == end_state.name:
            multiplier = 1
            res = 0
            i = 0
            while True:
                sigmas = [get_sigma(u_e) for u_e in get_avg_u_e_s(i + 1)] + [mu] * i
                addition = max(0, multiplier * prod(get_avg_u_e_s(i)) * prod(get_avg_u_i_s(i)) \
                           * get_convolution(time, sigmas))
                if addition > 1:
                    break
                res += addition
                i += 1
                multiplier *= mu_la_la
                if addition < res * frac_sign or min_i < i:
                    break
        elif INFECTED == start_state.name and EXPOSED == end_state.name:
            multiplier = lambda_i
            res = 0
            i = 0
            while True:
                sigmas = [get_sigma(u_e) for u_e in get_avg_u_e_s(i + 1)] + [mu] * (i + 1)
                addition = max(0, multiplier * prod(get_avg_u_e_s(i)) * prod(get_avg_u_i_s(i + 1)) \
                           * get_convolution(time, sigmas))
                if addition > 1:
                    break
                res += addition
                i += 1
                multiplier *= mu_la_la
                if addition < res * frac_sign or min_i < i:
                    break
        elif EXPOSED == start_state.name and INFECTED == end_state.name:
            multiplier = mu
            res = 0
            i = 0
            while True:
                sigmas = [get_sigma(u_e) for u_e in get_avg_u_e_s(i + 1)] + [mu] * (i + 1)
                addition = max(0, multiplier * prod(get_avg_u_e_s(i)) * prod(get_avg_u_i_s(i)) \
                           * get_convolution(time, sigmas))
                if addition > 1:
                    break
                res += addition
                i += 1
                multiplier *= mu_la_la
                if addition < res * frac_sign or min_i < i:
                    break
        elif EXPOSED == start_state.name and EXPOSED == end_state.name:
            res = np.exp(-mu * time)

            multiplier = mu * lambda_i
            i = 1
            while True:
                sigmas = [get_sigma(u_e) for u_e in get_avg_u_e_s(i - 1)] + [mu] * (i + 1)
                addition = max(0, multiplier * prod(get_avg_u_e_s(i - 1)) * prod(get_avg_u_i_s(i)) \
                           * get_convolution(time, sigmas))
                if addition > 1:
                    break
                res += addition
                i += 1
                multiplier *= mu_la_la
                if addition < res * frac_sign or min_i < i:
                    break
        return res

    def get_forest_log_likelihood(self, forest, num_trees, total_time_till_now):
        res = 0

        num_hidden_trees = num_trees - len(forest)
        if num_hidden_trees:
            res += num_hidden_trees * np.log(self.get_unsampled_p(total_time_till_now))

        def work(tree):
            return self.get_tree_log_likelihood(tree=tree)

        # if len(forest) > 1:
        #     with ThreadPoolExecutor(max_workers=len(forest)) as executor:
        #         vs = list(executor.map(work, forest))
        # else:
        #     vs = [work(forest[0])]
        vs = [work(tree) for tree in forest]
        vs = np.array(vs, dtype=np.float128)
        if np.any(pd.isna(vs)):
            return np.nan

        return res + np.sum(vs)

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
                                 for state in self.states], dtype=np.float128)
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


def prod(values):
    res = 1
    for v in values:
        res *= v
    return res