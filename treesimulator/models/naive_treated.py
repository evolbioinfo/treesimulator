import numpy as np

from tree.likelihood import get_unsampled_p, get_convolution, get_avg_unsampled_p
from treesimulator.models import Model, SAMPLING, TRANSMISSION, TRANSITION, State

NAIVE = 'n'
TREATED = 't'


class NaiveTreatedModel(Model):

    def num_params(self, type=None):
        if SAMPLING == type:
            return 2
        if TRANSMISSION == type:
            return 2
        if TRANSITION == type:
            return 1
        return 5

    def _get_states(self):
        treated_state = State(name=TREATED, index=1)
        naive_state = State(name=NAIVE, index=0, next_state=treated_state)

        naive_state.recipient = naive_state
        treated_state.recipient = naive_state

        return np.array([naive_state, treated_state])

    def params2rates(self, ps, sampled_pis=None, **kwargs):
        """
        Converts parameters into a rate array.

        :param ps: parameters, in the following order:
            treatment_elsewhere, lambda_n, lambda_t, treatment_first, treatment_second
        """
        if sampled_pis is None:
            mu_n, lambda_n, lambda_t, psi_n, psi_t = ps
            self.rates[0, 0] = mu_n
            self.rates[1, :] = lambda_n, lambda_t
            self.rates[2, :] = psi_n, psi_t
            self.rates[3, :] = self.get_sd()
        else:
            self.rates[1, :] = ps[:2]
            self.rates[2, :] = ps[2:]
            # pi_n = spi_n * psi_t / (spi_n * psi_t + spi_t * psi_n)
            pi_n = sampled_pis[0] * self.rates[2, 1] / sampled_pis.dot(self.rates[2, [1, 0]])
            pi_t = 1 - pi_n
            self.rates[-1, :] = [pi_n, pi_t]
            # mu_n = pi_t (avg_lambda - avg_psi + psi_t) / pi_n
            self.rates[0, 0] = pi_t / pi_n * (self.rates[1, :].dot(self.rates[-1, :])
                                              - self.rates[2, :].dot(self.rates[-1, :]) + self.rates[2, 1])

    def get_bounds(self, lb, ub, sampled_pis=None, **kwargs):
        return np.array([[lb, ub]] * (5 if sampled_pis is None else 4))

    def _get_sd_formulas(self):
        mu = self.rates[0, 0]
        l_n, l_t = self.rates[1, :]
        s_n, s_t = self.rates[2, :]

        p = [
            -s_t + s_n - l_n + l_t,
            l_n - s_n + mu + s_t,
            -mu
        ]
        roots = np.roots(p)
        roots = roots[np.where(np.isreal(roots) & (~np.isnan(roots)) & (roots >= 0) & (roots <= 1))]
        pi_t = roots[0]

        return np.array([1 - pi_t, pi_t])

    def map_state(self, state):
        state = str(state)
        if state.startswith('n'):
            return self.states[0]
        if state.startswith('t'):
            return self.states[1]
        return self.states

    def get_name(self):
        return 'naive-treated'

    def get_unsampled_p(self, time, root_state=None, simple=False, **kwargs):
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

        def __e_n(t):
            mu = self.rates[0, 0]
            l_n = self.rates[1, 0]
            s_n = self.rates[2, 0]
            e_n_only = get_unsampled_p(time, l_n, mu + s_n)
            sigma_n, sigma_t = self.rates[:-1, :].sum(axis=0)
            e_n_avg = get_avg_unsampled_p(0, t, l_n, mu + s_n)
            e_t_avg = self.get_avg_unsampled_p(0, t, self.states[1], simple=True)
            e_n_transm_n_only_mut_t = \
                mu * l_n * e_n_avg * e_t_avg * get_convolution(t, [sigma_n] * 2)
            e_mut_t_etc = mu * (1 - np.exp(-t * sigma_n)) / sigma_n * e_t_avg
            return e_n_only + e_n_transm_n_only_mut_t + e_mut_t_etc

        def __e_t(t):
            sigma_t = self.rates[:-1, 1].sum()
            e_t_only = np.exp(-t * sigma_t)
            e_n_avg = get_avg_unsampled_p(0, t, self.rates[1, 0], self.rates[[0, 2], 0].sum()) if simple \
                else self.get_avg_unsampled_p(0, t, self.states[0])
            l_t = self.rates[1, 1]
            e_t_transm_n_only = l_t * e_n_avg * get_convolution(t, [sigma_t] * 2)
            return e_t_only + e_t_transm_n_only

        if root_state is None:
            return self.rates[-1, :].dot([__e_n(time), __e_t(time)])
        if NAIVE == root_state.name:
            return __e_n(time)
        if TREATED == root_state.name:
            return __e_t(time)
        raise ValueError('Unknown root state {}'.format(root_state.name))
