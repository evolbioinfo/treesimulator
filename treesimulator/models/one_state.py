import numpy as np

from tree.likelihood import get_unsampled_p
from treesimulator.models import Model, SAMPLING, TRANSMISSION, TRANSITION, State

STATIONARY_DIST = np.array([1])

MU_RATES = np.array([0])

S = 's'


class OneStateModel(Model):

    def num_params(self, type=None):
        if SAMPLING == type:
            return 1
        if TRANSMISSION == type:
            return 1
        if TRANSITION == type:
            return 0
        return 2

    def _get_states(self):
        s = State(name=S, index=0)
        s.recipient = s
        return np.array([s])

    def params2rates(self, ps, sampled_pis=None, **kwargs):
        self.rates[1:3, 0] = ps
        self.rates[-1, 0] = 1

    def map_state(self, state):
        return self.states[0]

    def get_name(self):
        return 'one-state'

    def get_unsampled_p(self, time, root_state=None, **kwargs):
        """
        The probability for a tree with a transmission rate r and sampling rate s to evolve unsampled over time t:
        E(r, s, t) = (r + s) / (r + s e^{(r + s) t})

        :param time: time
        :type time: float
        :param root_state: (optional) root state, if None, the root state will be averages using equilibrium frequencies
        :type root_state: treesimulator.models.State

        :return: probability for a tree to evolve unsampled over time t
        :rtype: float
        """
        return get_unsampled_p(time, *self.rates[1:-1, 0])

    def get_avg_unsampled_p(self, time_start, time_end, root_state=None, **kwargs):
        """
        Calculates the probability for a tree to evolve unsampled, averaged over times between time_end and time_start,
        by dividing the integral of unsampled prob over that interval by the interval's length:
        integral(r + s)/(r + s exp((r + s) t)) dt = (t (r + s) - log(s e^(t (r + s)) + r))/r + constant.

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

        r, s = self.rates[1:-1, 0]

        def integral(t):
            res = (t * (r + s) - np.log(s * np.exp(t * (r + s)) + r)) / r
            if np.isnan(res):
                print(t, r, s)
            return res

        return (integral(time_end) - integral(time_start)) / (time_end - time_start)
