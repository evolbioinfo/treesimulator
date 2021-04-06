import logging
from abc import abstractmethod

import numpy as np

DIGITS_TO_KEEP = 3

SAMPLING = 'sampling'
TRANSMISSION = 'transmission'
TRANSITION = 'transition'


def get_avg_rate(rates, weights):
    """
    Returns average rate for given parameters.
    :param rates: np.array of rates
    :param weights: np.array of weights
    :return: float, average rate
    """
    return rates.dot(weights)


def get_error_fractions(estimated_rates, rates):
    return np.nan_to_num((estimated_rates - rates) / rates)


def describe_error(label, estimated_rates, expected_rates):
    separator = '\t' if len(estimated_rates.shape) == 1 or estimated_rates.shape[1] == 1 else '\n'
    logging.info('Estimated {} are:{}{}'.format(label, separator, np.round(estimated_rates, DIGITS_TO_KEEP)))

    if expected_rates is not None:
        error_fractions = get_error_fractions(estimated_rates, expected_rates)
        logging.info('The error fractions in percentage:{}{}'
                     .format(separator, np.round(100 * error_fractions)))
        logging.info('\n==========================\n')


class Model(object):
    def __init__(self, minus_avg_sigma=True, ps=None, *args, **kwargs):
        super(Model, self).__init__()
        self.__states = self._get_states()
        self.__ps = np.array(ps) if ps is not None else np.ones(len(self.states), dtype=float)
        self.__rates = np.zeros(shape=(4, len(self.states)), dtype=np.float)
        self.__sigmas = None
        self.__avg_sigma = None
        self.__minus_avg = minus_avg_sigma

    @property
    def ps(self):
        return self.__ps

    @property
    def states(self):
        return self.__states

    @property
    def rates(self):
        """
        Get rate array with states as columns,
            and transition, transmission, sampling rates and equilibrium frequencies as rows.

        :return rate array
        :rtype np.array
        """
        return self.__rates

    @property
    def sigmas(self):
        """
        Get rate sum array (minus average sigma if needed).

        :return rate sum array
        :rtype np.array
        """
        if self.__sigmas is None:
            self.__sigmas = np.sum(self.rates[:-1, :], axis=0)
            self.__avg_sigma = np.min(self.__sigmas) if self.__minus_avg else 0
            self.__sigmas -= self.__avg_sigma
        return self.__sigmas

    @property
    def avg_sigma(self):
        """
        Get average rate sum.

        :return average rate sum
        :rtype float
        """
        return self.__avg_sigma

    @abstractmethod
    def num_params(self, type=None):
        return 0

    @abstractmethod
    def _get_states(self):
        pass

    @abstractmethod
    def params2rates(self, ps, sampled_pis=None, **kwargs):
        """
        Converts parameters into a rate array.
        """
        pass

    def get_bounds(self, lb, ub, sampled_pis=None, **kwargs):
        return np.array([[lb, ub]] * self.num_params())

    def _get_sd_formulas(self):
        """
        Given the rates, finds (horizontal) stationary distributions via formulas.

        :return: np.array of stat dist values
        """
        n = len(self.states)
        return np.repeat(1 / n, n)

    def get_sd(self):
        """
        Given the transition rates and transmission rates,
        finds stationary distributions (horizontal).
        :param rates: np.array of rates.
        :return: np.array of stat distribution.
        """
        return self._get_sd_formulas()

    @abstractmethod
    def map_state(self, state):
        pass

    def get_interesting_states(self, mutation='any'):
        return list(state.name for state in self.states)

    @abstractmethod
    def get_name(self):
        pass

    def get_sampled_pis(self):
        return self.rates[-1, :] * self.rates[2, :] / self.rates[-1, :].dot(self.rates[2, :])

    def get_avg_rates(self):
        return self.rates[-1, :].dot(self.rates[1:3, :].transpose())

    def get_expected_sampled_pis(self):
        return self.rates[2, :] * self.rates[-1, :] / self.rates[-1, :].dot(self.rates[2, :])

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
        pass

    def get_avg_unsampled_p(self, time_start, time_end, root_state=None, **kwargs):
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
        return (self.get_unsampled_p(time_start, root_state, **kwargs)
                + self.get_unsampled_p(time_end, root_state, **kwargs)
                + 2 * self.get_unsampled_p((time_end + time_start) / 2, root_state, **kwargs)) / 4



class State:
    def __init__(self, name, index, next_state=None, recipient=None):
        self._name = name
        self._i = index
        self._next = next_state
        self._recipient = recipient

    @property
    def index(self):
        return self._i

    @property
    def name(self):
        return self._name

    @property
    def next_state(self):
        return self._next

    @next_state.setter
    def next_state(self, next_state):
        self._next = next_state

    @property
    def recipient(self):
        return self._recipient

    @recipient.setter
    def recipient(self, recipient):
        self._recipient = recipient

    def is_symmetric(self):
        return self == self.recipient

    def __str__(self):
        return self._name

    def __repr__(self):
        return self.name

    def __eq__(self, other):
        """Override the default equals behavior"""
        if isinstance(other, self.__class__):
            return self._name == other._name
        return NotImplemented

    def __ne__(self, other):
        """Define a non-equality test"""
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        return NotImplemented

    def __hash__(self):
        """Override the default hash behavior (that returns the id or the object)"""
        return hash(self._name)

    def get_rate_sum(self, rates):
        """
        Returns the sum of sampling, transmission and state transition rates for this state.
        :param rates: array of rates
        :return: float, rate sum
        """
        return rates[:-1, self.index].sum()
