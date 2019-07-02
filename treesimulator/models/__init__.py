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
        logging.info('Expected {} are:{}{}'.format(label, separator, np.round(expected_rates, DIGITS_TO_KEEP)))
        error_fractions = get_error_fractions(estimated_rates, expected_rates)
        logging.info('The error fractions are:{}{}'.format(separator, np.round(error_fractions, DIGITS_TO_KEEP)))
        logging.info('The error fractions in percentage:{}{}'
                     .format(separator, np.round(100 * error_fractions)))
        logging.info('\n==========================\n')


class Model(object):
    def __init__(self):
        super(Model, self).__init__()
        self._states = self._get_states()

    @property
    def states(self):
        return self._states

    @abstractmethod
    def num_params(self, type=None):
        return 0

    @abstractmethod
    def _get_states(self):
        pass

    @abstractmethod
    def params2rates(self, ps, rates=None, nested_rates=None, sampled_pis=None, **kwargs):
        """
        Convert a list of parameters into sampling, transition and transmission rate dictionaries
        :param ps: a list of parameter values, in the following order: drm_loss, drm_gain, treatment_elsewhere,
                  lambda_nr, lambda_ns, lambda_ts, lambda_tr, treatment_first, treatment_second
        :return: tuple of rate dictionaries: (change_rates, bifurcation_rates, sampling_rates)
        """
        pass

    def get_constraints(self, rates=None, nested_rates=None, sampled_pis=None):

        def check_constraints(ps):
            rs = self.params2rates(ps, rates, nested_rates, sampled_pis)
            if np.any(rs < 0):
                return -1
            if np.any(rs[3:] > 1):
                return -1
            return 1

        return {'type': 'ineq', 'fun': check_constraints},

    def get_bounds(self, lb, ub, nested_rates=None, sampled_pis=None, **kwargs):
        return np.array([[lb, ub]] * self.num_params())

    def _get_sd_formulas(self, rates):
        n = len(self.states)
        return np.repeat(1 / n, n)

    def get_sd(self, rates):
        """
        Given the transition rates and transmission rates,
        finds stationary distributions (horizontal).
        :param rates: np.array of rates.
        :return: np.array of stat distribution.
        """
        return self._get_sd_formulas(rates)

    @abstractmethod
    def map_state(self, state):
        pass

    def get_interesting_states(self, mutation='any'):
        return list(state.name for state in self.states)

    @abstractmethod
    def get_name(self):
        pass

    def get_nested_model(self):
        return None

    @staticmethod
    def get_sampled_pis(rates):
        return rates[-1, :] * rates[2, :] / rates[-1, :].dot(rates[2, :])

    @staticmethod
    def get_avg_rates(rates):
        return rates[-1, :].dot(rates[1:3, :].transpose())


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
