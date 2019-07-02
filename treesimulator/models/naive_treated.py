import logging

import numpy as np

from treesimulator.models import Model, SAMPLING, TRANSMISSION, TRANSITION, State
from treesimulator.models.one_state import OneStateModel

NAIVE = 'n'
TREATED = 't'


def avg_rates2nt_rates(avg_lambda, avg_psi, lambda_n, psi_n, sampled_pi_n):
    pi_n = avg_psi / psi_n * sampled_pi_n
    pi_t = 1 - pi_n

    lambda_t = (avg_lambda - pi_n * lambda_n) / pi_t
    psi_t = (avg_psi - pi_n * psi_n) / pi_t
    mu_n = pi_t / pi_n * (avg_lambda - avg_psi + psi_t)

    return mu_n, lambda_n, lambda_t, psi_n, psi_t, pi_n, pi_t


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

    def params2rates(self, ps, rates=None, nested_rates=None, sampled_pis=None, **kwargs):
        """
        Convert a list of parameters into sampling, transition and transmission rate dictionaries
        :param ps: a list of parameter values, in the following order: treatment_elsewhere,
                  lambda_n, lambda_t, treatment_first, treatment_second
        :return: tuple of rate dictionaries: (change_rates, bifurcation_rates, sampling_rates)
        """
        if rates is None:
            rates = np.zeros(shape=(4, 2), dtype=np.float)
        if sampled_pis is None:
            mu_n, lambda_n, lambda_t, psi_n, psi_t = ps
            rates[0, 0] = mu_n
            rates[1, :] = lambda_n, lambda_t
            rates[2, :] = psi_n, psi_t
            rates[3, :] = self.get_sd(rates)
            return rates
        if nested_rates is not None:
            avg_lambda, avg_psi = nested_rates[1:3, 0]
            lambda_n, psi_n = ps
            return NaiveTreatedModel.get_rates_from_other_rates(sampled_pis, avg_lambda, avg_psi, lambda_n, psi_n)
        else:
            rates[1, :] = ps[:2]
            rates[2, :] = ps[2:]
            # pi_n = spi_n * psi_t / (spi_n * psi_t + spi_t * psi_n)
            pi_n = sampled_pis[0] * rates[2, 1] / sampled_pis.dot(rates[2, [1, 0]])
            pi_t = 1 - pi_n
            rates[-1, :] = [pi_n, pi_t]
            # mu_n = pi_t (avg_lambda - avg_psi + psi_t) / pi_n
            rates[0, 0] = pi_t / pi_n * (rates[1, :].dot(rates[-1, :]) - rates[2, :].dot(rates[-1, :]) + rates[2, 1])
            return rates

    def get_bounds(self, lb, ub, nested_rates=None, sampled_pis=None, **kwargs):
        """
        Convert a list of parameters into sampling, transition and transmission rate dictionaries
        :param ps: a list of parameter values, in the following order: lambda_n, psi_n
        :return: tuple of rate dictionaries: (change_rates, bifurcation_rates, sampling_rates)
        """

        if sampled_pis is None:
            return np.array([[lb, ub]] * 5)

        if nested_rates is None:
            return np.array([[lb, ub]] * 4)

        avg_lambda, avg_psi = nested_rates[1:3, 0]
        sampled_pi_n = sampled_pis[0]

        bounds = np.array([[lb, ub]] * 2)

        # pi_n = avg_psi * sampled_pi_n / psi_n <= 1 => psi_n >= avg_psi * sampled_pi_n
        bounds[1, 0] = max(bounds[1, 0], avg_psi * sampled_pi_n)
        # psi_n >= avg_psi * sampled_pi_n max(psi) / (max(psi) - avg_psi * sampled_pi_t)
        bounds[1, 0] = max(bounds[1, 0],
                           avg_psi * sampled_pi_n * bounds[1, 1] / (bounds[1, 1] - avg_psi * sampled_pis[1]))

        # 0 <= lambda_t = (avg_lambda - pi_n * lambda_n) / pi_t
        # => lambda_n <= avg_lambda / pi_n = avg_lambda max(psi_n) / (avg_psi sampled_pi_n)
        bounds[0, 1] = min(bounds[0, 1], avg_lambda * bounds[1, 1] / (avg_psi * sampled_pi_n))

        return bounds

    def _get_sd_formulas(self, rates):
        """
        Given the transition rates and transmission rates,
        finds stationary distributions (horizontal) via formulas.
        :param rates: np.array of rates
        :return: list of dictionaries for possible stationary distributions in a form {state_i: value_i}
        """

        mu = rates[0, 0]
        l_n, l_t = rates[1, :]
        s_n, s_t = rates[2, :]

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
        return 'nt'

    def get_nested_model(self):
        return OneStateModel()

    @staticmethod
    def get_rates_from_avg_rates(pi_n_sampled, avg_lambda, avg_psi, lambda_n, kappa_n):
        rates = np.zeros(shape=(4, 2), dtype=np.float)

        pi_n = avg_lambda / (kappa_n + avg_lambda - avg_psi)
        pi_t = 1 - pi_n

        psi_n = avg_psi * pi_n_sampled / pi_n

        lambda_t = (avg_lambda - pi_n * lambda_n) / pi_t
        psi_t = (avg_psi - pi_n * psi_n) / pi_t

        mu_n = pi_t / pi_n * (avg_lambda - avg_psi + psi_t)

        rates[0, 0] = mu_n
        rates[1, :] = lambda_n, lambda_t
        rates[2, :] = psi_n, psi_t
        rates[3, :] = pi_n, pi_t

        return rates

    @staticmethod
    def kappa2psi(pi_n_sampled, avg_lambda, avg_psi, kappa_n):
        pi_n = avg_lambda / (kappa_n + avg_lambda - avg_psi)
        return avg_psi * pi_n_sampled / pi_n

    @staticmethod
    def get_rates_from_other_rates(sampled_pis, avg_lambda, avg_psi, lambda_n, psi_n, pi_n=None):
        rates = np.zeros(shape=(4, 2), dtype=np.float)

        pi_n = (sampled_pis[0] * avg_psi / psi_n) if pi_n is None else pi_n
        pi_t = 1 - pi_n

        mu_n = avg_lambda / pi_n - avg_lambda + avg_psi - psi_n

        lambda_t = (avg_lambda - pi_n * lambda_n) / pi_t
        psi_t = (avg_psi - pi_n * psi_n) / pi_t

        rates[0, 0] = mu_n
        rates[1, :] = lambda_n, lambda_t
        rates[2, :] = psi_n, psi_t
        rates[3, :] = pi_n, pi_t

        return rates

    @staticmethod
    def get_naive_rates(rates):
        expected_n_rates = rates[1:3, 0] * 1
        # if kappa:
        #     expected_n_rates[1] += rates[0, 0]
        return expected_n_rates
