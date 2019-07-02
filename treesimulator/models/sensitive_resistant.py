import numpy as np

from treesimulator.models import Model, SAMPLING, TRANSMISSION, TRANSITION, State
from treesimulator.models.one_state import OneStateModel

SENSITIVE = 's'
RESISTANT = 'r'


def avg_rates2sr_rates(avg_lambda, avg_psi, lambda_s, psi_s, sampled_pi_s):
    pi_s = avg_psi / psi_s * sampled_pi_s
    pi_r = 1 - pi_s

    lambda_r = (avg_lambda - pi_s * lambda_s) / pi_r
    psi_r = (avg_psi - pi_s * psi_s) / pi_r

    mu_n = pi_r / pi_s * (avg_lambda - avg_psi + psi_r)

    return mu_n, lambda_s, lambda_r, psi_s, psi_r, pi_s, pi_r


class SensitiveResistantModel(Model):

    def num_params(self, type=None):
        if SAMPLING == type:
            return 2
        if TRANSMISSION == type:
            return 2
        if TRANSITION == type:
            return 2
        return 6

    def _get_states(self):
        resistant_state = State(name=RESISTANT, index=1)
        sensitive_state = State(name=SENSITIVE, index=0, next_state=resistant_state)

        sensitive_state.recipient = sensitive_state
        resistant_state.recipient = resistant_state
        resistant_state.next_state = sensitive_state

        return np.array([sensitive_state, resistant_state])

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
            mu_s, lambda_s, lambda_r, psi_s, psi_r = ps
        else:
            if nested_rates is not None:
                avg_lambda, avg_psi = nested_rates[1:3, 0]
                lambda_s, psi_s = ps
            else:
                avg_lambda, avg_psi, lambda_s, psi_s = ps

            mu_s, lambda_s, lambda_r, psi_s, psi_r, pi_n, pi_t = \
                avg_rates2sr_rates(avg_lambda, avg_psi, lambda_s, psi_s, sampled_pis[0])

        rates[0, 0] = mu_s
        rates[1, :] = lambda_s, lambda_r
        rates[2, :] = psi_s, psi_r

        if sampled_pis is None:
            pi_n, pi_t = self.get_sd(rates)

        rates[3, :] = pi_n, pi_t

        return rates

    def get_bounds(self, lb, ub, nested_rates=None, sampled_pis=None, **kwargs):
        """
        Convert a list of parameters into sampling, transition and transmission rate dictionaries
        :param ps: a list of parameter values, in the following order: lambda_n, psi_n
        :return: tuple of rate dictionaries: (change_rates, bifurcation_rates, sampling_rates)
        """

        if nested_rates is None:
            return np.array([[lb, ub]] * 4)

        avg_lambda, avg_psi = nested_rates[1:3, 0]
        sampled_pi_n = sampled_pis[0]

        bounds = np.array([[lb, ub]] * 2)

        # pi_n = avg_psi * sampled_pi_n / psi_n <= 1 => psi_n >= avg_psi * sampled_pi_n
        bounds[1, 0] = max(bounds[1, 0], avg_psi * sampled_pi_n)
        # psi_n >= avg_psi * sampled_pi_n max(psi) / (max(psi) - avg_psi * sampled_pi_t)
        bounds[1, 0] = max(bounds[1, 0], avg_psi * sampled_pi_n * bounds[1, 1] / (bounds[1, 1] - avg_psi * sampled_pis[1]))

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
            -s_n + s_t - l_t + l_n,
            2 * l_t - l_n + mu + s_n - s_t,
            -l_t
        ]
        roots = np.roots(p)
        roots = roots[np.where(np.isreal(roots) & (~np.isnan(roots)) & (roots >= 0) & (roots <= 1))]
        pi_n = roots[0]

        return np.array([pi_n, 1 - pi_n])

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
    def kappa2psi(sampled_pis, avg_lambda, avg_psi, kappa_n, kappa_ns, kappa_s):
        pi_n = avg_lambda / (kappa_n + avg_lambda - avg_psi)
        pi_ns = (pi_n * kappa_n - sampled_pis[0] * avg_psi) / kappa_ns
        pi_nr = pi_n - pi_ns

        pi_s = ((1 - pi_nr) * avg_lambda + (pi_nr - sampled_pis[0]) * avg_psi) / (kappa_s + avg_lambda - avg_psi)
        return avg_psi * sampled_pis[1:3].sum() / pi_s
