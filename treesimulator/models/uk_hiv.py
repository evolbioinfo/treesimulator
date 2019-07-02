import numpy as np

from treesimulator.models import Model, SAMPLING, TRANSMISSION, TRANSITION, State
from treesimulator.models.naive_treated import NaiveTreatedModel, avg_rates2nt_rates

TR = 'tr'
NR = 'nr'
TS = 'ts'
NS = 'ns'

US = 'us'
UR = 'ur'


class UKHivModel(Model):

    def num_params(self, type=None):
        if SAMPLING == type:
            return 3
        if TRANSMISSION == type:
            return 3
        if TRANSITION == type:
            return 3
        return 9

    def _get_states(self):
        tr_state = State(name=TR, index=3)
        ts_state = State(name=TS, index=2, next_state=tr_state)
        ns_state = State(name=NS, index=1, next_state=ts_state)
        nr_state = State(name=NR, index=0, next_state=ns_state)

        nr_state.recipient = nr_state
        ns_state.recipient = ns_state
        tr_state.recipient = nr_state
        ts_state.recipient = ns_state

        return np.array([nr_state, ns_state, ts_state, tr_state])

    def params2rates(self, ps, rates=None, nested_rates=None, sampled_pis=None, **kwargs):
        """
        Get rate matrix from the parameter vector.
        :param ps: np.array([mu_nr, mu_ns, mu_ts, lambda_n, lambda_ts, lambda_tr, psi_n, psi_ts, psi_tr])
        :param rates: rate matrix 3 x 4 to be updated: rows are mutations, transmissions, sampling, columns are states: nr, ns, ts, tr
        :param kwargs:
        :return: rates
        """
        update_n_params = nested_rates is None or rates is None or sampled_pis is None

        if rates is None:
            rates = np.zeros(shape=(4, 4), dtype=np.float)

        if sampled_pis is None:
            rates[0, :] = np.hstack((ps[:3], [0]))
            rates[1, :] = np.hstack(([ps[3]], ps[3: 6])) if len(ps) <= 9 else ps[3: 7]
            rates[2, :] = np.hstack(([ps[6]], ps[6:])) if len(ps) <= 9 else ps[7:]
            rates[3, :] = self.get_sd(rates)
            return rates

        if nested_rates is None:
            avg_lambda, avg_psi, lambda_n, psi_n, lambda_tr, psi_tr = ps
            mu_n, lambda_n, lambda_t, psi_n, psi_t, pi_n, pi_t = avg_rates2nt_rates(avg_lambda, avg_psi, lambda_n,
                                                                                    psi_n, sampled_pis[:2].sum())
        else:
            lambda_tr, psi_tr = ps
            lambda_n = nested_rates[1, 0]
            psi_n = nested_rates[2, 0]
            pi_n, pi_t = nested_rates[3, :]
            avg_lambda, avg_psi = nested_rates[3, :].dot(nested_rates[1, :]), nested_rates[3, :].dot(nested_rates[2, :])
            mu_n = nested_rates[0, 0]

        if update_n_params:
            rates[1, :2] = [lambda_n, lambda_n]
            rates[2, :2] = [psi_n, psi_n]
            rates[3, :2] = sampled_pis[:2] * avg_psi / psi_n
            rates[0, 1] = pi_n / rates[3, 1] * mu_n

        pi_tr = avg_psi / psi_tr * sampled_pis[3]
        pi_ts = pi_t - pi_tr

        rates[3, 2:] = [pi_ts, pi_tr]
        rates[1, 2:] = [(avg_lambda - pi_n * lambda_n - pi_tr * lambda_tr) / pi_ts, lambda_tr]
        rates[2, 2:] = [(avg_psi - pi_n * psi_n - pi_tr * psi_tr) / pi_ts, psi_tr]
        rates[0, 0] = -avg_lambda + avg_psi + lambda_n - psi_n + pi_tr / rates[3, 0] * lambda_tr
        rates[0, 2] = pi_tr / pi_ts * (avg_lambda - avg_psi + psi_tr)

        return rates

    def get_bounds(self, lb, ub, nested_rates=None, sampled_pis=None, **kwargs):
        """
        Convert a list of parameters into sampling, transition and transmission rate dictionaries
        :param ps: a list of parameter values, in the following order: lambda_tr, psi_tr
        :return: tuple of rate dictionaries: (change_rates, bifurcation_rates, sampling_rates)
        """
        if nested_rates is None:
            return np.array([[lb, ub]] * 6)

        lambda_n = nested_rates[1, 0]
        pi_n, pi_t = nested_rates[3, :]

        avg_lambda, avg_psi = nested_rates[3, :].dot(nested_rates[1, :]), nested_rates[3, :].dot(nested_rates[2, :])
        sampled_pi_tr = sampled_pis[3]

        bounds = np.array([[lb, ub]] * 2)

        # psi_tr >= avg_psi * sampled_pi_tr / pi_t
        bounds[1, 0] = max(bounds[1, 0], avg_psi * sampled_pi_tr / pi_t)

        # lambda_tr <= (avg_lambda - pi_n lambda_n) max(psi_tr) / (avg_psi sampled_pi_tr)
        bounds[0, 1] = min(bounds[0, 1], (avg_lambda - pi_n * lambda_n) * bounds[1, 1] / (avg_psi * sampled_pi_tr))

        return bounds

    def _get_sd_formulas(self, rates):
        """
        Given the transition rates and transmission rates,
        finds stationary distributions (horizontal) via formulas.
        :return: np.array of stat dist values
        """
        mu_nr, mu_ns, mu_ts, _ = rates[0, :]
        l_nr, l_ns, l_ts, l_tr = rates[1, :]
        s_nr, s_ns, s_ts, s_tr = rates[2, :]

        a = -l_nr + mu_nr + s_nr
        b = -l_ns + mu_ns + s_ns
        c = mu_ts + s_ts
        d = s_tr

        p = [
            1,
            a + b + c + d,
            a * (b + c + d) + b * (c + d) + c * d - l_ts * mu_ns,
            a * (b * c + c * d + b * d) + b * c * d - l_ts * mu_ns * (a + d),
            a * b * c * d - mu_ns * l_ts * a * d - l_tr * mu_nr * mu_ns * mu_ts
        ]

        for L in sorted((np.real(L) for L in np.roots(p) if np.isreal(L) and not np.isnan(L)),
                        reverse=True):
            pi_ts = 1 / (1 + mu_ts / (L + d) + (L + c) / mu_ns + mu_ts * l_tr / ((L + d) * (L + a)))
            pi_tr = pi_ts * mu_ts / (L + d)
            pi_ns = pi_ts * (L + c) / mu_ns
            pi_nr = pi_tr * l_tr / (L + a)
            res = np.array([pi_nr, pi_ns, pi_ts, pi_tr])
            if np.all(res >= 0):
                return res
        return None

    def map_state(self, state):
        if state == NR:
            return self.states[0]
        if state == NS:
            return self.states[1]
        if state == TS:
            return self.states[2]
        if state == TR:
            return self.states[3]
        if state == US:
            return [self.states[1], self.states[2]]
        if state == UR:
            return [self.states[0], self.states[3]]
        return self.states

    def get_interesting_states(self, mutation='any'):
        return [TR, NR, UR] + ([NS, TS, US] if 'any' == mutation else [])

    def get_name(self):
        return 'uk'

    def get_nested_model(self):
        return NaiveTreatedModel()

    # def find_rates(self, avg_lambda, avg_psi, sampled_pis):
    #     # let's assume our x is [pi_nr, pi_ns, pi_ts, pi_ts, mu_..., ..., lambda_..., ..., psi_..., ...]
    #     pi_nr = Symbol("pi_nr", positive=True)
    #     pi_ns = Symbol("pi_ns", positive=True)
    #     pi_ts = Symbol("pi_ts", positive=True)
    #     pi_tr = Symbol("pi_tr", positive=True)
    #
    #     mu_nr = Symbol("mu_nr", positive=True)
    #     mu_ns = Symbol("mu_ns", positive=True)
    #     mu_ts = Symbol("mu_ts", positive=True)
    #
    #     la_n = Symbol("la_n", positive=True)
    #     la_ts = Symbol("la_ts", positive=True)
    #     la_tr = Symbol("la_tr", positive=True)
    #
    #     ps_n = Symbol("ps_nr", positive=True)
    #     ps_ts = Symbol("ps_ts", positive=True)
    #     ps_tr = Symbol("ps_tr", positive=True)
    #
    #     s_pi_nr, s_pi_ns, s_pi_ts, s_pi_tr = sampled_pis
    #
    #     avg_diff = avg_lambda - avg_psi
    #
    #     equations = [
    #         pi_nr + pi_ns + pi_ts + pi_tr - 1,
    #         #
    #         # pi_nr * ps_nr + pi_ns * ps_ns + pi_ts * ps_ts + pi_tr * ps_tr - avg_psi,
    #         # pi_nr * la_nr + pi_ns * la_ns + pi_ts * la_ts + pi_tr * la_tr - avg_lambda,
    #         #
    #         pi_nr * ps_n - s_pi_nr * avg_psi,
    #         pi_ns * ps_n - s_pi_ns * avg_psi,
    #         pi_ts * ps_ts - s_pi_ts * avg_psi,
    #         pi_tr * ps_tr - s_pi_tr * avg_psi,
    #         #
    #         pi_nr * (avg_diff + ps_n + mu_nr - la_n) - pi_tr * la_tr,
    #         pi_ns * (avg_diff + ps_n + mu_ns - la_n) - pi_nr * mu_nr - pi_ts * la_ts,
    #         pi_ts * (avg_diff + ps_ts + mu_ts) - pi_ns * mu_ns,
    #         pi_tr * (avg_diff + ps_tr) - pi_ts * mu_ts,
    #     ]
    #
    #     if self.num_params(SAMPLING) == 1:
    #         equations += [
    #             ps_n - ps_ts,
    #             ps_n - ps_tr
    #         ]
    #
    #     return solve(equations)

    @staticmethod
    def get_rates_from_avg_rates(sampled_pis, avg_lambda, avg_psi, lambda_n, kappa_n,
                                 lambda_s=None, kappa_s=None):
        rates = np.zeros(shape=(4, 4), dtype=np.float)

        lambda_nr, lambda_ns = lambda_n, lambda_n

        pi_n = avg_lambda / (kappa_n + avg_lambda - avg_psi)
        pi_t = 1 - pi_n

        psi_n = avg_psi * sampled_pis[:2].sum() / pi_n

        psi_nr, psi_ns = psi_n, psi_n

        lambda_t = (avg_lambda - pi_n * lambda_n) / pi_t
        psi_t = (avg_psi - pi_n * psi_n) / pi_t

        mu_n = pi_t / pi_n * (avg_lambda - avg_psi + psi_t)

        pi_nr = avg_psi / psi_nr * sampled_pis[0]
        pi_ns = pi_n - pi_nr

        mu_ns = pi_n / pi_ns * mu_n

        pi_s = ((1 - pi_nr) * avg_lambda + (pi_nr - sampled_pis[0]) * avg_psi) / (kappa_s + avg_lambda - avg_psi)
        pi_r = 1 - pi_s

        pi_ts = pi_s - pi_ns
        pi_tr = pi_r - pi_nr

        psi_ts = avg_psi * sampled_pis[2] / pi_ts
        psi_tr = (psi_t * pi_t - psi_ts * pi_ts) / pi_tr

        lambda_ts = (pi_s * lambda_s - pi_ns * lambda_ns) / pi_ts
        lambda_tr = (pi_t * lambda_t - pi_ts * lambda_ts) / pi_tr

        mu_nr = avg_psi - avg_lambda + lambda_n - psi_n + pi_tr / pi_nr * lambda_tr
        mu_ts = pi_tr / pi_ts * (-avg_psi + avg_lambda + psi_tr)

        rates[0, :3] = mu_nr, mu_ns, mu_ts
        rates[1, :] = lambda_nr, lambda_ns, lambda_ts, lambda_tr
        rates[2, :] = psi_nr, psi_ns, psi_ts, psi_tr
        rates[3, :] = pi_nr, pi_ns, pi_ts, pi_tr

        if np.all(rates[0, :-1] > 0) and np.all(rates[1:, :] > 0) and np.all(rates[-1, :] <= 1):
            return rates

        return None

    @staticmethod
    def get_rates_from_avg_one_state_rates(sampled_pis, avg_lambda, avg_psi, lambda_n, kappa_n,
                                           lambda_s=None, kappa_s=None, lambda_ns=None, kappa_ns=None):
        rates = np.zeros(shape=(4, 4), dtype=np.float)

        pi_n = avg_lambda / (kappa_n + avg_lambda - avg_psi)
        pi_ns = (pi_n * kappa_n - sampled_pis[0] * avg_psi) / kappa_ns
        pi_nr = pi_n - pi_ns

        psi_nr = avg_psi * sampled_pis[0] / pi_nr
        psi_ns = avg_psi * sampled_pis[1] / pi_ns

        lambda_nr = (pi_n * lambda_n - pi_ns * lambda_ns) / pi_nr

        mu_ns = kappa_ns - psi_ns

        pi_s = ((1 - pi_nr) * avg_lambda + (pi_nr - sampled_pis[0]) * avg_psi) / (kappa_s + avg_lambda - avg_psi)
        pi_r = 1 - pi_s
        pi_ts = pi_s - pi_ns
        pi_tr = pi_r - pi_nr

        psi_ts = avg_psi * sampled_pis[2] / pi_ts
        psi_tr = avg_psi * sampled_pis[3] / pi_tr

        lambda_ts = (pi_s * lambda_s - pi_ns * lambda_ns) / pi_ts
        lambda_tr = (avg_lambda - pi_s * lambda_s - pi_nr * lambda_nr) / pi_tr

        mu_ts = pi_tr / pi_ts * (avg_lambda - avg_psi + psi_tr)
        mu_nr = lambda_nr - psi_nr - avg_lambda + avg_psi + pi_tr / pi_nr * lambda_tr

        rates[0, :3] = mu_nr, mu_ns, mu_ts
        rates[1, :] = lambda_nr, lambda_ns, lambda_ts, lambda_tr
        rates[2, :] = psi_nr, psi_ns, psi_ts, psi_tr
        rates[3, :] = pi_nr, pi_ns, pi_ts, pi_tr

        # if np.all(rates[0, :-1] > 0) and np.all(rates[1:, :] > 0) and np.all(rates[-1, :] <= 1):
        #     return rates
        #
        # return None
        return rates

    @staticmethod
    def get_rates_from_other_rates(sampled_pis, avg_lambda, avg_psi, lambda_n, psi_n, lambda_s, psi_s, lambda_ns,
                                   psi_ns):
        rates = np.zeros(shape=(4, 4), dtype=np.float)

        pi_n = sampled_pis[:2].sum() * avg_psi / psi_n
        pi_ns = sampled_pis[1] * avg_psi / psi_ns
        pi_nr = pi_n - pi_ns
        pi_s = sampled_pis[1:3].sum() * avg_psi / psi_s
        pi_r = 1 - pi_s
        pi_tr = pi_r - pi_nr
        pi_ts = pi_s - pi_ns

        psi_nr = avg_psi * sampled_pis[0] / pi_nr
        psi_ts = avg_psi * sampled_pis[2] / pi_ts
        psi_tr = avg_psi * sampled_pis[3] / pi_tr

        lambda_nr = (pi_n * lambda_n - pi_ns * lambda_ns) / pi_nr
        lambda_ts = (pi_s * lambda_s - pi_ns * lambda_ns) / pi_ts
        lambda_tr = (avg_lambda - pi_s * lambda_s - pi_nr * lambda_nr) / pi_tr

        mu_s = (pi_tr * (avg_lambda - avg_psi) + sampled_pis[-1] * avg_psi) / pi_s
        mu_ts = mu_s * pi_s / pi_ts
        mu_n = avg_lambda / pi_n - avg_lambda + avg_psi - psi_n
        mu_ns = mu_n * pi_n / pi_ns
        mu_nr = lambda_nr - psi_nr - avg_lambda + avg_psi + pi_tr / pi_nr * lambda_tr

        rates[0, :3] = mu_nr, mu_ns, mu_ts
        rates[1, :] = lambda_nr, lambda_ns, lambda_ts, lambda_tr
        rates[2, :] = psi_nr, psi_ns, psi_ts, psi_tr
        rates[3, :] = pi_nr, pi_ns, pi_ts, pi_tr

        # if np.all(rates[0, :-1] > 0) and np.all(rates[1:, :] > 0) and np.all(rates[-1, :] <= 1):
        return rates

        # return None

    @staticmethod
    def kappa2psi(sampled_pis, avg_lambda, avg_psi, kappa_n, kappa_ns):
        pi_n = avg_lambda / (kappa_n + avg_lambda - avg_psi)
        pi_ns = (pi_n * kappa_n - sampled_pis[0] * avg_psi) / kappa_ns
        return avg_psi * sampled_pis[1] / pi_ns

    @staticmethod
    def get_naive_rates(rates):
        return rates[-1, :2].dot(rates[1:3, :2].transpose()) / rates[-1, :2].sum()

    @staticmethod
    def get_ns_rates(rates):
        return rates[1:3, 1]

    @staticmethod
    def get_nr_rates(rates):
        return rates[1:3, 0]

    @staticmethod
    def get_sensitive_rates(rates):
        return rates[-1, 1:3].dot(rates[1:3, 1:3].transpose()) / rates[-1, 1:3].sum()

    @staticmethod
    def get_resistant_rates(rates):
        return rates[-1, [0, 3]].dot(rates[1:3, [0, 3]].transpose()) / rates[-1, [0, 3]].sum()

    @staticmethod
    def get_naive_treated_rates(rates):
        expected_nt_rates = np.zeros(shape=(4, 2))
        expected_nt_rates[-1, 0] = rates[-1, :2].sum()
        expected_nt_rates[-1, 1] = rates[-1, 2:].sum()
        expected_nt_rates[1:3, 0] = rates[-1, :2].dot(rates[1:3, :2].transpose()) / expected_nt_rates[-1, 0]
        expected_nt_rates[1:3, 1] = rates[-1, 2:].dot(rates[1:3, 2:].transpose()) / expected_nt_rates[-1, 1]
        expected_nt_rates[0, 0] = rates[-1, 1] * rates[0, 1] / expected_nt_rates[-1, 0]
        return expected_nt_rates
