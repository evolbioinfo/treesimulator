import numpy as np

from treesimulator.models import Model, SAMPLING, TRANSMISSION, TRANSITION, State

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

    def params2rates(self, ps, sampled_pis=None, **kwargs):
        """
        Converts parameters into a rate array.

        :param ps: parameters, in the following order:
            reversion, treatment_elsewhere, drm,
            lambda_nr, lambda_ns, lambda_tr, lambda_ts,
            treatment_nr, treatment_ns, treatment_ts, treatment_tr
        """
        if sampled_pis is None:
            self.rates[0, :] = np.hstack((ps[:3], [0]))
            self.rates[1, :] = np.hstack(([ps[3]], ps[3: 6])) if len(ps) <= 9 else ps[3: 7]
            self.rates[2, :] = np.hstack(([ps[6]], ps[6:])) if len(ps) <= 9 else ps[7:]
            self.rates[3, :] = self.get_sd()
        else:
            # todo: find formulas
            raise ValueError('Not implemented')

    def get_bounds(self, lb, ub, sampled_pis=None, **kwargs):
        return np.array([[lb, ub]] * (9 if sampled_pis is None else 8))

    def _get_sd_formulas(self):
        mu_nr, mu_ns, mu_ts, _ = self.rates[0, :]
        l_nr, l_ns, l_ts, l_tr = self.rates[1, :]
        s_nr, s_ns, s_ts, s_tr = self.rates[2, :]

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
