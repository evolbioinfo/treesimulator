import numpy as np

from treesimulator.models import Model, SAMPLING, TRANSMISSION, TRANSITION, State

SENSITIVE = 's'
RESISTANT = 'r'


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

    def params2rates(self, ps, sampled_pis=None, **kwargs):
        """
        Converts parameters into a rate array.

        :param ps: parameters, in the following order:
            drm, reversion, lambda_s, lambda_r, treatment_s, treatment_r
        """
        if sampled_pis is None:
            self.rates[0, :] = ps[:2]
            self.rates[1, :] = ps[2:4]
            self.rates[2, :] = ps[4:]
            self.rates[3, :] = self.get_sd()
        else:
            self.rates[0, 0] = ps[0]
            self.rates[1, :] = ps[1:3]
            self.rates[2, :] = ps[3:]
            # pi_s = (psi_r / psi_s) / (spi_r / spi_s + psi_r / psi_s)
            pi_s = self.rates[2, 1] / self.rates[2, 0] / (sampled_pis[1] / sampled_pis[0] + self.rates[2, 1]
                                                          / self.rates[2, 0])
            pi_r = 1 - pi_s
            self.rates[-1, :] = [pi_s, pi_r]
            # mu_r = lambda_r - psi_r - avg_lambda + avg_psi + mu_s * pi_s / pi_r
            self.rates[0, 1] = max(self.rates[1, 1] - self.rates[2, 1] - self.rates[1, :].dot(self.rates[-1, :])
                                   + self.rates[2, :].dot(self.rates[-1, :])
                                   + self.rates[0, 0] * self.rates[-1, 0] / self.rates[-1, 1], 0)
            if self.rates[0, 1] == 0:
                self.rates[3, :] = self.get_sd()

    def get_bounds(self, lb, ub, sampled_pis=None, **kwargs):
        return np.array([[lb, ub]] * (6 if sampled_pis is None else 5))

    def _get_sd_formulas(self):
        mu_s, mu_r = self.rates[0, :]
        l_s, l_r = self.rates[1, :]
        s_s, s_r = self.rates[2, :]

        p = [
            -s_s + s_r - l_r + l_s,
            mu_r + mu_s + l_r - l_s + s_s - s_r,
            -mu_r
        ]
        roots = np.roots(p)
        roots = roots[np.where(np.isreal(roots) & (~np.isnan(roots)) & (roots >= 0) & (roots <= 1))]
        pi_s = roots[0]

        return np.array([pi_s, 1 - pi_s])

    def map_state(self, state):
        state = str(state)
        if state.endswith('s'):
            return self.states[0]
        if state.endswith('r'):
            return self.states[1]
        return self.states

    def get_name(self):
        return 'sensitive-resistant'
