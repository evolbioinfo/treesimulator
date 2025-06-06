import logging

import numpy as np
from scipy.optimize import least_squares
import sympy

EXPOSED = 'E'
INFECTIOUS = 'I'
SUPERSPREADER = 'S'
CONTACT = 'C'

SUPERSPREADING_FRACTION = f'f_{SUPERSPREADER}'

SS_TRANSMISSION_RATIO = f'X_{SUPERSPREADER}'

INCUBATION_FRACTION = f'f_{EXPOSED}'





class Model(object):

    def __init__(self, states=None,
                 transition_rates=None, transmission_rates=None, removal_rates=None, ps=None,
                 state_frequencies=None,
                 n_recipients=None,
                 *args, **kwargs):
        self.__pis = None
        self.__states = np.array(states)
        num_states = len(self.states)
        self.__ps = np.array(ps) if ps is not None else np.ones(num_states, dtype=float)
        self.check_ps()
        self.__n_recipients = np.array(n_recipients) if n_recipients is not None else np.ones(num_states, dtype=float)
        self.check_n_recipients()
        self.__transmission_rates = np.array(transmission_rates) if transmission_rates is not None \
            else np.zeros(shape=(num_states, num_states), dtype=float)
        self.check_transmission_rates()
        self.__transition_rates = np.array(transition_rates) if transition_rates is not None \
            else np.zeros(shape=(num_states, num_states), dtype=float)
        self.check_transition_rates()
        self.__removal_rates = np.array(removal_rates) if removal_rates is not None \
            else np.zeros(shape=num_states, dtype=float)
        self.check_removal_rates()
        if state_frequencies is not None:
            self.__pis = np.array(state_frequencies)
            self.check_frequencies()

    def clone(self):
        return Model(self.states, self.transition_rates, self.transmission_rates, self.removal_rates, self.ps)

    @property
    def ps(self):
        return self.__ps

    @ps.setter
    def ps(self, ps):
        self.__ps = ps

    @property
    def n_recipients(self):
        return self.__n_recipients

    @n_recipients.setter
    def n_recipients(self, n_recipients):
        self.__n_recipients = n_recipients
        self.check_n_recipients()

    @property
    def states(self):
        return self.__states

    def _state_frequencies_with_sympy(self):
        MU_IJ, LA_IJ, PSI_I = self.transition_rates, self.transmission_rates, self.removal_rates
        LA_I_ = LA_IJ.sum(axis=1)
        m = len(self.states)
        PI_I = np.array(sympy.symbols(' '.join(f"x_{k}" for k in range(m))))
        eqs = [sympy.Eq(PI_I.sum(), 1)]

        # pi_k dN = dN_k
        # dN = sum_i N_i sum_j la_ij dt - sum_i N_i psi_i dt = [sum_i (la_i_ - psi_i) pi_i] N dt,
        #   where la_i_ = sum_j la_ij;
        # dN_k = -N_k sum_j mu_kj dt + sum_j N_j mu_jk dt + sum_j N_j la_jk dt - N_k psi_k dt =
        #   = [-pi_k (sum_j mu_kj + psi_k) + sum_j pi_j (mu_jk + la_jk)] N dt

        dN_div_N_dt = PI_I.dot(LA_I_ - PSI_I)

        for k in range(m - 1):
            pi_k = PI_I[k]
            dN_k_div_N_dt = -pi_k * (MU_IJ[k, :].sum() + PSI_I[k]) + PI_I.dot(MU_IJ[:, k] + LA_IJ[:, k])
            eqs.append(sympy.Eq(pi_k * dN_div_N_dt - dN_k_div_N_dt, 0))

        solution = sympy.solve(eqs)
        if isinstance(solution, list):
            solution = next((s for s in solution if all((_ >= 0) and (_ <= 1) for _ in s.values())), None)
            if solution is None:
                raise ValueError('Sympy could not find a solution')
        self.state_frequencies = np.array([float(solution[PI_I[k]]) for k in range(m)])

    def _state_frequencies_with_least_squares(self):
        MU_IJ, LA_IJ, PSI_I = self.transition_rates, self.transmission_rates, self.removal_rates
        LA_I_ = LA_IJ.sum(axis=1)
        m = len(self.states)

        def func(PI_I):
            res = [PI_I.sum() - 1]

            # pi_k dN = dN_k
            # dN = sum_i N_i sum_j la_ij dt - sum_i N_i psi_i dt = [sum_i (la_i_ - psi_i) pi_i] N dt,
            #   where la_i_ = sum_j la_ij;
            # dN_k = -N_k sum_j mu_kj dt + sum_j N_j mu_jk dt + sum_j N_j la_jk dt - N_k psi_k dt =
            #   = [-pi_k (sum_j mu_kj + psi_k) + sum_j pi_j (mu_jk + la_jk)] N dt
            dN_div_N_dt = PI_I.dot(LA_I_ - PSI_I)

            for k in range(m - 1):
                pi_k = PI_I[k]
                dN_k_div_N_dt = -pi_k * (MU_IJ[k, :].sum() + PSI_I[k]) + PI_I.dot(MU_IJ[:, k] + LA_IJ[:, k])
                res.append(pi_k * dN_div_N_dt - dN_k_div_N_dt)
            return res

        self.state_frequencies = least_squares(func, x0=np.ones(m) / m, bounds=[0, 1]).x

    @property
    def state_frequencies(self):
        if self.__pis is None:
            if len(self.states) == 1:
                self.__pis = np.array([1.])
            else:
                try:
                    self._state_frequencies_with_sympy()
                except Exception as e1:
                    try:
                        self._state_frequencies_with_least_squares()
                    except Exception as e2:
                        logging.warning(f'Could not calculate the equilibrium frequencies due to {e1} and {e2}, '
                                        'setting them to equal ones!')
                        m = len(self.states)
                        self.__pis = np.ones(m) / m
        return self.__pis

    def check_frequencies(self):
        if np.any(np.round(self.__pis, 6) < 0):
            raise ValueError('Equilibrium frequencies cannot be negative')
        if np.any(np.round(self.__pis, 6) > 1):
            raise ValueError('Equilibrium frequencies cannot be greater than one')
        if np.round(self.__pis.sum(), 2) != 1:
            raise ValueError(f'Equilibrium frequencies must sum up to one but they sum up to {self.__pis.sum():g} instead')
        # Ensure they sum up to one perfectly
        self.__pis = np.maximum(self.__pis, 0)
        self.__pis /= self.__pis.sum()


    @state_frequencies.setter
    def state_frequencies(self, pis):
        self.__pis = np.array(pis)
        self.check_frequencies()

    @property
    def transition_rates(self):
        """
        Get transition rate matrix with states as columns and rows.

        :return rate array
        :rtype np.array
        """
        return self.__transition_rates

    @transition_rates.setter
    def transition_rates(self, rates):
        self.__transition_rates = rates
        self.check_transition_rates()

    @property
    def transmission_rates(self):
        """
        Get transmission rate matrix with states as columns and rows.

        :return rate array
        :rtype np.array
        """
        return self.__transmission_rates

    @transmission_rates.setter
    def transmission_rates(self, rates):
        self.__transmission_rates = rates
        self.check_transmission_rates()

    @property
    def removal_rates(self):
        """
        Get removal rate array with states as columns.

        :return rate array
        :rtype np.array
        """
        return self.__removal_rates

    @removal_rates.setter
    def removal_rates(self, rates):
        self.__removal_rates = rates
        self.check_removal_rates()

    def get_name(self):
        return 'MTBD'

    def check_transition_rates(self):
        n_states = len(self.states)
        if self.transition_rates.shape != (n_states, n_states):
            raise ValueError("Transition matrix shape is wrong, should be {}x{}.".format(n_states, n_states))
        if not np.all(self.transition_rates >= 0):
            raise ValueError("Transition rates cannot be negative")
        if np.any(np.diag(self.transition_rates) != 0):
            raise ValueError("Transition rates from the state to itself (i.e., diagonal transition rate matrix values) must be zero")

    def check_transmission_rates(self):
        n_states = len(self.states)
        if self.transmission_rates.shape != (n_states, n_states):
            raise ValueError("Transmission matrix shape is wrong, should be {}x{}.".format(n_states, n_states))
        if not np.all(self.transmission_rates >= 0):
            raise ValueError("Transmission rates cannot be negative")

    def check_removal_rates(self):
        n_states = len(self.states)
        if self.removal_rates.shape != (n_states,):
            raise ValueError("Removal rate vector length is wrong, should be {}.".format(n_states))
        if not np.all(self.removal_rates >= 0):
            raise ValueError("Removal rates cannot be negative")

    def check_ps(self):
        n_states = len(self.states)
        if self.ps.shape != (n_states,):
            raise ValueError("Sampling probability vector length is wrong, should be {}.".format(n_states))
        if not np.all(self.ps >= 0):
            raise ValueError("Sampling probabilities cannot be negative")
        if not np.all(self.ps <= 1):
            raise ValueError('Sampling probabilities cannot be greater than 1')

    def check_n_recipients(self):
        n_states = len(self.states)
        if self.n_recipients.shape != (n_states,):
            raise ValueError("Recipient number vector length is wrong, should be {}.".format(n_states))
        if not np.all(self.n_recipients >= 1):
            raise ValueError('The number of recipients cannot be below 1 '
                             '(put the transmission rate to zero to prevent transmission from a certain state)')

    def get_epidemiological_parameters(self):
        """Converts rate parameters to the epidemiological ones"""
        pis = self.state_frequencies
        transition_rates_per_state = self.transition_rates.sum(axis=1)
        transmission_rates_per_state = self.transmission_rates.sum(axis=1)
        Rs = np.ones(len(self.states), dtype=float) * np.inf
        irremovable_mask = self.removal_rates == 0
        Rs[~irremovable_mask] = \
            transmission_rates_per_state[~irremovable_mask] * self.n_recipients[~irremovable_mask] \
            / self.removal_rates[~irremovable_mask]
        Rs[(transmission_rates_per_state == 0) & irremovable_mask] = 1
        res = {}
        n_states = len(self.states)
        is_mult = np.any(self.n_recipients != 1)
        for i in range(n_states):
            state_i = self.states[i]
            if self.removal_rates[i]:
                res[f'R_{state_i}'] = Rs[i]
            res[f'd_{state_i}'] = 1 / (self.removal_rates[i] + transition_rates_per_state[i])
            if n_states > 1:
                res[f'pi_{state_i}'] = pis[i]
            if is_mult:
                res[f'n_recipients_{state_i}'] = self.n_recipients[i]
            for j in range(n_states):
                state_j = self.states[j]
                if i != j and self.transition_rates[i][j]:
                    res[f'mu_{state_i}{state_j}'] = self.transition_rates[i][j]
                if self.transmission_rates[i][j]:
                    res[f'la_{state_i}{state_j}'] = self.transmission_rates[i][j]
            if self.removal_rates[i]:
                res[f'psi_{state_i}'] = self.removal_rates[i]
                res[f'p_{state_i}'] = self.ps[i]
        return res


class BirthDeathExposedInfectiousModel(Model):

    def __init__(self, mu, la, psi, p=0.5, *args, **kwargs):
        """
        :param mu: transition E->I
        :param la: transmission from I to E
        :param psi: removal of I
        :param p: sampling probability of I
        """
        mus = np.zeros(shape=(2, 2), dtype=float)
        mus[0, 1] = mu
        las = np.zeros(shape=(2, 2), dtype=float)
        las[1, 0] = la
        psis = np.zeros(shape=2, dtype=float)
        psis[1] = psi

        Model.__init__(self, states=[EXPOSED, INFECTIOUS],
                       transition_rates=mus, transmission_rates=las, removal_rates=psis, ps=[0, p],
                       *args, **kwargs)

    @property
    def state_frequencies(self):
        mu = self.transition_rates[0, 1]
        la = self.transmission_rates[1, 0]
        psi = self.removal_rates[1]
        mu_plus_psi = mu + psi
        if la == psi:
            pi_i = mu / mu_plus_psi
        else:
            two_la_minus_psi = 2 * (la - psi)
            det = np.power(np.power(mu_plus_psi, 2) + 2 * mu * two_la_minus_psi, 1 / 2)
            pi_i = (det - mu_plus_psi) / two_la_minus_psi
        return np.array([1 - pi_i, pi_i])

    def get_name(self):
        return 'BDEI'

    def get_epidemiological_parameters(self):
        """Converts rate parameters to the epidemiological ones"""
        result = Model.get_epidemiological_parameters(self)
        result[INCUBATION_FRACTION] = result[f'd_{EXPOSED}'] / (result[f'd_{EXPOSED}'] + result[f'd_{INFECTIOUS}'])
        return result


class BirthDeathModel(Model):

    def __init__(self, la, psi, p=0.5, *args, **kwargs):
        """
        :param la: transmission
        :param psi: removal
        :param p: sampling probability
        """
        las = la * np.ones(shape=(1, 1), dtype=float)
        Model.__init__(self, states=[INFECTIOUS], transmission_rates=las, removal_rates=[psi], ps=[p], *args, **kwargs)

    def get_name(self):
        return 'BD'


class BirthDeathWithSuperSpreadingModel(Model):

    def __init__(self, la_nn, la_ns, la_sn, la_ss, psi, p=0.5, *args, **kwargs):
        """
        :param la_nn: transmission from I to I
        :param la_ns: transmission from I to S
        :param la_sn: transmission from S to I
        :param la_ss: transmission from S to S
        :param psi: removal
        :param p: sampling
        """
        las = np.zeros(shape=(2, 2), dtype=float)
        # Avoid division of zero by zero with these ifs
        s_ratio = la_ss / la_ns if la_ss != la_ns else 1
        n_ratio = la_sn / la_nn if la_sn != la_nn else 1
        if np.abs(s_ratio - n_ratio) > 1e-3:
            raise ValueError(
                'transmission ratio constraint is violated: la_ss / la_ns ({}) must be equal to la_sn / la_nn ({})'
                .format(s_ratio, n_ratio))
        las[0, 0] = la_nn
        las[0, 1] = la_ns
        las[1, 0] = la_sn
        las[1, 1] = la_ss
        psis = psi * np.ones(shape=2, dtype=float)
        Model.__init__(self, states=[INFECTIOUS, SUPERSPREADER],
                       transmission_rates=las, removal_rates=psis, ps=[p, p], *args, **kwargs)

    @property
    def state_frequencies(self):
        la_ss = self.transmission_rates[1, 1]
        la_sn = self.transmission_rates[1, 0]
        la_ns = self.transmission_rates[0, 1]
        la_nn = self.transmission_rates[0, 0]
        f = la_ss / (la_ss + la_sn) if (la_ss + la_sn > 0) else la_ns / (la_nn + la_ns)
        return np.array([1 - f, f])

    def get_name(self):
        return 'BDSS'

    def get_epidemiological_parameters(self):
        """Converts rate parameters to the epidemiological ones"""
        pis = self.state_frequencies
        # R0 could also be expressed
        # as (self.transmission_rates[0, 0] + self.transmission_rates[1, 1]) / self.removal_rates[0]
        # in cases with 1 recipient
        # (it gives the same result due to constraints)
        result = Model.get_epidemiological_parameters(self)
        # avoid division of zero by zero
        result[SS_TRANSMISSION_RATIO] = self.transmission_rates[1, 1] / self.transmission_rates[0, 1] \
            if self.transmission_rates[1, 1] != self.transmission_rates[0, 1] else 1
        result[SUPERSPREADING_FRACTION] = pis[-1]
        return result


class BirthDeathExposedInfectiousWithSuperSpreadingModel(Model):

    def __init__(self, mu_n, mu_s, la_n, la_s, psi, p=0.5, *args, **kwargs):
        """
        :param mu_n: transition from E to I
        :param mu_s: transition from E to S
        :param la_n: transmission from I to E
        :param la_s: transmission from S to E
        :param psi: removal
        :param p: sampling
        """
        mus = np.zeros(shape=(3, 3), dtype=float)
        mus[0, 1] = mu_n
        mus[0, 2] = mu_s

        las = np.zeros(shape=(3, 3), dtype=float)
        las[1, 0] = la_n
        las[2, 0] = la_s

        psis = psi * np.ones(shape=3, dtype=float)
        psis[0] = 0

        Model.__init__(self, states=[EXPOSED, INFECTIOUS, SUPERSPREADER],
                       transition_rates=mus, transmission_rates=las, removal_rates=psis, ps=[0, p, p],
                       *args, **kwargs)

    def get_name(self):
        return 'BDEISS'

    def get_epidemiological_parameters(self):
        """Converts rate parameters to the epidemiological ones"""
        pis = self.state_frequencies

        result = Model.get_epidemiological_parameters(self)
        result[SS_TRANSMISSION_RATIO] = self.transmission_rates[2, 0] / self.transmission_rates[1, 0]
        result[SUPERSPREADING_FRACTION] = pis[2] / (pis[1] + pis[2])
        result[INCUBATION_FRACTION] = result[f'd_{EXPOSED}'] / (result[f'd_{EXPOSED}'] + result[f'd_{INFECTIOUS}'])
        return result


class CTModel(Model):
    """
    Contact-tracing model adds two parameters:
        * upsilon -- the probability to notify a contact upon sampling
        * phi -- the removal rate after being notified
    """

    def __init__(self, model, phi=np.inf, upsilon=0.5, allow_irremovable_states=False, *args, **kwargs):
        """
        
        :param model: initial model, whose states will now be contact-traced
        :param phi: sampling rate after being notified
        :param upsilon: probability to notify a contact upon sampling
        :param allow_irremovable_states: if set to True and the initial model included "irremovable" states
            (i.e., whose removal rate was zero, e.g., E in the BDEI model), then even after notification 
            their removal rate will stay zero, and the corresponding individuals will become "removable" (at a rate phi)
            only once they change the state to a "removable" one (e.g., from E-notified to I-notified in BDEI-CT).
            By default, allow_irremovable_states is False and all states become "removable" once notified.
        :param args: 
        :param kwargs: 
        """
        transition_rates = np.pad(model.transition_rates, ((0, model.transition_rates.shape[0]),
                                                           (0, model.transition_rates.shape[1])),
                                  mode='constant', constant_values=0)
        transition_rates[model.transition_rates.shape[0]:, model.transition_rates.shape[1]:] = model.transition_rates

        transmission_rates = np.pad(model.transmission_rates, ((0, model.transmission_rates.shape[0]),
                                                               (0, model.transmission_rates.shape[1])),
                                    mode='constant', constant_values=0)
        transmission_rates[model.transmission_rates.shape[0]:, :model.transmission_rates.shape[1]] \
            = model.transmission_rates

        n_removal_rates = model.removal_rates.shape[0]
        removal_rates = np.pad(model.removal_rates, (0, n_removal_rates), mode='constant', constant_values=phi)
        rhos = np.pad(model.ps, (0, model.ps.shape[0]), mode='constant', constant_values=1)
        if allow_irremovable_states:
            # If there was no way to remove a certain state (e.g. E in BDEI),
            # then notification should not change its "irremovable" status
            mask = np.zeros(2 * n_removal_rates, dtype=bool)
            mask[n_removal_rates:] = (model.removal_rates == 0)
            removal_rates[mask] = 0
            rhos[mask] = 0
        self.__irremovable_states = allow_irremovable_states

        pis = None
        try:
            pis = self._ct_state_frequencies(transition_rates, transmission_rates, removal_rates, rhos, upsilon)
            pis /= pis.sum()
        except Exception as e:
            logging.warning(f'Could not calculate the CT equilibrium frequencies due to {e}')

        Model.__init__(self, states=[_ for _ in model.states] + [f'{_}-{CONTACT}' for _ in model.states],
                       transition_rates=transition_rates,
                       transmission_rates=transmission_rates,
                       removal_rates=removal_rates,
                       ps=rhos,
                       n_recipients=np.concatenate([model.n_recipients, model.n_recipients]),
                       state_frequencies=pis,
                       *args, **kwargs)
        self.__upsilon = upsilon
        self.__phi = phi
        self.check_upsilon()
        self.model = model

    def _ct_state_frequencies(self, MU_IJ, LA_IJ, PSI_I, RHO_I, upsilon):
        LA_I_ = LA_IJ.sum(axis=1)
        LA__J = LA_IJ.sum(axis=0)
        LA__J_plus_LA_I_ = (LA__J + LA_I_)
        MU_I_ = MU_IJ.sum(axis=1)
        EXIT_I = MU_I_ + PSI_I
        PSI_RHO_UPS_I = PSI_I * RHO_I * upsilon
        m = len(PSI_I)
        half_m = int(m / 2)  # non-notified

        notification = []
        for k in range(half_m):
            prob_psi_j_before_k = PSI_I / (EXIT_I + EXIT_I[k])
            prob_psi_j_before_k[(EXIT_I + EXIT_I[k]) == 0] = 0

            zero_mask = LA__J_plus_LA_I_ == 0
            prob_la_bw_j_and_k = np.zeros(m, dtype=float)
            prob_la_bw_j_and_k[~zero_mask] = (LA_IJ[k, ~zero_mask] + LA_IJ[~zero_mask, k]) / LA__J_plus_LA_I_[~zero_mask]


            notification.append(prob_la_bw_j_and_k * PSI_RHO_UPS_I * prob_psi_j_before_k) # * frac_unnotified_over_k
            # notification.append(PSI_RHO_UPS_I * prob_psi_j_before_k) # * pi_k

        def func(PI_I):
            res = [PI_I.sum() - 1]
            dN_div_N_dt = PI_I.dot(LA_I_ - PSI_I)

            for k in range(half_m):
                pi_k = PI_I[k]
                pi_k_c = PI_I[k + half_m]
                frac_unnotified_over_k = pi_k / (pi_k + pi_k_c)

                dN_k_div_N_dt = (-pi_k * (MU_IJ[k, :].sum() + PSI_I[k])
                                 + PI_I.dot(MU_IJ[:, k] + LA_IJ[:, k]
                                            - notification[k] * frac_unnotified_over_k))
                res.append(pi_k * dN_div_N_dt - dN_k_div_N_dt)
            for k_c in range(half_m, m - 1):
                pi_k_c = PI_I[k_c]
                k = k_c - half_m
                pi_k = PI_I[k]
                frac_unnotified_over_k = pi_k / (pi_k + pi_k_c)

                dN_k_c_div_N_dt = (-pi_k_c * (MU_IJ[k_c, :].sum() + PSI_I[k_c])
                                 + PI_I.dot(MU_IJ[:, k_c] + LA_IJ[:, k_c]
                                            + notification[k] * frac_unnotified_over_k))
                res.append(pi_k_c * dN_div_N_dt - dN_k_c_div_N_dt)
            return res

        res = least_squares(func, x0=np.ones(m) / m, bounds=[0, 1])
        return res.x


    def clone(self):
        return CTModel(model=self.model, phi=self.__phi, upsilon=self.__upsilon,
                       allow_irremovable_states=self.__irremovable_states)

    @property
    def upsilon(self):
        return self.__upsilon

    @property
    def phi(self):
        return self.__phi

    def get_name(self):
        return self.model.get_name() + '-CT'

    def check_upsilon(self):
        assert (0 <= self.__upsilon <= 1)

    def get_epidemiological_parameters(self):
        result = self.model.get_epidemiological_parameters()
        result['upsilon'] = self.upsilon

        for state_i in self.model.states:
            if f'pi_{state_i}' in result:
                del result[f'pi_{state_i}']
        for state_i, freq in zip(self.states, self.state_frequencies):
            result[f'pi_{state_i}'] = freq

        n_contact_states = len(self.model.states)
        for state_i, phi in zip(self.states[n_contact_states: ], self.removal_rates[n_contact_states:]):
            result[f'phi_{state_i}'] = phi
            if phi:
                result[f'd_{state_i}'] = 1 / phi
        return result
