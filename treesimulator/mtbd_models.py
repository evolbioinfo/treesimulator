import numpy as np
from scipy.optimize import fsolve

EXPOSED = 'e'
INFECTED = 'i'
SUPERSPREADER = 's'

SAMPLING = 'sampling'
TRANSMISSION = 'transmission'
TRANSITION = 'transition'


class Model(object):

    def __init__(self, states=None,
                 transition_rates=None, transmission_rates=None, removal_rates=None, ps=None, *args, **kwargs):
        super(Model, self).__init__()
        self.__states = np.array(states)
        num_states = len(self.states)
        self.__ps = np.array(ps) if ps is not None else np.ones(num_states, dtype=float)
        self.__transmission_rates = np.array(transmission_rates) if transmission_rates is not None \
            else np.zeros(shape=(num_states, num_states), dtype=np.float)
        self.__transition_rates = np.array(transition_rates) if transition_rates is not None \
            else np.zeros(shape=(num_states, num_states), dtype=np.float)
        self.__removal_rates = np.array(removal_rates) if removal_rates is not None \
            else np.zeros(shape=num_states, dtype=np.float)
        self.check_rates()

    def clone(self):
        return Model(self.states, self.transition_rates, self.transmission_rates, self.removal_rates, self.ps)

    @property
    def ps(self):
        return self.__ps

    @ps.setter
    def ps(self, ps):
        self.__ps = ps

    @property
    def states(self):
        return self.__states

    @property
    def state_frequencies(self):
        MU, LA, PSI = self.transition_rates, self.transmission_rates, self.removal_rates
        m = len(self.states)

        def func(PI):
            SIGMA = PI.dot(LA.sum(axis=1) - PSI)
            res = [PI.sum() - 1]
            for k in range(m - 1):
                pi_k = PI[k]
                res.append(pi_k * SIGMA + pi_k * (PSI[k] + MU[k, :].sum()) - PI.dot(MU[:, k] + LA[:, k]))
            return res

        PI = fsolve(func, np.ones(m) / m)
        if np.any(PI < 0) or np.any(PI > 1):
            return np.ones(m) / m
        return PI

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

    def get_name(self):
        return 'MTBD'

    def check_rates(self):
        n_states = len(self.states)
        assert (self.transition_rates.shape == (n_states, n_states))
        assert (self.transmission_rates.shape == (n_states, n_states))
        assert (self.removal_rates.shape == (n_states,))
        assert (self.ps.shape == (n_states,))
        assert (np.all(self.transition_rates >= 0))
        assert (np.all(self.transmission_rates >= 0))
        assert (np.all(self.removal_rates >= 0))
        assert (np.all(self.ps >= 0) and np.all(self.ps <= 1))

    def get_epidemiological_parameters(self):
        """Converts rate parameters to the epidemiological ones"""
        return {'transitions': self.transition_rates, 'transmissions': self.transmission_rates,
                'removals': self.removal_rates, 'sampling': self.ps}


class BirthDeathExposedInfectiousModel(Model):

    def __init__(self, mu, la, psi, p=0.5, *args, **kwargs):
        """
        :param mu: transition E->I
        :param la: transmission from I to E
        :param psi: removal of I
        :param p: sampling probability of I
        """
        mus = np.zeros(shape=(2, 2), dtype=np.float)
        mus[0, 1] = mu
        las = np.zeros(shape=(2, 2), dtype=np.float)
        las[1, 0] = la
        psis = np.zeros(shape=2, dtype=np.float)
        psis[1] = psi

        Model.__init__(self, states=[EXPOSED, INFECTED],
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
        return {'R0': self.transmission_rates[1, 0] / self.removal_rates[1],
                'infectious time': 1 / self.removal_rates[1],
                'incubation period': 1 / self.transition_rates[0, 1],
                'sampling probability': self.ps[1]}


class BirthDeathModel(Model):

    def __init__(self, la, psi, p=0.5, *args, **kwargs):
        """
        :param la: transmission
        :param psi: removal
        :param p: sampling probability
        """
        las = la * np.ones(shape=(1, 1), dtype=np.float)
        Model.__init__(self, states=[INFECTED], transmission_rates=las, removal_rates=[psi], ps=[p], *args, **kwargs)

    def get_name(self):
        return 'BD'

    def get_epidemiological_parameters(self):
        """Converts rate parameters to the epidemiological ones"""
        return {'R0': self.transmission_rates[0, 0] / self.removal_rates[0],
                'infectious time': 1 / self.removal_rates[0],
                'sampling probability': self.ps[0]}


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
        las = np.zeros(shape=(2, 2), dtype=np.float)
        s_ratio = la_ss / la_ns
        n_ratio = la_sn / la_nn
        if np.abs(s_ratio - n_ratio) > 1e-3:
            raise ValueError(
                'transmission ratio constraint is violated: la_ss / la_ns ({}) must be equal to la_sn / la_nn ({})'
                    .format(s_ratio, n_ratio))
        las[0, 0] = la_nn
        las[0, 1] = la_ns
        las[1, 0] = la_sn
        las[1, 1] = la_ss
        psis = psi * np.ones(shape=2, dtype=np.float)
        Model.__init__(self, states=[INFECTED, SUPERSPREADER],
                       transmission_rates=las, removal_rates=psis, ps=[p, p], *args, **kwargs)

    @property
    def state_frequencies(self):
        la_ss = self.transmission_rates[1, 1]
        la_sn = self.transmission_rates[1, 0]
        f = la_ss / (la_ss + la_sn)
        return np.array([1 - f, f])

    def get_name(self):
        return 'BDSS'

    def get_epidemiological_parameters(self):
        """Converts rate parameters to the epidemiological ones"""
        return {'R0': (self.transmission_rates[0, 0] + self.transmission_rates[1, 1]) / self.removal_rates[0],
                'infectious time': 1 / self.removal_rates[1],
                'sampling probability': self.ps[1],
                'superspreading transmission ratio': self.transmission_rates[1, 1] / self.transmission_rates[0, 1],
                'superspreading fraction':
                    self.transmission_rates[1, 1] / (self.transmission_rates[1, 1] + self.transmission_rates[1, 0])}


class PNModel(Model):
    def __init__(self, model, removal_rate=np.inf, pn=0.5, *args, **kwargs):
        Model.__init__(self, states=np.pad(model.states, (0, 1), mode='constant', constant_values='n'),
                       transition_rates=np.pad(model.transition_rates, ((0, 1), (0, 1)),
                                               mode='constant', constant_values=0),
                       transmission_rates=np.pad(model.transmission_rates, ((0, 1), (0, 1)),
                                                 mode='constant', constant_values=0),
                       removal_rates=np.pad(model.removal_rates, (0, 1), mode='constant', constant_values=removal_rate),
                       ps=np.pad(model.ps, (0, 1), mode='constant', constant_values=1),
                       state_frequencies=np.pad(model.state_frequencies, (0, 1), mode='constant', constant_values=0),
                       *args, **kwargs)
        self.__pn = pn
        self.check_p()
        self.model = model

    def clone(self):
        return PNModel(self.model, self.removal_rates[-1], self.__pn)

    @property
    def pn(self):
        return self.__pn

    @property
    def partner_removal_rate(self):
        """
        Get partner removal rate

        :return partner removal rate
        :rtype np.float
        """
        return self.removal_rates[-1]

    def get_name(self):
        return self.model.get_name() + '-PN'

    def check_p(self):
        assert (0 <= self.__pn <= 1)

    def get_epidemiological_parameters(self):
        res = self.model.get_epidemiological_parameters()
        res['notification probability'] = self.pn
        res['removal time after notification'] = 1 / self.removal_rates[-1]
        return res
