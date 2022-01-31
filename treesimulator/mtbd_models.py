import numpy as np

EXPOSED = 'e'
INFECTED = 'i'
SUPERSPREADER = 's'


SAMPLING = 'sampling'
TRANSMISSION = 'transmission'
TRANSITION = 'transition'


class Model(object):
    def __init__(self, states=None,
                 transition_rates=None, transmission_rates=None, removal_rates=None, ps=None,
                 state_frequencies=None, *args, **kwargs):
        super(Model, self).__init__()
        self.__states = np.array(states)
        num_states = len(self.states)
        self.__state_freqs = ((1 / num_states) * np.ones(num_states, dtype=float)) \
            if state_frequencies is None else np.array(state_frequencies)
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

    @property
    def states(self):
        return self.__states

    @property
    def state_frequencies(self):
        return self.__state_freqs

    @property
    def transition_rates(self):
        """
        Get transition rate matrix with states as columns and rows.

        :return rate array
        :rtype np.array
        """
        return self.__transition_rates

    @property
    def transmission_rates(self):
        """
        Get transmission rate matrix with states as columns and rows.

        :return rate array
        :rtype np.array
        """
        return self.__transmission_rates

    @property
    def removal_rates(self):
        """
        Get removal rate array with states as columns.

        :return rate array
        :rtype np.array
        """
        return self.__removal_rates

    def get_name(self):
        return 'MTBD'

    def check_rates(self):
        n_states = len(self.states)
        assert(self.transition_rates.shape == (n_states, n_states))
        assert(self.transmission_rates.shape == (n_states, n_states))
        assert(self.removal_rates.shape == (n_states,))
        assert(self.ps.shape == (n_states,))
        assert(np.all(self.transition_rates >= 0))
        assert(np.all(self.transmission_rates >= 0))
        assert(np.all(self.removal_rates >= 0))
        assert(np.all(self.ps >= 0) and np.all(self.ps <= 1))

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
                       transition_rates=mus, transmission_rates=las, removal_rates=psis, ps=[0, p], *args, **kwargs)

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
        if la_ss / la_ns != la_sn / la_nn:
            raise ValueError('transmission ratio constraint is violated: la_ss / la_ns must be equal to la_sn / la_nn')
        las[0, 0] = la_nn
        las[0, 1] = la_ns
        las[1, 0] = la_sn
        las[1, 1] = la_ss
        psis = psi * np.ones(shape=2, dtype=np.float)
        f = la_ss / (la_ss + la_sn)
        Model.__init__(self, states=[INFECTED, SUPERSPREADER],
                       transmission_rates=las, removal_rates=psis, ps=[p, p],
                       state_frequencies=[1 - f, f], *args, **kwargs)

    def get_name(self):
        return 'BDSS'

    def get_epidemiological_parameters(self):
        """Converts rate parameters to the epidemiological ones"""
        return {'R0': (self.transmission_rates[0, 0] + self.transition_rates[1, 1]) / self.removal_rates[0],
                'infectious time': 1 / self.removal_rates[1],
                'sampling probability': self.ps[1],
                'superspreading transmission ratio': self.transmission_rates[1, 1] / self.transmission_rates[0, 1],
                'superspreading fraction':
                    self.transmission_rates[1, 1] / (self.transmission_rates[1, 1] + self.transmission_rates[1, 0])}
