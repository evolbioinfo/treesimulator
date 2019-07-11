import unittest

import numpy as np
from numpy import random

from treesimulator.models import get_avg_rate
from treesimulator.models.naive_treated import NaiveTreatedModel
from treesimulator.models.uk_hiv import UKHivModel
from treesimulator.tree_generator import simulate_tree_gillespie


def get_simulated_sd(rates, model, time, n_obs=1e5):
    counts = np.zeros(len(model.states))
    while counts.sum() < n_obs:
        counts += simulate_tree_gillespie(model.states, rates, time, root_state=None, return_tree=False)[1]
    return counts / counts.sum()


class TestSD(unittest.TestCase):

    def test_same_rates_nt(self):
        model = NaiveTreatedModel()
        rates = np.repeat(5 * random.random(), model.num_params())
        rates = model.params2rates(rates)
        sd = rates[-1, :]
        avg_s_rt = get_avg_rate(rates[2, :], sd)
        time = 2 / avg_s_rt
        sd_sim = get_simulated_sd(rates, model, time)
        for s in model.states:
            diff = abs(sd[s.index] - sd_sim[s.index]) / sd[s.index]
            self.assertLess(diff, .2,
                            msg='Calculated eq. freq. for {} does not resemble simulated one: {} vs {}'
                            .format(s.name, sd, sd_sim))

    def test_same_rates_uk(self):
        model = UKHivModel()
        rates = np.repeat(5 * random.random(), model.num_params())
        rates = model.params2rates(rates)
        sd = rates[-1, :]
        avg_s_rt = get_avg_rate(rates[2, :], sd)
        time = 2 / avg_s_rt
        sd_sim = get_simulated_sd(rates, model, time)
        for s in model.states:
            diff = abs(sd[s.index] - sd_sim[s.index]) / sd[s.index]
            self.assertLess(diff, .25,
                            msg='Calculated eq. freq. for {} does not resemble simulated one: {} vs {}'
                            .format(s.name, sd[s.index], sd_sim[s.index]))

    def test_random_rates_nt(self):
        model = NaiveTreatedModel()
        rates = np.maximum(5 * np.random.rand(model.num_params()), .01)
        rates = model.params2rates(rates)
        # let transmission rates be at most 3 times more than sampling rates
        rates[1, :] = np.minimum(rates[2, :] * 3, rates[1, :])
        sd = model.get_sd(rates)
        avg_s_rt = get_avg_rate(rates[2, :], sd)
        time = 2 / avg_s_rt
        sd_sim = get_simulated_sd(rates, model, time)
        for s in model.states:
            diff = abs(sd[s.index] - sd_sim[s.index]) / sd[s.index]
            self.assertLess(diff, .2,
                            msg='Calculated eq. freq. for {} does not resemble simulated one: {} vs {}'
                            .format(s.name, sd[s.index], sd_sim[s.index]))

    def test_random_rates_uk(self):
        model = UKHivModel()
        rates = np.maximum(5 * np.random.rand(model.num_params()), .01)
        rates = model.params2rates(rates)
        # let transmission rates be at most 3 times more than sampling rates
        rates[1, :] = np.minimum(rates[2, :] * 3, rates[1, :])
        sd = model.get_sd(rates)
        avg_s_rt = get_avg_rate(rates[2, :], sd)
        time = 2 / avg_s_rt
        sd_sim = get_simulated_sd(rates, model, time)
        for s in model.states:
            diff = abs(sd[s.index] - sd_sim[s.index]) / sd[s.index]
            self.assertLess(diff, .25,
                            msg='Calculated eq. freq. for {} does not resemble simulated one: {} vs {}'
                            .format(s.name, sd, sd_sim))

    def test_reasonable_rates_nt(self):
        model = NaiveTreatedModel()
        rates = model.params2rates([1 / 2,
                                    2, 1 / 5,
                                    3 / 4, 3 / 4])
        sd = rates[-1, :]
        avg_s_rt = get_avg_rate(rates[2, :], sd)
        time = 2 / avg_s_rt
        sd_sim = get_simulated_sd(rates, model, time)
        for s in model.states:
            diff = abs(sd[s.index] - sd_sim[s.index]) / sd[s.index]
            self.assertLess(diff, .2,
                            msg='Calculated eq. freq. for {} does not resemble simulated one: {} vs {}'
                            .format(s.name, sd, sd_sim))

    def test_reasonable_rates_uk(self):
        model = UKHivModel()
        rates = model.params2rates([1 / 6, 1 / 3, 1 / 5,
                                    1.5, 2, 1 / 5, 1 / 2,
                                    1 / 2, 1 / 2, 1 / 7, 1 / 5])
        sd = rates[-1, :]
        avg_s_rt = get_avg_rate(rates[2, :], sd)
        time = 2 / avg_s_rt
        sd_sim = get_simulated_sd(rates, model, time)
        for s in model.states:
            diff = abs(sd[s.index] - sd_sim[s.index]) / sd[s.index]
            self.assertLess(diff, .25,
                            msg='Calculated eq. freq. for {} does not resemble simulated one: {} vs {}'
                            .format(s.name, sd, sd_sim))
