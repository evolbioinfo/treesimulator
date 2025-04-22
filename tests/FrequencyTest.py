import unittest

import numpy as np

from treesimulator.mtbd_models import BirthDeathExposedInfectiousModel, BirthDeathWithSuperSpreadingModel, \
    BirthDeathModel, Model, EXPOSED, INFECTIOUS, SUPERSPREADER



class FrequencyTest(unittest.TestCase):

    def test_bd(self):
        model = BirthDeathModel(la=np.random.random() * 10, psi=np.random.random() * 10, p=np.random.random())
        self.assertAlmostEqual(1, model.state_frequencies[0], places=6)

    def test_bdei_R_is_one(self):
        mu = np.random.random() * 10
        psi = np.random.random() * 10
        la = psi
        model = BirthDeathExposedInfectiousModel(mu=mu, la=la, psi=psi, p=np.random.random())
        pi_I = mu / (mu + psi)
        self.assertAlmostEqual(1 - pi_I, model.state_frequencies[0], places=6)
        self.assertAlmostEqual(pi_I, model.state_frequencies[1], places=6)

    def test_bdei_R_not_one(self):
        mu = np.random.random() * 10
        psi = np.random.random() * 10
        la = np.random.random() * 10
        model = BirthDeathExposedInfectiousModel(mu=mu, la=la, psi=psi, p=np.random.random())
        pi_I = (-(mu + psi) + np.sqrt(np.power(mu - psi, 2) + 4 * mu * la)) / (2 * (la - psi))
        self.assertAlmostEqual(1 - pi_I, model.state_frequencies[0], places=6)
        self.assertAlmostEqual(pi_I, model.state_frequencies[1], places=6)

    def test_bdss(self):
        psi = np.random.random() * 10
        la = np.random.random() * 10
        f_ss = np.random.random()
        x_ss = 1 + np.random.random() * 24
        model = BirthDeathWithSuperSpreadingModel(la_sn=x_ss * la * (1 - f_ss), la_ss=x_ss * la * f_ss,
                                                  la_nn=(1 - f_ss) * la, la_ns=f_ss * la,
                                                  psi=psi, p=np.random.random())
        self.assertEqual(1 - f_ss, model.state_frequencies[0])
        self.assertEqual(f_ss, model.state_frequencies[1])

    def test_mtbd_bdei_R_is_one(self):
        mu = np.random.random() * 10
        psi = np.random.random() * 10
        la = psi
        rho = np.random.random()
        model = Model(states=[EXPOSED, INFECTIOUS],
                      transition_rates=np.array([[0, mu], [0, 0]]),
                      transmission_rates=np.array([[0, 0], [la, 0]]),
                      removal_rates=np.array([0, psi]),
                      ps=np.array([0, rho]))
        bdei_model = BirthDeathExposedInfectiousModel(mu=mu, la=la, psi=psi, p=rho)
        print(mu, la, psi, rho)
        self.assertAlmostEqual(bdei_model.state_frequencies[0], model.state_frequencies[0], places=6)
        self.assertAlmostEqual(bdei_model.state_frequencies[1], model.state_frequencies[1], places=6)

    def test_mtbd_bdei_R_not_one(self):
        mu = np.random.random() * 10
        psi = np.random.random() * 10
        la = np.random.random() * 10
        rho = np.random.random()
        model = Model(states=[EXPOSED, INFECTIOUS],
                      transition_rates=np.array([[0, mu], [0, 0]]),
                      transmission_rates=np.array([[0, 0], [la, 0]]),
                      removal_rates=np.array([0, psi]),
                      ps=np.array([0, rho]))
        bdei_model = BirthDeathExposedInfectiousModel(mu=mu, la=la, psi=psi, p=rho)
        bdei_pis = bdei_model.state_frequencies
        pis = model.state_frequencies
        self.assertAlmostEqual(bdei_pis[0], pis[0], places=6)
        self.assertAlmostEqual(bdei_pis[1], pis[1], places=6)

    def test_mtbd_bdss(self):
        psi = np.random.random() * 10
        la = np.random.random() * 10
        f_ss = np.random.random()
        x_ss = 1 + np.random.random() * 24
        rho = np.random.random()
        la_nn = (1 - f_ss) * la
        la_ns = f_ss * la
        la_sn = x_ss * la * (1 - f_ss)
        la_ss = x_ss * la * f_ss
        bdss_model = BirthDeathWithSuperSpreadingModel(la_sn=la_sn, la_ss=la_ss,
                                                       la_nn=la_nn, la_ns=la_ns,
                                                       psi=psi, p=rho)
        model = Model(states=[INFECTIOUS, SUPERSPREADER],
                      transmission_rates=np.array([[la_nn, la_ns], [la_sn, la_ss]]),
                      removal_rates=np.array([psi, psi]),
                      ps=np.array([rho, rho]))
        bdss_pis = bdss_model.state_frequencies
        pis = model.state_frequencies
        self.assertAlmostEqual(bdss_pis[0], pis[0], places=6)
        self.assertAlmostEqual(bdss_pis[1], pis[1], places=6)
