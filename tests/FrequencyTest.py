import unittest

import numpy as np

from treesimulator.mtbd_models import BirthDeathExposedInfectiousModel, BirthDeathWithSuperSpreadingModel, \
    BirthDeathModel, Model, EXPOSED, INFECTIOUS, SUPERSPREADER, BirthDeathExposedInfectiousWithSuperSpreadingModel, \
    CTModel


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
        self.assertAlmostEqual(np.float64(1 - pi_I), model.state_frequencies[0], places=6)
        self.assertAlmostEqual(np.float64(pi_I), model.state_frequencies[1], places=6)

    def test_bdss(self):
        psi = np.random.random() * 10
        la = np.random.random() * 10
        f_ss = np.random.random()
        x_ss = 1 + np.random.random() * 24
        model = BirthDeathWithSuperSpreadingModel(la_sn=x_ss * la * (1 - f_ss), la_ss=x_ss * la * f_ss,
                                                  la_nn=(1 - f_ss) * la, la_ns=f_ss * la,
                                                  psi=psi, p=np.random.random())
        self.assertAlmostEqual(np.float64(1 - f_ss), model.state_frequencies[0], places=6)
        self.assertAlmostEqual(np.float64(f_ss), model.state_frequencies[1], places=6)

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


    def test_bdeiss_mu_inf_vs_bdss(self):
        mu = 1e8
        psi = np.random.random() * 10
        la = np.random.random() * 10
        f_ss = np.random.random()
        x_ss = 1 + np.random.random() * 24
        rho = np.random.random()
        mu_i = (1 - f_ss) * mu
        mu_s = f_ss * mu
        bdeiss_model = BirthDeathExposedInfectiousWithSuperSpreadingModel(mu_n=mu_i, mu_s=mu_s, la_n=la, la_s=x_ss * la, \
                                                                          psi=psi, p=rho)
        bdss_model = BirthDeathWithSuperSpreadingModel(la_nn=(1 - f_ss) * la, la_ns=f_ss * la, \
                                                       la_sn=(1 - f_ss) * x_ss * la, la_ss=f_ss * x_ss * la, \
                                                       psi=psi, rho=rho)
        bdeiss_pis = bdeiss_model.state_frequencies[1:]
        pis = bdss_model.state_frequencies
        print(bdeiss_pis, pis)
        self.assertAlmostEqual(bdeiss_pis[0], pis[0], places=6)
        self.assertAlmostEqual(bdeiss_pis[1], pis[1], places=6)


    def test_bdss_f_ss_0_vs_bd(self):
        psi = np.random.random() * 10
        la = np.random.random() * 10
        f_ss = 0
        x_ss = 1
        rho = np.random.random()
        bdss_model = BirthDeathWithSuperSpreadingModel(la_nn=(1 - f_ss) * la, la_ns=f_ss * la, \
                                                       la_sn=(1 - f_ss) * x_ss * la, la_ss=f_ss * x_ss * la, \
                                                       psi=psi, rho=rho)
        bdss_pis = bdss_model.state_frequencies
        print(bdss_pis)
        self.assertAlmostEqual(bdss_pis[0], 1, places=6)

    def test_bdei_mu_inf_vs_bd(self):
        psi = np.random.random() * 10
        la = np.random.random() * 10
        rho = np.random.random()
        bdei_model = BirthDeathExposedInfectiousModel(mu=1e8, la=la, psi=psi, rho=rho)
        bdei_pis = bdei_model.state_frequencies
        print(bdei_pis)
        self.assertAlmostEqual(bdei_pis[1], 1, places=6)


    def test_bdeiss_f_ss_0_vs_bdei(self):
        mu = np.random.random() * 10
        psi = np.random.random() * 10
        la = np.random.random() * 10
        f_ss = 0
        x_ss = 1
        rho = np.random.random()
        mu_i = (1 - f_ss) * mu
        mu_s = f_ss * mu
        bdeiss_model = BirthDeathExposedInfectiousWithSuperSpreadingModel(mu_n=mu_i, mu_s=mu_s, la_n=la, la_s=x_ss * la, \
                                                                          psi=psi, p=rho)
        bdei_model = BirthDeathExposedInfectiousModel(mu=mu, la=la, psi=psi, rho=rho)
        bdeiss_pis = bdeiss_model.state_frequencies
        pis = bdei_model.state_frequencies
        print(bdeiss_pis, pis)
        self.assertAlmostEqual(bdeiss_pis[0], pis[0], places=6)
        self.assertAlmostEqual(bdeiss_pis[1], pis[1], places=6)


    def test_bdeissct_ups_0_vs_bdeiss(self):
        mu = np.random.random() * 10
        psi = np.random.random() * 10
        la = np.random.random() * 10
        f_ss = np.random.random()
        x_ss = 1 + np.random.random() * 24
        x_c = 1 + np.random.random() * 250
        rho = np.random.random()
        mu_i = (1 - f_ss) * mu
        mu_s = f_ss * mu
        bdeiss_model = BirthDeathExposedInfectiousWithSuperSpreadingModel(mu_n=mu_i, mu_s=mu_s, la_n=la, la_s=x_ss * la, \
                                                                          psi=psi, p=rho)
        bdeissct_model = CTModel(bdeiss_model, phi=psi * x_c, upsilon=0)
        bdeiss_pis = bdeiss_model.state_frequencies
        pis = bdeissct_model.state_frequencies
        print(bdeiss_pis, pis)
        self.assertAlmostEqual(bdeiss_pis[0], pis[0], places=3)
        self.assertAlmostEqual(bdeiss_pis[1], pis[1], places=3)
        self.assertAlmostEqual(bdeiss_pis[2], pis[2], places=3)
