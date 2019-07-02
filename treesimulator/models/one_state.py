import logging

import numpy as np
import pandas as pd
from scipy import integrate
from scipy.optimize import minimize

from treesimulator import TIME_TILL_NOW
from treesimulator.models import Model, SAMPLING, TRANSMISSION, TRANSITION, State

STATIONARY_DIST = np.array([1])

MU_RATES = np.array([0])

S = 's'


class OneStateModel(Model):

    def __init__(self):
        super(OneStateModel, self).__init__()

    def num_params(self, type=None):
        if SAMPLING == type:
            return 1
        if TRANSMISSION == type:
            return 1
        if TRANSITION == type:
            return 0
        return 2

    def _get_states(self):
        s = State(name=S, index=0)
        s.recipient = s
        return np.array([s])

    def params2rates(self, ps, rates=None, nested_rates=None, sampled_pis=None, **kwargs):
        if rates is None:
            rates = np.zeros(shape=(4, 1), dtype=np.float)
        rates[1:3, 0] = ps
        rates[-1, 0] = 1
        return rates

    def map_state(self, state):
        return self.states[0]

    def get_name(self):
        return 'one-state'

    @staticmethod
    def get_error(ps, num_tips, forest, sampled_fraction):
        lambda_rate = ps[0]
        psi_rate = lambda_rate * sampled_fraction
        time_till_now = next((getattr(tree, TIME_TILL_NOW) + tree.dist / 2) for tree in forest if tree is not None)
        exponent = np.exp((lambda_rate + psi_rate) * time_till_now)

        error = num_tips / len(forest) * (exponent - 1) / (1 + sampled_fraction * exponent) * np.log(
            1 - sampled_fraction) - np.log(1 + sampled_fraction) + np.log(1 + sampled_fraction * exponent)
        return np.abs(error)

    @staticmethod
    def get_too_pdf(num_lineages_now, lambda_rate, psi_rate, t):
        if lambda_rate == psi_rate:
            return np.log(num_lineages_now) + np.log(t) * (num_lineages_now - 1) - np.log(1 + t) * (
                        num_lineages_now + 1)
        return np.log(num_lineages_now) + np.log(lambda_rate) * num_lineages_now + np.log(lambda_rate - psi_rate) * 2 \
               + np.log(1 - np.exp(-(lambda_rate - psi_rate) * t)) * (num_lineages_now - 1) \
               - (lambda_rate - psi_rate) * t \
               - np.log(lambda_rate - psi_rate * np.exp(-(lambda_rate - psi_rate) * t)) * (num_lineages_now + 1)

    @staticmethod
    def get_time_of_origin_pdf(num_lineages_now, lambda_rate, psi_rate, time_till_now):

        def get_v(ps):
            if np.any(pd.isnull(ps)):
                return np.nan
            t = ps[0]

            res = OneStateModel.get_too_pdf(num_lineages_now, lambda_rate, psi_rate, t)
            return np.inf if pd.isnull(res) else -res

        bounds = np.array([np.array([time_till_now, time_till_now + 100])])

        for _ in range(5):
            logging.debug('Bounds are:\n {}\n'.format(bounds))

            vs = np.around(np.random.uniform(bounds[:, 0], bounds[:, 1]), 2)

            fres = minimize(get_v, x0=vs, method='SLSQP', bounds=bounds, options={'eps': 0.001})
            if fres.success and not np.any(np.isnan(fres.x)):
                return fres.x[0]
        return None

    @staticmethod
    def get_time_of_origin_p(num_lineages_now, lambda_rate, psi_rate, time_till_now):

        def get_v(ps):
            if np.any(pd.isnull(ps)):
                return np.nan
            t = ps[0]

            # time_proportion = (t - time_till_now) / t
            # res = np.abs(OneStateModel.get_num_lineages(num_lineages_now, time_proportion, rho) - 1)
            # logging.info('TOO: {}, {} -> {}'.format(t, rho, res))
            # return np.inf if pd.isnull(res) else res

            res = np.log(lambda_rate) * (num_lineages_now - 1) + np.log(lambda_rate - psi_rate) * 2 \
                  + np.log(1 - np.exp(-(lambda_rate - psi_rate) * t)) * (num_lineages_now - 1) \
                  - (lambda_rate - psi_rate) * t \
                  - np.log(lambda_rate - psi_rate * np.exp(-(lambda_rate - psi_rate) * t)) * (num_lineages_now + 1)

            # res = np.power(lambda_rate, n - 1) * np.power(lambda_rate - psi_rate, 2) \
            #        * np.power(1 - np.exp(-(lambda_rate - psi_rate) * t), n - 1) \
            #       * (np.exp(-(lambda_rate - psi_rate) * t)
            #          / np.power(lambda_rate - psi_rate * np.exp(-(lambda_rate - psi_rate) * t), n + 1))
            return np.inf if pd.isnull(res) else -res

        bounds = np.array([np.array([time_till_now, time_till_now + 100])])

        for _ in range(5):
            logging.debug('Bounds are:\n {}\n'.format(bounds))

            vs = np.around(np.random.uniform(bounds[:, 0], bounds[:, 1]), 2)

            fres = minimize(get_v, x0=vs, method='SLSQP', bounds=bounds, options={'eps': 0.001})
            if fres.success and not np.any(np.isnan(fres.x)):
                logging.debug('Estimated TOO as {} via P ({})'.format(fres.x, np.exp(-fres.fun)))
                return fres.x[0]
        return None

    @staticmethod
    def get_num_trees(total_num_tips, num_sampled_tips, lambda_rate, psi_rate, time_till_now, num_visible_trees):

        def get_v_tips_per_tree(ps):
            if np.any(pd.isnull(ps)):
                return np.nan
            num_lineages_now = ps[0]
            res = OneStateModel.get_too_pdf(num_lineages_now, lambda_rate, psi_rate, time_till_now)
            return np.inf if pd.isnull(res) else -res

        bounds = np.array([np.array([num_sampled_tips / num_visible_trees, total_num_tips / num_visible_trees])])

        num_tips_per_visible_tree = None
        for _ in range(5):
            logging.debug('Bounds are:\n {}\n'.format(bounds))

            vs = np.around(np.random.uniform(bounds[:, 0], bounds[:, 1]), 2)

            fres = minimize(get_v_tips_per_tree, x0=vs, method='SLSQP', bounds=bounds, options={'eps': .01})
            if fres.success and not np.any(np.isnan(fres.x)):
                num_tips_per_visible_tree = fres.x[0]

        if num_tips_per_visible_tree is None:
            return None

        logging.info('\tEstimated {} tips per visible tree'.format(num_tips_per_visible_tree))

        num_hidden_tips_in_hidden_trees = total_num_tips - num_tips_per_visible_tree * num_visible_trees

        def get_v_trees(ps):
            if np.any(pd.isnull(ps)):
                return np.nan
            num_lineages_now = num_hidden_tips_in_hidden_trees / (ps[0] - num_visible_trees)
            res = OneStateModel.get_too_pdf(num_lineages_now, lambda_rate, psi_rate, time_till_now)
            return np.inf if pd.isnull(res) else -res

        bounds = np.array([np.array([num_visible_trees, num_visible_trees + num_hidden_tips_in_hidden_trees])])

        num_trees = None
        for _ in range(5):
            logging.debug('Bounds are:\n {}\n'.format(bounds))

            vs = np.around(np.random.uniform(bounds[:, 0], bounds[:, 1]), 2)

            fres = minimize(get_v_trees, x0=vs, method='SLSQP', bounds=bounds, options={'eps': 1})
            if fres.success and not np.any(np.isnan(fres.x)):
                num_trees = fres.x[0]

        if num_trees is None:
            return None

        return num_trees

    @staticmethod
    def get_num_trees_too(num_lineages_now, lambda_rate, psi_rate, time_till_now):
        t = OneStateModel.get_time_of_origin_pdf(num_lineages_now, lambda_rate, psi_rate, time_till_now)
        # t = OneStateModel.get_time_of_origin_p(num_lineages_now, lambda_rate, psi_rate, time_till_now)
        # logging.info('Estimated TOO as {} via p'.format(t))
        if pd.isnull(t):
            return None
        rho = psi_rate / lambda_rate
        time_proportion = (t - time_till_now) / t
        return np.round(OneStateModel.get_num_lineages(num_lineages_now, time_proportion, rho))

    @staticmethod
    def get_num_lineages(num_lineages_now, time_proportion, rho):

        def get_v(t):
            log_res = np.log(np.exp(-(2 - time_proportion) * t) - np.exp(-2 * t)) \
                      - np.log(1 - rho * np.exp(-(1 - time_proportion) * t)) \
                      + np.log(1 - np.exp(-t)) * (num_lineages_now - 2) \
                      - np.log(1 - rho * np.exp(-t)) * (num_lineages_now + 1)
            return np.power(np.e, log_res)

        return 1 + num_lineages_now * (num_lineages_now - 1) * np.power(1 - rho, 2) \
               * integrate.quad(get_v, 0, np.inf)[0]
