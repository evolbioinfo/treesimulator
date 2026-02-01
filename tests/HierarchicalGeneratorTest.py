import logging
import unittest

import numpy as np

from generator import simulate_epidemic_gillespie, reconstruct_tree, annotate_tree_with_time
from treesimulator.hierarchical_generator import generate_hierarchically
from treesimulator.mtbd_models import BirthDeathModel
from treesimulator import TIME
from bdct.bd_model import infer

logging.getLogger().handlers = []
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s | %(levelname)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")

class HierarchicalGeneratorTest(unittest.TestCase):

    def test_bd(self):
        R = 1 + np.random.random() * 9
        psi = 0.5 + np.random.random()
        rho = 0.1 + 0.9 * np.random.random()
        la = R * psi
        model = BirthDeathModel(la=la, psi=psi, p=rho)
        N = 10000

        logging.info(f'R={R}, la={la}, psi={psi}, rho={rho}')
        tree =  generate_hierarchically([model], max_sampled=N,
                                        state_frequencies=None, root_state=None)

        # sampled_ids = []
        # while len(sampled_ids) < N:
        #     id2parent_id, id2node_time, sampled_ids, id2state, max_time, observed_nums = \
        #         simulate_epidemic_gillespie([model], max_sampled=N)
        # tree2 = reconstruct_tree(id2parent_id, id2node_time, sampled_ids, id2state, max_time)
        # annotate_tree_with_time(tree2)
        # logging.info('Generated a tree with {} sampled tips'.format(len(tree2)))
        #
        # _, [[la_min2, la_max2], [psi_min2, psi_max2], [_, _]] = \
        #     infer([tree2], max(getattr(_, TIME) for _ in tree2), la=None, psi=None, p=rho, ci=True, num_attemps=1)

        _, [[la_min, la_max], [psi_min, psi_max], [_, _]] = \
            infer([tree], max(getattr(_, TIME) for _ in tree), la=None, psi=None, p=rho, ci=True, num_attemps=1)

        # self.assertAlmostEqual(la_est2, la, places=2)
        # self.assertAlmostEqual(psi_est2, psi, places=2)
        # self.assertLess(la_min2, la)
        # self.assertLess(la, la_max2)
        # self.assertLess(psi_min2, psi)
        # self.assertLess(psi, psi_max2)

        # self.assertAlmostEqual(la_est, la, places=2)
        # self.assertAlmostEqual(psi_est, psi, places=2)
        self.assertLess(la_min, la)
        self.assertLess(la, la_max)
        self.assertLess(psi_min, psi)
        self.assertLess(psi, psi_max)