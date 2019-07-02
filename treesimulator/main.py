import logging

from pandas.tests.extension.numpy_.test_numpy_nested import np

from treesimulator import STATE
from treesimulator.models.naive_treated import NaiveTreatedModel
from treesimulator.tree_generator import simulate_tree_gillespie

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Simulates a tree.")
    parser.add_argument('--n_tips', default=200, type=int, help="desired number of simulated tips "
                                                                "(set to infinity to on time only)")
    parser.add_argument('--time', default=np.inf, type=float,
                        help="max time over which the simulation is to be done "
                             "(set to infinity to condition on the number of tips only)")
    parser.add_argument('--log', default=logging.INFO, type=str,
                        help="amount of log information to display: set to logging.DEBUG to see more")
    parser.add_argument('--rates', default=[1 / 2, 1, 1 / 5, 1 / 2, 1 / 8], type=float, nargs=5,
                        help='parameter values to be used for simulation: '
                             'state_change_rate, transmission_rate_naive, transmission_rate_treated,'
                             ' sampling_rate_naive, sampling_rate_treated')
    parser.add_argument('--nwk', type=str, required=True,
                        help='path to the file where to save the simulated tree')
    params = parser.parse_args()
    logging.basicConfig(level=params.log, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    model = NaiveTreatedModel()
    rates = model.params2rates(params.rates)

    logging.info('Rates are:\n{}\n'.format(rates))

    tree = simulate_tree_gillespie(model.states, rates, max_time=params.time, max_sampled=params.n_tips,
                                   state_feature=STATE)
    logging.info('Simulated a tree with {} sampled tips, under constraint of max {} tips and over time {}.'
                 .format(len(tree) if tree else 0, params.n_tips, params.time))
    if tree:
        tree.write(features=[STATE], outfile=params.nwk)
