import logging

import numpy as np

from treesimulator import save_forest, save_log
from treesimulator.generator import generate
from treesimulator.mtbd_models import Model


def main():
    """
    Entry point for tree/forest generation with a generic MTBD model with command-line arguments.
    :return: void
    """
    import argparse

    parser = \
        argparse.ArgumentParser(description="Simulates a tree (or a forest of trees) for given MTBD model parameters. "
                                            "If a simulation leads to less than --min_tips tips, it is repeated.")
    parser.add_argument('--min_tips', required=True, type=int,
                        help="desired minimal bound on the total number of simulated leaves. "
                             "For a tree simulation, "
                             "if --min_tips and --max_tips are equal, exactly that number of tips will be simulated. "
                             "If --min_tips is less than --max_tips, "
                             "a value randomly drawn between one and another will be simulated.")
    parser.add_argument('--max_tips', required=True, type=int,
                        help="desired maximal bound on the total number of simulated leaves"
                             "For a tree simulation, "
                             "if --min_tips and --max_tips are equal, exactly that number of tips will be simulated. "
                             "If --min_tips is less than --max_tips, "
                             "a value randomly drawn between one and another will be simulated.")
    parser.add_argument('--T', required=False, default=np.inf, type=float,
                        help="Total simulation time. If specified, a forest will be simulated instead of one tree. "
                             "The trees in this forest will be simulated during the given time, "
                             "till the --min_tips number is reached. If after simulating the last tree, "
                             "the forest exceeds the --max_tips number, the process will be restarted.")
    parser.add_argument('--states', nargs='+', type=str, help="model states")
    parser.add_argument('--transition_rates', nargs='+', type=float,
                        help="transition rate matrix, row after row, in the same order as model states, "
                             "e.g. ig a model has 2 states given as --states E I,"
                             "then here we expect E->E E->I I->E I->I.")
    parser.add_argument('--transmission_rates', nargs='+', type=float,
                        help="transmission rate matrix, row after row, in the same order as model states, "
                             "e.g. ig a model has 2 states given as --states E I,"
                             "then here we expect E->E E->I I->E I->I.")
    parser.add_argument('--removal_rates', nargs='+', type=float,
                        help="removal rate array, in the same order as model states, "
                             "e.g. ig a model has 2 states given as --states E I,"
                             "then here we expect removal(E) removal(I).")
    parser.add_argument('--sampling_probabilities', nargs='+', type=float,
                        help="sampling probability array, in the same order as model states, "
                             "e.g. ig a model has 2 states given as --states E I,"
                             "then here we expect p(E) p(I).")
    parser.add_argument('--log', required=True, default=None, type=str, help="output log file")
    parser.add_argument('--nwk', required=True, default=None, type=str, help="output tree or forest file")
    params = parser.parse_args()
    logging.getLogger().handlers = []
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    logging.info(
        'MTBD model parameters are:\n\ttransition_rates=={}\n\ttransmission_rates={}\n\ttransmission_rates={}\n\tps={}'
            .format(params.transition_rates, params.transmission_rates,
                    params.removal_rates, params.sampling_probabilities))
    logging.info('Total time T={}'.format(params.T))

    n_states = len(params.states)
    model = Model(states=params.states,
                  transmission_rates=np.array(params.transmission_rates).reshape((n_states, n_states)),
                  transition_rates=np.array(params.transition_rates).reshape((n_states, n_states)),
                  removal_rates=params.removal_rates, ps=params.sampling_probabilities)

    forest, (total_tips, u, T) = generate(model, params.min_tips, params.max_tips, T=params.T, state_frequencies=[1])

    save_forest(forest, params.nwk)
    save_log(model, total_tips, T, u, params.log)


if '__main__' == __name__:
    main()
