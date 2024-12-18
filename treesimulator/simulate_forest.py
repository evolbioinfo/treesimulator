import logging

import numpy as np

from treesimulator import save_forest, save_log, save_ltt
from treesimulator.generator import generate, observed_ltt
from treesimulator.mtbd_models import Model, CTModel


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
                             "e.g. if a model has 2 states given as --states E I,"
                             "then here we expect E->E E->I I->E I->I.")
    parser.add_argument('--transmission_rates', nargs='+', type=float,
                        help="transmission rate matrix, row after row, in the same order as model states, "
                             "e.g. if a model has 2 states given as --states E I,"
                             "then here we expect E->E E->I I->E I->I.")
    parser.add_argument('--removal_rates', nargs='+', type=float,
                        help="removal rate array, in the same order as model states, "
                             "e.g. if a model has 2 states given as --states E I,"
                             "then here we expect removal(E) removal(I).")
    parser.add_argument('--sampling_probabilities', nargs='+', type=float,
                        help="sampling probability array, in the same order as model states, "
                             "e.g. if a model has 2 states given as --states E I,"
                             "then here we expect p(E) p(I).")
    parser.add_argument('--upsilon', required=False, default=0, type=float, help='notification probability')
    parser.add_argument('--phi', required=False, default=0, type=float, help='notified removal rate')
    parser.add_argument('--max_notified_contacts', required=False, default=1, type=int,
                        help='maximum number of notified contracts per person')
    parser.add_argument('--avg_recipients', nargs='*', type=float,
                        help='average number of recipients per transmission '
                             'for each donor state (in the same order as the model states). '
                             'By default, only one-to-one transmissions are allowed, '
                             'but if larger numbers are given then one-to-many transmissions become possible.')
    parser.add_argument('--log', required=True, type=str, help="output log file")
    parser.add_argument('--nwk', required=True, type=str, help="output tree or forest file")
    parser.add_argument('--ltt', required=False, default=None, type=str, help="output LTT file")
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="print information on the progress of the tree generation (to console)")
    params = parser.parse_args()
    logging.getLogger().handlers = []
    logging.basicConfig(level=logging.DEBUG if params.verbose else logging.INFO,
                        format='%(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    n_states = len(params.states)
    transition_rates = np.array(params.transition_rates).reshape((n_states, n_states))
    transmission_rates = np.array(params.transmission_rates).reshape((n_states, n_states))
    if not params.avg_recipients:
        params.avg_recipients = [1] * n_states

    is_mult = np.any(np.array(params.avg_recipients) != 1)
    logging.info(
        'MTBD{} model parameters are:\n'
        '\ttransition_rates=\n{}\n'
        '\ttransmission_rates=\n{}\n'
        '\tremoval_rates={}\n'
        '\tsampling probabilities={}{}'
        .format('-MULT' if is_mult else '',
                transition_rates, transmission_rates, params.removal_rates,
                params.sampling_probabilities,
                '\n\tavg_recipient_numbers={}'.format(params.avg_recipients) if is_mult else ''))
    model = Model(states=params.states,
                  transmission_rates=transmission_rates,
                  transition_rates=transition_rates,
                  removal_rates=params.removal_rates, ps=params.sampling_probabilities)
    if params.upsilon and params.upsilon > 0:
        logging.info('PN parameters are:\n\tphi={}\n\tupsilon={}'.format(params.phi, params.upsilon))
        model = CTModel(model=model, upsilon=params.upsilon, phi=params.phi)

    if params.T < np.inf:
        logging.info('Total time T={}'.format(params.T))

    forest, (total_tips, u, T), ltt = generate(model, params.min_tips, params.max_tips, T=params.T,
                                               max_notified_contacts=params.max_notified_contacts)

    save_forest(forest, params.nwk)
    save_log(model, total_tips, T, u, params.log)
    if params.ltt:
        save_ltt(ltt, observed_ltt(forest, T), params.ltt)


if '__main__' == __name__:
    main()
