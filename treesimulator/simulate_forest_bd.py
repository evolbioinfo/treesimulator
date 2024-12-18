import logging

import numpy as np

from treesimulator import save_forest, save_log, save_ltt
from treesimulator.generator import generate, observed_ltt
from treesimulator.mtbd_models import BirthDeathModel, CTModel


def main():
    """
    Entry point for tree/forest generation with the BD model with command-line arguments.
    :return: void
    """
    import argparse

    parser = \
        argparse.ArgumentParser(description="Simulates a tree (or a forest of trees) for given BD model parameters. "
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
    parser.add_argument('--la', required=True, type=float, help="transmission rate")
    parser.add_argument('--psi', required=True, type=float, help="removal rate")
    parser.add_argument('--p', required=True, type=float, help='sampling probability')
    parser.add_argument('--upsilon', required=False, default=0, type=float, help='notification probability')
    parser.add_argument('--max_notified_contacts', required=False, default=1, type=int,
                        help='maximum number of notified contacts per person')
    parser.add_argument('--avg_recipients', required=False, default=1, type=float,
                        help='average number of recipients per transmission. '
                             'By default one (one-to-one transmission), '
                             'but if a larger number is given then one-to-many transmissions become possible.')
    parser.add_argument('--phi', required=False, default=0, type=float, help='notified removal rate')
    parser.add_argument('--log', required=True, type=str, help="output log file")
    parser.add_argument('--nwk', required=True, type=str, help="output tree or forest file")
    parser.add_argument('--ltt', required=False, default=None, type=str, help="output LTT file")
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help="describe generation process")
    params = parser.parse_args()
    logging.getLogger().handlers = []
    logging.basicConfig(level=logging.DEBUG if params.verbose else logging.INFO,
                        format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    is_mult = params.avg_recipients != 1
    logging.info('BD{} model parameters are:\n\tlambda={}\n\tpsi={}\n\tp={}{}'
                 .format('-MULT' if is_mult else '',
                         params.la, params.psi, params.p,
                         '\n\tr={}'.format(params.avg_recipients) if is_mult else ''))
    model = BirthDeathModel(p=params.p, la=params.la, psi=params.psi, n_recipients=[params.avg_recipients])
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
