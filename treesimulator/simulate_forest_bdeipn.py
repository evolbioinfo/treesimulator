import logging

import numpy as np

from treesimulator import save_forest, save_log, save_ltt
from treesimulator.generator import generate, observed_ltt
from treesimulator.mtbd_models import BirthDeathExposedInfectiousModel, PNModel


def main():
    """
    Entry point for tree/forest generation with the BDEI model with command-line arguments.
    :return: void
    """
    import argparse

    parser = \
        argparse.ArgumentParser(description="Simulates a tree (or a forest of trees) for given BDEI-PN model parameters. "
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
    parser.add_argument('--mu', required=True, type=float, help="becoming infectious rate")
    parser.add_argument('--la', required=True, type=float, help="transmission rate")
    parser.add_argument('--psi', required=True, type=float, help="removal rate")
    parser.add_argument('--p', required=True, type=float, help='sampling probability')
    parser.add_argument('--pn', required=True, type=float, help='notification probability')
    parser.add_argument('--partner_psi', required=True, type=float, help='partner removal rate')
    parser.add_argument('--log', required=True, type=str, help="output log file")
    parser.add_argument('--nwk', required=True, type=str, help="output tree or forest file")
    parser.add_argument('--ltt', required=False, default=None, type=str, help="output LTT file")
    params = parser.parse_args()
    logging.getLogger().handlers = []
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    logging.info('BDEI model parameters are:\n\tmu={}\n\tlambda={}\n\tpsi={}\n\tp={}'
                 .format(params.mu, params.la, params.psi, params.p))
    logging.info('Total time T={}'.format(params.T))

    model = PNModel(model=BirthDeathExposedInfectiousModel(p=params.p, mu=params.mu, la=params.la, psi=params.psi),
                    pn=params.pn, removal_rate=params.partner_psi)

    forest, (total_tips, u, T), ltt = generate(model, params.min_tips, params.max_tips, T=params.T)

    save_forest(forest, params.nwk)
    save_log(model, total_tips, T, u, params.log)
    if params.ltt:
        save_ltt(ltt, observed_ltt(forest, T), params.ltt)


if '__main__' == __name__:
    main()
