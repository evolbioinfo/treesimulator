import logging

import numpy as np

from treesimulator import save_forest, save_log, save_ltt
from treesimulator.generator import generate, observed_ltt
from treesimulator.mtbd_models import BirthDeathWithSuperSpreadingModel, PNModel


def main():
    """
    Entry point for tree/forest generation with the BDSS model with command-line arguments.
    :return: void
    """
    import argparse

    parser = \
        argparse.ArgumentParser(description="Simulates a tree (or a forest of trees) for given BDSS model parameters. "
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
    parser.add_argument('--la_ss', required=True, type=float, help="superspreader-to-superspreader transmission rate")
    parser.add_argument('--la_sn', required=True, type=float, help="superspreader-to-normalspreader transmission rate")
    parser.add_argument('--la_nn', required=True, type=float, help="normalspreader-to-normalspreader transmission rate")
    parser.add_argument('--la_ns', required=True, type=float, help="normalspreader-to-superspreader transmission rate")
    parser.add_argument('--psi', required=True, type=float, help="removal rate")
    parser.add_argument('--p', required=True, type=float, help='sampling probability')
    parser.add_argument('--upsilon', required=False, default=0, type=float, help='notification probability')
    parser.add_argument('--phi', required=False, default=0, type=float, help='partner removal rate')
    parser.add_argument('--max_notified_partners', required=False, default=1, type=int,
                        help='maximum number of notified partners per person')
    parser.add_argument('--log', required=True, type=str, help="output log file")
    parser.add_argument('--nwk', required=True, type=str, help="output tree or forest file")
    parser.add_argument('--ltt', required=False, default=None, type=str, help="output LTT file")
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help="describe generation process")
    params = parser.parse_args()
    logging.getLogger().handlers = []
    logging.basicConfig(level=logging.DEBUG if params.verbose else logging.INFO,
                        format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    logging.info('BDSS model parameters are:\n\tla_ss={}\n\tla_sn={}\n\tla_ns={}\n\tla_nn={}\n\tpsi={}\n\tp={}'
                 .format(params.la_ss, params.la_sn, params.la_ns, params.la_nn, params.psi, params.p))
    logging.info('Total time T={}'.format(params.T))

    model = BirthDeathWithSuperSpreadingModel(p=params.p,
                                              la_ss=params.la_ss, la_sn=params.la_sn,
                                              la_ns=params.la_ns, la_nn=params.la_nn, psi=params.psi)
    if params.upsilon and params.upsilon > 0:
        logging.info('PN model parameters are:\n\tphi={}\n\tupsilon={}'.format(params.phi, params.upsilon))
        model = PNModel(model=model, upsilon=params.upsilon, partner_removal_rate=params.phi)

    forest, (total_tips, u, T), ltt = generate(model, params.min_tips, params.max_tips, T=params.T,
                                               max_notified_partners=params.max_notified_partners)

    save_forest(forest, params.nwk)
    save_log(model, total_tips, T, u, params.log)
    if params.ltt:
        save_ltt(ltt, observed_ltt(forest, T), params.ltt)


if '__main__' == __name__:
    main()
