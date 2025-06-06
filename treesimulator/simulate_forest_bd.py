import logging
import numpy as np
from treesimulator import save_forest, save_log, save_ltt
from treesimulator.generator import generate, observed_ltt
from treesimulator.mtbd_models import BirthDeathModel, CTModel


def main():
    """
    Entry point for tree/forest generation with the BD-CT-Skyline model with command-line arguments.

    For skyline models, the first model (models[0]) always starts at time 0, and the time points list
    specifies when to switch to each subsequent model.
    :return: void
    """
    import argparse

    parser = \
        argparse.ArgumentParser(description="Simulates a tree (or a forest of trees) for given BD-CT-Skyline model parameters. "
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


    parser.add_argument('--la', required=True, nargs='+', type=float,
                        help="List of transmission rates (one per skyline interval).")
    parser.add_argument('--psi', required=True, nargs='+', type=float,
                        help="List of removal rates (one per skyline interval).")
    parser.add_argument('--p', required=True, nargs='+', type=float,
                        help="List of sampling probabilities (one per skyline interval).")
    parser.add_argument('--skyline_times', nargs='*', type=float,
                        help="List of time points specifying when to switch from model i to model i+1 in the Skyline."
                             "Must be sorted in ascending order and contain one less elements "
                             "than the number of models in the Skyline."
                             "The first model always starts at time 0.")

    # Contact tracing parameters
    parser.add_argument('--upsilon', nargs='*', type=float,
                        help='List of notification probabilities (one per skyline interval). Omit for no contact tracing.')
    parser.add_argument('--phi', nargs='*', type=float,
                        help='List of notified removal rates (one per skyline interval). Only used if upsilon is specified.')
    parser.add_argument('--max_notified_contacts', required=False, default=1, type=int,
                        help='Maximum notified contacts')
    parser.add_argument('--allow_irremovable_states', action='store_true', default=False,
                        help='If specified and the initial model had a zero removal rate (psi=0), '
                             'then even after notification the removal rate will stay zero (phi=0) '
                             'and the given value of phi ignored.')

    parser.add_argument('--avg_recipients', required=False, default=1, type=float,
                        help='average number of recipients per transmission. '
                             'By default one (one-to-one transmission), '
                             'but if a larger number is given then one-to-many transmissions become possible.')

    parser.add_argument('--log', required=False, default=None, type=str, help="output file to log BD model parameters")
    parser.add_argument('--nwk', required=True, type=str, help="output tree or forest file")
    parser.add_argument('--ltt', required=False, default=None, type=str, help="output LTT file")
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help="describe generation process")
    parser.add_argument('-s', '--seed', type=int, default=None, help='random seed for reproducibility')

    params = parser.parse_args()
    logging.getLogger().handlers = []
    logging.basicConfig(level=logging.DEBUG if params.verbose else logging.INFO,
                        format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    # Check that all parameter arrays have the same length
    n_models = len(params.la)
    if n_models != len(params.psi) or n_models != len(params.p):
        raise ValueError("All parameter lists (la, psi, p) must have the same length")
    is_ct = params.upsilon
    if is_ct:
        if n_models != len(params.upsilon) or n_models != len(params.phi):
            raise ValueError("Contact-tracing parameter lists must have the same length "
                             "as the other parameter lists (la, psi, p)")

    if n_models > 1 and (not params.skyline_times or len(params.skyline_times) != n_models - 1):
        raise ValueError(f'One should specify {n_models - 1} skyline times for {n_models}, '
                         f'got {len(params.skyline_times)} instead.')

    if params.T < np.inf:
        logging.info('Total time T={}'.format(params.T))

    is_mult = params.avg_recipients != 1

    models = []
    for i in range(n_models):
        logging.info('{}BD{}{} model parameters are:\n\tlambda={}\n\tpsi={}\n\tp={}{}{}'
                     .format('For time interval {}-{},\n'.format(0 if i == 0 else params.skyline_times[i - 1],
                                                                 params.skyline_times[i] if i < (n_models - 1) else '...')
                                                          if n_models > 1 else '',
                             '-MULT' if is_mult else '',
                             '-CT' if is_ct else '',
                             params.la[i], params.psi[i], params.p[i],
                             '\n\tavg_recipient_number={}'.format(params.avg_recipients) if is_mult else '',
                             '\n\tphi={}\n\tupsilon={}'.format(params.phi[i], params.upsilon[i]) if is_ct else ''))
        model = BirthDeathModel(p=params.p[i], la=params.la[i], psi=params.psi[i], n_recipients=[params.avg_recipients])
        if is_ct:
            model = CTModel(model=model, upsilon=params.upsilon[i], phi=params.phi[i],
                            allow_irremovable_states=params.allow_irremovable_states)
        models.append(model)

    forest, (total_tips, u, T, observed_frequencies), ltt = \
        generate(models, skyline_times=params.skyline_times, T=params.T,
                 min_tips=params.min_tips, max_tips=params.max_tips, max_notified_contacts=params.max_notified_contacts,
                 random_seed=params.seed)

    save_forest(forest, params.nwk)

    if params.log:
        save_log(models, params.skyline_times, total_tips, T, u, params.log,
                 kappa=params.max_notified_contacts, observed_frequencies=observed_frequencies)
    if params.ltt:
        save_ltt(ltt, observed_ltt(forest, T), params.ltt)

if '__main__' == __name__:
    main()