import logging

import numpy as np

from treesimulator import save_forest, save_log, save_ltt
from treesimulator.generator import generate, observed_ltt
from treesimulator.mtbd_models import Model, CTModel


def main():
    """
    Entry point for tree/forest generation with a generic MTBD-CT-Skyline model with command-line arguments.


    For skyline models, the first model (models[0]) always starts at time 0, and the time points list
    specifies when to switch to each subsequent model.
    :return: void
    """
    import argparse

    parser = \
        argparse.ArgumentParser(description="Simulates a tree (or a forest of trees) "
                                            "for given MTBD-CT-Skyline model parameters. "
                                            "If a simulation leads to less than --min_tips tips, it is repeated.")

    # Common parameters
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

    # States parameter
    parser.add_argument('--states', required=True, nargs='+', type=str, help="model states")

    parser.add_argument('--skyline_times', nargs='+', type=float,
                        help="List of time points specifying when to switch from model i to model i+1 in the Skyline."
                             "Must be sorted in ascending order and contain one less elements "
                             "than the number of models in the Skyline."
                             "The first model always starts at time 0.")

    parser.add_argument('--transition_rates', nargs='+', type=float,
                        help="transition rate matrix, row after row, in the same order as model states, "
                             "e.g. if a model has 2 states given as --states E I,"
                             "then here we expect E->E E->I I->E I->I."
                             "All the transition rate from the state to itself  (e.g., E->E and I->I) must be zero."
                             "If the Skyline is used, transition rate matrix for the model at the tree root "
                             "is followed by the transition rate matrix of the next model, etc. "
                             "For instance, for two models (m1 and m2) with the states E and I, we expect:"
                             "E->E_m1 E->I_m1 I->E_m1 I->I_m1 E->E_m2 E->I_m2 I->E_m2 I->I_m2.")
    parser.add_argument('--transmission_rates', nargs='+', type=float,
                        help="transmission rate matrix, row after row, in the same order as model states, "
                             "e.g. if a model has 2 states given as --states E I,"
                             "then here we expect E->E E->I I->E I->I."
                             "If the Skyline is used, transmission rate matrix for the model at the tree root "
                             "is followed by the transmission rate matrix of the next model, etc. "
                             "For instance, for two models (m1 and m2) with the states E and I, we expect:"
                             "E->E_m1 E->I_m1 I->E_m1 I->I_m1 E->E_m2 E->I_m2 I->E_m2 I->I_m2.")
    parser.add_argument('--removal_rates', nargs='+', type=float,
                        help="removal rate array, in the same order as model states, "
                             "e.g. if a model has 2 states given as --states E I,"
                             "then here we expect removal(E) removal(I)."
                             "If the Skyline is used, removal rates for the model at the tree root "
                             "are followed by the removal rates of the next model, etc. "
                             "For instance, for two models (m1 and m2) with the states E and I, we expect:"
                             "removal(E)_m1 removal(I)_m1 removal(E)_m2 removal(I)_m2.")
    parser.add_argument('--sampling_probabilities', nargs='+', type=float,
                        help="sampling probability array, in the same order as model states, "
                             "e.g. if a model has 2 states given as --states E I,"
                             "then here we expect p(E) p(I)."
                             "If the Skyline is used, sampling probabilities for the model at the tree root "
                             "are followed by the sampling probabilities of the next model, etc. "
                             "For instance, for two models (m1 and m2) with the states E and I, we expect:"
                             "p(E)_m1 p(I)_m1 p(E)_m2 p(I)_m2.")


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
                        help='Maximum notified contacts.')
    parser.add_argument('--allow_irremovable_states', action='store_true', default=False,
                        help='If specified and the initial model included "irremovable" states '
                             '(i.e., whose removal rate was zero, e.g., E in the BDEI model), '
                             'then even after notification their removal rate will stay zero, '
                             'and the corresponding individuals will become "removable" (at a rate phi) '
                             'only once they change the state to a "removable" one '
                             '(e.g., from E-C to I-C in BDEI-CT).')

    parser.add_argument('--avg_recipients', nargs='*', type=float,
                        help='average number of recipients per transmission '
                             'for each donor state (in the same order as the model states). '
                             'By default, only one-to-one transmissions are allowed, '
                             'but if larger numbers are given then one-to-many transmissions become possible.')

    parser.add_argument('--root_state', type=str, default=None,
                        help='state at the beginning of the root branch(es). '
                             'By default, one of the model states will be chosen randomly (and independently for each tree) '
                             'based on the state equilibrium frequencies.')


    parser.add_argument('--log', required=False, default=None, type=str, help="output file to log model parameters")
    parser.add_argument('--nwk', required=True, type=str, help="output tree or forest file")
    parser.add_argument('--ltt', required=False, default=None, type=str, help="output LTT file")
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help="describe generation process")
    parser.add_argument('-s', '--seed', type=int, default=None, help='random seed for reproducibility')

    params = parser.parse_args()
    logging.getLogger().handlers = []
    logging.basicConfig(level=logging.DEBUG if params.verbose else logging.INFO,
                        format='%(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    n_states = len(params.states)
    if params.skyline_times and len(params.skyline_times) > 0:
        n_models = len(params.skyline_times) + 1
    else:
        n_models = 1

    try:
        transition_rates = np.array(params.transition_rates).reshape((n_states, n_states, n_models))
        transmission_rates = np.array(params.transmission_rates).reshape((n_states, n_states, n_models))
        removal_rates = np.array(params.removal_rates).reshape((n_states, n_models))
        sampling_probabilities = np.array(params.sampling_probabilities).reshape((n_states, n_models))
    except:
        raise ValueError(f'Got {n_models - 1} skyline times, and {n_states} states, '
                         f'hence expecting {n_models} model(s) with {n_states} x {n_states} x {n_models} transition rates, '
                         f'{n_states} x {n_states} x {n_models} transmission rates, '
                         f'{n_states} x {n_models} removal rates, '
                         f'and {n_states} x {n_models} sampling probabilities. '
                         f'However (at least some of) your inputs do not correspond to these dimensions, '
                         f'please check them.')


    is_ct = params.upsilon
    if is_ct:
        if n_models != len(params.upsilon) or n_models != len(params.phi):
            raise ValueError("Contact-tracing parameter lists must have the same length "
                             f"as the number of models in the skyline ({n_models})")

    if not params.avg_recipients:
        params.avg_recipients = [1] * n_states
    is_mult = np.any(np.array(params.avg_recipients) != 1)

    if params.T < np.inf:
        logging.info('Total time T={}'.format(params.T))

    models = []
    for i in range(n_models):
        logging.info('{}MTBD{}{} model parameters are:\n'
                     '\tstates={}\n'
                     '\ttransition_rates=\n{}\n'
                     '\ttransmission_rates=\n{}\n'
                     '\tremoval_rates={}\n'
                     '\tsampling probabilities={}'
                     '{}{}'
                     .format('For time interval {}-{},\n'.format(0 if i == 0 else params.skyline_times[i - 1],
                                                                 params.skyline_times[i] if i < (
                                                                         n_models - 1) else '...')
                             if n_models > 1 else '',
                             '-MULT' if is_mult else '',
                             '-CT' if is_ct else '',
                             params.states,
                             transition_rates[i, :, :],
                             transmission_rates[i, :, :],
                             removal_rates[i, :],
                             sampling_probabilities[i, :],
                             '\n\tavg_recipient_numbers={}'.format(params.avg_recipients) if is_mult else '',
                             '\n\tphi={}\n\tupsilon={}'.format(params.phi[i], params.upsilon[i]) if is_ct else ''))
        model = Model(states=params.states,
                      transmission_rates=transmission_rates[i, :, :],
                      transition_rates=transition_rates[i, :, :],
                      removal_rates=removal_rates[i, :], ps=sampling_probabilities[i, :],
                      n_recipients=params.avg_recipients)
        if is_ct:
            model = CTModel(model=model, upsilon=params.upsilon[i], phi=params.phi[i],
                            allow_irremovable_states=params.allow_irremovable_states)
        models.append(model)

    forest, (total_tips, u, T, observed_frequencies), ltt = \
        generate(models, skyline_times=params.skyline_times, T=params.T,
                 min_tips=params.min_tips, max_tips=params.max_tips, max_notified_contacts=params.max_notified_contacts,
                 root_state=params.root_state, random_seed=params.seed)

    save_forest(forest, params.nwk)

    if params.log:
        save_log(models, params.skyline_times, total_tips, T, u, params.log, kappa=params.max_notified_contacts,
                 observed_frequencies=observed_frequencies)
    if params.ltt:
        save_ltt(ltt, observed_ltt(forest, T), params.ltt)


if '__main__' == __name__:
    main()