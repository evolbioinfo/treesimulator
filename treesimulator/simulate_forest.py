import logging

import numpy as np

from treesimulator import save_forest, save_log, save_ltt
from treesimulator.generator import generate, observed_ltt
from treesimulator.mtbd_models import Model, CTModel


def main():
    """
    Entry point for tree/forest generation with a generic MTBD model, supporting both
    standard and skyline (time-varying) models.

    For skyline models, the first model (models[0]) always starts at time 0, and the time points list
    specifies when to switch to each subsequent model.
    :return: void
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Simulates a tree (or a forest of trees) for given MTBD model parameters. "
                    "Supports both standard models and skyline (time-varying) models. "
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

    # Modified parameter approach for skyline models
    parser.add_argument('--t', nargs='+', type=float,
                        help="time points for skyline models, specifying when to switch from model i to model i+1. "
                             "The number of models will be determined as len(t) + 1 since model[0] always starts at time 0.")

    # New approach: Specify matrices for each time point using a flattened format
    # Format: transition_matrices[timepoint_idx][row][col]
    parser.add_argument('--transition_matrices', required=True, nargs='+', type=float,
                        help="Flattened transition rate matrices for each time point. "
                             "For n states and m time points, provide m*n*n values.")

    parser.add_argument('--transmission_matrices', required=True, nargs='+', type=float,
                        help="Flattened transmission rate matrices for each time point. "
                             "For n states and m time points, provide m*n*n values.")

    parser.add_argument('--removal_vectors', required=True, nargs='+', type=float,
                        help="Flattened removal rate vectors for each time point. "
                             "For n states and m time points, provide m*n values.")

    parser.add_argument('--sampling_vectors', required=True, nargs='+', type=float,
                        help="Flattened sampling probability vectors for each time point. "
                             "For n states and m time points, provide m*n values.")

    # Contact tracing parameters
    parser.add_argument('--upsilon', nargs='+', type=float, default=[0.0],
                        help='List of notification probabilities for each interval. Use 0 for no contact tracing.')
    parser.add_argument('--phi', nargs='+', type=float, default=[0.0],
                        help='List of notified removal rates for each interval. Used with contact tracing.')
    parser.add_argument('--max_notified_contacts', required=False, default=1, type=int,
                        help='maximum number of notified contracts per person')

    # Averaging parameters
    parser.add_argument('--avg_recipients', nargs='*', type=float,
                        help='average number of recipients per transmission '
                             'for each donor state (in the same order as the model states). '
                             'By default, only one-to-one transmissions are allowed, '
                             'but if larger numbers are given then one-to-many transmissions become possible.')

    # Output parameters
    parser.add_argument('--log', required=True, type=str, help="output log file")
    parser.add_argument('--nwk', required=True, type=str, help="output tree or forest file")
    parser.add_argument('--ltt', required=False, default=None, type=str, help="output LTT file")
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="print information on the progress of the tree generation (to console)")

    params = parser.parse_args()
    logging.getLogger().handlers = []
    logging.basicConfig(level=logging.DEBUG if params.verbose else logging.INFO,
                        format='%(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    # Validation and processing
    n_states = len(params.states)

    # Determine number of timepoints from time points provided
    n_timepoints = 1 if params.t is None else len(params.t) + 1

    # Validate time points
    if params.t is not None:
        # Verify that times are sorted
        for i in range(len(params.t) - 1):
            if params.t[i] >= params.t[i + 1]:
                raise ValueError(
                    f"Time points must be in ascending order. Found t[{i}] = {params.t[i]} >= t[{i + 1}] = {params.t[i + 1]}")

        logging.info(
            f"Creating a skyline model with {n_timepoints} intervals at time points: [0.0, {', '.join(map(str, params.t))}]")

    # Reshape the flattened matrices/vectors into proper structures
    transition_matrices = []
    transmission_matrices = []
    removal_vectors = []
    sampling_vectors = []

    # Expected lengths
    expected_transition_len = n_timepoints * n_states * n_states
    expected_removal_len = n_timepoints * n_states

    # Validate lengths
    if len(params.transition_matrices) != expected_transition_len:
        raise ValueError(
            f"Expected {expected_transition_len} values for transition matrices, got {len(params.transition_matrices)}")
    if len(params.transmission_matrices) != expected_transition_len:
        raise ValueError(
            f"Expected {expected_transition_len} values for transmission matrices, got {len(params.transmission_matrices)}")
    if len(params.removal_vectors) != expected_removal_len:
        raise ValueError(
            f"Expected {expected_removal_len} values for removal vectors, got {len(params.removal_vectors)}")
    if len(params.sampling_vectors) != expected_removal_len:
        raise ValueError(
            f"Expected {expected_removal_len} values for sampling vectors, got {len(params.sampling_vectors)}")

    # Reshape the inputs
    for t in range(n_timepoints):
        # Extract and reshape transition matrix for this time point
        start_idx = t * n_states * n_states
        end_idx = start_idx + n_states * n_states
        trans_matrix = np.array(params.transition_matrices[start_idx:end_idx]).reshape((n_states, n_states))
        transition_matrices.append(trans_matrix)

        # Extract and reshape transmission matrix for this time point
        trans_matrix = np.array(params.transmission_matrices[start_idx:end_idx]).reshape((n_states, n_states))
        transmission_matrices.append(trans_matrix)

        # Extract removal vector for this time point
        start_idx = t * n_states
        end_idx = start_idx + n_states
        removal_vectors.append(np.array(params.removal_vectors[start_idx:end_idx]))

        # Extract sampling vector for this time point
        sampling_vectors.append(np.array(params.sampling_vectors[start_idx:end_idx]))

    # Handle default values for avg_recipients
    if not params.avg_recipients:
        params.avg_recipients = [1] * n_states

    is_mult = np.any(np.array(params.avg_recipients) != 1)

    # Create models from these structured parameters
    models = []
    for t in range(n_timepoints):
        # Create model for this time point
        model = Model(
            states=params.states,
            transition_rates=transition_matrices[t],
            transmission_rates=transmission_matrices[t],
            removal_rates=removal_vectors[t],
            ps=sampling_vectors[t],
            n_recipients=params.avg_recipients
        )

        # Apply contact tracing if specified
        upsilon_value = params.upsilon[t] if t < len(params.upsilon) else params.upsilon[-1]
        phi_value = params.phi[t] if t < len(params.phi) else params.phi[-1]

        if upsilon_value > 0:
            model = CTModel(model=model, upsilon=upsilon_value, phi=phi_value)
            logging.info(f'Added contact tracing with upsilon={upsilon_value}, phi={phi_value} for time point {t}')

        models.append(model)

        # Display time point info
        if t == 0:
            time_info = "starting at time 0"
        else:
            time_info = f"for times >= {params.t[t - 1]}"

        if t < n_timepoints - 1:
            time_info += f" and < {params.t[t]}"

        logging.info(f"Model {t + 1} {time_info}:")
        logging.info(f'\ttransition_rates=\n{transition_matrices[t]}')
        logging.info(f'\ttransmission_rates=\n{transmission_matrices[t]}')
        logging.info(f'\tremoval_rates={removal_vectors[t]}')
        logging.info(f'\tsampling probabilities={sampling_vectors[t]}')
        if is_mult:
            logging.info(f'\tavg_recipient_numbers={params.avg_recipients}')

    if params.T < np.inf:
        logging.info('Total time T={}'.format(params.T))

    # Use the generic generator to handle both single and skyline models
    forest, (total_tips, u, T), ltt = generate(
        models,
        params.min_tips,
        params.max_tips,
        T=params.T,
        skyline_times=params.t if n_timepoints > 1 else None,
        max_notified_contacts=params.max_notified_contacts
    )

    save_forest(forest, params.nwk)
    save_log(models[0], total_tips, T, u, params.log)
    if params.ltt:
        save_ltt(ltt, observed_ltt(forest, T), params.ltt)


if '__main__' == __name__:
    main()