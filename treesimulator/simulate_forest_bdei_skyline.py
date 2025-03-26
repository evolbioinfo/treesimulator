import logging
import numpy as np
from treesimulator import save_forest, save_log, save_ltt
from treesimulator.generator_skyline import generate, observed_ltt
from treesimulator.mtbd_models import BirthDeathExposedInfectiousModel, CTModel


def main():
    """
    Entry point for tree/forest generation using the BDEI-Skyline model with command-line arguments.
    Now implemented using a list-based approach rather than a dedicated skyline class.
    Supports time-varying contact tracing parameters.

    For skyline models, the first model (models[0]) always starts at time 0, and the time points list
    specifies when to switch to each subsequent model.
    :return: void
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Simulates a tree (or a forest of trees) with the BDEI-Skyline model using a list-based approach.")

    parser.add_argument('--min_tips', required=True, type=int, help="Minimum number of simulated leaves.")
    parser.add_argument('--max_tips', required=True, type=int, help="Maximum number of simulated leaves.")
    parser.add_argument('--T', required=False, default=np.inf, type=float, help="Total simulation time.")
    parser.add_argument('--mu', required=True, nargs='+', type=float,
                        help="List of E->I transition rates for each interval.")
    parser.add_argument('--la', required=True, nargs='+', type=float,
                        help="List of transmission rates for each interval.")
    parser.add_argument('--psi', required=True, nargs='+', type=float,
                        help="List of removal rates for each interval.")
    parser.add_argument('--p', required=True, nargs='+', type=float,
                        help="List of sampling probabilities for each interval.")
    parser.add_argument('--t', nargs='+', type=float,
                        help="List of time points specifying when to switch from model i to model i+1. "
                             "Length should be one less than the number of parameter sets, since model[0] "
                             "always starts at time 0. Must be in ascending order.")

    # Contact tracing parameters - now as lists for time-varying support
    parser.add_argument('--upsilon', nargs='+', type=float, default=[0.0],
                        help='List of notification probabilities for each interval. Use 0 for no contact tracing.')
    parser.add_argument('--phi', nargs='+', type=float, default=[0.0],
                        help='List of notified removal rates for each interval. Used with contact tracing.')
    parser.add_argument('--max_notified_contacts', required=False, default=1, type=int,
                        help='Maximum notified contacts')

    parser.add_argument('--avg_recipients', required=False, default=1, type=float,
                        help='Average number of recipients per transmission.')
    parser.add_argument('--log', required=True, type=str, help="Output log file")
    parser.add_argument('--nwk', required=True, type=str, help="Output tree or forest file")
    parser.add_argument('--ltt', required=False, default=None, type=str, help="Output LTT file")
    parser.add_argument('-v', '--verbose', action='store_true', help="Verbose output")

    params = parser.parse_args()

    # Set up logging AFTER parsing arguments (moved from before parser definition)
    logging.getLogger().handlers = []
    logging.basicConfig(level=logging.DEBUG if params.verbose else logging.INFO,
                        format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    try:
        # If upsilon/phi are single values, replicate them for all time points
        if len(params.upsilon) == 1 and len(params.mu) > 1:
            params.upsilon = params.upsilon * len(params.mu)
        if len(params.phi) == 1 and len(params.mu) > 1:
            params.phi = params.phi * len(params.mu)

        # Validate all parameters have the same length
        if len(params.mu) != len(params.la) or len(params.mu) != len(params.psi) or \
                len(params.mu) != len(params.p) or len(params.upsilon) != len(params.mu) or \
                len(params.phi) != len(params.mu):
            raise ValueError("All parameter arrays (mu, la, psi, p, upsilon, phi) must have the same length")

        # For skyline model, time points should be one less than parameter count
        if len(params.mu) > 1:  # If we have multiple parameter sets, it's a skyline model
            # Check if time points are provided
            if params.t is None:
                raise ValueError("Time points (--t) must be provided when using multiple parameter sets")

            # Check number of time points
            if len(params.t) != len(params.mu) - 1:
                raise ValueError(
                    f"For skyline models, the number of time points must be one less than the number of parameter sets "
                    f"since model[0] always starts at time 0. Got {len(params.mu)} parameter sets and {len(params.t)} time points.")

            # Verify that times are sorted
            for i in range(len(params.t) - 1):
                if params.t[i] >= params.t[i + 1]:
                    raise ValueError(
                        f"Time points must be in ascending order. Found t[{i}] = {params.t[i]} >= t[{i + 1}] = {params.t[i + 1]}")

        # Log the configuration
        is_mult = params.avg_recipients != 1
        mult_str = '-MULT' if is_mult else ''
        logging.info(f'BDEI{mult_str}-Skyline parameters are:')
        logging.info(f'mu values: {params.mu}')
        logging.info(f'lambda values: {params.la}')
        logging.info(f'psi values: {params.psi}')
        logging.info(f'p values: {params.p}')

        if len(params.mu) > 1:
            # Include time 0 in the display of time points
            display_times = ['0.0'] + [str(t) for t in params.t]
            logging.info(f'Time points: [{", ".join(display_times)}]')

        logging.info(f'upsilon values: {params.upsilon}')
        logging.info(f'phi values: {params.phi}')
        if is_mult:
            logging.info(f'Average recipients: {params.avg_recipients}')

        # Create a list of BDEI models
        models = []
        for i in range(len(params.mu)):
            # Display time interval info
            if i == 0:
                time_info = "starting at time 0"
            else:
                time_info = f"for times >= {params.t[i - 1]}"

            if i < len(params.mu) - 1:
                time_info += f" and < {params.t[i]}"

            model_name = f'BDEI{i + 1}'
            logging.info(
                f'Creating model {model_name} {time_info} with mu={params.mu[i]}, la={params.la[i]}, psi={params.psi[i]}, p={params.p[i]}')

            # Create a BDEI model with the parameters for this time interval
            model = BirthDeathExposedInfectiousModel(
                mu=params.mu[i],
                la=params.la[i],
                psi=params.psi[i],
                p=params.p[i],
                n_recipients=[params.avg_recipients, params.avg_recipients]
            )

            # Apply contact tracing if specified for this time interval
            if params.upsilon[i] > 0:
                model = CTModel(model=model, upsilon=params.upsilon[i], phi=params.phi[i])
                logging.info(
                    f'Added contact tracing to model {model_name} with upsilon={params.upsilon[i]}, phi={params.phi[i]}')

            models.append(model)

        if params.T < np.inf:
            logging.info(f'Total time T={params.T}')

        # Generate forest using the skyline model approach (list of models)
        forest, (total_tips, u, max_time), ltt = generate(
            models,
            min_tips=params.min_tips,
            max_tips=params.max_tips,
            T=params.T,
            skyline_times=params.t if len(params.mu) > 1 else None,  # Only pass time points for skyline models
            max_notified_contacts=params.max_notified_contacts
        )

        # Save outputs
        save_forest(forest, params.nwk)
        # For logging, use the first model (without the skyline parameter which isn't supported)
        save_log(models[0], total_tips, max_time, u, params.log)
        if params.ltt:
            save_ltt(ltt, observed_ltt(forest, max_time), params.ltt)

        logging.info("Simulation completed successfully")

    except RuntimeError as e:
        logging.error(f"Runtime error during simulation: {e}")
    except ValueError as e:
        logging.error(f"Value error during simulation: {e}")
    except Exception as e:
        logging.error(f"Unexpected error: {e}")


if '__main__' == __name__:
    main()