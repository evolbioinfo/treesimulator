import logging
import numpy as np
from treesimulator import save_forest, save_log, save_ltt
from treesimulator.generator import generate, observed_ltt
from treesimulator.mtbd_models import BirthDeathWithSuperSpreadingModel, CTModel


def main():
    """
    Entry point for tree/forest generation using the BDSS-Skyline model with command-line arguments.
    Now implemented using a list-based approach rather than a dedicated skyline class.
    Supports time-varying contact tracing parameters.

    For skyline models, the first model (models[0]) always starts at time 0, and the time points list
    specifies when to switch to each subsequent model.
    :return: void
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Simulates a tree (or a forest of trees) with the BDSS-Skyline model using a list-based approach.")

    parser.add_argument('--min_tips', required=True, type=int, help="Minimum number of simulated leaves.")
    parser.add_argument('--max_tips', required=True, type=int, help="Maximum number of simulated leaves.")
    parser.add_argument('--T', required=False, default=np.inf, type=float, help="Total simulation time.")
    parser.add_argument('--la_nn', required=True, nargs='+', type=float,
                        help="List of transmission rates from normal to normal for each interval.")
    parser.add_argument('--la_ns', required=True, nargs='+', type=float,
                        help="List of transmission rates from normal to super for each interval.")
    parser.add_argument('--la_sn', required=True, nargs='+', type=float,
                        help="List of transmission rates from super to normal for each interval.")
    parser.add_argument('--la_ss', required=True, nargs='+', type=float,
                        help="List of transmission rates from super to super for each interval.")
    parser.add_argument('--psi', required=True, nargs='+', type=float,
                        help="List of removal rates for each interval.")
    parser.add_argument('--p', required=True, nargs='+', type=float,
                        help="List of sampling probabilities for normal spreaders for each interval.")
    parser.add_argument('--p_s', required=True, nargs='+', type=float,
                        help="List of sampling probabilities for superspreaders for each interval.")
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

    parser.add_argument('--avg_recipients', nargs=2, default=[1, 1], type=float,
                        help='Average number of recipients per transmission for each donor state (normal spreader, superspreader).')
    parser.add_argument('--log', required=True, type=str, help="Output log file")
    parser.add_argument('--nwk', required=True, type=str, help="Output tree or forest file")
    parser.add_argument('--ltt', required=False, default=None, type=str, help="Output LTT file")
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help="Verbose output")

    params = parser.parse_args()
    logging.getLogger().handlers = []
    logging.basicConfig(level=logging.DEBUG if params.verbose else logging.INFO,
                        format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    try:
        # If upsilon/phi are single values, replicate them for all time points
        if len(params.upsilon) == 1 and len(params.la_nn) > 1:
            params.upsilon = params.upsilon * len(params.la_nn)
        if len(params.phi) == 1 and len(params.la_nn) > 1:
            params.phi = params.phi * len(params.la_nn)

        # Check if all parameter arrays have the same length
        param_lengths = [
            len(params.la_nn), len(params.la_ns), len(params.la_sn), len(params.la_ss),
            len(params.psi), len(params.p), len(params.p_s),
            len(params.upsilon), len(params.phi)
        ]

        if len(set(param_lengths)) > 1:
            raise ValueError("All parameter arrays must have the same length!")

        # For skyline model, time points should be one less than parameter count
        if len(params.la_nn) > 1:  # If we have multiple parameter sets, it's a skyline model
            # Check if time points are provided
            if params.t is None:
                raise ValueError("Time points (--t) must be provided when using multiple parameter sets")

            # Check number of time points
            if len(params.t) != len(params.la_nn) - 1:
                raise ValueError(
                    f"For skyline models, the number of time points must be one less than the number of parameter sets "
                    f"since model[0] always starts at time 0. Got {len(params.la_nn)} parameter sets and {len(params.t)} time points.")

            # Verify that times are sorted
            for i in range(len(params.t) - 1):
                if params.t[i] >= params.t[i + 1]:
                    raise ValueError(
                        f"Time points must be in ascending order. Found t[{i}] = {params.t[i]} >= t[{i + 1}] = {params.t[i + 1]}")

        # Log the configuration
        is_mult = np.any(np.array(params.avg_recipients) != 1)
        mult_str = '-MULT' if is_mult else ''
        logging.info(f'BDSS{mult_str}-Skyline parameters are:')
        logging.info(f'la_nn values: {params.la_nn}')
        logging.info(f'la_ns values: {params.la_ns}')
        logging.info(f'la_sn values: {params.la_sn}')
        logging.info(f'la_ss values: {params.la_ss}')
        logging.info(f'psi values: {params.psi}')
        logging.info(f'p values: {params.p}')
        logging.info(f'p_s values: {params.p_s}')

        if len(params.la_nn) > 1:
            # Include time 0 in the display of time points
            display_times = ['0.0'] + [str(t) for t in params.t]
            logging.info(f'Time points: [{", ".join(display_times)}]')

        logging.info(f'upsilon values: {params.upsilon}')
        logging.info(f'phi values: {params.phi}')
        if is_mult:
            logging.info(f'Average recipients: normal={params.avg_recipients[0]}, super={params.avg_recipients[1]}')

        # Create a list of BDSS models
        models = []
        for i in range(len(params.la_nn)):
            # Display time interval info
            if i == 0:
                time_info = "starting at time 0"
            else:
                time_info = f"for times >= {params.t[i - 1]}"

            if i < len(params.la_nn) - 1:
                time_info += f" and < {params.t[i]}"

            model_name = f'BDSS{i + 1}'

            # Check transmission ratio constraint for this interval
            if params.la_ns[i] == 0 or params.la_nn[i] == 0:
                # If either denominator is zero, check if all are consistent with zero
                if (params.la_ns[i] == 0 and params.la_ss[i] == 0) or (params.la_nn[i] == 0 and params.la_sn[i] == 0):
                    pass  # This is fine, both ratios are effectively 0/0 which we'll treat as equal
                else:
                    raise ValueError(
                        f'Transmission ratio constraint cannot be verified for interval {i + 1}: Cannot divide by zero '
                        f'(la_ns={params.la_ns[i]}, la_nn={params.la_nn[i]})'
                    )
            else:
                s_ratio = params.la_ss[i] / params.la_ns[i]
                n_ratio = params.la_sn[i] / params.la_nn[i]
                if np.abs(s_ratio - n_ratio) > 1e-3:
                    raise ValueError(
                        f'Transmission ratio constraint is violated for interval {i + 1}: '
                        f'la_ss / la_ns ({s_ratio}) must be equal to la_sn / la_nn ({n_ratio})'
                    )

            logging.info(f'Creating model {model_name} {time_info} with:')
            logging.info(
                f'  la_nn={params.la_nn[i]}, la_ns={params.la_ns[i]}, la_sn={params.la_sn[i]}, la_ss={params.la_ss[i]}')
            logging.info(f'  psi={params.psi[i]}, p={params.p[i]}, p_s={params.p_s[i]}')

            # Create a BDSS model with the parameters for this time interval
            model = BirthDeathWithSuperSpreadingModel(
                la_nn=params.la_nn[i],
                la_ns=params.la_ns[i],
                la_sn=params.la_sn[i],
                la_ss=params.la_ss[i],
                psi=params.psi[i],
                p=params.p[i],
                p_s=params.p_s[i],
                n_recipients=params.avg_recipients
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
        forest, (total_tips, u, T), ltt = generate(
            models,
            min_tips=params.min_tips,
            max_tips=params.max_tips,
            T=params.T,
            skyline_times=params.t if len(params.la_nn) > 1 else None,  # Only pass time points for skyline models
            max_notified_contacts=params.max_notified_contacts
        )

        # Save outputs
        save_forest(forest, params.nwk)
        save_log(models[0], total_tips, T, u, params.log)
        if params.ltt:
            save_ltt(ltt, observed_ltt(forest, T), params.ltt)

        logging.info("Simulation completed successfully")

    except RuntimeError as e:
        logging.error(f"Runtime error during simulation: {e}")
    except ValueError as e:
        logging.error(f"Value error during simulation: {e}")
    except Exception as e:
        logging.error(f"Unexpected error: {e}")


if '__main__' == __name__:
    main()