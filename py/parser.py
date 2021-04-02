#!/usr/bin/env python3

import argparse
from argparse import HelpFormatter, ArgumentDefaultsHelpFormatter
from operator import attrgetter


# https://stackoverflow.com/questions/
# 12268602/sort-argparse-help-alphabetically
class SortingHelpFormatter(ArgumentDefaultsHelpFormatter, HelpFormatter):
    def add_arguments(self, actions):
        actions = sorted(actions, key=attrgetter('option_strings'))
        super(SortingHelpFormatter, self).add_arguments(actions)


def global_parser(sys_argv):
    ap = argparse.ArgumentParser(formatter_class=SortingHelpFormatter)

    ap.add_argument(
        '--force', dest='force', default=False, action='store_true',
        help='Overrides failsafes for e.g. overwriting directories.'
    )
    ap.add_argument(
        '--purge', dest='purge', default=False,
        action='store_true',
        help='Deletes all stored data in the cache directory.'
    )
    ap.add_argument(
        '-i', '--input', dest='config_location', type=str,
        default='config.yaml',
        help="Input file location."
    )

    subparsers = ap.add_subparsers(
        help='Choices for various priming, execution and post-processing '
        'protocols.', dest='protocol'
    )

    # (1) ---------------------------------------------------------------------
    _ = subparsers.add_parser(
        "prime", formatter_class=SortingHelpFormatter,
        description='Prime the computation for submission by creating the '
        'appropriate directories and writing the SLURM submit file.'
    )

    # required_prime = prime_sp.add_argument_group(
    #     'required'
    # )

    # required_prime.add_argument(
    #     '-T', '--timesteps', dest='timesteps', type=int, required=True,
    #     help='The log10 number of timesteps of the standard simulation/'
    #     'approximate maximum simulation time for the Gillespie simulation. '
    #     'Note that in the case of Gillespie, the simulation is guaranteed to '
    #     'reach at least this value.'
    # )
    # required_prime.add_argument(
    #     '-S', '--n-sim', dest='nsim', type=int, required=True,
    #     help='Total number of tracers to simulate. Each tracer will have its '
    #     'own randomly initialized energy landscape.'
    # )
    # required_prime.add_argument(
    #     '-N', '--n-spin', dest='nspin', type=int, required=True,
    #     help='Total number of tracers to simulate. Each tracer will have its '
    #     'own randomly initialized energy landscape.'
    # )
    # required_prime.add_argument(
    #     '-D', '--dynamics', dest='dynamics', type=int, choices=[0, 1],
    #     required=True,
    #     help='Sets the dynamics of the simulation. Use 0 for standard '
    #     'dynamics and 1 for Gillespie.'
    # )
    # required_prime.add_argument(
    #     '-L', '--landscape', dest='landscape', type=int, choices=[0, 1],
    #     required=True,
    #     help='Sets the energy landscape of the simulation. Use 0 for EREM and '
    #     '1 for REM.'
    # )
    # required_prime.add_argument(
    #     '-B', '--beta', dest='beta', type=float, required=True,
    #     help='Sets the inverse temperature for the simulation.'
    # )

    # prime_sp.add_argument(
    #     '--beta-critical', dest='beta_critical', type=float, default=None,
    #     help='Sets the critical inverse temperature for the simulation. If '
    #     'None, then this will default to 1.0 for EREM and 1.1778 for REM. '
    #     'Note that the latter is ~sqrt(2 log 2).'
    # )

    # prime_sp.add_argument(
    #     '--dw', dest='dw', type=float, default=0.5,
    #     help='Sets the parameter in the aging functions.'
    # )

    # prime_sp.add_argument(
    #     '--loopN', dest='loopN', default=False, action='store_true',
    #     help='Uses the loop over N dynamics.'
    # )

    # prime_sp.add_argument(
    #     '--local', dest='local', default=False, action='store_true',
    #     help='Saves a bash script to the working directory instead of a SLURM '
    #     'script to the cache.'
    # )

    # (2) ---------------------------------------------------------------------
    execute_sp = subparsers.add_parser(
        "execute", formatter_class=SortingHelpFormatter,
        description='Runs all primed jobs.'
    )

    execute_sp.add_argument(
        '--run1', dest='run_one', default=False,
        action='store_true',
        help='If flagged, will run/submit only one primed job.'
    )

    execute_sp.add_argument(
        '--local', dest='local', default=False,
        action='store_true',
        help='Runs only a single primed job locally, used for debugging.'
    )

    # (3) ---------------------------------------------------------------------
    subparsers.add_parser(
        "eval", formatter_class=SortingHelpFormatter,
        description='Evaluates all results in the cache, saving them to final.'
    )

    # Quick post processing on the value for beta_critical
    args = ap.parse_args(sys_argv)

    return args
