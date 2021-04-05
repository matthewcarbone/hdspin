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
    eval_parser = subparsers.add_parser(
        "eval", formatter_class=SortingHelpFormatter,
        description='Evaluates all results in the cache, saving them to final.'
    )

    eval_parser.add_argument(
        '--directory', dest='specified_directory', default=None,
        help="Will only evaluate the directory matching this string in the "
        "cache"
    )

    # Quick post processing on the value for beta_critical
    args = ap.parse_args(sys_argv)

    return args
