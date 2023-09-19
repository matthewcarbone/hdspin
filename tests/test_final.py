import argparse
from argparse import HelpFormatter, ArgumentDefaultsHelpFormatter
import json
from operator import attrgetter
from pathlib import Path
import sys

import numpy as np


class SortingHelpFormatter(ArgumentDefaultsHelpFormatter, HelpFormatter):
    def add_arguments(self, actions):
        actions = sorted(actions, key=attrgetter('option_strings'))
        super(SortingHelpFormatter, self).add_arguments(actions)


def add_check_2_args(ap):
    """Adds the QM9 parser options."""

    ap.add_argument(
        '-d1', '--directory1', dest='directory1', type=str,
        required=True, help='First directory of results.'
    )
    ap.add_argument(
        '-d2', '--directory2', dest='directory2', type=str,
        required=True, help='Second directory of results.'
    )
    ap.add_argument(
        '-s', '--scale', dest="scale", type=float, default=1.0,
        help='Sets the standard deviation scale for comparisons.'
    )
    ap.add_argument(
        '-t', '--threshold', dest="threshold", type=float, default=0.95,
        help='Threshold for comparisons.'
    )
    ap.add_argument(
        '-f', '--files', dest="files", type=str, nargs="+", default=["energy.txt"],
        help='Tiles to run.'
    )
    ap.add_argument(
        "--take-last-prop", dest="take_last_prop", type=float, default=1.0,
        help='Take the last proportion of the simulation average, '
        'allowing one to disregard short-time effects.'
    )


def global_parser(sys_argv):
    ap = argparse.ArgumentParser(formatter_class=SortingHelpFormatter)
    # ap.add_argument(
    #     '-d', '--directory', dest='directory', type=str, required=True,
    #     help='The location of the directory containing the results. Must '
    #     'contain a config.json file, and final and grids directories.'
    # )
    
    subparsers = ap.add_subparsers(
        help='Core testing options.', dest='test_type'
    )

    compare_subparser = subparsers.add_parser(
        "compare", formatter_class=SortingHelpFormatter,
        description='For comparing two sets of results.'
    )
    add_check_2_args(compare_subparser)

    return ap.parse_args(sys_argv)


def read_json(path):
    with open(path, 'r') as infile:
        dat = json.load(infile)
    return dat


def check_both_similar(x1, x2, scale=1.0):
    """Takes two arrays. Each array should have shape (N, >=3), where the three
    are the x value, y value and standard deviation of the y value (standard
    error is also sometimes in the next slot). This function compares these
    arrays. It asserts that for arbitrary xi,

    xi[:, 1] - xi[:, 2] * scale < xj[:, 1]

    and

    xj[:, 1] < xi[:, 1] + xi[:, 2] * scale

    Essentially checking whether the curves fall within +/- scale standard
    deviations of each other.
    
    Parameters
    ----------
    x1 : np.ndarray
    x2 : np.ndarray
    scale : float, optional
    
    Returns
    -------
    float
        The proportion of True conditions. Ideally, this is always 1.0.
    """

    cond1 = (x1[:, 1] - x1[:, 2] * scale) < x2[:, 1]
    cond2 = x2[:, 1] < (x1[:, 1] + x1[:, 2] * scale)
    cond3 = (x2[:, 1] - x2[:, 2] * scale) < x1[:, 1]
    cond4 = x1[:, 1] < (x2[:, 1] + x2[:, 2] * scale)
    return (cond1 & cond2 & cond3 & cond4).mean()


def run_all_check_both_similar(args):
    """Executes all check_both_similar tests."""

    print(">>> Running comparison tests <<<")

    # Define some paths
    df1 = Path(args.directory1) / "final"
    df2 = Path(args.directory2) / "final"

    # Load results for standard things
    problems = []
    for file in args.files:

        path1 = df1 / file
        if not path1.exists():
            print(f"⚠️ Path 1: {path1} DNE, continuing...")

        path2 = df2 / file
        if not path2.exists():
            print(f"⚠️ Path 2: {path2} DNE, continuing...")

        x1 = np.loadtxt(path1)
        x2 = np.loadtxt(path2)
        L = int(x1.shape[0] * args.take_last_prop)
        x1 = x1[-L:, :]
        x2 = x2[-L:, :]

        similarity = check_both_similar(x1, x2, scale=args.scale)

        if similarity < args.threshold:
            message = f"❌ {file} had sim={similarity:.03f}<{args.threshold}"
            problems.append(message)
        else:
            message = f"✅ {file} had sim={similarity:.03f}>={args.threshold}"
        print(message)

    return problems
    

if __name__ == '__main__':
    args = global_parser(sys.argv[1:])

    problems = []

    if args.test_type == "compare":
        problems.extend(run_all_check_both_similar(args))
    else:
        # Load configs
        config1 = read_json(Path(args.directory1) / "config.json")
        config2 = read_json(Path(args.directory2) / "config.json")

        # Get beta
        beta = config1["beta"]
        assert beta == config2["beta"]

    if len(problems) > 0:
        print("-" * 80)
        for problem in problems:
            print(problem)
        raise RuntimeError("There were problems. See above.")
    
    print("Checks all successful.")
