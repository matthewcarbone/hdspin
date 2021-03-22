#!/usr/bin/env python3

__author__ = "Matthew R. Carbone & Marco Baity-Jesi"
__maintainer__ = "Matthew Carbone"
__email__ = "x94carbone@gmail.com"
__status__ = "Prototype"


import sys

from py.parser import global_parser
from py.eval import Evaluator
from py.executors import Primer, Executor, cleanup

# Main method -----------------------------------------------------------------

if __name__ == '__main__':

    args = global_parser(sys.argv[1:])

    if args.purge:
        cleanup()

    elif args.protocol == 'prime':
        primer = Primer(args.config_location)
        primer.prime()

    elif args.protocol == 'execute':
        executor = Executor(args.config_location, args.run_one, args.local)

    elif args.protocol == 'eval':
        ev = Evaluator(args)
        ev.eval_traj()
        ev.eval_psi_config()
        ev.eval_aging_config()
        ev.eval_psi_basin()
        ev.eval_aging_basin()

    else:
        raise RuntimeError("Unknown error")
