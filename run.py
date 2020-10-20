#!/usr/bin/env python3

__author__ = "Matthew R. Carbone & Marco Baity-Jesi"
__maintainer__ = "Matthew Carbone"
__email__ = "x94carbone@gmail.com"
__status__ = "Prototype"

import os
import shutil
import sys

from py.parser import global_parser
from py import utils as u
from py.eval import Evaluator

# Main method -----------------------------------------------------------------

if __name__ == '__main__':

    args = global_parser(sys.argv[1:])

    if args.cleanup:
        will_delete = args.cache
        user_input = 'no'
        if args.force:
            user_input = 'yes'
        else:
            user_input = input(
                f"Delete tree {will_delete}, runfiles and job data? [yes]/no\n"
            )
        if user_input == 'yes':
            print("Deleting...")
            if os.path.exists(will_delete):
                shutil.rmtree(will_delete)
            if os.path.exists("local_submit.sh"):
                os.remove("local_submit.sh")
            if os.path.exists("job_data"):
                shutil.rmtree("job_data")
            if os.path.exists("main.out"):
                os.remove("main.out")
            print("Done")
        else:
            print("Exiting")

    elif args.protocol == 'prime':
        base_dir, max_iter = u.make_directory_and_configs(args)

        if args.local:
            u.write_bash_script(args, base_dir, max_iter)
        else:
            u.write_SLURM_script(args, base_dir, max_iter)

    elif args.protocol == 'execute':
        if not os.path.exists("exe/main.out"):
            raise RuntimeError("Run Make before using run.py execute")
        u.submit(args)

    elif args.protocol == 'eval':
        ev = Evaluator(args)
        ev.eval_traj()
        ev.eval_psi_config()
        ev.eval_aging_config()
        ev.eval_psi_basin()
        ev.eval_aging_basin()

    else:
        raise RuntimeError("Unknown error")
