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

# Main method -----------------------------------------------------------------

if __name__ == '__main__':

    args = global_parser(sys.argv[1:])

    if args.cleanup:
        will_delete = args.cache
        user_input = 'no'
        if args.force:
            user_input = 'yes'
        else:
            user_input = input(f"Delete tree {will_delete}? [yes]/no\n")
        if user_input == 'yes':
            print("Deleting...")
            shutil.rmtree(will_delete)
            print("Done")
        else:
            print("Exiting")

    elif args.protocol == 'prime':
        base_dir, max_iter = u.make_directory_and_configs(args)
        u.write_SLURM_script(args, base_dir, max_iter)

    elif args.protocol == 'execute':
        if not os.path.exists("main.out"):
            raise RuntimeError("Run Make before using run.py execute")
        u.submit(args)

    else:
        raise RuntimeError("Unknown error")
