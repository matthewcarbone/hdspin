#!/usr/bin/env python3

__author__ = "Matthew R. Carbone & Marco Baity-Jesi"
__maintainer__ = "Matthew Carbone"
__email__ = "x94carbone@gmail.com"
__status__ = "Prototype"

import numpy as np

import os
import shlex
import subprocess
import yaml
import warnings


MAX_LONG = 9223372036854775807


def listdir_fp(d):
    return [os.path.join(d, f) for f in os.listdir(d)]


def run_command(command, silent=True):
    """https://www.endpoint.com/blog/2015/01/28/
    getting-realtime-output-using-python"""

    process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)

    while True:
        output = process.stdout.readline()
        if output == b'' and process.poll() is not None:
            break
        if output and not silent:
            print(output.strip().decode())

    rc = process.poll()
    return rc


def submit(args):
    """Submits all jobs in the cache to the job controller."""

    jobs = os.listdir(args.cache)
    print("Preparing to submit jobs:")
    for j in jobs:
        print(j)

    if not args.force:
        user_input = input("Continue? [yes]/no\n")
        if user_input != 'yes':
            print("Exiting")
            return

    print("Submitting jobs!")

    for j in jobs:
        script = os.path.join(args.cache, j, "scripts/submit.sh")
        run_command(f'mv {script} .')
        run_command("sbatch submit.sh")
        run_command(f'mv submit.sh {args.cache}/{j}/scripts')
