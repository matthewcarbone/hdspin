#!/usr/bin/env python3

__author__ = "Matthew R. Carbone & Marco Baity-Jesi"
__maintainer__ = "Matthew Carbone"
__email__ = "x94carbone@gmail.com"
__status__ = "Prototype"


import os
import shlex
import subprocess


def listdir_fp(d):
    return [os.path.join(d, f) for f in os.listdir(d)]


def get_cache():
    cache = os.environ.get("HDSPIN_CACHE_DIR")
    assert cache is not None
    assert isinstance(cache, str)
    return cache


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
