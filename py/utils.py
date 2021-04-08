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


def get_general_filename(
    base_fname, inherent_structure, energetic_threshold,
    extra_text=None, extension=".txt"
):
    """

    Parameters
    ----------
    base_fname : str
        The base filename, such as "_psi_basin".
    inherent_structure : bool
    energetic_threshold : bool, optional
        If None, will ignore appending any string corresponding to if this
        is the energetic threshold or not.
    extra_text : str, optional
    extension : str, optional
        The etension to append at the end of the final string (the default
        is ".txt").

    Returns
    -------
    str
        The full file name.
    """

    f = base_fname

    if energetic_threshold is not None:
        if energetic_threshold:
            f += "_E"
        else:
            f += "_S"

    if int(inherent_structure) == 1:
        f += "_IS"
    elif int(inherent_structure) == 2:
        f += "_proxy_IS"

    if extra_text is not None:
        assert isinstance(extra_text, str)
        f += extra_text

    f += extension

    return f
