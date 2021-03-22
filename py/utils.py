#!/usr/bin/env python3

__author__ = "Matthew R. Carbone & Marco Baity-Jesi"
__maintainer__ = "Matthew Carbone"
__email__ = "x94carbone@gmail.com"
__status__ = "Prototype"


import os


def get_cache():
    cache = os.environ.get("HDSPIN_CACHE_DIR")
    assert cache is not None
    assert isinstance(cache, str)
    return cache
