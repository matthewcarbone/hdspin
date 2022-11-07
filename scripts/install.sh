#!/bin/bash

python -m pip install --upgrade wheel setuptools
pip install --verbose ".[test]"
