#!/usr/bin/env python3

from copy import copy
import json
from pathlib import Path
from shutil import copyfile
import sys


def save_json(d, path):
    with open(path, 'w') as outfile:
        json.dump(d, outfile, indent=4, sort_keys=True)


def read_json(path):
    with open(path, 'r') as infile:
        dat = json.load(infile)
    return dat


if __name__ == '__main__':

    d = read_json("input.json")
    original_d = copy(d)
    nvals = d.pop("N_spins")
    betas = d.pop("beta")
    dynamics = d.pop("dynamics_or_threshold")

    path = Path(sys.argv[1])

    for beta in betas:
        for N in nvals:

            beta = round(beta, 4)

            fname = path / Path(f"{N:02}_{beta:.04f}")
            fname.mkdir(exist_ok=True, parents=True)

            save_json(original_d, fname / Path("base_input.json"))

            copyfile("submit.sh", str(fname / Path("submit.sh")))
            copyfile("calc.sbatch.sh", str(fname / Path("calc.sbatch.sh")))
            copyfile("eval_1.sbatch.sh", str(fname / Path("eval_1.sbatch.sh")))

            _d = copy(d)

            if dynamics in ["gillespie", "standard"]:
                _d["dynamics"] = dynamics
            else:
                assert isinstance(dynamics, float)
                if beta < dynamics:
                    _d["dynamics"] = "standard"
                else:
                    _d["dynamics"] = "gillespie"

            _d["N_spins"] = N
            _d["beta"] = beta

            save_json(_d, fname / Path("input.json"))
