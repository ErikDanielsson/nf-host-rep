#! /usr/bin/env python
import json
import csv
import sys

"""
This file adds the model parameters specified in the file `param_fn` to
an existing phyjson file. It should be a csv file with a header with names:

- mu
- beta
- lambda01
- lambda10
- lambda12
- lambda21

and single row with the corresponding parameter values we should simulates from
"""
phyjson_in_fn = sys.argv[1]
params_fn = sys.argv[2]
phyjson_out_fn = sys.argv[3]

with open(phyjson_in_fn) as fh:
    phyjson = json.load(fh)
print(phyjson)
with open(params_fn) as fh:
    reader = csv.reader(fh, delimiter="\t")
    header_line = map(str, next(reader))
    params_line = map(float, next(reader))
    params = dict(zip(header_line, params_line))
print(params)

phyjson.update(params)
print(phyjson)

with open(phyjson_out_fn, "w") as fh:
    json.dump(phyjson, fh, indent=2)
