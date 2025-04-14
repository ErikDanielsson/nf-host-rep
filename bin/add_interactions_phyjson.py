#! /usr/bin/env python
import json
import csv
import sys
import re

phyjson_in_fn = sys.argv[1]
interactions_csv = sys.argv[2]
phyjson_out_fn = sys.argv[3]

with open(interactions_csv) as fh:
    reader = csv.reader(fh, delimiter="\t")
    # The first line is the header
    _ = next(reader)
    interactions_and_label = [(i[0], i[1:]) for i in reader]
    pattern = re.compile(r"Symbiont(\d+)")
    # Sort the interactions one last time...
    sorted_interactions = sorted(
        interactions_and_label, key=lambda x: int(pattern.match(x[0]).groups()[0])
    )
    # Flatten the interactions in row major order
    interactions = [
        int(i) for _, parasite_ints in sorted_interactions for i in parasite_ints
    ]

with open(phyjson_in_fn) as fh:
    phyjson = json.load(fh)
    phyjson["interactions"] = interactions

with open(phyjson_out_fn, "w") as fh:
    json.dump(phyjson, fh, indent=2)
