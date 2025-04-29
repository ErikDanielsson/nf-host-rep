#! /usr/bin/env python
import json
import sys
import csv

interaction_sim_fn = sys.argv[1]
csv_out_fn = sys.argv[2]

with open(interaction_sim_fn) as fh:
    # Since we are using incremental printing we will read the first line of
    # the file and then parse it. Since the TreePPL simulation script samples
    # directly from the true distribution we don't need to add any burnin
    json_line = next(fh)
    parsed = json.loads(json_line)
    data = parsed["__data__"]
    data.pop("tree")
    interactions = data["interactions"]
    parsed_interactions = [
        (i["__data__"]["label"], i["__data__"]["repertoire"]) for i in interactions
    ]
    sorted_interactions = sorted(parsed_interactions)

with open(csv_out_fn, "w") as fh:
    # Write the interactions to csv so that we can subsequently write them to Nexus in R...
    nhosts = len(parsed_interactions[0][1])
    writer = csv.writer(fh, delimiter="\t")
    header = [""] + [f"Host{i+1}" for i in range(nhosts)]
    writer.writerow(header)
    for label, parasite_ints in sorted_interactions:
        writer.writerow([f"Symbiont{label}"] + parasite_ints)
