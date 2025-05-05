#! /usr/bin/env python
import json
import sys
import csv

interaction_sim_fn = sys.argv[1]
host_name_to_inds = sys.argv[2]
symbiont_name_to_inds = sys.argv[3]
csv_out_fn = sys.argv[4]

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

inds_to_host_name = {}
with open(host_name_to_inds) as fh:
    reader = csv.reader(fh, delimiter="\t")
    for row in reader:
        name = row[0]
        index = int(row[1])
        inds_to_host_name[index] = name

inds_to_symbiont_name = {}
with open(symbiont_name_to_inds) as fh:
    reader = csv.reader(fh, delimiter="\t")
    for row in reader:
        name = row[0]
        index = int(row[1])
        inds_to_symbiont_name[index] = name

with open(csv_out_fn, "w") as fh:
    # Write the interactions to csv so that we can subsequently write them to Nexus in R...
    nhosts = len(parsed_interactions[0][1])
    writer = csv.writer(fh, delimiter="\t")
    header = [""] + [inds_to_host_name[i] for i in range(1, nhosts + 1)]
    writer.writerow(header)
    for label, parasite_ints in sorted_interactions:
        writer.writerow([inds_to_symbiont_name[label]] + parasite_ints)
