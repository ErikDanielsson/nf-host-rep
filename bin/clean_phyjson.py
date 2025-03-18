#!/usr/bin/env python

import sys
import json


def mod_symbiont_tree(json_obj):
    """
    Convert any age field in the json output to float
    """

    def mod_symbiont_tree_rec(ndict):
        """
        Recursive helper
        """
        if "__data__" not in ndict:
            return
        else:
            node = ndict["__data__"]
            node["age"] = float(node["age"])
            if "right" in node:
                mod_symbiont_tree_rec(node["right"])
            if "left" in node:
                mod_symbiont_tree_rec(node["left"])

    symbiont_tree_obj = json_obj["symbiont_tree"]
    new_obj = mod_symbiont_tree_rec(symbiont_tree_obj)
    return new_obj


def mod_host_distances(json_obj):
    """
    Convert any host distance in the json output to float
    """
    json_obj["host_distances"] = [float(d) for d in json_obj["host_distances"]]


# Read the input and output filenames
json_fn = sys.argv[1]
new_json_fn = sys.argv[2]

# Load the original json
with open(json_fn) as fh:
    json_obj = json.load(fh)

# Convert age fields
mod_symbiont_tree(json_obj)

# Convert host distance fields
mod_host_distances(json_obj)

# Write the correct json file
with open(new_json_fn, "w") as fh:
    json.dump(json_obj, fh, indent=4)
