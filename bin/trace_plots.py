#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import seaborn
import pandas as pd
import numpy as np
import json
from pathlib import Path


def revbayes_outfile_name(iter):
    return f"out.{iter}.log"


def treeppl_outfile_name(iter):
    return f"output.{iter}.json"


def read_revbayes(rev_out_files):
    parsed_rev_out_files = []
    for outfile in rev_out_files:
        parsed_rev_out_files.append(pd.read_csv(outfile, sep="\t"))
    return parsed_rev_out_files


def read_treeppl(tppl_out_files):
    parsed_tppl_out_files = []
    for outfile in tppl_out_files:
        with open(outfile) as fh:
            parsed_tppl_out_files.append(json.load(fh))
    return parsed_tppl_out_files


def proc_treeppl(tppl_out_files):
    samples_per_file = []
    for file in tppl_out_files:
        samples_per_file.append([s["__data__"] for s in file["samples"]])
    return samples_per_file


def main():
    rev_out_files = []
    tppl_out_files = []
    nsims = (len(sys.argv) - 1) // 2
    print(sys.argv)
    for i in range(1, nsims + 1):
        if not sys.argv[i].endswith(".log"):
            raise Exception(f"{sys.argv[i]} is not a RB output file. {nsims} {i}")
        rev_out_files.append(sys.argv[i])
    for i in range(nsims + 1, 2 * nsims + 1):
        if not sys.argv[i].endswith(".json"):
            raise Exception(f"{sys.argv[i]} is not a TPPL output file")
        tppl_out_files.append(sys.argv[i])

    rev_out_files = read_revbayes(rev_out_files)
    tppl_out_files = read_treeppl(tppl_out_files)
    tppl_samples = proc_treeppl(tppl_out_files)

    fig, axs = plt.subplots(2, nsims)
    for i, (rev_out_file, tppl_sample) in enumerate(zip(rev_out_files, tppl_samples)):
        mu_sample = [s["mu"] for s in tppl_sample]
        axs[0, i].plot(rev_out_file["Iteration"], rev_out_file["clock_host"])
        axs[1, i].plot(range(len(mu_sample)), mu_sample)
    plt.savefig("trace_plot.png")


main()
