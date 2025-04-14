#! /usr/bin/env python
import csv
import sys


def rev_row_to_row(row, writer):
    name = row[0]
    if name.startswith("Symbiont"):
        writer.writerow([name] + [int(c) for c in row[1]])


in_fn = sys.argv[1]
out_fn = sys.argv[2]
with open(in_fn) as in_fh, open(out_fn, "w", newline="") as out_fh:
    reader = csv.reader(in_fh)
    writer = csv.writer(out_fh, delimiter=",")
    initial_row = next(reader)
    nchars = len(initial_row[1])
    writer.writerow([""] + [f"Host{i+1}" for i in range(nchars)])
    rev_row_to_row(initial_row, writer)
    for r in reader:
        rev_row_to_row(r, writer)
