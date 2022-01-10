import numpy as np
import math
import sys
from functools import cmp_to_key

def read_simplices():
    f = open(input_file, "r")
    simplices = []
    for line in f:
        if line.startswith("#s"):
            dim = int(line[3:].split("-")[0].strip())
            ind = int(line[3:].split("-")[1].strip())
            simplices.append([dim, ind])
    return simplices

def read_mat():
    f = open(mat_file, "r")
    cols = []
    for line in f:
        col = []
        if line.strip() != "":
            for tok in line.strip().split(" "):
                dim = int(tok.split("-")[0])
                ind = int(tok.split("-")[1])
                col.append([dim, ind])
        cols.append(col)
    return cols

mat_file = "./mat-r.txt"
input_file = "./output/4pts_visu/4pts_pc_d3tMr0fu00.txt"
f = open("../thesis/gen_latex/4pts_mat_r.txt", "w")

simplices = read_simplices()
cols = read_mat()

bounds = []
for c in cols:
    m = 1
    for e in c:
        ind = simplices.index(e)
        m = max(m, ind)
    bounds.append(m)


use_bounds = False

for r in range(len(simplices)):
    for c in range(len(simplices)):
        col = cols[c]
        rowele = simplices[r]
        if rowele in col:
            f.write(" 1 ")
        elif (r < bounds[c]) or (not use_bounds):
            f.write(" \\cdot ")
        else:
            f.write(" ")
        if c < len(simplices) - 1:
            f.write(" & ")
        elif r < len(simplices) - 1:
            f.write(" \\\\\n")
