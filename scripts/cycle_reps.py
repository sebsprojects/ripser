# Parse .init from cmd line arg
# Check if all cycles really are cycles
# Visualize data set (relative part in different color)
# Visualize cycles (plus relative edges)

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from matplotlib import colors as mcol
from matplotlib.lines import Line2D
import numpy as np
import math
import sys
import configparser

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
PI = math.pi
tex_fonts = {
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 11,
    "font.size": 11,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    'text.latex.preamble': r'\usepackage{amsfonts}'
}
#plt.rcParams.update(tex_fonts)
doc_w = 418.25372 / 72.27
figw = doc_w * 0.5

fig, ax = plt.subplots(figsize=(figw, figw))
ax.set_aspect('equal')

if len(sys.argv) < 2:
    print("error: missing config file-path")

config_path = sys.argv[1]
cp = configparser.ConfigParser()
cp.read(config_path)

input_path = cp["ripser"]["file_path"]
input_type = cp["ripser"]["input_type"]
output_path = cp["ripser"]["output_path"]
relative_subcomplex_string = cp["ripser"]["relative_subcomplex"]
relative_subcomplex = []

if not relative_subcomplex_string == "":
    rel_toks = relative_subcomplex_string.split(",")
    for tok in rel_toks:
        start = 0
        end = -1
        toks = tok.split("-")
        if len(toks) == 1:
            start = int(toks[0])
            end = int(toks[0])
        else:
            if len(toks[0]) > 0:
                start = int(toks[0])
            if len(toks[1]) > 0:
                end = int(toks[1])
        relative_subcomplex.append([start, end])


def read_samples(input_path, dims=[0,1]):
    input_file = open(input_path, "r")
    data = [[], []]
    for (i, line) in enumerate(input_file):
        toks = line.split(",")
        data[0].append(float(toks[dims[0]]))
        data[1].append(float(toks[dims[1]]))
    for (i, el) in enumerate(relative_subcomplex):
        relative_subcomplex[i][1] = len(data[0]) - 1 if el[1] < 0 else el[1]
    rxy = [[],[]]
    xy = [[],[]]
    for i, dat in enumerate(zip(data[0], data[1])):
        is_rel = False
        for j in relative_subcomplex:
            if i >= j[0] and i <= j[1]:
                is_rel = True
                break
        if is_rel:
            rxy[0].append(dat[0])
            rxy[1].append(dat[1])
        else:
            xy[0].append(dat[0])
            xy[1].append(dat[1])
    return [data, xy, rxy]

def read_reps(input_path, output_path):
    reps_filename = output_path + "/" + input_path.split("/")[-1] + "_reps.txt"
    reps_file = open(reps_filename, "r");
    reps = []
    for line in reps_file:
        if line.startswith("#"):
            interval = list(map(float, line[1:].split()))
            reps.append({ "birth": interval[0], "death": interval[1], "rep": []})
        elif not line == "\n":
            reps[-1]["rep"].append(list(map(int, line.split())))
    return reps

def plot_samples(data_abs, data_rel):
    #ax.set_xlim(-1.2, 1.2)
    #ax.set_ylim(-1.2, 1.2)
    ax.scatter(data_abs[0], data_abs[1], s=7, c=colors[0])
    ax.scatter(data_rel[0], data_rel[1], s=7, c=colors[1])

def plot_reps(reps, data_all, selection):
    for s in selection:
        lines = []
        for r in reps[s]["rep"]:
            lines.append([(data_all[0][r[0]], data_all[1][r[0]]), (data_all[0][r[1]], data_all[1][r[1]])])
        lc = mc.LineCollection(lines)
        ax.add_collection(lc)

def verify_reps(reps):
    for rep in reps:
        parity = {}
        for l in rep["rep"]:
            if l[0] not in parity:
                parity[l[0]] = 0
            if l[1] not in parity:
                parity[l[1]] = 0
            parity[l[0]] = parity[l[0]] + 1
            parity[l[1]] = parity[l[1]] + 1
        for (key, val) in parity.items():
            if val != 2:
                pass
                #print("ERROR", rep, key, val)

data_all, data_abs, data_rel = read_samples(input_path)
reps = read_reps(input_path, output_path)

plot_samples(data_abs, data_rel)
plot_reps(reps, data_all, [0])
plt.show()
#verify_reps(reps)
#plot_samples(data_abs, data_rel)
