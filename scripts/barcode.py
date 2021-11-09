import matplotlib
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from matplotlib import colors as mcol
from matplotlib.lines import Line2D
import numpy as np
import math
import sys

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

def read_reps(output_file):
    reps_file = open(output_file, "r");
    reps = []
    current_rep = { "dim": -1, "birth": -1, "death": -1, "rep": []}
    for line in reps_file:
        if line.startswith("#"):
            if current_rep["dim"] >= 0:
                reps.append(current_rep)
            current_rep = { "dim": -1, "birth": -1, "death": -1, "rep": []}
            toks = list(line[1:].split())
            current_rep["dim"] = int(toks[0].split('=')[1])
            current_rep["birth"] = float(toks[-2])
            current_rep["death"] = float(toks[-1])
        elif not line == "\n":
            current_rep["rep"].append(list(map(int, line.split())))
    return reps

def plot_barcode(reps_1, reps_2):
    bar_lines = []
    bar_colors = []
    h = 0.1
    for r in reps_1:
        if r["dim"] == 0:
            continue
        color_index = 1 if r["dim"] == 0 else 0
        bar_lines.append([(r["birth"],h), (r["death"],h)])
        bar_colors.append(np.array(mcol.to_rgba(colors[color_index], 1)))
        h += 0.02
    h = 0.09
    for r in reps_2:
        if r["dim"] == 0:
            continue
        color_index = 0 if r["dim"] == 0 else 1
        bar_lines.append([(r["birth"],h), (r["death"],h)])
        bar_colors.append(np.array(mcol.to_rgba(colors[color_index], 1)))
        h += 0.02
    lc = mc.LineCollection(bar_lines, colors=bar_colors, linewidths=1.5)
    ax.add_collection(lc)

reps_abs = read_reps("output/rings100.pc_reps_thresh1d959776_relnone.txt")
reps_rel = read_reps("output/rings100.pc_reps_thresh1d959776_rel0-49.txt")

ax.set_xlim((0,1))
ax.set_ylim((0,0.5))

plot_barcode(reps_abs, reps_rel)

plt.show()
