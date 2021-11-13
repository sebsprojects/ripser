import matplotlib
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from matplotlib import colors as mcol
from matplotlib.lines import Line2D
import numpy as np
import math
import sys
from functools import cmp_to_key

if len(sys.argv) < 1:
    print("error: missing config file-path")
    sys.exit()

# PARAMETER
# -----------------------------------------------------------------------------
dims = [0,1,2]
doc_w = 418.25372 / 72
figw = doc_w * 1.2  # FACTOR 1 TOO SMALL FOR SOME STUPID REASON
input_file = sys.argv[1]


# MATPLOTLIB SETUP
# -----------------------------------------------------------------------------

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
plt.rcParams.update(tex_fonts)
fig, ax = plt.subplots(figsize=(figw, figw))


# FUNCTION DEFS
# -----------------------------------------------------------------------------

def comp_bar(a, b):
    if a[0] < b[0] or (a[0] == b[0] and a[1] < b[1]):
        return -1
    elif a == b:
        return 0
    else:
        return 1

def read_barcode(simplex_lookup):
    f = open(input_file, "r");
    intervals = []
    max_bound = 0
    offs = 1
    for line in f:
        if line.startswith("#b"):
            dim = int(line[2])
            if dim + 1 > len(intervals):
                intervals.append([dim,[]])
            toks = [t.strip() for t in line[3:].split("-")]
            sint = toks[1].replace("[","").replace(")", "").split(",")
            start_index = int(sint[0])
            iint = [(simplex_lookup.index([dim, start_index])) + offs]
            if sint[1].strip() != "":
                end_index = int(sint[1])
                iint.append(simplex_lookup.index([dim+1, end_index]) + offs)
                max_bound = max(max_bound, iint[1])
            max_bound = max(max_bound, iint[0])
            intervals[dim][1].append(iint)
    for dim, b in enumerate(intervals):
        for i, interval in enumerate(b[1]):
            if len(interval) == 1:
                intervals[dim][1][i].append(max_bound + 1)
        intervals[dim][1] = sorted(intervals[dim][1], key=cmp_to_key(comp_bar))
    return (max_bound, intervals)

def read_simplices():
    f = open(input_file, "r")
    simplex_lookup = []
    for line in f:
        if line.startswith("#s"):
            dim = int(line[2])
            index = int(line[3:].split("-")[0].strip())
            simplex_lookup.append([dim, index])
    return simplex_lookup

def plot_barcode(max_bound, intervals):
    h_inc = (max_bound + 1) * 0.02
    h_skip = h_inc * 1.5
    h = 0
    ax.set_xlim(1, max_bound+1)
    ax.set_xticks(range(1, max_bound+1))
    ax.xaxis.grid(True, linestyle="dotted")
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for bc in intervals:
        dim = bc[0]
        bars = bc[1]
        h += h_skip
        for b in bars:
            ax.plot(b, [h, h], color=colors[dim])
            if b != bars[-1]:
                h += h_inc
    h += h_skip
    ax.set_ylim(0,h)
    ax.set_aspect("equal")
    #ax2 = ax.twiny()
    #ax2.set_xlim(ax.get_xlim())
    #ax2.set_xticks(ax.get_xticks())
    #ax2.spines["right"].set_visible(False)
    #ax2.spines["left"].set_visible(False)
    fig.set_figheight(h * 0.5)


# MAIN
# -----------------------------------------------------------------------------

simplex_lookup = read_simplices()
max_bound, intervals = read_barcode(simplex_lookup)
plot_barcode(max_bound, intervals)
plt.savefig("../thesis/img/test.pdf", format='pdf', bbox_inches="tight", pad_inches=0.05)
#plt.show()
