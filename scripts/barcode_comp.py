import matplotlib
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from matplotlib import colors as mcol
from matplotlib.lines import Line2D
import numpy as np
import math
import sys
from functools import cmp_to_key
from matplotlib.ticker import FormatStrFormatter


if len(sys.argv) < 1:
    print("error: missing config file-path")
    sys.exit()


# PARAMETER
# -----------------------------------------------------------------------------
dims = [0,1,2]
doc_w = 418.25372 / 72
figw = doc_w * 1.2  # FACTOR 1 TOO SMALL FOR SOME STUPID REASON


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


# FUNCTION DEFS
# -----------------------------------------------------------------------------

def comp_bar(a, b):
    if a[0] < b[0] or (a[0] == b[0] and a[1] < b[1]):
        return -1
    elif a == b:
        return 0
    else:
        return 1

def read_barcode(input_file):
    f = open(input_file, "r");
    intervals = []
    max_bound = 0
    bounds = []
    for line in f:
        if line.startswith("#b"):
            dim = int(line[2])
            if dim + 1 > len(intervals):
                intervals.append([dim,[]])
            toks = [t.strip() for t in line[3:].split("-")]
            sint = toks[0].replace("[","").replace(")", "").split(",")
            if sint[0] == sint[1]:
                continue
            iint = [float(sint[0])]
            bounds.append(float(sint[0]))
            if sint[1].strip() != "":
                iint.append(float(sint[1]))
                bounds.append(float(sint[1]))
                max_bound = max(max_bound, iint[1])
            max_bound = max(max_bound, iint[0])
            intervals[dim][1].append(iint)
    for dim, b in enumerate(intervals):
        for i, interval in enumerate(b[1]):
            if len(interval) == 1:
                intervals[dim][1][i].append(max_bound * 1.1)
        intervals[dim][1] = sorted(intervals[dim][1], key=cmp_to_key(comp_bar))
    return (sorted(list(set(bounds))), intervals)

def plot_barcode(ax, bounds, intervals):
    h_inc = (bounds[-1] * 1.1) * 0.02
    h_skip = h_inc * 1.5
    h = 0
    ax.set_xlim(0, bounds[-1] * 1.1)
    #ax.set_xticks(bounds)
    #ax.xaxis.grid(True, linestyle="dotted")
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3g'))
    for bc in intervals:
        dim = bc[0]
        bars = bc[1]
        if len(bars) == 0:
            continue
        h += h_skip
        for b in bars:
            #print(b, h)
            ax.plot(b, [h, h], color=colors[dim])
            if b != bars[-1]:
                h += h_inc
    h += h_skip
    ax.set_ylim(0,h)
    #ax2 = ax.twiny()
    #ax2.set_xlim(ax.get_xlim())
    #ax2.set_xticks(bounds)
    #ax2.spines["right"].set_visible(False)
    #ax2.spines["left"].set_visible(False)
    #ax.set_aspect("equal")
    #fig.set_figwidth(figw)


# MAIN
# -----------------------------------------------------------------------------

n = len(sys.argv) - 1
fig, axs = plt.subplots(n, 1, figsize=(figw, figw))

for i in range(n):
    input_file = sys.argv[i + 1]
    bounds, intervals = read_barcode(input_file)
    plot_barcode((axs if n == 1 else (axs.flat[i])), bounds, intervals)
#plt.savefig("../thesis/img/test1.pdf", format='pdf', bbox_inches="tight", pad_inches=0.05)
plt.savefig("diam_rings100_abs_0-49_0-74.pdf", format='pdf', bbox_inches="tight", pad_inches=0.05)
#plt.show()
