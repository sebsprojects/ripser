import matplotlib
import matplotlib.pyplot as plt
import matplotlib.figure as mplfig
import matplotlib.axes as mplax
from matplotlib import collections as mc
from matplotlib import colors as mcol
from matplotlib.lines import Line2D
import matplotlib.patches as pat
import numpy as np
import math
import sys
from functools import cmp_to_key

if len(sys.argv) < 2:
    print("error: missing config file-path")
    sys.exit()

# PARAMETER
# -----------------------------------------------------------------------------


# MATPLOTLIB SETUP
# -----------------------------------------------------------------------------

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
PI = math.pi
tex_fonts = {
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 9,
    "font.size": 9,
    "legend.fontsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    'text.latex.preamble': r'\usepackage{amsfonts}\usepackage{bm}'
}
plt.rcParams.update(tex_fonts)


# FUNCTION DEFS
# -----------------------------------------------------------------------------

def comp_bar(a, b):
    if a[0] < b[0] or (a[0] == b[0] and a[1] > b[1]):
        return -1
    elif a == b:
        return 0
    else:
        return 1

def read_barcode(input_file, simplex_lookup):
    f = open(input_file, "r");
    intervals = []
    max_bound = 0
    offs = 1
    for line in f:
        if line.startswith("#b"):
            dim = int(line[2])
            if dim + 1 > len(intervals):
                intervals.append([dim,[]])
            toks = [t.strip() for t in line[3:].split(";")]
            sint = toks[1].replace("[","").replace(")", "").split(",")
            start_index = int(sint[0].split("-")[0])
            iint = [(simplex_lookup.index([dim, start_index])) + offs]
            if sint[1].strip() != "":
                end_index = int(sint[1].split("-")[0])
                iint.append(simplex_lookup.index([dim+1, end_index]) + offs)
                max_bound = max(max_bound, iint[1])
            max_bound = max(max_bound, iint[0])
            intervals[dim][1].append(iint)
    for dim, b in enumerate(intervals):
        intervals[dim][1] = sorted(intervals[dim][1], key=cmp_to_key(comp_bar))
    return (max_bound, intervals)

def read_simplices(input_file, exclude_rel=False):
    f = open(input_file, "r")
    simplex_lookup = []
    rel = range(8)
    for line in f:
        if line.startswith("#s"):
            dim = int(line[2])
            index = int(line[3:].split("-")[0].strip())
            #simplices = [int(s.strip()) for s in line[3:].split("-")[1].split("'")]
            simplex_lookup.append([dim, index])
    return simplex_lookup

def init_ax(ax, max_bound):
    ax.set_xlim(0, max_bound)
    ax.set_ylim(0,1)
    ax.set_xticks(range(0, max_bound+1))
    ax.xaxis.grid(True, linestyle="dotted")
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)


# MAIN
# -----------------------------------------------------------------------------

input_file = sys.argv[1]
simplex_lookup = read_simplices(input_file)
max_bound, intervals = read_barcode(input_file, simplex_lookup)

docw = 418.25372 / 72
figw = docw

dim_count = len(intervals)
bar_count = 0
for bc in intervals:
    bar_count += len(bc[1])

subplotxmargin = docw * 0.02
subplotymargin = docw * 0.04

barh = 0.006 * docw
barskip = barh * 1.5
dimskip = barh * 4

h = bar_count * (barh + barskip)
h += (dim_count + 1) * dimskip
h -= dim_count * barskip

yticks = []
ylabels = []
y = 0
for bc in intervals:
    y += dimskip
    tl = len(bc[1]) * barh + (len(bc[1]) - 1) * barskip
    y += tl * 0.5
    yticks.append(y / h)
    y += tl * 0.5
    ylabels.append(r'$\bm{p=' + str(bc[0]) + "}$")

x = 2 * subplotymargin
y = subplotymargin
w = figw - 2 * subplotymargin - subplotxmargin
figh = h + 2 * subplotymargin

fig = mplfig.Figure(figsize=(figw, figh))
ax_dummy = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])

init_ax(ax_dummy, max_bound)
init_ax(ax, max_bound)
ax_dummy.xaxis.set_ticks_position("top")
ax.set_yticks(yticks, ylabels)
for i, l in enumerate(ax.get_yticklabels()):
    l.set_color(colors[i])
ax.yaxis.set_tick_params(length=0)


y = 0
for bc in intervals:
    dim = bc[0]
    bars = bc[1]
    y += dimskip
    for b in bars:
        r = 0.25
        essential = False
        if len(b) == 1:
            essential = True
            b.append(max_bound)
        rect = pat.Rectangle((b[0], y / h), b[1] - b[0] - r * 0.5, barh / h, linewidth=0, facecolor=colors[dim])
        print(w / max_bound)
        c1 = pat.Ellipse((b[0], (y + 0.5 * barh) / h), r, r * (w / 16), facecolor=colors[dim])
        c2 = pat.Ellipse((b[1], (y + 0.5 * barh) / h), r, r * (w / 16), fill=essential, facecolor=colors[dim], ec=colors[dim], clip_on=False, lw=1.5)
        ax.add_patch(rect)
        ax.add_patch(c1)
        ax.add_patch(c2)
        y += barh
        if b != bars[-1]:
            y += barskip


fig.savefig("test.pdf", format='pdf')
fig.savefig("../thesis/img/01-barcode-4pts-index-abshom.pdf", format='pdf')
