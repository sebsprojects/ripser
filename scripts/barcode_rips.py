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
            toks = [t.strip() for t in line[3:].split(";")]
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
                intervals[dim][1][i].append(max_bound * 1.05)
        intervals[dim][1] = sorted(intervals[dim][1], key=cmp_to_key(comp_bar))
    return (sorted(list(set(bounds))), intervals)

def init_ax(ax, bounds, xlabels=[]):
    ax.set_xlim(0, bounds[-1] * 1.05)
    ax.set_ylim(0,1)
    xlabels = xlabels if xlabels != [] else [str(round(b, 1)).format(1) for b in bounds]
    ax.set_xticks(bounds, xlabels)
    ax.xaxis.grid(True, linestyle="dotted")
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)


# MAIN
# -----------------------------------------------------------------------------

input_file = sys.argv[1]
bounds, intervals = read_barcode(input_file)

docw = 418.25372 / 72
figw = docw

dim_count = len([bc for bc in intervals if bc[1] != []])
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
    if bc[1] == []:
        continue
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

xlabels = [ r'$0$', "", r'$\sqrt{17}$', r'6', r'$\sqrt{41}$']
xlabelss = [ r'$0$', r'$4$', "", r'6', r'$\sqrt{41}$']
init_ax(ax_dummy, bounds, xlabels)
init_ax(ax, bounds, xlabelss)
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
        rect = pat.Rectangle((b[0], y / h), b[1] - b[0], barh / h, linewidth=0, facecolor=colors[dim])
        ax.add_patch(rect)
        y += barh
        if b != bars[-1]:
            y += barskip


fig.savefig("test.pdf", format='pdf')
fig.savefig("../thesis/img/01-barcode-4pts-rips-abshom.pdf", format='pdf')
