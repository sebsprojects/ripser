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

from matplotlib.ticker import AutoMinorLocator

# MATPLOTLIB SETUP
# -----------------------------------------------------------------------------

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
PI = math.pi
tex_fonts = {
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 8,
    "axes.titlesize": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    'text.latex.preamble': r'\usepackage{amsfonts}\usepackage{bm}\usepackage{amsmath}'
}
plt.rcParams.update(tex_fonts)

seed = 42


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
    reps = []
    max_bound = 0
    bounds = []
    for line in f:
        if line.startswith("#b"):
            dim = int(line[2])
            if dim + 1 > len(intervals):
                intervals.append([dim,[]])
                reps.append([dim,[]])
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
            #simps = toks[2].split("::")
            #verts = [[int(v.strip()) for v in simp.split("'")] for simp in simps]
            #reps[dim][1].append(verts)
    for dim, b in enumerate(intervals):
        for i, interval in enumerate(b[1]):
            if len(interval) == 1:
                intervals[dim][1][i].append(max_bound * 1.05)
        #intervals[dim][1] = sorted(intervals[dim][1], key=cmp_to_key(comp_bar))
    return (sorted(list(set(bounds))), intervals)

def init_ax(ax):
    ax.set_yticks([])
    ax.set_xlim([0.15,0.65])
    ax.set_yticklabels([])
    ax.xaxis.grid(True, linestyle="dotted", which="both", alpha=0.5)
    #ax.spines["left"].set_visible(False)
    #ax.spines["right"].set_visible(False)
    ax.tick_params(axis='both', which='major', pad=1)
    ax.minorticks_on()
    minor_locator = AutoMinorLocator(2)
    ax.xaxis.set_minor_locator(minor_locator)

def plot_barcode(ax, intervals, intervals_comp):
    y = barskip
    dim = 1
    t = 0.005
    for b in intervals:
        col = colors[0]
        for bcomp in intervals_comp:
            if b == bcomp:
                col = colors[1]
                break
            #if b[0] == bcomp[0]:
            #    col = "#FFBE7D"
            #if abs(b[0] - bcomp[0]) < t and abs(b[1] - bcomp[1]) < t:
            #    col = colors[2]
            #elif b[1] == bcomp[1]:
            #    col = colors[3]
        rect = pat.Rectangle((b[0], y / h), b[1] - b[0], barh / h, linewidth=0, facecolor=col)
        ax.add_patch(rect)
        y += barh
        y += barskip


# MAIN
# -----------------------------------------------------------------------------

datasetname = "random100"
rel = 69
r = (rel + 1) / 100

#input_file_1 = "./output/barcode_comp/random_point_cloud_100_4__d2tMr1fu00.txt"
#input_file_2 = "./output/barcode_comp/random_point_cloud_100_4__d2tMr1fu00_rel0-" + str(rel) + ".txt"
bounds_1, intervals_1 = read_barcode(input_file_1)
bounds_2, intervals_2 = read_barcode(input_file_2)

docw = 418.25372 / 72
figw = docw

dim_count = 1
bar_count = 0
#for bc in intervals:
#    bar_count += len(bc[1])
bar_count = max(len(intervals_1[1][1]), len(intervals_2[1][1]))

margin = 0.05 * docw

barh = 0.0025 * docw
barskip = barh * 1.5
dimskip = barh * 4

h = bar_count * (barh + barskip) + barskip
#h += (dim_count + 1) * dimskip
#h -= dim_count * barskip

w = figw * 0.5 - margin
figh = h + 2 * margin
x = margin * 0.5
y = margin * 0.6

fig = mplfig.Figure(figsize=(figw, figh))
ax_dummy = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
ax_dummy.xaxis.set_ticks_position("top")
init_ax(ax)
init_ax(ax_dummy)

plot_barcode(ax, intervals_1[1][1], intervals_2[1][1])

x = margin * 1.5 + w
y = margin * 0.6
ax_dummy = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
ax_dummy.xaxis.set_ticks_position("top")
init_ax(ax)
init_ax(ax_dummy)

plot_barcode(ax, intervals_2[1][1], intervals_1[1][1])

ax = fig.add_axes([ margin / figw, (figh - 1.5 * margin) / figh, (figw - 2 * margin) / figw, margin / figh])
ax.text(0.5, 1, r'$\underline{~\text{Setup: ' + str(datasetname) + r'},~p=1,~r=0~\text{vs.}~r='+str(r) + r'~}$', ha="center")
ax.axis("off")

fig.savefig("test.pdf", format='pdf')
fig.savefig("../thesis/img/02-barcodecomp-" + datasetname + "-rel" + str(rel) + ".pdf", format='pdf')
