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

def read_dataset():
    f = open(input_file, "r")
    points = []
    for line in f:
        if line.find("input path") != -1:
            path = line.split(":")[1].strip()
            d = open(path, "r")
            for l in d:
                if l.strip() == "":
                    continue
                else:
                    points.append([float(p.strip()) for p in l.split(",")])
    return points


def read_simplices(input_file):
    f = open(input_file, "r")
    simplices = []
    for line in f:
        if line.startswith("#s"):
            dim = int(line[2])
            index = int(line[3:].split("-")[0].strip())
            vertices = [int(v) for v in line[3:].split("-")[1].strip().split("'")]
            simplices.append([dim, index, vertices])
    return simplices

def init_ax(ax):
    pass

def transform_point(p, sx, sy, tx, ty):
    return [p[0] * sx + tx, p[1] * sy + ty]

def plot_rep(ax, rep, col, sx, sy, tx, ty,omit_dim2=False):
    if rep == []:
        return
    for rep_simplex in rep[1]:
        for s in simplex_lookup:
            if s[0] == 0:
                ps = [transform_point(points[v], sx, sy, tx, ty) for v in s[2]]
                c = col if s[0] == rep[0] and s[2][0] in rep[1] else "lightgrey"
                zorder = 3 if rep[0] == 0 else 0
                plot_point(ax, ps[0], c, zorder)
            if s[0] == 1:
                ps = [transform_point(points[v], sx, sy, tx, ty) for v in s[2]]
                if rep[0] == 1 and s[1] in rep[1]:
                    c = col
                    zorder = 2
                else:
                    c = "lightgrey"
                    zorder = -1
                # Draw only the outer outline for dim=2 and dim=3
                if rep[0] >= 2:
                    if s[1] == 5 or s[1] == 0:
                        zorder = -2
                    else:
                        zorder = 3
                plot_line(ax, ps, c, zorder)
    # MANUAL DIM=2
    if not omit_dim2:
        tps = [transform_point(points[v], sx, sy, tx, ty) for v in [1,2,3]]
        plot_trig(ax, tps, "lightgrey", -3, 0.5)
        tps = [transform_point(points[v], sx, sy, tx, ty) for v in [0,2,3]]
        plot_trig(ax, tps, "lightgrey", -3, 0.5)

def plot_point(ax, p, c, zorder):
    ax.plot([p[0]], [p[1]], "o", color=c, markersize=3, zorder=0)

def plot_line(ax, ps, c, zorder):
    cs = [c, c]
    lc = mc.LineCollection([ps], colors=cs, zorder=zorder)
    ax.add_collection(lc)

def plot_trig(ax, ps, c, zorder, alpha=0.5):
    trig = plt.Polygon(ps, zorder=zorder, alpha=alpha, fill=True, facecolor=c, lw=0)
    ax.add_patch(trig)

# MAIN
# -----------------------------------------------------------------------------

input_file = sys.argv[1]
simplex_lookup = read_simplices(input_file)
points = read_dataset()

z0 = [[0, [0,1,2,3]], [0, [2,0]], [0, [1]], [0, [0]]]
z12 = [[1, [4,5]], [1, [5]], [1, [0]], [2, [2]]]


reps = [z0, z12]

docw = 418.25372 / 72
rep_h = docw * 0.07
rep_w = rep_h * (3 / 2)
rep_skip = rep_h * 0.5
rep_skipx = rep_h
dim_skip = rep_h * 1

w = 4 * rep_w + 3 * rep_skipx + 2 * rep_skip
h = 2 * rep_h * 1 + dim_skip + 1 * rep_skip


subplotymargin = docw * 0.00
x = subplotymargin
y = subplotymargin
figw = w + 2 * subplotymargin
figh = h + 2 * subplotymargin

print(figw / docw)

fig = mplfig.Figure(figsize=(figw, figh))
ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
ax.set_xlim([0,1])
ax.set_ylim([0,1])
ax.set_xticks([])
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
init_ax(ax)

#sx = (1 / w) * (h / 4)
#sy = (1 / h) * (h / 4)

sy = (rep_h / 4) / h
sx = sy * (h / w)

x = 3 * sx + rep_skip / w
y = 2 * sy + rep_skip / h

for rs in reps:
    for rep in rs:
        omit = (rs == reps[-1]) and (rep == rs[-1])
        plot_rep(ax, rep, colors[0], sx, sy, x, y, omit)
        if rep != rs[-1]:
            x += (rep_w + rep_skipx) / w
    if rs != reps[-1]:
        y += (rep_h + rep_skip) / h
        x = 3 * sx + rep_skip / w

rep = reps[1][-1]
for s in simplex_lookup:
    if s[0] == rep[0] and s[1] == rep[1][0]:
        tps = [transform_point(points[v], sx, sy, x, y) for v in s[2]]
        plot_trig(ax, tps, colors[0], -3, 0.8)
        tps = [transform_point(points[v], sx, sy, x, y) for v in [3,1,2]]
        plot_trig(ax, tps, "lightgrey", 0)

fig.savefig("test.pdf", format='pdf')
fig.savefig("../thesis/img/01-cohom-reps.pdf", format='pdf')
