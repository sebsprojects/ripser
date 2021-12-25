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

def plot_rep(ax, rep, sx, sy, tx, ty):
    for rep_simplex in rep[1]:
        for s in simplex_lookup:
            if s[0] == 0:
                ps = [transform_point(points[v], sx, sy, tx, ty) for v in s[2]]
                c = colors[1] if s[0] == rep[0] and s[2][0] in rep[1] else "lightgrey"
                zorder = 3 if rep[0] == 0 else 0
                plot_point(ax, ps[0], c, zorder)
            if s[0] == 1:
                ps = [transform_point(points[v], sx, sy, tx, ty) for v in s[2]]
                if rep[0] == 1 and s[1] in rep[1]:
                    c = colors[1]
                    zorder = 2
                else:
                    c = "lightgrey"
                    zorder = -1
                plot_line(ax, ps, c, zorder)
            if s[0] == 2:
                ps = [transform_point(points[v], sx, sy, tx, ty) for v in s[2]]
                if rep[0] == 2 and s[1] in rep[1]:
                    c = colors[1]
                    zorder = -2
                else:
                    c = "lightgrey"
                    zorder = -3
                plot_trig(ax, ps, c, zorder)


def plot_point(ax, p, c, zorder):
    ax.plot([p[0]], [p[1]], "o", color=c, markersize=3, zorder=0)

def plot_line(ax, ps, c, zorder):
    cs = [c, c]
    lc = mc.LineCollection([ps], colors=cs, zorder=zorder)
    ax.add_collection(lc)

def plot_trig(ax, ps, c, zorder):
    trig = plt.Polygon(ps, color=c, zorder=zorder)
    ax.add_patch(trig)

# MAIN
# -----------------------------------------------------------------------------

input_file = sys.argv[1]
simplex_lookup = read_simplices(input_file)
points = read_dataset()

print(simplex_lookup)

z0 = [[0, [3]], [0, [2]], [0, [1]], [0, [0]]]
z1 = [[1, [4,2,3,1]], [1, [5,3,1]], [1, [0,2,1]]]
z2 = [[2, [0,1,2,3]]]

b0 = [[0, []], [0, [2,3]], [0, [1,2]], [0, [0, 2]]]
b1 = [[1, [4,2,3,1]], [1, [5,4,2]], [1, [0,4,3]]]
b2 = [[2, [0,1,2,3]]]


reps = [[z0, b0], [z1, b1], [z2, b2]]

docw = 418.25372 / 72
figw = docw
h = docw

rep_h = docw * 0.05
rep_skip = rep_h * 0.75
dim_skip = rep_h * 1

subplotxmargin = docw * 0.02
subplotymargin = docw * 0.04
x = 2 * subplotymargin
y = subplotymargin
w = figw - 2 * subplotymargin - subplotxmargin
figh = h + 2 * subplotymargin

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

y = 2 * sy + rep_skip / h
rep_w = rep_h * (3 / 2)
x = 3 * sx + rep_skip / w

for [zs, bs] in reps:
    for rep in zs:
        plot_rep(ax, rep, sx, sy, x, y)
        x += (rep_w + rep_skip) / w
    y += (rep_h + rep_skip) / h
    x = 3 * sx + rep_skip / w
    for rep in bs:
        plot_rep(ax, rep, sx, sy, x, y)
        x += (rep_w + rep_skip) / w
    y += (rep_h + dim_skip) / h
    x = 3 * sx + rep_skip / w

fig.savefig("test.pdf", format='pdf')
#fig.savefig("../thesis/img/01-barcode-4pts-index-abshom.pdf", format='pdf')
