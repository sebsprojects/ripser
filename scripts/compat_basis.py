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
    "font.size": 10,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
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
    ax.plot([p[0]], [p[1]], "o", color=c, markersize=4, zorder=0)

def plot_line(ax, ps, c, zorder,lw=1.5):
    cs = [c, c]
    lc = mc.LineCollection([ps], colors=cs, zorder=zorder,linewidths=lw)
    ax.add_collection(lc)

def plot_trig(ax, ps, c, zorder, alpha=0.5):
    trig = plt.Polygon(ps, color=c, zorder=zorder, alpha=alpha, fill=True, facecolor=c, lw=0)
    ax.add_patch(trig)

# MAIN
# -----------------------------------------------------------------------------

input_file = sys.argv[1]
simplex_lookup = read_simplices(input_file)
points = read_dataset()

print(simplex_lookup)

z = [[0, [3]], [0, [2]],   [0, [1]],   [0, [0]]]
b = [[],       [0, [2,3]], [0, [1,2]], [0, [0,2]]]

z += [[1, [1]], [1, [1,3]], [1, [2]], [1, [2,3,1,4]], [1, [3,1,5]], [1, [2,1,0]]]
b += [[],       [],         [],       [1, [4,2,3,1]], [1, [5,4,2]], [1, [0,4,3]]]

z += [[2, [3]], [2, [3,2]], [2, [1]], [2, [1,2,3,0]]]
b += [[],       [],         [],       [2, [1,2,3,0]]]

z += [[3, [0]]]
b += [[]]



docw = 418.25372 / 72
rep_h = docw * 0.07
rep_skip = rep_h * 0.3
rep_skipx = rep_h * 1
dim_skip = rep_h * 1

w = 5 * rep_h * (3 / 2) + 4 * rep_skipx + 2 * rep_skip
h = 6 * rep_h + 2 * dim_skip + 5 * rep_skip + 0.6 * dim_skip + rep_h * 0.6

subplotymargin = docw * 0.01

x = subplotymargin
y = subplotymargin

figw = w + 2 * subplotymargin
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

grey = "dimgrey"
alpha_grey = 0.8

y = (h - rep_h) / h
rep_w = rep_h * (3 / 2)
x = 3 * sx + rep_skip / w
for i in range(3):
    for j in range(5):
        ind = i * 5 + j
        col = grey if ind != 0 and b[ind] == [] else colors[0]
        od2 = ind >= 10
        plot_rep(ax, z[ind], col, sx, sy, x, y - (rep_h + rep_skip) / h, od2)
        if ind == 10:
            ps = [transform_point(points[v], sx, sy, x, y - (rep_h + rep_skip) / h) for v in simplex_lookup[10][2]]
            plot_trig(ax, ps, grey, -3, alpha_grey)
            ps = [transform_point(points[v], sx, sy, x, y - (rep_h + rep_skip) / h) for v in simplex_lookup[11][2]]
            plot_trig(ax, ps, "lightgrey", -3, 0.5)
        if ind == 11:
            #ps = [transform_point(points[v], sx, sy, x, y - (rep_h + rep_skip) / h) for v in simplex_lookup[8][2]]
            #plot_line(ax, ps, "lightgrey", 3)
            ps = [transform_point(points[v], sx, sy, x, y - (rep_h + rep_skip) / h) for v in simplex_lookup[10][2]]
            plot_trig(ax, ps, grey, -3, alpha_grey)
            ps = [transform_point(points[v], sx, sy, x, y - (rep_h + rep_skip) / h) for v in simplex_lookup[11][2]]
            plot_trig(ax, ps, grey, -3, alpha_grey)
        if ind == 12:
            ps = [transform_point(points[v], sx, sy, x, y - (rep_h + rep_skip) / h) for v in simplex_lookup[13][2]]
            plot_trig(ax, ps, "lightgrey", -3, 0.5)
            ps = [transform_point(points[v], sx, sy, x, y - (rep_h + rep_skip) / h) for v in simplex_lookup[12][2]]
            plot_trig(ax, ps, grey, -3, alpha_grey)
        if ind == 13:
            ps = [transform_point(points[v], sx, sy, x, y - (rep_h + rep_skip) / h) for v in simplex_lookup[10][2]]
            plot_trig(ax, ps, colors[0], -3, 0.8)
            ps = [transform_point(points[v], sx, sy, x, y - (rep_h + rep_skip) / h) for v in simplex_lookup[11][2]]
            plot_trig(ax, ps, colors[0], -3, 0.8)
            ps = [transform_point(points[v], sx, sy, x, y) for v in simplex_lookup[10][2]]
            plot_trig(ax, ps, colors[1], -3, 0.7)
            ps = [transform_point(points[v], sx, sy, x, y) for v in simplex_lookup[11][2]]
            plot_trig(ax, ps, colors[1], -3, 0.7)
        if ind == 14:
            ps = [transform_point(points[v], sx, sy, x, y - (rep_h + rep_skip) / h) for v in simplex_lookup[10][2]]
            plot_trig(ax, ps, grey, -3, alpha_grey)
            ps = [transform_point(points[v], sx, sy, x, y - (rep_h + rep_skip) / h) for v in simplex_lookup[11][2]]
            plot_trig(ax, ps, grey, -3, alpha_grey)
        plot_rep(ax, b[ind], colors[1], sx, sy, x, y, od2)
        #ax.text(x, y + rep_h * 0.8 / h, r'$K_{' + str(ind + 1) + "}$", ha="center")
        ax.text(x, y + rep_h * 0.8 / h, r'$i=' + str(ind + 1) + "$", ha="center")
        x += (rep_w + rep_skipx) / w
    y -= (rep_h * 2 + rep_skip * 1) / h
    x = 3 * sx + rep_skip / w
    if i < 2:
        plot_line(ax, [[0,y], [w, y]], "black", 1, 0.75)
        y -= dim_skip * 1.3 / h


x += (rep_w + rep_skipx) / w
ps = [transform_point(points[v], sx * 0.25, sy * 0.25, x - 0.01, y + 0.01) for v in simplex_lookup[10][2]]
plot_trig(ax, ps, grey, -3, alpha_grey)
ps = [transform_point(points[v], sx * 0.29, sy * 0.29, x + 0.01, y + 0.008) for v in simplex_lookup[11][2]]
plot_trig(ax, ps, grey, -3, alpha_grey)

x += (rep_w + rep_skipx) * 2 / w
ps = [transform_point(points[v], sx * 0.25, sy * 0.25, x - 0.05, y + 0.01) for v in simplex_lookup[10][2]]
plot_trig(ax, ps, colors[0], -3, 0.8)
ps = [transform_point(points[v], sx * 0.29, sy * 0.29, x - 0.025, y + 0.008) for v in simplex_lookup[11][2]]
plot_trig(ax, ps, colors[0], -3, 0.8)
ps = [transform_point(points[v], sx * 0.29, sy * 0.29, x + 0.025, y + 0.008) for v in simplex_lookup[13][2]]
plot_trig(ax, ps, colors[0], -3, 0.8)
ps = [transform_point(points[v], sx * 0.29, sy * 0.29, x + 0.05, y + 0.01) for v in simplex_lookup[12][2]]
plot_trig(ax, ps, colors[0], -3, 0.8)

x += (rep_w + rep_skipx) * 1 / w
ps = [transform_point(points[v], sx * 0.25, sy * 0.25, x, y + 0.01) for v in simplex_lookup[10][2]]
plot_trig(ax, ps, grey, -3, alpha_grey)
ps = [transform_point(points[v], sx * 0.25, sy * 0.25, x - 0.001, y + 0.01) for v in simplex_lookup[11][2]]
plot_trig(ax, ps, grey, -3, alpha_grey)


fig.savefig("test.pdf", format='pdf')
fig.savefig("../thesis/img/01-compat-basis.pdf", format='pdf')
