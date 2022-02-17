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

def comp_bar(a, b):
    if a[0] > b[0] or (a[0] == b[0] and a[1] > b[1]):
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
            #if sint[0] == sint[1]:
            #    continue
            iint = [float(sint[0])]
            bounds.append(float(sint[0]))
            if sint[1].strip() != "":
                iint.append(float(sint[1]))
                bounds.append(float(sint[1]))
                max_bound = max(max_bound, iint[1])
            max_bound = max(max_bound, iint[0])
            intervals[dim][1].append(iint)
            simps = toks[2].split("::")
            verts = [[int(v.strip()) for v in simp.split("'")] for simp in simps]
            reps[dim][1].append(verts)
    for dim, b in enumerate(intervals):
        for i, interval in enumerate(b[1]):
            if len(interval) == 1:
                intervals[dim][1][i].append(max_bound * 1.05)
        #intervals[dim][1] = sorted(intervals[dim][1], key=cmp_to_key(comp_bar))
    return (sorted(list(set(bounds))), intervals, reps)

def init_ax(ax):
    pass

def transform_point(p, sx, sy, tx, ty):
    return [p[0] * sx + tx, p[1] * sy + ty]

def plot_rep(ax, rep, col, sx, sy, tx, ty):
    if rep == []:
        return
    for rep_simplex in rep[1]:
        ps = [transform_point(points[v], sx, sy, tx, ty) for v in rep_simplex]
        c = col
        if rep[0] == 0:
            plot_point(ax, ps[0], c, 1)
        if rep[0] == 1:
            plot_line(ax, ps, c, 2)
        if rep[0] == 2:
            plot_trig(ax, ps, c, -3)
        if rep[0] ==3:
            plot_trig(ax, [p for p in [ps[2],ps[0],ps[3],ps[1]]], c, -3)

def plot_outline(ax, sx, sy, tx, ty, incl_interior):
    col = "lightgrey"
    for i, p in enumerate(points):
        plot_point(ax, transform_point(p, sx, sy, tx, ty), col, -1)
        for j, q in enumerate(points):
            if j >= i:
                continue
            if i == 3 and j == 2 and (not incl_interior):
                continue
            if i == 1 and j == 0 and (not incl_interior):
                continue
            ps = [transform_point(k, sx, sy, tx, ty) for k in [p,q]]
            plot_line(ax, ps, col, -2)
    if incl_interior:
        ps = [transform_point(points[k], sx, sy, tx, ty) for k in [0,1,2]]
        plot_trig(ax, ps, col, -2)
        ps = [transform_point(points[k], sx, sy, tx, ty) for k in [0,1,3]]
        plot_trig(ax, ps, col, -2)

def plot_point(ax, p, c, zorder):
    ax.plot([p[0]], [p[1]], "o", color=c, markersize=3.5, zorder=0)

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
points = read_dataset()
bounds, intervals, reps = read_barcode(input_file)

docw = 418.25372 / 72
rep_h = docw * 0.09
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
y = 2 * sy + rep_skip / h + (rep_h + rep_skip) / h

c = 0
for rs in reps:
    for rep in rs[1]:
        plot_outline(ax, sx, sy, x, y, (c < 4))
        plot_rep(ax, [rs[0], rep], colors[0], sx, sy, x, y)
        c+=1
        x += (rep_w + rep_skipx) / w
        if c == 4:
            y = y - (rep_h + rep_skip) / h
            x = 3 * sx + rep_skip / w

fig.savefig("test.pdf", format='pdf')
#fig.savefig("../thesis/img/01-cohom-reps.pdf", format='pdf')
