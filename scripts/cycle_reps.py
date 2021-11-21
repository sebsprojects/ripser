import matplotlib
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from matplotlib import colors as mcol
from matplotlib.lines import Line2D
import numpy as np
import math
import sys
from functools import cmp_to_key

if len(sys.argv) < 2:
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

def read_intervals(input_file, target_dim):
    f = open(input_file, "r");
    intervals = []
    max_bound = 0
    offs = 1
    for line in f:
        if line.startswith("#b"):
            dim = int(line[2])
            if dim != target_dim:
                continue
            toks = [t.strip() for t in line[3:].split("-")]
            intv = toks[0].strip()
            simps = toks[2].split("::")
            verts = [[int(v.strip()) for v in simp.split("'")] for simp in simps]
            intervals.append([intv, verts])
    return intervals

def read_vertices(input_file):
    f = open(input_file, "r")
    vertex_lookup = []
    rel = range(8)
    for line in f:
        if line.startswith("#s0"):
            index = int(line[3:].split("-")[0].strip())
            vertex_lookup.append(index)
    return vertex_lookup

def read_dataset():
    f = open(input_file, "r")
    points = []
    for line in f:
        if line.find("dataset at file path") != -1:
            path = line.split(":")[1].strip()
            d = open(path, "r")
            for l in d:
                if l.strip() == "":
                    continue
                else:
                    points.append([float(p.strip()) for p in l.split(",")])
    return points

def plot_rep(ax, interval, vertices, points):
    simps = interval[1]
    for simp in simps:
        ps = [points[vertices[i]] for i in simp]
        print(ps)
        if len(ps) == 1:
            plot_point(ax, ps, colors[1])
        if len(ps) == 2:
            plot_line(ax, ps, colors[1])
        if len(ps) == 3:
            plot_trig(ax, ps, colors[1])

def plot_point(ax, p, c):
    ax.plot([p[0]], [p[1]], "o", color=c)

def plot_line(ax, ps, c):
    cs = [c, c]
    lc = mc.LineCollection([ps], colors=cs)
    ax.add_collection(lc)

def plot_trig(ax, ps, c):
    trig = plt.Polygon(ps, color=c)
    ax.add_patch(trig)

def init_ax(ax, xlim, ylim, t):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(r'$I=' + t + "$")

# MAIN
# -----------------------------------------------------------------------------

fig, ax = plt.subplots(1, 1, figsize=(figw, figw))
input_file = sys.argv[1]

intervals = read_intervals(input_file, 1)
vertices = read_vertices(input_file)
points = read_dataset()
i = 1
init_ax(ax, [-4,4], [-4,4], intervals[i][0])
plot_rep(ax, intervals[i], vertices, points)

plt.savefig("test.pdf", format='pdf', bbox_inches="tight", pad_inches=0.05)
#plt.show()
