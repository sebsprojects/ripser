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
    "axes.labelsize": 9,
    "font.size": 9,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    'text.latex.preamble': r'\usepackage{amsfonts}'
}
plt.rcParams.update(tex_fonts)


# FUNCTION DEFS
# -----------------------------------------------------------------------------

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

def read_simplices():
    f = open(input_file, "r")
    simplices = []
    for line in f:
        if line.startswith("#s"):
            vertices = [int(v) for v in line[3:].split("-")[1].strip().split("'")]
            simplices.append(vertices)
    return simplices

def init_ax(ax, i):
    ax.set_aspect("equal")
    ax.set_xlim([-1.2, 1.2])
    ax.set_ylim([-1.2, 1.2])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(r'$K_{' + str(i) + "}$")


def plot_filtration(axs, simplices, points):
    for i in range(len(simplices)):
        init_ax(axs[i], i + 1)
        for j in range(i + 1):
            dim = len(simplices[j])
            c = colors[1] if i == j else colors[0]
            if dim == 1:
                plot_point(axs[i], points[simplices[j][0]], c)
            elif dim == 2:
                ps = [points[simplices[j][k]] for k in range(2)]
                plot_line(axs[i], ps, c)
            else:
                for m in range(dim):
                    a = m
                    b = m + 3
                    ps = [points[simplices[j][k % dim]] for k in range(a,b)]
                    c = colors[1] if i == j else "lightskyblue"
                    plot_trig(axs[i], ps, c)


def plot_point(ax, p, c):
    ax.plot([p[0]], [p[1]], "o", color=c)

def plot_line(ax, ps, c):
    cs = [c, c]
    lc = mc.LineCollection([ps], colors=cs)
    ax.add_collection(lc)

def plot_trig(ax, ps, c):
    trig = plt.Polygon(ps, color=c)
    ax.add_patch(trig)

# MAIN
# -----------------------------------------------------------------------------


read_dataset()
simplices = read_simplices()
points = read_dataset()

print(simplices)

num_w = 7
num_h = 5
fig, axs = plt.subplots(num_h, num_w, figsize=(figw, figw * 0.9))
plot_filtration(axs.flat, simplices, points)


plt.savefig("../thesis/img/test2.pdf", format='pdf', bbox_inches="tight", pad_inches=0.05)
#plt.show()
