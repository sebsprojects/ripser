import matplotlib
import matplotlib.figure as mplfig
import matplotlib.axes as mplax
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


# MATPLOTLIB SETUP
# -----------------------------------------------------------------------------

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
PI = math.pi
tex_fonts = {
    "text.usetex": True,
    "font.family": "serif",
#    "axes.labelsize": 9,
    "font.size": 9,
#    "legend.fontsize": 9,
#    "xtick.labelsize": 9,
#    "ytick.labelsize": 9,
    'text.latex.preamble': r'\usepackage{amsfonts}'
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

def read_simplices():
    f = open(input_file, "r")
    simplices = []
    diams = []
    for line in f:
        if line.startswith("#s"):
            vertices = [int(v) for v in line[3:].split("-")[1].strip().split("'")]
            diam = float(line[3:].split("-")[2].strip())
            if diam not in diams:
                diams.append(diam)
            simplices.append([diam, vertices])
    return diams, simplices

def get_limits(points, yfac=1.0):
    x_lim = [points[0][0], points[0][0]]
    y_lim = [points[0][1], points[0][1]]
    for p in points[1:]:
        x_lim[0] = min(p[0], x_lim[0])
        x_lim[1] = max(p[0], x_lim[1])
        y_lim[0] = min(p[1], y_lim[0])
        y_lim[1] = max(p[1], y_lim[1])
    x_len = x_lim[1] - x_lim[0]
    y_len = y_lim[1] - y_lim[0]
    x_mid = x_lim[0] + 0.5 * x_len
    y_mid = y_lim[0] + 0.5 * y_len
    offs = 0.65 * max(x_len, y_len)
    return [x_mid - offs, x_mid + offs], [y_mid - yfac * offs, y_mid + yfac * offs]

def init_ax(ax, xlim, ylim, s):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(s, pad=5)

def plot_filtration(axs, simplices, diams, points, title_list):
    xlim, ylim = get_limits(points)
    for i in range(len(diams)):
        init_ax(axs[i], xlim, ylim, title_list[i])
        for simp in simplices:
            if simp[0] > diams[i]:
                break
            dim = len(simp[1])
            diam = simp[0]
            c = colors[1] if diam == diams[i] else colors[0]
            if dim == 1:
                plot_point(axs[i], points[simp[1][0]], c)
            elif dim == 2:
                ps = [points[simp[1][k]] for k in range(2)]
                plot_line(axs[i], ps, c)


def plot_point(ax, p, c):
    ax.plot([p[0]], [p[1]], "o", color=c, zorder=3, markersize=5)

def plot_line(ax, ps, c):
    cs = [c, c]
    lc = mc.LineCollection([ps], colors=cs)
    ax.add_collection(lc)

def plot_trig(ax, ps, c):
    trig = plt.Polygon(ps, color=c, alpha=0.5)
    ax.add_patch(trig)

# MAIN
# -----------------------------------------------------------------------------

input_file = sys.argv[1]
points = read_dataset()
small_lims = get_limits(points)

print(points)

fpoints = []
for [x,y] in points:
    fpoints.append([x - 0, y])
    fpoints.append([x + 5,-y])

big_lims = get_limits(fpoints, 0.5)
print(big_lims)

docw = 418.25372 / 72

subw = docw * 0.14
marginw = docw * 0.02
marginh = docw * 0.043

figw = 6 * subw + 5 * marginw
figh = subw + marginh + marginw

fig = mplfig.Figure(figsize=(figw, figh))

print(figw / docw)

title_list = ["K", "L", "M"]

axs = []
y = figh - marginh - subw
x = marginw
ax = fig.add_axes([x / figw, y / figh, subw / figw, subw / figh])
init_ax(ax, small_lims[0], small_lims[1], r'$L \subseteq K$')
axs.append(ax)
x += subw + marginw
ax = fig.add_axes([x / figw, y / figh, subw * 2 / figw, subw / figh])
init_ax(ax, big_lims[0], big_lims[1], r'$N \subseteq M$')
axs.append(ax)
x += 2 * subw + marginw
ax = fig.add_axes([x / figw, y / figh, subw / figw, subw / figh])
init_ax(ax, small_lims[0], small_lims[1], r'$K \setminus L$')
axs.append(ax)
x += subw + marginw
ax = fig.add_axes([x / figw, y / figh, subw * 2 / figw, subw / figh])
init_ax(ax, big_lims[0], big_lims[1], r'$M \setminus N$')
axs.append(ax)

print(fpoints)

for p in points:
    if p[0] < 0:
        c = "lightgrey"
    else:
        c = colors[0]
    plot_point(axs[0], p, c)
for i in range(len(points)):
    for j in range(i):
        if points[i][0] < 0 or points[j][0] < 0:
            c = "lightgrey"
        else:
            c = colors[0]
        plot_line(axs[0], [points[i], points[j]], c)
plot_trig(axs[0], points[0:3], "lightgrey")
plot_trig(axs[0], [points[1], points[3], points[0]], "lightgrey")

for p in points:
    if p[0] > 0:
        c = "lightgrey"
    else:
        c = colors[0]
    plot_point(axs[2], p, c)
for i in range(len(points)):
    for j in range(i):
        if points[i][0] > 0 and points[j][0] > 0:
            c = "lightgrey"
        else:
            c = colors[0]
        plot_line(axs[2], [points[i], points[j]], c)
plot_trig(axs[2], points[0:3], colors[0])
plot_trig(axs[2], [points[1], points[3], points[0]], colors[0])


def norm2(p1, p2):
    return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)

for p in fpoints:
    if p[0] < -0.5:
        c = "lightgrey"
    else:
        c = colors[0]
    plot_point(axs[1], p, c)
for i in range(len(fpoints)):
    for j in range(i):
        if norm2(fpoints[i], fpoints[j]) > 6.5:
            continue
        if fpoints[i][0] < -0.5 or fpoints[j][0] < -0.5:
            c = "lightgrey"
        else:
            c = colors[0]
        plot_line(axs[1], [fpoints[i], fpoints[j]], c)
plot_trig(axs[1], [fpoints[0], fpoints[2], fpoints[4]], "lightgrey")
plot_trig(axs[1], [fpoints[0], fpoints[2], fpoints[6]], "lightgrey")
plot_trig(axs[1], [fpoints[1], fpoints[3], fpoints[5]], colors[0])
plot_trig(axs[1], [fpoints[1], fpoints[3], fpoints[7]], colors[0])

for p in fpoints:
    if p[0] > -1:
        c = "lightgrey"
    else:
        c = colors[0]
    plot_point(axs[3], p, c)
for i in range(len(fpoints)):
    for j in range(i):
        if norm2(fpoints[i], fpoints[j]) > 6.5:
            continue
        if fpoints[i][0] > -1 and fpoints[j][0] > -1:
            c = "lightgrey"
        else:
            c = colors[0]
        plot_line(axs[3], [fpoints[i], fpoints[j]], c)
plot_trig(axs[3], [fpoints[0], fpoints[2], fpoints[4]], colors[0])
plot_trig(axs[3], [fpoints[0], fpoints[2], fpoints[6]], colors[0])
plot_trig(axs[3], [fpoints[1], fpoints[3], fpoints[5]], "lightgrey")
plot_trig(axs[3], [fpoints[1], fpoints[3], fpoints[7]], "lightgrey")

fig.savefig("filt.pdf", format='pdf')
fig.savefig("../thesis/img/01-relative.pdf", format='pdf')
