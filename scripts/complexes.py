import matplotlib
import matplotlib.figure as mplfig
import matplotlib.axes as mplax
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from matplotlib import colors as mcol
from matplotlib.lines import Line2D
import mpl_toolkits.mplot3d as a3
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
    d = open(input_file, "r")
    points = []
    for l in d:
        if l.strip() == "":
            continue
        else:
            points.append([float(p.strip()) for p in l.split(",")])
    return points

def get_limits(points):
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
    return [x_mid - offs, x_mid + offs], [y_mid - offs, y_mid + offs]

def init_ax(ax, xlim, ylim, s):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(r'' + s, pad=5)

def plot_point(ax, p, c):
    ax.plot([p[0]], [p[1]], "o", color=c, zorder=3, markersize=5)

def plot_line(ax, ps, c):
    cs = [c, c]
    lc = mc.LineCollection([ps], colors=cs)
    ax.add_collection(lc)

def plot_trig(ax, ps, c):
    trig = plt.Polygon(ps, color=c, alpha=0.4)
    ax.add_patch(trig)

# MAIN
# -----------------------------------------------------------------------------

input_file = sys.argv[1]
points = read_dataset()
xlim, ylim = get_limits(points)

docw = 418.25372 / 72
subw = docw * 0.14
marginw = docw * 0.02
marginh = docw * 0.043

num_w = 6
num_h = 1

figw = (subw + marginw) * num_w + marginw
figh = (subw + marginh) * num_h + marginw
fig = mplfig.Figure(figsize=(figw, figh))
print(figw / docw)
title_list = ["$\{~\}$", "$X$", "$L$", "$K$", "$\\Delta X$", "$\\widehat{\\Delta X}$"]
axs = []
for j in range(num_h):
    y = figh - marginh - subw - j * (subw + marginh)
    for i in range(num_w):
        x = marginw + i * (subw + marginw)
        ax = fig.add_axes([x / figw, y / figh, subw / figw, subw / figh])
        init_ax(ax, xlim, ylim, title_list[i])
        axs.append(ax)

for p in points:
    for i in range(1,num_w - 1):
        plot_point(axs[i], p, colors[0])

plot_line(axs[2], [points[0], points[2]], colors[0])
plot_line(axs[2], [points[1], points[2]], colors[0])
plot_line(axs[2], [points[0], points[1]], colors[0])
plot_line(axs[2], [points[3], points[0]], colors[0])


plot_line(axs[3], [points[0], points[2]], colors[0])
plot_line(axs[3], [points[1], points[2]], colors[0])
plot_line(axs[3], [points[0], points[1]], colors[0])
plot_line(axs[3], [points[3], points[0]], colors[0])
plot_line(axs[3], [points[1], points[3]], colors[0])
plot_trig(axs[3], [points[0], points[2], points[1]], colors[0])

plot_line(axs[4], [points[0], points[2]], colors[0])
plot_line(axs[4], [points[1], points[2]], colors[0])
plot_line(axs[4], [points[0], points[1]], colors[0])
plot_line(axs[4], [points[3], points[0]], colors[0])
plot_line(axs[4], [points[1], points[3]], colors[0])
plot_line(axs[4], [points[2], points[3]], colors[0])
plot_trig(axs[4], [points[0], points[2], points[1]], colors[0])
plot_trig(axs[4], [points[1], points[3], points[0]], colors[0])

# 3D PLOT
ax3 = fig.add_axes([x / figw + 0.01, y / figh - 0.01, subw / figw - 0.02, subw / figh - 0.02], projection="3d")
ax3.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax3.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax3.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax3.set_axis_off()
xyz = np.array([[ 0, 1, 0], [ 1, 1, 0], [ 0, 1, 1], [0, 0, 0]])
tri = np.array([[ 0, 1, 2], [ 0, 1, 3], [ 0, 2, 3], [1, 2, 3]])
vts = xyz[tri, : ]
tri = a3.art3d.Poly3DCollection(vts, linewidths=1, edgecolors='tab:blue')
tri.set_alpha(0.25)
ax3.scatter3D(xyz[:,0], xyz[:,1], xyz[:,2], depthshade=False, s=15)
ax3.add_collection3d(tri)


fig.savefig("complexes.pdf", format='pdf')
fig.savefig("../thesis/img/00-complexes.pdf", format='pdf')
