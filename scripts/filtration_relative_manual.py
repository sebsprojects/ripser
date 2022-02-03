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
    for line in f:
        if line.startswith("#s"):
            dim = int(line[2])
            index = int(line[3:].split("-")[1].strip())
            vertices = [int(v) for v in line[3:].split("-")[2].strip().split("'")]
            simplices.append([dim, index, vertices])
    return simplices

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

def init_ax(ax, xlim, ylim, i):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(r'$K_{' + str(i) + "}$", pad=5)

def plot_point(ax, p, c):
    ax.plot([p[0]], [p[1]], "o", color=c, zorder=3, markersize=5)

def plot_line(ax, ps, c):
    cs = [c, c]
    lc = mc.LineCollection([ps], colors=cs)
    ax.add_collection(lc)

def plot_trig(ax, ps, c):
    if c == colors[1]:
        alpha=0.5
    else:
        alpha=0.4
    trig = plt.Polygon(ps, color=c, alpha=alpha)
    ax.add_patch(trig)

# MAIN
# -----------------------------------------------------------------------------

input_file = sys.argv[1]
simplices = read_simplices()
points = read_dataset()

docw = 418.25372 / 72
subw = docw * 0.14
marginw = docw * 0.02
marginh = docw * 0.043

num_w = 6
num_h = 1

figw = (subw + marginw) * num_w - subw * 0.5 + marginw
figh = (subw + marginh) * num_h + marginw
fig = mplfig.Figure(figsize=(figw, figh))
print(figw / docw)
xlim, ylim = get_limits(points)
axs = []
for j in range(num_h):
    y = figh - marginh - subw - j * (subw + marginh)
    for i in range(num_w):
        if j * num_w +i > 15:
            break
        x = marginw + i * (subw + marginw)
        ax = fig.add_axes([x / figw, y / figh, subw / figw, subw / figh])
        init_ax(ax, xlim, ylim, j * (num_w) + i)
        axs.append(ax)

print(xlim, ylim)

axs[5].set_title("")
axs[5].set_axis_off()
axs[5].text(-2, 0.25, r"\Large\textbf{.\,.\,.}", horizontalalignment="center", verticalalignment="center")


gray = "lightgrey"
hcol = colors[1]

plot_point(axs[1], points[simplices[0][2][0]], gray)
plot_point(axs[1], points[simplices[1][2][0]], hcol) # 1

plot_point(axs[2], points[simplices[0][2][0]], gray)
plot_point(axs[2], points[simplices[1][2][0]], colors[0])
plot_point(axs[2], points[simplices[2][2][0]], hcol) # 1

plot_point(axs[3], points[simplices[0][2][0]], gray)
plot_point(axs[3], points[simplices[1][2][0]], colors[0])
plot_point(axs[3], points[simplices[2][2][0]], colors[0])
plot_point(axs[3], points[simplices[3][2][0]], gray)
plot_line(axs[3], [points[simplices[4][2][k]] for k in range(2)], hcol) # 1

plot_point(axs[4], points[simplices[0][2][0]], gray)
plot_point(axs[4], points[simplices[1][2][0]], colors[0])
plot_point(axs[4], points[simplices[2][2][0]], colors[0])
plot_point(axs[4], points[simplices[3][2][0]], gray)
plot_line(axs[4], [points[simplices[4][2][k]] for k in range(2)], colors[0])
plot_line(axs[4], [points[simplices[5][2][k]] for k in range(2)], gray)
plot_line(axs[4], [points[simplices[6][2][k]] for k in range(2)], hcol) # 1

plot_point(axs[4], points[simplices[0][2][0]], gray)
plot_point(axs[4], points[simplices[1][2][0]], colors[0])
plot_point(axs[4], points[simplices[2][2][0]], colors[0])
plot_point(axs[4], points[simplices[3][2][0]], gray)
plot_line(axs[4], [points[simplices[4][2][k]] for k in range(2)], colors[0])
plot_line(axs[4], [points[simplices[5][2][k]] for k in range(2)], gray)
plot_line(axs[4], [points[simplices[6][2][k]] for k in range(2)], hcol) # 1

fig.savefig("filt.pdf", format='pdf')
fig.savefig("../thesis/img/01-filt-relative.pdf", format='pdf')
