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

seed = 42

def gen_intro_ex_2(n):
    np.random.seed(seed)
    r = 1
    samples = []
    randfact = 0.15
    for t in range(n * 2):
        rx = np.random.rand() * randfact
        ry = np.random.rand() * randfact
        x = r * math.sin(t * 2 * PI / (n * 2)) + (rx - randfact * 0.5)
        y = r * math.cos(t * 2 * PI / (n * 2)) + (ry - randfact * 0.5)
        samples.append((x,y))
    for t in range(n):
        rx = np.random.rand() * randfact
        ry = np.random.rand() * randfact
        x = r * 0.5 * math.sin(t * 2 * PI / n) + (rx - randfact * 0.5)
        y = r * 0.5 * math.cos(t * 2 * PI / n) + (ry - randfact * 0.5)
        samples.append((x,y))
    return samples

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
            simps = toks[2].split("::")
            verts = [[int(v.strip()) for v in simp.split("'")] for simp in simps]
            reps[dim][1].append(verts)
    for dim, b in enumerate(intervals):
        for i, interval in enumerate(b[1]):
            if len(interval) == 1:
                intervals[dim][1][i].append(max_bound * 1.05)
        #intervals[dim][1] = sorted(intervals[dim][1], key=cmp_to_key(comp_bar))
    return (sorted(list(set(bounds))), intervals, reps)

def init_ax(ax, bounds, xticks=[], xticklabels=[]):
    ax.set_xlim(bounds)
    ax.set_ylim(0,1)
    if xticks != []:
        ax.set_xticks(xticks, xticklabels)
        ax.xaxis.grid(True, linestyle="dotted")
    ax.set_yticks([])
    #ax.spines["left"].set_visible(False)
    #ax.spines["right"].set_visible(False)


def init_plotax(ax):
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])


def plot_rep(ax, samples, rep, num):
    if num > 1:
        x, y = zip(*samples[30:])
        ax.scatter(x, y, s=1, zorder=2)
        x, y = zip(*samples[:30])
        ax.scatter(x, y, c="lightgrey", s=1, zorder=2)
    elif num == 0:
        x, y = zip(*samples[30:])
        ax.scatter(x, y, c="lightgrey", s=1, zorder=2)
        x, y = zip(*samples[:30])
        ax.scatter(x, y, c=colors[1], s=1, zorder=2)
    elif num == 1:
        x, y = zip(*[x for i, x in enumerate(samples[30:]) if i % 2 == 0])
        ax.scatter(x, y, c=colors[3], s=1, zorder=5)
        x, y = zip(*samples[:30])
        ax.scatter(x, y, c="lightgrey", s=1, zorder=5)
        x, y = zip(*[x for i, x in enumerate(samples[30:]) if i % 2 == 1])
        ax.scatter(x, y, c=colors[0], s=1, zorder=5)
    lines = []
    cs = []
    for i in range(0, len(samples)):
        for j in range(0, i):
            a = samples[i]
            b = samples[j]
            if [j,i] in rep:
                lines.append([a,b])
                if num == 0:
                    cs.append(colors[1])
                elif num ==1:
                    col1 = colors[0]
                    col2 = colors[3]
                    d = 1.5
                    l1 = Line2D([a[0],b[0]],[a[1],b[1]], lw=1, c=col1, linestyle="dotted")
                    l1.set_dashes([d, d])
                    l2 = Line2D([a[0],b[0]],[a[1],b[1]], lw=1, c=col2, linestyle="dotted")
                    l2.set_dashes([0, d, d, 0])
                    ax.add_artist(l1)
                    ax.add_artist(l2)
                else:
                    cs.append(colors[0])
    lc = mc.LineCollection(lines, colors=cs, linewidths=1, zorder=3)
    ax.add_collection(lc)


# MAIN
# -----------------------------------------------------------------------------

input_file = sys.argv[1]
bounds, intervals, reps = read_barcode(input_file)
data = gen_intro_ex_2(15)


print(reps[1])

docw = 418.25372 / 72
figw = docw

dim_count = 1
bar_count = 0
for bc in intervals:
    bar_count += len(bc[1])

subplotxmargin = docw * 0.02
marginw = subplotxmargin
marginh = docw * 0.04

repw = (docw - 7 * subplotxmargin) / 6
barh = 0.006 * docw
barskip = barh * 1.5
dimskip = barh * 4

bar_count = len(intervals[1][1])
h = (bar_count - 1) * (barh + barskip) + barskip
h_top = 2 * (barh + barskip) + barskip

w = figw - 2 * subplotxmargin
figh = h + 2 * repw + 2 * h_top + 7 * marginw
x = subplotxmargin
y = figh - subplotxmargin * 2 - h_top

fig = mplfig.Figure(figsize=(figw, figh))


ax = fig.add_axes([ x / figw, y / figh, w / figw, h_top / figh])
ax.xaxis.set_ticks_position("top")
ax.set_yticks([])
ax.set_xlim([0,1.8])
ax.tick_params(axis='both', which='major', pad=1)

hh = barskip
b = intervals[1][1][0]
rect = pat.Rectangle((b[0], hh / h_top), b[1] - b[0], barh / h_top, linewidth=0, facecolor=colors[1])
ax.add_patch(rect)

y = y - h_top - 0.5 * marginw

ax = fig.add_axes([ x / figw, y / figh, w / figw, h_top / figh])
ax.set_yticks([])
ax.set_xticks([])
ax.set_xlim([0,1.8])

b = intervals[1][1][1]
rect = pat.Rectangle((b[0], hh / h_top), b[1] - b[0], barh / h_top, linewidth=0, facecolor=colors[3])
ax.add_patch(rect)

y = y - h - 0.5 * marginw

ax_dummy = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
init_ax(ax_dummy, [0,1.8])
ax_dummy.tick_params(axis='both', which='major', pad=1)

#ax.set_yticks(yticks, ylabels)
#for i, l in enumerate(ax.get_yticklabels()):
#    l.set_color(colors[i])
#ax.yaxis.set_tick_params(length=0)

hh = barskip
for i, b in enumerate(intervals[1][1][1:]):
    rect = pat.Rectangle((b[0], hh / h), b[1] - b[0], barh / h, linewidth=0, facecolor=colors[0])
    ax_dummy.add_patch(rect)
    hh += barh
    if b != intervals[1][1][-1]:
        hh += barskip

y = y - h - marginw * 0.5

for i in range(12):
    ax = fig.add_axes([x / figw, y / figh, repw / figw, repw / figh])
    init_plotax(ax)
    plot_rep(ax, data, reps[1][1][i], i)
    x += repw + marginw
    if i == 5:
        x = marginw
        y = y - repw - marginh * 0.5




fig.savefig("test.pdf", format='pdf')
fig.savefig("../thesis/img/02-relrep-relcomplex.pdf", format='pdf')
