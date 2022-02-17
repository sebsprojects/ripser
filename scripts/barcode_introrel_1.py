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
        intervals[dim][1] = sorted(intervals[dim][1], key=cmp_to_key(comp_bar))
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

def plot_rep(ax, samples, rep):
    x, y = zip(*samples[0:15])
    ax.scatter(x, y, s=1, zorder=2)
    x, y = zip(*samples[15:])
    ax.scatter(x, y, c="lightgrey", s=1, zorder=2)
    lines = []
    cs = []
    for i in range(0, len(samples)):
        for j in range(0, i):
            a = samples[i]
            b = samples[j]
            if [j,i] in rep:
                lines.append([a,b])
    lc = mc.LineCollection(lines, linewidths=1, zorder=3)
    ax.add_collection(lc)


# MAIN
# -----------------------------------------------------------------------------

input_file = sys.argv[1]
bounds, intervals, reps = read_barcode(input_file)

print(reps[1])

docw = 418.25372 / 72
figw = docw

dim_count = 1
bar_count = 0
for bc in intervals:
    bar_count += len(bc[1])

subplotxmargin = docw * 0.02

repw = (docw - 7 * subplotxmargin) / 6
barh = 0.006 * docw
barskip = barh * 1.5
dimskip = barh * 4

yticks = []
ylabels = []

#y = 0
#bc = intervals[1]
#y += dimskip
#tl = len(bc[1]) * barh + (len(bc[1]) - 1) * barskip
#y += tl * 0.5
#yticks.append(y / h)
#y += tl * 0.5
#ylabels.append(r'$\bm{p=' + str(bc[0]) + "}$")

x = subplotxmargin
y = subplotxmargin * 2
w = figw - 3 * subplotxmargin - repw
figh = repw + 4 * subplotxmargin
h = repw

fig = mplfig.Figure(figsize=(figw, figh))
ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
ax_dummy = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])

init_ax(ax, [0,0.9], [])
init_ax(ax_dummy, [0,0.9], intervals[1][1][0], [r'$\approx 0.32$', r'$\approx 0.85$'])
ax_dummy.xaxis.set_ticks_position("top")

#ax.set_yticks(yticks, ylabels)
#for i, l in enumerate(ax.get_yticklabels()):
#    l.set_color(colors[i])
#ax.yaxis.set_tick_params(length=0)

y = 0
bc = intervals[1]
dim = bc[0]
bars = bc[1]
y += dimskip
for b in bars:
    rect = pat.Rectangle((b[0], y / h), b[1] - b[0], barh / h, linewidth=0, facecolor=colors[0])
ax_dummy.add_patch(rect)
y += barh
if b != bars[-1]:
    y += barskip

ax = fig.add_axes([ (w + 2*subplotxmargin) / figw, (subplotxmargin * 2) / figh, repw / figw, repw / figh])
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlim(-1.2, 1.2)
ax.set_ylim(-1.2, 1.2)
ax.set_aspect('equal')

plot_rep(ax, gen_intro_ex_2(15)[30:] + gen_intro_ex_2(15)[0:30], reps[1][1][0])

fig.savefig("test.pdf", format='pdf')
fig.savefig("../thesis/img/02-relrep-subcomplex.pdf", format='pdf')
