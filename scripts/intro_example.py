import matplotlib
import matplotlib.figure as mplfig
import matplotlib.axes as mplax
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from matplotlib import collections as mc
from matplotlib import colors as mcol
from matplotlib.lines import Line2D
import numpy as np
import math
import sys
from functools import cmp_to_key

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
PI = math.pi
tex_fonts = {
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 11,
    "font.size": 9,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    'text.latex.preamble': r'\usepackage{amsfonts}'
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

def write_dataset(samples, filename):
    f = open(filename + ".pc", "w")
    for p in samples:
        f.write(str(p[0]) + "," + str(p[1]) + "\n")
    f.close()

    f = open(filename + ".dm", "w")
    for i in range(1, len(samples)):
        for j in range(len(samples)):
            if j < i:
                a = samples
                dist = math.sqrt(math.pow(a[i][0] - a[j][0], 2) +
                                 math.pow(a[i][1] - a[j][1], 2))
                f.write(str(dist))
                if j != i - 1:
                    f.write(",")
        if i != len(samples) - 1:
            f.write("\n")
    f.close()

def diff2(a, b):
    return math.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)

def plot_edges(ax, samples, thresh):
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(r'$t=' + str(thresh) + "$", pad=3)
    x, y = zip(*samples)
    ax.scatter(x, y, s=7, zorder=2)
    lines = []
    cs = []
    for i in range(0, len(samples)):
        for j in range(0, i):
            a = samples[i]
            b = samples[j]
            dist_ab = diff2(a,b)
            if dist_ab <= thresh:
                if True:
                    cs.append(colors[0])
                else:
                    cs.append(colors[0])
                lines.append([a, b])
            for k in range(0,j):
                c = samples[k]
                dist_abc = max(dist_ab, diff2(a,c))
                dist = max(dist_abc, diff2(b,c))
                if dist <= thresh:
                    trig = plt.Polygon([a,b,c], alpha=0.7, zorder=0)
                    ax.add_patch(trig)
    lc = mc.LineCollection(lines, colors=cs, linewidths=1, zorder=1)
    ax.add_collection(lc)

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

def init_barax(ax):
    ax.set_yticks([])
    ax.set_xlim(0, 0.9)
    #ax.spines["left"].set_visible(False)
    #ax.spines["right"].set_visible(False)

input_file = sys.argv[1]
bounds, intervals, reps = read_barcode(input_file)

print(reps[1])

docw = 418.25372 / 72
subw = docw * 0.225
marginw = docw * 0.02
marginh = docw * 0.043
nx = 4
ny = 2

plot_skip = docw * 0.05

barh = 0.006 * docw
barskip = barh * 1.5
bar_count = len(intervals[1][1])
h = bar_count * (barh + barskip) + barskip

figw = nx * subw + (nx+1) * marginw
figh = ny * subw + (ny+1) * marginh + h + plot_skip

fig = mplfig.Figure(figsize=(figw, figh))

print(figw / docw)

data = gen_intro_ex_2(15)

axs = []
y = figh - marginh - subw
x = marginw

threshs = [0, 0.319, 0.369, 0.4149, 0.501, 0.545, 0.592, 0.851]
reps_to_plot = [[], [0], [0, 1], [], [], [], [], []]

for i in range(0, nx * ny):
    ax = fig.add_axes([x / figw, y / figh, subw / figw, subw / figh])
#   init_ax(ax, small_lims[0], small_lims[1], str(n))
    axs.append(ax)
    plot_edges(ax, data, threshs[i])
    x += subw + marginw
    if i == 3:
        x = marginw
        y = y - subw - marginh

x = marginw
y = y - h - plot_skip
subw = figw - 2 * marginw

barax = fig.add_axes([x / figw, y / figh, subw / figw, h / figh])
barax_dummy = fig.add_axes([ x / figw, y / figh, subw / figw, h / figh])
init_barax(barax_dummy)
init_barax(barax)

barax_dummy.xaxis.set_ticks_position("top")
#barax_dummy.set_xticks([round(t,2) for t in threshs])
barax_dummy.set_xticks(threshs)
barax_dummy.set_xticklabels([])
barax_dummy.xaxis.grid(True, linestyle="dotted")

y = barskip
for b in intervals[1][1]:
    rect = pat.Rectangle((b[0], y / h), b[1] - b[0], barh / h, linewidth=0, facecolor=colors[0])
    barax_dummy.add_patch(rect)
    y += barh
    if b != intervals[1][1][-1]:
        y += barskip

fig.savefig("intro.pdf", format='pdf')
fig.savefig("../thesis/img/00-intro-example.pdf", format='pdf')

#plot_samples(gen_intro_ex_2(25))
#plot_edges(gen_intro_ex_2(10))
#write_dataset(gen_intro_ex_2(10), "../datasets/intro_ex2_10")
