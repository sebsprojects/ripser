from os import listdir
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.figure as mplfig
import matplotlib.axes as mplax
from matplotlib import colors as mcol
import math

#base = "./rel_test/output/random100"
#name = "random_point_cloud_100"
#n = 100

#base = "./rel_test/output/sphere192"
#name = "sphere_3_192"
#n = 192

#base = "./rel_test/output/o3_256"
#name = "o3_1024"
#n = 256

base = "./rel_test/output/covid"
name = "covid_landdistmat"
n = 5000
dim = 1

tex_fonts = {
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 9,
    "font.size": 7,
    "legend.fontsize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    'text.latex.preamble': r'\usepackage{amsfonts}\usepackage{bm}'
}
plt.rcParams.update(tex_fonts)

files = listdir(base)
rel_ends = []
recs = []
for f in files:
    toks = f.split("_")
    rel_end = 0
    cl_app = False
    for t in toks:
        if "rel" in t:
            rel_end = int(t.split("-")[1]) + 1
        if "rrcl-app" in t:
            cl_app = True
    rel_ends.append(rel_end)
    if name in f and cl_app:
        rec = { 'rel': rel_end, 'app': cl_app, 'red': [], "add": [], "pop": [], "push": [], "time": [], "class": [], "zero": [], "appa": [], "cl": [] }
        m = open(base + "/" + f, "r")
        for line in m:
            if "red_count" in line:
                rec["red"].append(int(line.split("=")[1]))
            if "class_count" in line:
                rec["class"].append(int(line.split("=")[1]))
            if "zero_count" in line:
                rec["zero"].append(int(line.split("=")[1]))
            if "push_count" in line:
                rec["push"].append(int(line.split("=")[1]))
            if "pop_count" in line:
                rec["pop"].append(int(line.split("=")[1]))
            if "add_count" in line:
                rec["add"].append(int(line.split("=")[1]))
            if "time" in line:
                rec["time"].append(float(line.split("=")[1]))
            if "app_count" in line:
                rec["appa"].append(float(line.split("=")[1]))
            if "cl_count" in line:
                rec["cl"].append(float(line.split("=")[1]))
        recs.append(rec)
        m.close()

docw = 418.25372 / 72
figw = docw


margin = 0.05 * docw
w = docw * 0.24 - margin
h = w
figh = h * dim + 2.5 * margin

x = margin
y = margin * dim + h * (dim - 1)

fig = mplfig.Figure(figsize=(figw, figh))

rel_ends = sorted(rel_ends)
xs = [ re / n for re in rel_ends ] + [ 1.0 ]

def collect_metric(label, dim, relative=True):
    ys = []
    max_val = 0
    for re in rel_ends:
        for r in recs:
            if r["rel"] == re:
                if r["app"]:
                    val = r[label][dim]
                    ys.append(val)
                    max_val = max(max_val, val)
    ys.append(0)
    max_below = 1.0
    if relative and ys[0] > 0:
        for i, y in enumerate(reversed(ys)):
            if y / ys[0] < 1.0:
                max_below = xs[-i - 1]
            else:
                break
        return max_val / ys[0], max_below, [y / ys[0] for y in ys]
    else:
        return max_val, ys

def init_ax(ax):
    ax.set_xlim([0,1])
    ax.set_xticks([i / 10 for i in range(11)])
    ax.set_xticklabels(["" if i % 2 == 1 else i / 10 for i in range(11)])
    ax.xaxis.grid(True, which="Major", linestyle="dotted")
    ax.yaxis.grid(True, which="Major", linestyle="dotted")

def plot_guides(ax, mb):
    col = "tab:red"
    st = "dashed"
    lw = 0.6
    ax.plot(xs, [1 for i in xs], color=col, lw=lw, linestyle=st)
    if mb < 1:
        ax.axvline(mb, color=col, lw=lw, linestyle=st)

ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
init_ax(ax)
y_lim, mb, ys = collect_metric("red", 1)
ax.set_ylim([0, y_lim * 1.1])
ax.plot(xs, ys)
ax.set_title(r'columns', pad=4)
plot_guides(ax, mb)

x += w + margin
ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
init_ax(ax)
y_lim, mb, ys = collect_metric("add", 1)
ax.set_ylim([0, y_lim * 1.1])
ax.plot(xs, ys)
ax.set_title(r'additions', pad=4)
plot_guides(ax, mb)

x += w + margin
ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
init_ax(ax)
y_lim, mb, ys = collect_metric("push", 1)
ax.set_ylim([0, y_lim * 1.1])
ax.plot(xs, ys)
ax.set_title(r'pushes', pad=4)
plot_guides(ax, mb)

x += w + margin
ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
init_ax(ax)
y_lim, mb, ys = collect_metric("pop", 1)
ax.set_title(r'pops', pad=4)
ax.set_ylim([0, y_lim * 1.1])
ax.plot(xs, ys)
plot_guides(ax, mb)

if dim == 2:
    x = margin
    y = margin * 1

    ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
    init_ax(ax)
    y_lim, mb, ys = collect_metric("red", 2)
    ax.set_ylim([0, y_lim * 1.1])
    ax.plot(xs, ys)
    plot_guides(ax, mb)

    x += w + margin
    ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
    init_ax(ax)
    y_lim, mb, ys = collect_metric("add", 2)
    ax.set_ylim([0, y_lim * 1.1])
    ax.plot(xs, ys)
    plot_guides(ax, mb)

    x += w + margin
    ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
    init_ax(ax)
    y_lim, mb, ys = collect_metric("push", 2)
    ax.set_ylim([0, y_lim * 1.1])
    ax.plot(xs, ys)
    plot_guides(ax, mb)

    x += w + margin
    ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
    init_ax(ax)
    y_lim, mb, ys = collect_metric("pop", 2)
    ax.set_ylim([0, y_lim * 1.1])
    ax.plot(xs, ys)
    plot_guides(ax, mb)


fig.savefig("test.pdf", format='pdf')
#fig.savefig("../thesis/img/02-rel-random100-perf.pdf", format='pdf')
