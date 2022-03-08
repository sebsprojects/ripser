from os import listdir
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.figure as mplfig
import matplotlib.axes as mplax
from matplotlib import colors as mcol
import math

from matplotlib.ticker import AutoMinorLocator


base = "./rel_test/output/random100"
name = "random_point_cloud_100"
n = 100

#base = "./rel_test/output/sphere192"
#name = "sphere_3_192"
#n = 192

#base = "./rel_test/output/o3_256"
#name = "o3_1024"
#n = 256

#base = "./rel_test/output/covid_10000"
#name = "covid_landdistmat"
#n = 10000

tex_fonts = {
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 8,
    "font.size": 7,
    "legend.fontsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    'text.latex.preamble': r'\usepackage{amsfonts}\usepackage{bm}'
}
plt.rcParams.update(tex_fonts)

files = listdir(base)
rel_ends = []
recs = []
for f in files:
    toks = f.split("_")
    rel_end = 0
    for t in toks:
        if "rel" in t:
            rel_end = int(t.split("-")[1]) + 1
    rel_ends.append(rel_end)
    if name in f:
        rec = { 'rel_end': rel_end,
                "col_count": [],
                'red_count': [],
                "zero_count": [],
                "clearing_count": [],
                "app_count": [],
                "add_count": [],
                "add_simplex_count": [],
                "push_count": [],
                "pop_count": [],
                "app_facet_count": [],
                "app_cofacet_count": [],
                "time": [],
                "mem": []
        }
        m = open(base + "/" + f, "r")
        for line in m:
            if "col_count" in line:
                rec["col_count"].append(int(line.split("=")[1]))
            if "red_count" in line:
                rec["red_count"].append(int(line.split("=")[1]))
            if "zero_count" in line:
                rec["zero_count"].append(int(line.split("=")[1]))
            if "clearing_count" in line:
                rec["clearing_count"].append(int(line.split("=")[1]))
            if "app_count" in line:
                rec["app_count"].append(int(line.split("=")[1]))
            if "add_count" in line:
                rec["add_count"].append(int(line.split("=")[1]))
            if "add_simplex_count" in line:
                rec["add_simplex_count"].append(int(line.split("=")[1]))
            if "push_count" in line:
                rec["push_count"].append(int(line.split("=")[1]))
            if "pop_count" in line:
                rec["pop_count"].append(int(line.split("=")[1]))
            if "app_facet_count" in line:
                rec["app_facet_count"].append(int(line.split("=")[1]))
            if "app_cofacet_count" in line:
                rec["app_cofacet_count"].append(int(line.split("=")[1]))
            if "time" in line:
                rec["time"].append(float(line.split("=")[1]))
            if "mem" in line:
                rec["mem"].append(float(line.split("=")[1]))
        recs.append(rec)
        m.close()

docw = 418.25372 / 72
figw = docw * 0.75
margin = 0.05 * docw
w = docw * 0.24 - margin
h = w
figh = h * 2 + 2.5 * margin

x = margin
y = margin * 2 + h

fig = mplfig.Figure(figsize=(figw, figh))

rel_ends = sorted(rel_ends)
xs = [ re / n for re in rel_ends ] + [ 1.0 ]

def collect_metric(label, dim):
    ys = []
    max_val = 0
    for re in rel_ends:
        for r in recs:
            if r["rel_end"] == re:
                val = r[label][dim]
                ys.append(val)
                max_val = max(max_val, val)
    ys.append(0)
    max_below = 1.0
    if ys[0] > 0:
        for i, y in enumerate(reversed(ys)):
            if y / ys[0] < 1.0:
                max_below = xs[-i - 1]
            else:
                break
        return max_val / ys[0], max_below, [y / ys[0] for y in ys]
    else:
        return max_val, 0, ys

def collect_relcol_metric(label, dim):
    ys = []
    max_val = 0
    for re in rel_ends:
        for r in recs:
            if r["rel_end"] == re:
                val = r[label][dim]
                relcol = r["col_count"][dim]
                print(val, relcol)
                ys.append(val / relcol)
                max_val = max(max_val, val / relcol)
    ys.append(0)
    max_below = 1.0
    return max_val, 0, ys

def collect_corrected_col1_metric():
    ys = []
    max_val = 0
    for re in rel_ends:
        for r in recs:
            if r["rel_end"] == re:
                val = r["red_count"][dim]
                cor = n - r["col_count"][0]
                print(val, cor)
                ys.append(val - cor)
    return max_val / ys[0], [y / ys[0] for y in ys] + [0]

def init_ax(ax):
    ax.set_xlim([0,1])
    ax.set_xticks([i * 2 / 10 for i in range(6)])
    ax.set_xticklabels([i * 2 / 10 for i in range(6)])
    ax.minorticks_on()
    minor_locator = AutoMinorLocator(2)
    ax.xaxis.set_minor_locator(minor_locator)
    ax.yaxis.set_tick_params(which='minor', bottom=False)
    ax.xaxis.grid(True, which="Both", linestyle="dotted")
    ax.yaxis.grid(True, which="Major", linestyle="dotted")
    ax.tick_params(axis='both', which='major', pad=1)

def plot_guides(ax, mb, include_hori=True):
    col = "tab:red"
    st = "solid"
    lw = 0.6
    if include_hori:
        ax.plot(xs, [1 for i in xs], color="black", lw=lw, linestyle=st)
    #if mb < 1:
    #    ax.axvline(mb, color=col, lw=lw, linestyle=st)

def plot_metric(ax, metrics, dim, title, labels, relcol=False):
    ylim = 0
    for i, m in enumerate(metrics):
        if relcol:
            l, mb, ys = collect_relcol_metric(m, dim)
        else:
            l, mb, ys = collect_metric(m, dim)
        label = ""
        if len(labels) > i:
            label=labels[i]
        ax.plot(xs, ys, label=label)
        plot_guides(ax, mb)
        ylim = max(ylim, l)
    ax.legend(prop={"size": 6}, borderpad=0.3, handlelength=1)
    ax.set_title(title, pad=4)
    ax.set_ylim([0, ylim * 1.05])

dim = 1
ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
init_ax(ax)
plot_metric(ax, ["time", "mem"], dim, "runtime", ["time", "vmSize"])
x += w + margin

ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
init_ax(ax)
plot_metric(ax, ["col_count", "red_count"], dim, "column count", ["total", "to reduce"])
#ylim, ys = collect_corrected_col1_metric()
#ax.plot(xs, ys, c="tab:red")
x += w + margin

ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
init_ax(ax)
plot_metric(ax, ["clearing_count", "app_count"], dim, "clearing, app count", ["clearing", "apparent"])
x += w + margin

x = margin
y -= h + margin * 1.3

ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
init_ax(ax)
plot_metric(ax, ["add_count", "add_simplex_count"], dim, "addition count", ["add", "add_simplex"])
x += w + margin

#ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
#init_ax(ax)
#plot_metric(ax, ["zero_count", "app_count"], dim, "rel zero, app count", ["zero", "apparent"], True)
#x += w + margin

ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
init_ax(ax)
plot_metric(ax, ["push_count","pop_count"], dim, "push, pop count", ["push", "pop"])
x += w + margin

ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
init_ax(ax)
plot_metric(ax, ["app_facet_count", "app_cofacet_count"], dim, "app enum count", ["facet", "cofacet"])
x += w + margin

#plot_metric(x, y, , dim)
#x += w + margin

fig.savefig("test.pdf", format='pdf')
#fig.savefig("../thesis/img/02-rel-random100-perf.pdf", format='pdf')
