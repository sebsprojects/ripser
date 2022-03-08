from os import listdir
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.figure as mplfig
import matplotlib.axes as mplax
from matplotlib import colors as mcol
import math

from matplotlib.ticker import AutoMinorLocator


datasetname = "random100"
base = "./rep_test/output"
n = 100

tex_fonts = {
    "text.usetex": True,
    "font.family": "serif",
    "axes.titlesize": 7,
    "axes.labelsize": 7,
    "font.size": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    'text.latex.preamble': r'\usepackage{amsfonts}\usepackage{bm}\usepackage{amsmath}'
}
plt.rcParams.update(tex_fonts)

files = listdir(base)

def gather_recs(nam):
    recs = []
    for f in files:
        if "repinverse" in f:
            t = "rep_inv"
        elif "rephom" in f:
            t = "rep_hom"
        if nam in f:
            rec = { "type": t,
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
    return recs

docw = 418.25372 / 72
figw = docw
figh = 0.4 * docw

margin = docw * 0.05

w = figw * 0.33 - margin
h = figh - margin * 1.4
x = margin * 0.75
y = margin * 0.6

fig = mplfig.Figure(figsize=(figw, figh))

num = 4
barw = ((num)/10) * 0.011 * (docw / w)
cols = [ "tab:blue", "tab:orange", "#a0cbe8", "#ffbe7d" ]
metrics = ["time", "add_simplex_count", "push_count", "pop_count"]
labels = [ r'time', "adds", "push", "pop" ]
barlabels = [ "rep_hom", "rep_inv" ]


def collect_metric(metric, dim):
    vals = []
    for rec in recs:
        vals.append(rec[metric][dim])
    return [ v / vals[0] for v in vals]

def get_xs(p):
    bw = barw * 1.1
    return [ p - 0.5 * bw, p + 0.5 * bw ]

def init_ax(ax):
    ax.set_ylim([0,25])
    ax.set_xlim([0,0.4])
    ax.tick_params(axis='both', which='major', pad=1)
    ax.yaxis.grid(True, which="both", linestyle="dotted", alpha=0.5)
    ax.plot([ x / 10 for x in range(11)], [ 1 for x in range(11)], c="black", lw=0.5)
    ax.minorticks_on()
    minor_locator = AutoMinorLocator(2)
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.yaxis.set_minor_locator(minor_locator)

recs = gather_recs("random")
rrecs = []
for i, u in enumerate(recs):
    for r in recs:
        if i == 0 and "rep_hom" == r["type"]:
            rrecs.append(r)
            break
        if i == 1 and "rep_inv" == r["type"]:
            rrecs.append(r)
            break
recs = rrecs

print([r["type"] for r in recs])
print(recs)

ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
init_ax(ax)
p = 0.05
xticks = []
dim = 1
for m in metrics:
    rects = ax.bar(get_xs(p), collect_metric(m, dim), barw, color=cols)
    xticks.append(p)
    p = p + 0.1
for i, r in enumerate(rects.patches):
    r.set_label(barlabels[i])
ax.set_xticks(xticks)
ax.set_xticklabels(labels)
ax.set_title(r'random100, $p = 1$', pad=5)
ax.legend(prop={"size": 6}, borderpad=0.3, handlelength=1, loc="upper right")

ax = fig.add_axes([ (margin * 1.75 + w) / figw, y / figh, w / figw, h / figh])
init_ax(ax)
p = 0.05
xticks = []
dim = 2
for m in metrics:
    rects = ax.bar(get_xs(p), collect_metric(m, dim), barw, color=cols)
    xticks.append(p)
    p = p + 0.1
for i, r in enumerate(rects.patches):
    r.set_label(barlabels[i])

ax.set_xticks(xticks)
ax.set_xticklabels(labels)
ax.set_title(r'random100, $p = 2$', pad=5)
ax.legend(prop={"size": 6}, borderpad=0.3, handlelength=1, loc="upper right")



recs = gather_recs("covid")
rrecs = []
for i, u in enumerate(recs):
    for r in recs:
        if i == 0 and "rep_hom" == r["type"]:
            rrecs.append(r)
            break
        if i == 1 and "rep_inv" == r["type"]:
            rrecs.append(r)
            break
recs = rrecs

print([r["type"] for r in recs])
print(recs)

x = margin * 2.75 + 2 * w
ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
init_ax(ax)
p = 0.05
xticks = []
dim = 1
for m in metrics:
    rects = ax.bar(get_xs(p), collect_metric(m, dim), barw, color=cols)
    xticks.append(p)
    p = p + 0.1
for i, r in enumerate(rects.patches):
    r.set_label(barlabels[i])

ax.set_xticks(xticks)
ax.set_xticklabels(labels)
ax.set_title(r'covid10000, $p = 1$', pad=5)
ax.legend(prop={"size": 6}, borderpad=0.3, handlelength=1, loc="upper right")

fig.savefig("test.pdf", format='pdf')
fig.savefig("../thesis/img/02-rep-hominv-comp-1.pdf", format='pdf')
