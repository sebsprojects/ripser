import matplotlib
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from matplotlib import colors as mcol
from matplotlib.lines import Line2D
import numpy as np
import math
import sys
from functools import cmp_to_key

if len(sys.argv) < 2:
    print("error: missing config file-path")
    sys.exit()

# PARAMETER
# -----------------------------------------------------------------------------
dims = [0,1,2]
doc_w = 418.25372 / 72
figw = doc_w * 1.2  # FACTOR 1 TOO SMALL FOR SOME STUPID REASON


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

def findex_comp_bar(a, b):
    if a[3][0] < b[3][0] or (a[3][0] == b[3][0] and a[3][1] < b[3][1]):
        return -1
    elif a[3] == b[3]:
        return 0
    else:
        return 1

def diam_comp_bar(a, b):
    if a[1][0] < b[1][0] or (a[1][0] == b[1][0] and a[1][1] < b[1][1]):
        return -1
    elif a[1] == b[1]:
        return 0
    else:
        return 1

def parse_simplex(s, dim, simplex_lookup, rel):
    toks = s.split("-")
    index = int(toks[0].strip())
    filt_index = []
    if len(simplex_lookup) > 0:
        filt_index = simplex_lookup.index([dim, index])
    if len(toks) > 1:
        vertices = [int(t) for t in toks[1].split("'")]
        is_rel = 1
        for v in vertices:
            if not v in rel:
                is_rel = 0
                break
        return [index, filt_index, vertices, is_rel]
    return [index, filt_index, [], -1]


def parse_bar(line, simplex_lookup=[], rel=[]):
    dim = int(line[2])
    toks = [t.strip() for t in line[3:].split(";")]
    # diams
    diam_toks = toks[0].replace("[","").replace(")", "").split(",")
    diams = [float(diam_toks[0].strip())]
    if diam_toks[1].strip() != "":
        diams.append(float(diam_toks[1].strip()))
    # indices
    index_toks = toks[1].replace("[","").replace(")", "").split(",")
    index_end_string = index_toks[1]
    start_simplex = parse_simplex(index_toks[0], dim, simplex_lookup, rel)
    indices = [start_simplex[0]]
    filt_indices = [start_simplex[1]]
    is_rel = start_simplex[3]
    if(index_toks[1].strip() != ""):
        end_simplex = parse_simplex(index_toks[1], dim + 1, simplex_lookup, rel)
        indices.append(end_simplex[0])
        filt_indices.append(end_simplex[1])
    return [dim, diams, indices, filt_indices, is_rel]


def read_barcode(input_file, simplex_lookup = [], rel = []):
    f = open(input_file, "r");
    intervals = []
    max_diam = 0
    offs = 1
    bars = []
    max_diam = 0
    max_index = 0
    max_findex = 0
    for line in f:
        if line.startswith("#b"):
            bar = parse_bar(line, simplex_lookup, rel)
            dim = bar[0]
            if len(rel) != 0 and bar[4] != 0:
                continue
            max_diam = max(max_diam, bar[1][0])
            max_index = max(max_index, bar[2][0])
            #max_findex = max(max_index, bar[3][0])
            if len(bar[1]) > 1:
                max_diam = max(max_diam, bar[1][1])
                max_index = max(max_index, bar[2][1])
                # max_findex = max(max_findex, bar[3][1])
            if len(bars) - 1 < dim:
                bars.append([])
            bars[dim].append(bar)
    for dim, bs in enumerate(bars):
        for b in bs:
            if len(b[1]) == 1:
                b[1].append(max_diam)
                b[2].append(max_index)
                b[3].append(max_findex)
        bars[dim] = sorted(bs, key=cmp_to_key(diam_comp_bar))
    return [max_diam, max_index, max_findex], bars


#TODO: Only works for a single interval
def read_rel(input_file):
    f = open(input_file, "r")
    for line in f:
        if line.find("relative") != -1:
            s = line.split(":")[1].strip()
            istart = int(s.split("-")[0])
            iend = int(s.split("-")[1])
            return range(istart,iend + 1)

#TODO: Adjust for absolute_subcomplex that does not start at 0
def read_simplices(input_file):
    f = open(input_file, "r")
    simplex_lookup = []
    for line in f:
        if line.startswith("#s"):
            dim = int(line[2])
            index = int(line[3:].split("-")[0].strip())
            simplex_lookup.append([dim, index])
    return simplex_lookup

def plot_barcode(ax, max_bound, bars, xbars1 = [], xbars2 = []):
    h_inc = 1
    h_skip = h_inc * 1.5
    h = 0
    ax.set_xlim(1, max_bound+1)
    #ax.set_xticks(range(1, max_bound+1))
    #ax.xaxis.grid(True, linestyle="dotted")
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for dim, bc in enumerate(bars):
        if dim == 0:
            continue
        h += h_skip * 3
        eps = max_bound * 0.01
        for bar in bc:
            c = "black"
            if xbars1 != []:
                for xbar in xbars1[dim]:
                    if bar[1] == xbar[1]:
                        c = "tab:blue"
                        break
                    #elif (abs(bar[1][0] - xbar[1][0]) < eps) and (abs(bar[1][1] - xbar[1][1]) < eps):
                    elif bar[1][0] == xbar[1][0] and abs(bar[1][1] - xbar[1][1]) < eps:
                        c = "tab:blue"
                        break
            if xbars2 != []:
                for xbar in xbars2[dim]:
                    if bar[1] == xbar[1]:
                        c = "tab:orange"
                        break
                    elif bar[1][0] == xbar[1][0] and abs(bar[1][1] - xbar[1][1]) < eps:
                        c = "tab:orange"
                        break
            ax.plot(bar[1], [h, h], color=c, linewidth=1)
            if bar != bc[-1]:
                h += h_inc
    h += h_skip
    ax.set_ylim(0, h)


# MAIN
# -----------------------------------------------------------------------------

fig, axs = plt.subplots(3, 1, figsize=(figw, 2 * figw))
big_file = sys.argv[1]
small_file = sys.argv[2]
rel_file = sys.argv[3]

big_bounds, big_bars = read_barcode(big_file)
small_bounds, small_bars = read_barcode(small_file)
rel_bounds, rel_bars = read_barcode(rel_file)

plot_barcode(axs.flat[0], big_bounds[0], small_bars, [], [])
plot_barcode(axs.flat[1], big_bounds[0], big_bars, [], [])
plot_barcode(axs.flat[2], rel_bounds[0], rel_bars, [], [])

max_b = max(big_bounds[0], small_bounds[0])
max_b = max(max_b, rel_bounds[0])
axs.flat[0].set_xlim(0, max_b)
axs.flat[0].set_title(r'\texttt{abs [0,249]}')
axs.flat[1].set_xlim(0, max_b)
axs.flat[1].set_title(r'\texttt{abs [0,499]}')
axs.flat[2].set_xlim(0, max_b)
axs.flat[2].set_title(r'\texttt{abs [0,499] rel [0,249]}')

plt.savefig("visu/covid_landdistmat_spike_2_dim1.pdf", format='pdf', bbox_inches="tight", pad_inches=0.05)
#plt.show()
