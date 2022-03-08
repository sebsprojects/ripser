import matplotlib
import matplotlib.pyplot as plt
import matplotlib.figure as mplfig
import matplotlib.axes as mplax
from matplotlib import colors as mcol
import math

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

def count(n, d):
    c = n
    r = 1
    for i in range(1, d + 1):
        c *= (n - i)
        r *= (i+1)
    return c / r

print(count(100,1))

docw = 418.25372 / 72
figw = docw * 0.75

margin = docw * 0.05
w = docw * 0.35 - 2 * margin
h = w
figh = h + 2 * margin

x = margin * 1.5
y = margin * 1.3
fig = mplfig.Figure(figsize=(figw, figh))

rs = [ i / 100 for i in range(0,101) ]
n = 100


ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
ys = [ count(r*n, 1) - count(r * n, 0) for r in rs]
ax.plot(rs, ys)

fig.savefig("test.pdf", format='pdf')
#fig.savefig("../thesis/img/02-rel-counts.pdf", format='pdf')
