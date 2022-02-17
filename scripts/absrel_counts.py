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
ax.set_xlabel(r'fraction of cut off points')
ax.set_ylabel(r'fraction of required columns')
ax.tick_params(axis='both', which='major', pad=1)
for d in range(3):
    line, = ax.plot(rs,[(count(n*(1-r), d)) / count(n,d) for r in rs])
    line.set_label(r'$p=' + str(d) + "$")
ax.set_xlim([0,1])
ax.set_xticks([i/5 for i in range(0,6)])
ax.set_ylim([0,1])
ax.xaxis.grid(True, which="Major", linestyle="dotted")
ax.yaxis.grid(True, which="Major", linestyle="dotted")
ax.legend(prop={"size": 7})

x = 4.5 * margin + w
ax = fig.add_axes([ x / figw, y / figh, w / figw, h / figh])
ax.set_xlabel(r'fraction of relative points')
ax.set_ylabel(r'fraction of required columns')
ax.set_xlim([0,1])
ax.set_ylim([0,1])
ax.set_xticks([i/5 for i in range(0,6)])
ax.set_aspect("equal")
ax.tick_params(axis='both', which='major', pad=1)
ax.xaxis.grid(True, which="Major", linestyle="dotted")
ax.yaxis.grid(True, which="Major", linestyle="dotted")
#init_ax(ax)
for d in range(3):
    ys = [(count(n,d) - count(r*n, d))/count(n,d) for r in rs]
    line, = ax.plot(rs, ys)
    line.set_label(r'$p=' + str(d) + "$")
#line, = ax.plot(y,r11)
#line.set_label(r'$r=0.5$')
#line, = ax.plot(y,r12)
#line.set_label(r'$r=0.25$')
#line, = ax.plot(y,r13)
#line.set_label(r'$r=0.1$')
ax.legend(prop={"size": 7})


#xx = 4 * margin + w * 2
#ax = fig.add_axes([ xx / figw, margin / figh, w / figw, h / figh])
#ax.set_xlim([0,1])
#ax.set_ylim([0,1])
#ax.set_xticks([i/5 for i in range(6)])
#ax.set_yticks([i/5 for i in range(6)])
#ax.set_xlabel(r'ratio of data points (absolute / total)')
#ax.set_ylabel(r'ratio of simplices')
#ax.xaxis.grid(True, which="Major", linestyle="dotted")
#ax.yaxis.grid(True, which="Major", linestyle="dotted")
#line, = ax.plot(x,d0,c=colors[3])
#line.set_label(r'$d=0$')
#line, = ax.plot(x,d1,c=colors[4])
#line.set_label(r'$d=1$')
#line, = ax.plot(x,d2,c=colors[5])
#line.set_label(r'$d=2$')
#ax.legend(prop={"size": 7})

fig.savefig("test.pdf", format='pdf')
fig.savefig("../thesis/img/02-rel-counts.pdf", format='pdf')
