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

def abs_count_d0(n):
    return n

def rel_count_d0(n, s):
    return n*s

def rel_abs_frac_d0(n, s):
    return rel_count_d0(n,s) / abs_count_d0(n)

def abs_count_d1(n):
    return n*(n-1)/2

def rel_count_d1(n, s):
    return abs_count_d1(n*s) + (n*s)*(n-n*s)

def rel_abs_frac_d1(n, s):
    return rel_count_d1(n,s) / abs_count_d1(n)

def abs_count_d2(n):
    return n*(n-1)*(n-2)/6

def rel_count_d2(n, s):
    return abs_count_d2(n*s) + (n*s)*(n-n*s)*(n-2)/2

def rel_abs_frac_d2(n,s):
    return rel_count_d2(n,s) / abs_count_d2(n)


n = 10000

x = []
d0 = []
d1 = []
d2 = []


res = 101
for i in range(res):
    x.append(i / res)
    d0.append(rel_abs_frac_d0(n, i / res))
    d1.append(rel_abs_frac_d1(n, i / res))
    d2.append(rel_abs_frac_d2(n, i / res))

y = []
a1 = []
r11 = []
r12 = []
r13 = []

a2 = []
r21 = []
r22 = []
r23 = []

for i in range(res):
    n = i * 100
    y.append(n)
    a1.append(abs_count_d1(n))
    r11.append(rel_count_d1(n, 0.5))
    r12.append(rel_count_d1(n, 0.25))
    r13.append(rel_count_d1(n, 0.1))
    a2.append(abs_count_d2(n))
    r21.append(rel_count_d2(n, 0.5))
    r22.append(rel_count_d2(n, 0.25))
    r23.append(rel_count_d2(n, 0.1))

docw = 418.25372 / 72
figw = docw
figh = docw * 0.33

margin = docw * 0.05
w = figw * 0.33 - 1.5 * margin
h = w

fig = mplfig.Figure(figsize=(figw, figh))
ax = fig.add_axes([ margin / figw, margin / figh, w / figw, h / figh])
#ax.set_xlabel(r'number of data points')
#ax.set_ylabel(r'number of 1-simplices')
ax.set_xlim([0,10000])
ax.set_ylim([0,50000000])
ax.xaxis.grid(True, which="Major", linestyle="dotted")
ax.yaxis.grid(True, which="Major", linestyle="dotted")
#init_ax(ax)
line, = ax.plot(y,a1)
line.set_label(r'$r=1$')
line, = ax.plot(y,r11)
line.set_label(r'$r=0.5$')
line, = ax.plot(y,r12)
line.set_label(r'$r=0.25$')
line, = ax.plot(y,r13)
line.set_label(r'$r=0.1$')
ax.legend(prop={"size": 7})

xx = 2.5 * margin + w
ax = fig.add_axes([ xx / figw, margin / figh, w / figw, h / figh])
line, = ax.plot(y,a2)
line.set_label(r'$r=1$')
line, = ax.plot(y,r21)
line.set_label(r'$r=0.5$')
line, = ax.plot(y,r22)
line.set_label(r'$r=0.25$')
line, = ax.plot(y,r23)
line.set_label(r'$r=0.1$')
ax.set_xlim([0,10000])
ax.set_ylim([0,200000000000])
ax.xaxis.grid(True, which="Major", linestyle="dotted")
ax.yaxis.grid(True, which="Major", linestyle="dotted")
ax.legend(prop={"size": 7})

xx = 4 * margin + w * 2
ax = fig.add_axes([ xx / figw, margin / figh, w / figw, h / figh])
ax.set_xlim([0,1])
ax.set_ylim([0,1])
ax.set_xticks([i/5 for i in range(6)])
ax.set_yticks([i/5 for i in range(6)])
#ax.set_xlabel(r'ratio of data points (absolute / total)')
#ax.set_ylabel(r'ratio of simplices')
ax.xaxis.grid(True, which="Major", linestyle="dotted")
ax.yaxis.grid(True, which="Major", linestyle="dotted")
line, = ax.plot(x,d0,c=colors[3])
line.set_label(r'$d=0$')
line, = ax.plot(x,d1,c=colors[4])
line.set_label(r'$d=1$')
line, = ax.plot(x,d2,c=colors[5])
line.set_label(r'$d=2$')
ax.legend(prop={"size": 7})

fig.savefig("test.pdf", format='pdf')
fig.savefig("../thesis/img/02-rel-counts.pdf", format='pdf')
