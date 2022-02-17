import sys

def read_rr(input_file):
    f = open(input_file, "r")
    rrs = []
    for line in f:
        toks = line.split(";")
        rr = [int(t.strip()) for t in toks]
        dim = rr[0]
        if len(rrs) <= dim:
            rrs.append([])
        rrs[dim].append(rr)
    return rrs

def gather_rr(rrs, lines):
    for dim, rr in enumerate(rrs):
        total_reds = len(rr)
        total_adds = 0
        total_pushes = 0
        total_pops = 0
        for r in rr:
            total_adds += r[3]
            total_pops += r[4]
            total_pushes += r[5]
        lines[dim * 4 + 0] += " & " + str(total_reds)
        lines[dim * 4 + 1] += " & " + str(total_adds)
        lines[dim * 4 + 2] += " & " + str(total_pushes)
        lines[dim * 4 + 3] += " & " + str(total_pops)
        if dim == 3:
            break
    return lines

pref = "./output/rr/random_point_cloud_50_16_"

input_files = []
input_files.append("_d3tMr1fu00_rr_hom_clearing.txt")
input_files.append("_d3tMr1fu00_rr_hom_clearing_opt.txt")
input_files.append("_d2tMr1fu00_rr_cohom_clearing.txt")
input_files.append("_d2tMr1fu00_rr_cohom_clearing_opt.txt")

labels = []
labels.append(r'\texttt{cl}')
labels.append(r'\texttt{cl-opt}')
labels.append(r'\texttt{co-cl}')
labels.append(r'\texttt{co-cl-opt}')

formatstring = "l | r r r r"
vlabels = [r'\texttt{\#red}', r'\texttt{\#add}',  r'\texttt{\#psh}', r'\texttt{\#pop}']

table = r'\begin{tabular}{' + formatstring + "}\n"
table += "  & " + labels[0] + " & " + labels[1] + " & " + labels[2] + " & " + labels[3] + r'\\' + "\n"
table += r'  \hline' + "\n"

lines = [vlabels[i] for i in range(4)]
lines += ([vlabels[i] for i in range(4)])
lines += ([vlabels[i] for i in range(4)])
lines += ([vlabels[i] for i in range(4)])

for i, f in enumerate(input_files):
    rrs = read_rr(pref + f)
    lines = gather_rr(rrs, lines)

for i, l in enumerate(lines):
    table += l + r' \\' + "\n"
    if i == 3 or i == 7 or i == 11:
        table += r'\hline' + "\n"

table += r'\end{tabular}' + "\n"
print(table)

f = open("../thesis/gen_latex/02-rep-hom-table-v.txt", "w")
f.write(table)
f.close()
