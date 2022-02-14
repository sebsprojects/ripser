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

def gather_rr(rrs, label):
    d = " & "
    table_string = label
    for dim, rr in enumerate(rrs):
        total_reds = len(rr)
        total_adds = 0
        total_pushes = 0
        total_pops = 0
        for r in rr:
            total_adds += r[3]
            total_pops += r[4]
            total_pushes += r[5]
        s = d + str(total_reds) + d + str(total_adds)
        s += d + str(total_pushes) + d + str(total_pops)
        table_string += s
        if dim == 2:
            break
    return table_string + r' \\'

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

formatstring = "l | r r r r | r r r r | r r r r"
header_dim = r'\texttt{\#red} & \texttt{\#add} & \texttt{\#psh} & \texttt{\#pop}'

table = r'\begin{tabular}{' + formatstring + "}\n"
table += "  version & " + header_dim + " & " + header_dim + r'\\' + "\n"
table += r'  \hline' + "\n"

for i, f in enumerate(input_files):
    rrs = read_rr(pref + f)
    tl = gather_rr(rrs, labels[i])
    table += "  " + tl + "\n"
table += r'\end{tabular}' + "\n"

f = open("../thesis/gen_latex/02-rep-hom-table.txt", "w")
f.write(table)
f.close()
