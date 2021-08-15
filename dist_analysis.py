with open("dist_circ100.txt") as f:
    lines = f.readlines()

dist = {}

for l in lines:
    toks = l.split(":")
    if len(toks) < 2:
        continue
    if toks[0] == "v":
        tup = toks[1].split(",")
        a = int(tup[0])
        b = int(tup[1])
        key = "v-" + str(a) + "-" + str(b)
        if not key in dist:
            dist[key] = 0
        dist[key] = 1 + dist[key]
    if toks[0] == "r":
        tup = toks[1].split(",")
        a = int(tup[0])
        b = int(tup[1])
        key = "r-" + str(a) + "-" + str(b)
        if not key in dist:
            dist[key] = 0
        dist[key] = 1 + dist[key]

freq = {}

print("Total searches: " + str(len(lines)) + " with " + str(len(dist)) + " distinct key combinations")

for (key, val) in dist.items():
    if not val in freq:
        freq[val] = 0
    freq[val] = freq[val] + 1

for (key, val) in sorted(freq.items()):
    print("# of searches = " + str(key) + " :: " + str(val))
