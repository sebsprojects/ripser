import subprocess

rs = ["ripsero3_cohom_rel", "ripsero3_cohom_rel_optimized"]
rel_configs = ["relative_subcomplex="]
for rel in [[0,i*10+9] for i in range(0, 9)]:
    rel_configs.append("relative_subcomplex=" + str(rel[0]) + "-" + str(rel[1]))

for r in rs:
    for c in rel_configs:
        config = open("rel_test/random100_part.ini", "r")
        config_dest = open("rel_test/random100.ini", "w")
        for l in config:
            config_dest.write(l)
        config_dest.write(c)
        config.close()
        config_dest.close()
        output = subprocess.Popen(["./" + r, "rel_test/random100.ini"])
