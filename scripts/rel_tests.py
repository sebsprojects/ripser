import subprocess
import math

#config_file = "rel_test/random100"
#n = 100
#config_file = "rel_test/sphere192"
#n = 192
#config_file = "rel_test/o3_512"
#n = 512
#config_file = "rel_test/o3_256"
#n = 256
config_file = "rel_test/covid_landdistmat"
n = 2500
#config_file = "rel_test/random50"
#n = 50
#config_file = "rel_test/dragon1000"
#n = 1000
#config_file = "rel_test/h3n2"
#n = 2723
#config_file = "rel_test/hiv"
#n = 1088

rs = ["ripsero3_cohom_rel_optimized" ]
rel_configs = ["relative_subcomplex="]

#xs = [ 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99 ]
xs = [ x / 20 for x in range(1, 20)]

for rel in xs:
    rel_configs.append("relative_subcomplex=0-" + str(int(math.floor(n * rel)) - 1))

for r in rs:
    for c in rel_configs:
        config = open(config_file + "_part.ini", "r")
        config_dest = open(config_file + ".ini", "w")
        for l in config:
            config_dest.write(l)
        config_dest.write(c)
        config.close()
        config_dest.close()
        p = subprocess.Popen(["./" + r, config_file + ".ini"])
        retval = p.wait()
