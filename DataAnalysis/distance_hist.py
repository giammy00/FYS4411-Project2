import numpy as np
import matplotlib.pyplot as plt
from path import OUTPUT_DIR
import os
def import_sample(N_part_list, threadnums):
    for N_particles in N_part_list:
        distance_hist = np.empty(0)
        for thread in threadnums:
            filename="sampledDistances_" + str(N_particles) + "_" + str(thread) + ".bin"
            x = np.fromfile( os.path.join(OUTPUT_DIR, "OPTIMIZE_RUNS", "lbfgs_2p_4nodes" , filename), dtype=float)
            distance_hist= np.append(distance_hist, x)

    return distance_hist
Np = [2,]
threads = np.arange(8)
bin_edges = np.arange(0,4.0, 0.05)
rij = import_sample(Np, threads)
plt.hist(rij, bins=bin_edges)
plt.grid(visible=True)
plt.ylabel("Counts", fontsize=12)
plt.xlabel("$r_{ij}$", fontsize=12)
plt.savefig("distanceshist.pdf", bbox_inches='tight')
plt.show()