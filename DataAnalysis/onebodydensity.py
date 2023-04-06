import numpy as np
import matplotlib.pyplot as plt
import os
OUTPUT_DIR = "/home/giammi/Desktop/FYS4411-Project1/Outputs"

def get_histogram(N, NX=100, NY=100, NZ=100):
    HIST_TOT = np.zeros((NX,NY,NZ), dtype=np.uint32)
    threadnums= np.arange(8)
    for thread in threadnums:
        fname = os.path.join(OUTPUT_DIR, "histogram"+str(N)+"_"+str(thread)+".bin")
        hist = np.fromfile(fname, dtype=np.uint32).reshape((NX, NY, NZ))
        HIST_TOT+=hist

    return HIST_TOT
#get histogram of position along xy plane, expect this to be symmetrical
HIST_TOT = get_histogram(10)
hist_xy = np.sum(HIST_TOT, axis=2)
#attempt fit to a gaussian:
log_hist = np.log(hist_xy)
#need the distance from each slot to the center
# Define the size of the matrix and the boundaries
size = 100
x_min, x_max, y_min, y_max = -1.5, 1.5, -1.5, 1.5
# Create a 2D array with x and y coordinates of each slot
x_coords = np.linspace(x_min, x_max, size)
y_coords = np.linspace(y_min, y_max, size)
xx, yy = np.meshgrid(x_coords, y_coords)
# Compute the distance of each slot from the center of the matrix
r2 = (xx)**2 + (yy)**2
plt.plot(r2.reshape(-1), log_hist.reshape(-1), linestyle='none', marker='.')
plt.show()
#get histogram along z axis:
#hist_z = np.sum(HIST_TOT, axis=(0,1))

