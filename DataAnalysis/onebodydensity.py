import numpy as np
import matplotlib.pyplot as plt
import os
OUTPUT_DIR = "/home/giammi/Desktop/FYS4411-Project1/Outputs"
# read histogram from binary file
histfilename = os.path.join(OUTPUT_DIR, "histogram10_5.bin")
NX = 100
NY= 100
NZ= 100
hist = np.fromfile(histfilename, dtype=np.uint32).reshape((NX, NY, NZ))
hist_2d = np.sum(hist, axis=0)
hist_norm = (hist_2d - np.min(hist_2d)) / (np.max(hist_2d) - np.min(hist_2d))
plt.imshow(hist_2d)
plt.show()