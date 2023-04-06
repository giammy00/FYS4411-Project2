import numpy as np
import matplotlib.pyplot as plt
import os
from path import OUTPUT_DIR
def get_histogram(N, NX=100, NY=100, NZ=100):
    '''get normalized histogram for experiment with N particles.
    maximum element of output is 1'''
    HIST_TOT = np.zeros((NX,NY,NZ), dtype=np.uint32)
    threadnums= np.arange(8)
    for thread in threadnums:
        fname = os.path.join(OUTPUT_DIR, "histogram"+str(N)+"_"+str(thread)+".bin")
        hist = np.fromfile(fname, dtype=np.uint32).reshape((NX, NY, NZ))
        HIST_TOT+=hist
    
    #HIST_TOT=HIST_TOT/np.max(HIST_TOT)
    return HIST_TOT

####PLOT TWO HISTOGRAMS 
def plot_hists():
    ticks = np.linspace(0, 100, 6)
    labels = np.linspace(-1.5, 1.5, 6)

    histogram = get_histogram(N=100)
    fig, ax = plt.subplots(1,2,figsize=(9,3.44))
    hist_xy = np.sum(histogram, axis=2)
    hist_xz = np.sum(histogram, axis=1)
    img =[]
    img.append( ax[0].imshow(hist_xy) ) 
    ax[0].set_title("x-y plane")
    img.append(ax[1].imshow(hist_xz) ) 
    ax[1].set_title("x-z plane")

    for ax_cur, imgcur  in zip(ax, img):
        ax_cur.set_xticks(ticks)
        ax_cur.set_yticks(ticks)
        ax_cur.set_xticklabels(['{:.1f}'.format(label) for label in labels])
        ax_cur.set_yticklabels(['{:.1f}'.format(label) for label in labels[::-1]])
        fig.colorbar(imgcur, ax=ax_cur, format='%.0e')
    plt.savefig("onebody_xy_xz_100.pdf", bbox_inches='tight')
    plt.show()

def gaussian_fit(n):
    
    hist_xy = np.sum(get_histogram(n), axis=2)
    #attempt fit to a gaussian:
    hist_xy_flat = hist_xy.reshape(-1)
    nnz_idx = np.nonzero(hist_xy_flat>np.exp(7))
    log_hist = np.log(hist_xy_flat[nnz_idx])
    #need the distance from each slot to the center
    # Define the size of the matrix and the boundaries
    size = 100
    x_min, x_max, y_min, y_max = -1.5, 1.5, -1.5, 1.5
    # Create a 2D array with x and y coordinates of each slot
    x_coords = np.linspace(x_min, x_max, size)
    y_coords = np.linspace(y_min, y_max, size)
    xx, yy = np.meshgrid(x_coords, y_coords)
    # Compute the distance of each slot from the origin 0,0
    r2 = (xx)**2 + (yy)**2
    r2 = r2.reshape(-1)
    r2 = r2[nnz_idx]
    #plt.plot(r2, log_hist, linestyle='none', marker='.')
    coeffs, var = np.polyfit(r2, log_hist, deg=1, cov=True )
    print(f"$\\alpha = {coeffs[0]/2} \\pm  {np.sqrt(var[0,0])/2}$ ")
    #plt.show()
    #get histogram along z axis:
#hist_z = np.sum(HIST_TOT, axis=(0,1))

if __name__ == "__main__":
    # for n in [10,50,100]:
    #     gaussian_fit(n)
    plot_hists()