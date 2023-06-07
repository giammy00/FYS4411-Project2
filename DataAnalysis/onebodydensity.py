import numpy as np
import matplotlib.pyplot as plt
import os
from path import OUTPUT_DIR
def get_histogram(N, folder="./Outputs", NX=100, NY=100, NZ=100):
    '''get unnormalized histogram for experiment with N particles.
    args:
        folder: must be the folder containing the histogram data
        NX, NY, NZ define the shape of the histogram.
    '''
    hist_tot_a = np.zeros((NX,NY,NZ), dtype=np.uint32)
    hist_tot_b = np.zeros((NX,NY,NZ), dtype=np.uint32)

    threadnums= np.arange(8)
    for thread in threadnums:
        fname_a = os.path.join(folder, "histogram1_"+str(N)+"_"+str(thread)+".bin")
        hist = np.fromfile(fname_a, dtype=np.uint32).reshape((NX, NY, NZ))
        hist_tot_a+=hist
        
        fname_b = os.path.join(folder, "histogram2_"+str(N)+"_"+str(thread)+".bin")
        hist = np.fromfile(fname_b, dtype=np.uint32).reshape((NX, NY, NZ))
        hist_tot_b+=hist
    
    #HIST_TOT=HIST_TOT/np.max(HIST_TOT)
    return [hist_tot_a, hist_tot_b]

####PLOT TWO HISTOGRAMS 
def plot_hists_all():
    ticks = np.linspace(0, 100, 4)
    labels = np.linspace(-1.5, 1.5, 4)
    fig, ax = plt.subplots(2,5,figsize=(9,4))
    for mode_idx, mode in enumerate(['lbfgs', 'gd']):

        for idx, NH in enumerate([2,4,6,8,10]):
            foldername = os.path.join(OUTPUT_DIR,"OPTIMIZE_RUNS" , mode+"_2p_"+str(NH)+"nodes")
            histograms = get_histogram(N=2, folder=foldername)
            overlapped_hist = histograms[0]+histograms[1]
            hist_xy = np.sum(overlapped_hist, axis=2)
            hist_xy = np.flip(hist_xy, axis=0)
            hist_xy = np.transpose(hist_xy)
            ax[mode_idx][idx].imshow(hist_xy)
            ax[mode_idx][idx].set_xticks(ticks)
            ax[mode_idx][idx].set_yticks(ticks)
            if idx==0:
                ax[mode_idx][idx].set_yticklabels(['{:.1f}'.format(label) for label in labels[::-1]])
            else:
                ax[mode_idx][idx].set_yticklabels([])

            if mode_idx==0:
                ax[mode_idx][idx].set_title(f"$N_H={NH}$", fontsize=14)
                ax[mode_idx][idx].set_xticklabels([])
            else:
                ax[mode_idx][idx].set_xticklabels(['{:.1f}'.format(label) for label in labels])

            #fig.colorbar(imgcur, ax=ax[idx], format='%.0e')
    ax[0][0].set_ylabel("L-BFGS", fontsize=14)
    ax[1][0].set_ylabel("MOMENTUM GD", fontsize=14)
    plt.savefig("onebody_all.pdf", bbox_inches='tight')
    plt.show()

def plot_hist():
    ticks = np.linspace(0, 100, 4)
    labels = np.linspace(-1.5, 1.5, 4)
    fig, ax = plt.subplots(1,1,figsize=(9,4))
    histograms = get_histogram(N=2)
    overlapped_hist = histograms[0]+histograms[1]
    hist_xy = np.sum(overlapped_hist, axis=2)
    hist_xy = np.flip(hist_xy, axis=0)
    hist_xy = np.transpose(hist_xy)
    ax.imshow(hist_xy)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_yticklabels(['{:.1f}'.format(label) for label in labels[::-1]])
    ax.set_yticklabels([])
    ax.set_title(f"$N_H={2}$", fontsize=14)
    ax.set_xticklabels([])
    ax.set_xticklabels(['{:.1f}'.format(label) for label in labels])

    #fig.colorbar(imgcur, ax=ax[idx], format='%.0e')
    plt.savefig("onebody_one.pdf", bbox_inches='tight')
    plt.show()

def plot_hist_seeds(seeds=["31415", "92653", "58979", "32384"]):
    ticks = np.linspace(0, 100, 4)
    labels = np.linspace(-1.5, 1.5, 4)
    fig, ax = plt.subplots(1,4,figsize=(9,4))
    print(ax)
    for i in range(0,len(seeds)):
        #fig, ax = plt.subplots(2,2,i)
        seed=seeds[i]
        folder="./Outputs/"+seed
        histograms = get_histogram(N=2, folder=folder)
        ax[i].set_title("seed="+seed, fontsize=14)
        overlapped_hist = histograms[0]+histograms[1]
        hist_xy = np.sum(overlapped_hist, axis=2)
        hist_xy = np.flip(hist_xy, axis=0)
        hist_xy = np.transpose(hist_xy)
        ax[i].imshow(hist_xy)
        ax[i].set_xticks(ticks)
        ax[i].set_yticks(ticks)
        if i==0:
            ax[i].set_yticklabels(['{:.1f}'.format(label) for label in labels[::-1]])
        else:
            ax[i].set_yticklabels([])
        #ax[i].set_yticklabels(['{:.1f}'.format(label) for label in labels[::-1]])
        #ax[i].set_yticklabels([])
        #ax.set_title(f"$N_H={2}$", fontsize=14)
        ax[i].set_xticklabels([])
        ax[i].set_xticklabels(['{:.1f}'.format(label) for label in labels])

    #fig.colorbar(imgcur, ax=ax[idx], format='%.0e')
    plt.savefig("onebody_seeds.pdf", bbox_inches='tight')
    plt.show()

def plot_hist_seeds_indiv(seeds=["31415", "92653", "58979", "32384"]):
    ticks = np.linspace(0, 100, 4)
    labels = np.linspace(-1.5, 1.5, 4)
    fig, ax = plt.subplots(3,4,figsize=(9,4))
    print(ax)
    for i in range(0,len(seeds)):
        #fig, ax = plt.subplots(2,2,i)
        seed=seeds[i]
        folder="./Outputs/"+seed
        histograms = get_histogram(N=2, folder=folder)
        ax[0][i].set_title("seed="+seed, fontsize=14)
        #overlapped_hist = histograms[0]+histograms[1]
        hist_xy = np.sum(histograms[0], axis=2)
        hist_xy = np.flip(hist_xy, axis=0)
        hist_xy = np.transpose(hist_xy)
        ax[0][i].imshow(hist_xy)
        ax[0][i].set_xticks([])
        ax[0][i].set_yticks(ticks)
        if i==0:
            ax[0][i].set_yticklabels(['{:.1f}'.format(label) for label in labels[::-1]])
        else:
            ax[0][i].set_yticklabels([])
        #ax[i].set_yticklabels(['{:.1f}'.format(label) for label in labels[::-1]])
        #ax[i].set_yticklabels([])
        #ax.set_title(f"$N_H={2}$", fontsize=14)
        #ax[0][i].set_xticklabels([])
        #ax[0][i].set_xticklabels(['{:.1f}'.format(label) for label in labels])

        hist_xy = np.sum(histograms[1], axis=2)
        hist_xy = np.flip(hist_xy, axis=0)
        hist_xy = np.transpose(hist_xy)
        ax[1][i].imshow(hist_xy)
        ax[1][i].set_xticks([])
        ax[1][i].set_yticks(ticks)
        if i==0:
            ax[1][i].set_yticklabels(['{:.1f}'.format(label) for label in labels[::-1]])
        else:
            ax[1][i].set_yticklabels([])
        #ax[i].set_yticklabels(['{:.1f}'.format(label) for label in labels[::-1]])
        #ax[i].set_yticklabels([])
        #ax.set_title(f"$N_H={2}$", fontsize=14)
        #ax[1][i].set_xticklabels([])
        #ax[1][i].set_xticklabels(['{:.1f}'.format(label) for label in labels])
        
        overlapped_hist = histograms[0]+histograms[1]
        hist_xy = np.sum(overlapped_hist, axis=2)
        hist_xy = np.flip(hist_xy, axis=0)
        hist_xy = np.transpose(hist_xy)
        ax[2][i].imshow(hist_xy)
        ax[2][i].set_xticks(ticks)
        ax[2][i].set_yticks(ticks)
        if i==0:
            ax[2][i].set_yticklabels(['{:.1f}'.format(label) for label in labels[::-1]])
        else:
            ax[2][i].set_yticklabels([])
        #ax[i].set_yticklabels(['{:.1f}'.format(label) for label in labels[::-1]])
        #ax[i].set_yticklabels([])
        #ax.set_title(f"$N_H={2}$", fontsize=14)
        ax[2][i].set_xticklabels([])
        ax[2][i].set_xticklabels(['{:.1f}'.format(label) for label in labels])

    #fig.colorbar(imgcur, ax=ax[idx], format='%.0e')
    plt.savefig("onebody_seeds_individuals.pdf", bbox_inches='tight')
    plt.show()

def gaussian_fit(n, alphaopt, betaopt):

    #get histograms along xy plane and z direction
    HIST = get_histogram(n)
    hist_xy = np.sum(HIST , axis=2)
    hist_z = np.sum(HIST, axis=(0,1))
    #flatten to perform linear fit
    hist_xy_flat = hist_xy.reshape(-1)
    #avoid log(0) and log of small numbers
    nnz_idx = np.nonzero(hist_xy_flat>np.exp(7))
    nnz_idx_z = np.nonzero(hist_z>np.exp(5))
    log_hist = np.log(hist_xy_flat[nnz_idx])
    log_hist_z = np.log(hist_z[nnz_idx_z])
    #need the distance from each slot to the center
    size = 100
    min, max = -1.5, 1.5
    # Create a 2D array with x and y coordinates of each slot
    x_coords =  np.linspace(min, max, size)
    y_coords =  np.linspace(min, max, size)
    zz =        np.linspace(min, max, size)
    xx, yy =    np.meshgrid(x_coords, y_coords)
    # Compute the distance of each slot from the origin 0,0
    r2 = (xx)**2 + (yy)**2
    r2 = r2.reshape(-1)
    r2 = r2[nnz_idx]
    z2 = zz[nnz_idx_z]**2
    
    #plt.plot(z2, log_hist_z, linestyle='none', marker='.')
    #fit the gaussian
    coeffs, var = np.polyfit(r2, log_hist, deg=1, cov=True )
    alpha = -coeffs[0]/2
    coeff_z, varz = np.polyfit(z2, log_hist_z, deg=1, cov=True)
    beta = -coeff_z[0]/(2*alpha)
    varbeta = np.sqrt(varz[0,0])/(2*alpha)
    print(f"${n}$ & $\{alpha:.3f}$ &  \
          ${alphaopt:.5f}$ & \
          ${beta:.4f} $  &\
          ${betaopt:.4f}$  & \
          ${alpha*beta:.4f}$ & \
          ${alphaopt*betaopt:.4f}$  \\\\ ")
    #plt.show()
    #get histogram along z axis:
#hist_z = np.sum(HIST_TOT, axis=(0,1))

if __name__ == "__main__":
    plot_hist_seeds_indiv()