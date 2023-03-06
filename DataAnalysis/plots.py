import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
plt.rcParams['text.usetex'] = True
#then define a function which plots only the data as function of alpha for a given dimension and returns an axes object 
def single_plot_E_alpha(data, ax=None, Ndim=3):
    ''' 
    Plot MonteCarlo energy estimates for various number of particles, with errobars corresponding to standard deviation.
    args:
        data: pd.DataFrame containing data, imported from output.txt.
        Ndim: plot only results corresponding to Ndim dimensions, should be 1, 2 or 3.
        ax: plt.Axes where I want to plot the results
    
    NB: assume there is only one parameter in the wavefunction (for now), don't know what happens when more parameters are given.'''

    current_data_dim = data[data['n_dimensions']==Ndim]

    nrs_particles = pd.unique(current_data_dim['#n_particles'])
    print(nrs_particles)
    for Nparticles in nrs_particles:
        print(Nparticles)
        #select only data related to current number of particles.
        current_data = current_data_dim[current_data_dim['#n_particles']==Nparticles]

        alphas = current_data['params']
        energy_pts = current_data['E']
        std_pts = np.sqrt(current_data['var'])
        ax.errorbar(alphas, energy_pts, std_pts, label=f"$N={Nparticles}$", linestyle='None', marker=".", capsize=5)
    ax.set_xlabel("$\\alpha$")
    ax.set_ylabel("trial energy")
    ax.legend()
    ax.grid()
    #ax.set_yscale('log')
    ax.set_title(f"{Ndim} dimensions")
    return

def plot_E_alpha_Ndim(data, dims=[3]):
    '''Make three subplots in a single figure, each plot for a different number of dimensions.'''
    fig, axs = plt.subplots(1,3)
    for i in dims:
        single_plot_E_alpha(data, ax=axs[i-1], Ndim=i) 

    return fig, axs
if __name__ == '__main__':
    # import output.txt from simulation as a Pandas dataframe
    data = pd.read_csv("Outputs/output.txt", sep="\t")
    dims = [1,2,3]
    fig, axs = plot_E_alpha_Ndim(data, dims=dims)
    plt.show()