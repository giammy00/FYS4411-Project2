# Variational Monte Carlo solver for a Bose-Einstein condensate
![Languages](https://img.shields.io/badge/Languages-C%2B%2B%20%7C%20Python-red)
![Size](https://img.shields.io/github/repo-size/giammy00/FYS4411-Project2)

This repository contains the code for the second project of the course [FYS4411 (spring 2023)](https://github.com/CompPhysics/ComputationalPhysics2).
This was forked from this template [repository](https://github.com/DavidSvejda2507/FYS4411-Project1.git), whose contributors are thankfully acknowledged. 
 
### Requirements
The code is written in C++. To compile it, you need to have the [nlopt library](https://nlopt.readthedocs.io/en/latest/) and [OpenMP](https://www.openmp.org/) installed. The data analysis is written in python.

### Compiling and running the project
Provided Cmake is installed, compile by running the script `compile_project` via
```bash
./compile_project
```
This produces an executable called `vmc` in the top directory. It should be run as:
```bash
./vmc Input/rbm_input.txt
```
where `Input/rbm_input.txt` is the file containing the parameters of the simulation.
There are two ways of running the code. 

1. If the input file contains the line `calculateGradients=1`, the code will compute the gradients and run an optimization algorithm for the parameters of the wave function. By default this is L-BFGS, but it can be set to be gradient descent with momentum by providing an extra command-line argument `GD`:
    ```bash
        ./vmc Input/simulation_input.txt GD
    ```
    Progress of the optimization is written to `Outputs/output.txt`, and the optimal parameters are saved in binary format to `Outputs/optimizedRBMParams.bin`.

2. If the input file contains the line `calculateGradients=0`, the code will __NOT__ compute the gradients. This is meant to spare the computational costs when the wave function is already optimized. It will import the parameters stored by the optimization run, and output a lot of files in the `Outputs` folder, which contain the sampled energies, distances, and the one-body histograms (in binary format). __Note that the output of this run can be very large in terms of memory and number of files.__

### Data Analysis
The python scripts for the data analysis are contained in the `DataAnalysis` folder, they assume that a simulation has been run without calculating gradients.

* `onebodydensity.py` computes and plots the histograms and the corresponding gaussian fit.
* `blocking.py` contains the blocking function provided by [Marius Jonsson](https://github.com/computative/block) and uses it to analyse the simulation data.
* `distance_hist.py` analyses the sampled relative distances to produce a histogram. 

### Contributors

 [Carl Petter Duehdal](https://github.com/carlpd) 

[Gianmarco Puleo](https://github.com/giammy00) 

[Ivan Scolan](https://github.com/ivanscolan)