# Variational Monte Carlo solver for a Bose-Einstein condensate

This repository contains the code for the first project of the course [FYS4411 (spring 2023)](https://github.com/CompPhysics/ComputationalPhysics2).
This was forked from this template [repository](https://github.com/mortele/variational-monte-carlo-fys4411.git), whose authors are thankfully acknowledged. 
 

### Compiling and running the project
Provided Cmake is installed, compile by running the script `compile_project` via
```bash
./compile_project
```
This produces an executable called `vmc` in the top directory. It should be run as:
```bash
./vmc Input/simulation_input.txt
```
where `Input/simulation_input.txt` is the file containing the parameters of the simulation.
There are two ways of running the code. 

1. If the input file contains the line `calculateGradients=1`, the code will compute the gradients and run an optimization algorithm for the parameters of the wave function. By default this is LBFGS, but it can be set to be gradient descent with momentum by providing an extra command-line argument `GD`:
    ```bash
        ./vmc Input/simulation_input.txt GD
    ```
    Progress of the optimization is written to `Outputs/output.txt`.

2. If the input file contains the line `calculateGradients=0`, the code will __NOT__ compute the gradients. This is meant to spare the computational costs when the wave function is already optimized. It will output a lot of files in the `Outputs` folder, which contain the sampled energies and the one-body histograms (in binary format).

## Authors
[Gianmarco Puleo](https://github.com/giammy00) <br>
[David Svejda](https://github.com/DavidSvejda2507)<br>
