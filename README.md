# Ising

This package implements Markov-chain-Monte-Carlo methods (MCMC) to study the equilibrium thermodynamic behavior of the Ising model on a square lattice. Currently it implements the classical Metropolis algorithm as well as Wolff's cluster algorithm.

You'll find a source file with all the simulation routines and a notebook that ties them together into parametric studies, produces nice plots and so on.  
Also included is a notebook for visualizing the dynamics See __"Live visualization"__ below.

It is mainly intended for students of statistical physics learning about phase transitions and critical phenomena.

For now the simulation is restricted to a square lattice with constant couplings and external fields. It should however be quite easy to extend the code to cover more general systems.

__If you spot an error or are missing a feature, please feel free to open an issue or a pull request!__

## Installation

### Compatibility
Julia >= 0.6 is required for the simulation and data analysis portion of the notebook.

At the moment the live-visualization requires Julia 0.6 and is not (yet) compatible with Julia 0.7.

### Setup

Clone the files to your computer
```bash
  git clone https://github.com/skleinbo/Ising.git
```

The simulation routines assembled in `src/mcmc.jl` use base Julia only.

If you want to follow along the analysis in the accompanying notebook, a few packages need to be installed.  
Run `src/build.jl` once to install any missing packages.

## Usage

Look around the source code to see which methods are available. Here's a minimal example to get you going:

```julia
julia> include("src/mcmc.jl");

julia> L=128; # _Linear_ system size

help> run_metropolis
[...]

julia> run_metropolis(L, 1/2.26, 0.; Tmax=25*10^3*L^2,sample_interval=10*L^2,sweep=10^3*L^2)
5-element Array{Float64,1}:
  -253.3488
 64566.7776
     0.002934521484375
     9.098881483078003e-6
     1.0010325963705214e-10
```

Admittedly the notebook is not that simple to understand if you hadn't had much exposure to Julia and DataFrames. Write your own routines!

### Live Visualization
<center><img src="img/window.png" width=300></img></center>

A few more packages are need to give us interactivity and render a scene. They will be installed automatically once you load `src/visualization.jl`.

If you execute the cells at the end of the notebook in order, two things should happen. You should be presented with three sliders inside the notebook that allow you to control temperature, field strength, and simulation speed. A window should open in which the system is presented as a checker board.
