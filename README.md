# Ising

This package implements Markov-chain-Monte-Carlo methods (MCMC) to study the equilibrium thermodynamic behavior of the Ising model on a square lattice. Currently it implements the classical Metropolis algorithm as well as Wolff's cluster algorithm.

You will find a source file with all the simulation routines and a notebook (`extras/MCMC_Ising.ipynb`) that ties them together into parametric studies, produces nice plots and so on.

To watch the algorithms in action, a visualization app based on [Makie](http://makie.juliaplots.org/stable/index.html) is provided. ![visualization](img/window.png)

For now the simulation is restricted to a square lattice with constant couplings and external fields. It should however be fairly straightforward to extend the code to cover more general systems.

__If you spot an error or miss a feature, please feel free to open an issue or a pull request!__

## Compatibility and Installation

Julia >= 1.3 is required.

Clone the repository to your computer
`git clone https://github.com/skleinbo/Ising.git`
or download the latest release.

The simulation routines assembled in `src/mcmc.jl` use base Julia only. [TODO: Wrap them in a module]

The interactive application requires [Makie](https://github.com/JuliaPlots/Makie.jl) and a few other packages (see `Project.toml`). Start Julia from the root of the cloned directory with `julia --project`. The first time you use the app, install the dependencies from the REPL with `] instantiate`. Then start the app via `include("bin/main.jl")`.

The analysis notebook utilizes various additional packages which are not installed by default. Install them once with
`]add Plots LaTeXStrings StatsPlots GLM DataFrames DataFramesMeta CSV Distributions`


## Usage
Load the MCMC methods with `include(src/mcmc.jl)`.

Help is available
```julia
help?> run_metropolis
search: run_metropolis _run_metropolis!

  run_metropolis(L, beta, h;Tmax=1, sweep=0, sample_interval=1)

  Sets up a random state and runs the Metropolis algorithm for a given set of parameters. Samples in defined intervals along the Markov-Chain.

  An initial thermal sweep to go to equilibrium may be specified.

  Returns an array of averaged observables: [E, E^2, m, m^2, m^4] with total energy E and magnetisation per spin m.

  Arguments
  ≡≡≡≡≡≡≡≡≡≡≡

    •    L::Integer: Linear system size

    •    beta::Float: Inverse temperature

    •    h::Float: External field

    •    Tmax::Integer: Number of steps

    •    sweep::Integer: Length of the initial sweep

    •    sample_interval::Integer: sample interval

  Example
  ≡≡≡≡≡≡≡≡≡

  julia> run_metropolis(50, 0., 0.;Tmax=50*10^3*50^2,sample_interval=10*50^2,sweep=10^3*50^2)
  5-element Array{Float64,1}:
    -39.9664      # <E>
   1753.5168      # <E^2>
      0.010382079999999908 # <m>
      0.00011572633600000021 # <m^2>
      1.6924325969920386e-8 # <m^4>
```

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

Admittedly the notebook may not be simple to understand if you hadn't had much exposure to Julia and DataFrames. And even if, it's not well documentes. __Write your own routines!__

### Live Visualization
`include("Visual_Ising.jl")` opens a window with sliders and buttons to adjust
parameters, a depiction of the current state, as well as energy and magnetization over time. Zoom in
and out with the mouse wheel.
"Speed" determines how many sweeps the program does per frame, i.e a value of `1/2` means
`1/2*L^2` attempts to flip a spin/cluster (Metropolis/Wolff) each frame (default: 60fps; change with `BASE_FPS[]=$fps`).

__Warning:__ Choosing this value to high will freeze the application, in particular when the system size is large too.

__Controls:__
  * To change the system size (default 128x128): `L[]=_linear system size_`
  * Adjust frame rate: `BASE_FPS[] = _fps_`
  * To change colors of the two states, set `color1` and `color2`
  * Press <kbd>p</kbd> to run/pause the simulation
  * Reset the state by pressing <kbd>r</kbd>
  * Speed up/slow down with <kbd>=</kbd> / <kbd>-</kbd> (QWERTY layout)
  * Reset the view with <kbd>c</kbd>
