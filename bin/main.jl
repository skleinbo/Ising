# import Pkg; Pkg.instantiate();

using GeometryBasics, Colors
import Observables
import Observables: on, off, notify

import GLMakie
import GLMakie: GLFW

include("../src/mcmc.jl");
include("../src/visualization.jl")

### GUI

fig = Figure(resolution=(768,1024))

light = AmbientLight(colorant"white")

state_plot = LScene(fig[1,1:2], show_axis=false, scenekw=(raw=true, light=[light], camera = cam2d!, autolimitaspect=1.0))
sub_plot_energy = Axis(fig[3,1], title="Energy")
sub_plot_mag = Axis(fig[3,2], title="Magnetization")

sub_controls = SliderGrid(fig[2,1], 
    (label="temp",  range=0f0:0.005f0:5f0),
    (label="field", range=-5f0:0.1f0:5f0),
    (label="speed", range=0.0:0.05:1.0, :startvalue=>0.5)
)
T = sub_controls.sliders[1].value
h = sub_controls.sliders[2].value
mult = sub_controls.sliders[3].value

rowsize!(fig.layout, 3, Relative(1/6))

color1 = colorant"black"
color2 = colorant"cornflowerblue"


### Logic

algorithm = Observable(:metropolis) # either :metropolis or :wolff

# Timey-Wimey
frame_node = Observable(0)
run_signal = Observable(false)
BASE_FPS = Observable(60)
BUFFER_LENGTH = 1024 # no. of stored observations

L = Observable(128) # linear system size

pos_and_marker = lift(get_primitives, L)
positions = lift(x->reshape(x[1], length(x[1])), pos_and_marker)
square = lift(x->x[2], pos_and_marker)
config0 = lift(frustratedConfiguration, L)
cluster = lift(L->zeros(Bool, L,L), L)

color_node = lift(L->[color1 for j in 1:L^2], L)

size_change = lift(L) do L
    run_signal[] = false

    update_cam!(state_plot.scene, FRect(-1,-1,2,2))
    display(fig);

    nothing
end

color_update = lift(config0) do config0
    # TODO: recalculating all colors is wasteful
    # -> return flipped spins -> adjust color
    color_node[] = reshape(color_gen(config0, color1, color2),L[]^2)
end


# Plot of the lattice
meshscatter!(state_plot, positions; color=color_node,
    markersize=1, marker=square,
    # scene=(camera=cam2d!, raw=true, shading=false),
    )

## Events

toggle_run!() = run_signal[] = !run_signal[]
toggle_algo!() = algorithm[] = algorithm[] == :metropolis ? :wolff : :metropolis

function set_title!()
    title_string = "Ising - $(titlecase(string(algorithm[])))"
    if mult[] == 0.0 || !run_signal[]
        title_string *= " - PAUSED"
    end
    GLFW.SetWindowTitle(fig.scene.current_screens[1].glscreen, title_string)
end

# Keys:
#   - r: reset state and plots
#   - c: reset camera
#   - +/- (QWERTY): adjust simulation speed
#   - s: short sweep: 1 time
#   - l: long sweep 1_000 times
#   - a: switch algorithm Metropolis/Wolff
on(fig.scene.events.keyboardbutton) do buttons
    if ispressed(fig.scene, Keyboard.p)
        # if the renderloop ended for some reason
        # it is restarted here
        if istaskdone(render_task)
            global render_task = renderloop()
        end
        toggle_run!()
    elseif ispressed(fig.scene, Keyboard.c)
        update_cam!(state_plot.scene, FRect(-1,-1,2,2))
    elseif ispressed(fig.scene, Keyboard.r)
        config0[] = frustratedConfiguration(L[]);
        cluster[] = zeros(Bool, L[],L[]);
        empty!(time_buffer[])
        empty!(E_buffer[])
        empty!(M_buffer[])
        notify(config0)
    elseif ispressed(fig.scene, Keyboard.s)
        run_signal[] = false
        _sweep!(config0, L[]^2)
        _measure!()
    elseif ispressed(fig.scene, Keyboard.l)
        run_signal[] = false
        _sweep!(config0, 1_000*L[]^2)
    elseif ispressed(fig.scene, Keyboard.equal)
        mult[] = min(mult[]+0.05, 1.0)
    elseif ispressed(fig.scene, Keyboard.minus)
        mult[] = max(mult[]-0.05, 0.0)
    elseif ispressed(fig.scene, Keyboard.a)
        toggle_algo!()
    end
    set_title!()
end
on(fig.scene.events.window_area) do wa
    update_cam!(state_plot.scene, FRect(-1,-1,2,2))
end

## Physical observables
time_buffer = Observable([0.0])
E_buffer = Observable([Point2(0.0,0.0)]) # pairs of (time, value)
M_buffer = Observable([Point2(0.0,0.0)])


E_plot = lines!(sub_plot_energy, E_buffer)
ylims!(sub_plot_energy, (-2.1, 0.1))
M_plot = lines!(sub_plot_mag, M_buffer)
ylims!(sub_plot_mag, (-1.1, 1.1))

# recalculate the plot limits when observables change.
# MakieLayout doesn't  adjust them automatically when
# new points get pushed to the plot.
plot_limits = Observables.on(E_buffer) do a
    trange = (E_buffer[][1][1], E_buffer[][end][1]+mult[])
    xlims!(sub_plot_energy, trange)
    xlims!(sub_plot_mag, trange)
    nothing
end

## Render- and simulation loop

function _sweep!(config0, n)
    if algorithm[] == :metropolis
        metropolis_sweep!(config0[], n, 1/T[], h[])
    elseif algorithm[] == :wolff
        wolff_sweep!(config0[], cluster[], ceil(Int, n/L[]^2), 1/T[], h[])
    end
    notify(config0)
    nothing
end

function _measure!()
    t = if length(time_buffer[])>0
            time_buffer[][end]
        else
             0.0
        end + mult[]
    push!(time_buffer[], t)
    push!(E_buffer[], Point2(t, H(config0[], 1.0, h[])/L[]^2))
    push!(M_buffer[], Point2(t, m(config0[])))
    if length(E_buffer[]) > BUFFER_LENGTH
        deleteat!(E_buffer[], 1)
        deleteat!(M_buffer[], 1)
        deleteat!(time_buffer[], 1)
    end
    Observables.notify(E_buffer)
    Observables.notify(M_buffer)
    nothing
end

function renderloop()
    @async begin
        while fig.scene.events.window_open[]
            t1 = time()
            if run_signal[] && mult[] > 0.0
                _sweep!(config0, round(Int, mult[] * L[]^2))
                # Measure observables
                # and push them to the relevant buffers.
                # Prune buffers if necessary.
                _measure!()
            end
            t2 = time()
            if (t2-t1) < 1/BASE_FPS[]
                sleep(1/BASE_FPS[] - (t2-t1))
            end
        end
        :ended
    end
end

render_task = renderloop()

display(fig)
