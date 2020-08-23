using GeometryBasics, Colors
import Observables
import Observables: AbstractObservable, on, off, async_latest

import GLMakie
using AbstractPlotting
using AbstractPlotting.MakieLayout


include("../src/mcmc.jl");
include("../src/visualization.jl")

### GUI

scene, layout = layoutscene(resolution=(768,1024))

state_plot = layout[1,1:2] = LScene(scene, raw=true, camera = cam2d!)
sub_plot_energy = layout[3,1] = LAxis(scene, title="Energy")
sub_plot_mag = layout[3,2] = LAxis(scene, title="Magnetization")

g_T = labelslider!(scene, "temp", 0f0:0.005f0:5f0)
T = g_T.slider.value
g_h = labelslider!(scene, "field", -5f0:0.1f0:5f0)
h = g_h.slider.value
g_mult = labelslider!(scene, "speed", 0.0:0.05:1.0, sliderkw=Dict(:startvalue=>0.5))
mult = g_mult.slider.value
sub_controls = layout[2,:] = GridLayout()
sub_controls[1,1] = g_T.layout
sub_controls[1,2] = g_h.layout
sub_controls[2,1:2] = g_mult.layout

rowsize!(layout, 3, Relative(1/4))

color1 = colorant"black"
color2 = colorant"cornflowerblue"


### Logic

algorithm = Node(:metropolis) # either :metropolis or :wolff

# Timey-Wimey
frame_node = Node(0)
run_signal = Node(false)
BASE_FPS = Node(60)
BUFFER_LENGTH = 1024 # no. of observations stored



L = Node(128) # linear system size

pos_and_marker = lift(get_primitives, L)
positions = lift(x->reshape(x[1], length(x[1])), pos_and_marker)
square = lift(x->x[2], pos_and_marker)
config0 = lift(frustratedConfiguration, L)
cluster = lift(L->zeros(Bool, L,L), L)

color_node = lift(L->[color1 for j in 1:L^2], L)

size_change = lift(L) do L
    run_signal[] = false

    update_cam!(state_plot.scene, FRect(-1,-1,2,2))
    display(scene);

    nothing
end

meshscatter!(state_plot, positions; color=color_node,
    markersize=1, camera=cam2d!, marker=square,
    raw=true, shading=false
    )
state_plot.scene[end][:light] = Vec{3,Float32}[[1.0, 1.0, 1.0], [0.1, 0.1, 0.1], [0.9, 0.9, 0.9], [0.0, 0.0, -20.0]]

## Events

# Keys:
#   - r: reset state and plots
#   - c: reset camera
#   - +/- (QWERTY): adjust simulation speed
on(scene.events.keyboardbuttons) do buttons
    if ispressed(scene, Keyboard.p)
        run_signal[] = !run_signal[]
    elseif ispressed(scene, Keyboard.c)
        update_cam!(state_plot.scene, FRect(-1,-1,2,2))
    elseif ispressed(scene, Keyboard.r)
        config0[] = frustratedConfiguration(L[]);
        cluster[] = zeros(Bool, L[],L[]);
        empty!(time_buffer[])
        empty!(E_buffer[])
        empty!(M_buffer[])
    elseif ispressed(scene, Keyboard.equal)
        mult[] = min(mult[]+0.05, 1.0)
    elseif ispressed(scene, Keyboard.minus)
        mult[] = max(mult[]-0.05, 0.0)
    end
end
on(scene.events.window_area) do wa
    update_cam!(state_plot.scene, FRect(-1,-1,2,2))
end

## Physical observables
time_buffer = Node([0.0])
E_buffer = Node([Point2(0.0,0.0)]) # pairs of (time, value)
M_buffer = Node([Point2(0.0,0.0)])

E_scatter = lines!(sub_plot_energy, E_buffer)
M_scatter = lines!(sub_plot_mag, M_buffer)

# recalculate the plot limits when observables change.
# MakieLayout doesn't  adjust them automatically when
# new points get pushed to the plot.
plot_limits = Observables.onany(E_buffer, M_buffer) do a,b
    MakieLayout.autolimits!(sub_plot_mag)
    MakieLayout.autolimits!(sub_plot_energy)
    nothing
end

## Render- and simulation loop

render_task = @async begin
    while scene.events.window_open[]
        t1 = time()
        if run_signal[] && mult[] > 0.0
            if algorithm[] == :metropolis
                metropolis_sweep!(config0[], round(Int, mult[] * L[]^2), 1/T[], h[])
            elseif algorithm[] == :wolff
                wolff_sweep!(round(Int, mult[] * L[]^2), config0[], cluster[], 1/T[], h[])
            end
            color_node[] = reshape(color_gen(config0[], color1, color2),L[]^2)
            # Measure observables
            # and push them to the relevant buffers.
            # Prune buffers if necessary.
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
            Observables.notify!(E_buffer)
            Observables.notify!(M_buffer)
        end
        t2 = time()
        if (t2-t1) < 1/BASE_FPS[]
            sleep(1/BASE_FPS[] - (t2-t1))
        end
    end
    :ended
end

display(scene)
@__DIR__#
