using Makie
using GeometryTypes, AbstractPlotting, Colors
import Observables
using Observables: AbstractObservable, on, off, async_latest
using Reactive

const BASE_FPS = 60
### GUI

s1, g_L = AbstractPlotting.textslider([2^n for n in 2:9], "L")
s2, g_mult = AbstractPlotting.textslider(-4.0:0.1:+2.0, "multiplier [log]")

s3, g_h = AbstractPlotting.textslider(-5f0:0.1f0:5f0, start=0f0,"field")
s4, g_T = AbstractPlotting.textslider(0f0:0.01f0:10f0, "temperature")

b_1 = AbstractPlotting.button(Theme(raw = true, camera = campixel!),"Pause/Play" )
b1_click = on(b_1[end][:clicks]) do c
    push!(run_signal, !run_signal.value)
end

### Logic
include("src/mcmc.jl");
include("src/visualization.jl")

# Timey-Wimey
frame_node = Node(0)
run_signal = Signal(true)
fps_sig = Reactive.fpswhen(run_signal,BASE_FPS);
tic_sig = map(fps_sig) do tic
    frame_node[] = (frame_node[]+1)%60
end


L = async_latest(g_L,1) #
g_L[] = 16


positions, square = get_primitives(L[])
config0 = frustratedConfiguration(L[]);
cluster = zeros(Bool, L[],L[]);

color_node = Node([RGB(0.,0.,0.) for j in 1:L[]^2])

init_map = lift(L) do L
    push!(run_signal, false)

    global positions,square = get_primitives(L)
    global config0 = frustratedConfiguration(L);
    global cluster = zeros(Bool, L,L);

    if @isdefined state_plot
        @info("Trying to close color signal")
        plt = state_plot.plots[1]
        plt.attributes[:color] |> close
    end
    global color_node[] = reshape(color_gen(config0,0f0),L^2)

    if @isdefined color_signal
        try
            empty!(frame_node.listeners)
            color_signal = nothing
        catch err
            @warn("Could not remove color_signal from frame_node")
        end
    end
    global color_signal = on(frame_node) do tic
            sweep!(config0, round(Int,10^g_mult.val * L^2), 1/g_T.val,g_h.val)
            color_node[] = reshape(color_gen(config0,0f0),L^2)
            nothing
    end
    global state_plot = Makie.meshscatter(reshape(positions,L^2); color=color_node,
        markersize=1,aspect_ratio=1,marker=GLNormalMesh(square), axis_type=axis3d!, camera=cam3d!,
        size=(800,800), raw=true
        )
    global scene = vbox(hbox(s1,s2,s3,s4,b_1),state_plot, parent=Scene(resolution=(1024,1024)))
    adjust_cam!(state_plot)
    display(scene);
    # push!(run_signal, true)
    nothing
end
#
