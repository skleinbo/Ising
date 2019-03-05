using Makie
using GeometryTypes, AbstractPlotting, Colors
import Observables
using Observables: AbstractObservable, on, off, async_latest

include("src/mcmc.jl");
include("src/visualization.jl")

### GUI
color1 = colorant"black"
color2 = colorant"cornflowerblue"
# s1, g_L = AbstractPlotting.textslider(16:16:512, "L")
s2, g_mult = AbstractPlotting.textslider(-3.0:0.1:+2.0, start=-2.0, "speed [log]")

s3, g_h = AbstractPlotting.textslider(-5f0:0.1f0:5f0, start=0f0,"field")
s4, g_T = AbstractPlotting.textslider(0f0:0.005f0:5f0, "temperature")

b1 = AbstractPlotting.button(Theme(raw = true, camera = campixel!),"Pause/Play" )
b2 = AbstractPlotting.button(Theme(raw = true, camera = campixel!),"Metropolis" )

b1_click = on(b1[end][:clicks]) do c
    run_signal[] = !run_signal[]
end

b2_click = on(b2[end][:clicks]) do c
    if algorithm[] == :metropolis
        algorithm[] = :wolff
        b2[end].input_args[end][] = "Wolff"
    else
        algorithm[] = :metropolis
        b2[end].input_args[end][] = "Metropolis"
    end
end

### Logic

algorithm = Node(:metropolis)

# Timey-Wimey
frame_node = Node(0)
run_signal = Node(false)
BASE_FPS = Node(60)


g_L = Node(128)
# L = async_latest(g_L,1) #


positions, square = get_primitives(g_L[])
config0 = frustratedConfiguration(g_L[]);
cluster = zeros(Bool, g_L[],g_L[]);

color_node = Node([RGB(0.,0.,0.) for j in 1:g_L[]^2])

render_map = lift(run_signal) do rs
    global simtask = @async begin
        if !rs
            @info "Render task ended."
            return nothing
        end
        while run_signal[]
            t1 = time()
            if algorithm[] == :metropolis
                metropolis_sweep!(config0, round(Int,10^g_mult[] * g_L[]^2), 1/g_T.val,g_h.val)
            elseif algorithm[] == :wolff
                wolff_sweep!(round(Int,10^g_mult[] * g_L[]^2), config0, cluster, 1/g_T[], g_h[])
            end
            color_node[] = reshape(color_gen(config0, color1, color2),g_L[]^2)
            t2 = time()
            if (t2-t1) < 1/BASE_FPS[]
                sleep(1/BASE_FPS[] - (t2-t1))
            end
            nothing
        end
    end
end

init_map = lift(g_L) do L
    run_signal[] = false

    global positions,square = get_primitives(L)
    global config0 = frustratedConfiguration(L);
    global cluster = zeros(Bool, L,L);

    empty!(color_node.listeners)
    color_node[] = reshape(color_gen(config0,color1, color2),L^2)

    global state_plot = Makie.meshscatter(reshape(positions,L^2); color=color_node,
        markersize=1,marker=GLNormalMesh(square), axis_type=axis2d!, camera=cam2d!,
        raw=true, shading=false
        )
    # state_plot.attributes[:padding][] = [0f0, 0f0, 0f0]
    # cam = cameracontrols(state_plot)
    # cam.area[] = HyperRectangle(-1.05f0, -1.05f0, 2.1f0,2.1f0)
    state_plot.theme.attributes[:backgroundcolor][] = :gray
    update_cam!(state_plot, FRect(-1,-1,2,2))
    state_plot[end][:light] = Vec{3,Float32}[[1.0, 1.0, 1.0], [0.1, 0.1, 0.1], [0.9, 0.9, 0.9], [0.0, 0.0, -20.0]]

    global scene = vbox(hbox(s2,s3,s4,b2,b1),state_plot, parent=Scene(windowtitle="Ising"))
    on(scene.events.keyboardbuttons) do buttons
        if ispressed(scene, Keyboard.p)
            run_signal[] = !run_signal[]
        elseif ispressed(scene, Keyboard.c)
            update_cam!(state_plot, FRect(-1,-1,2,2))
        elseif ispressed(scene, Keyboard.r)
            config0 = frustratedConfiguration(g_L[]);
            cluster = zeros(Bool, g_L[],g_L[]);
        elseif ispressed(scene, Keyboard.equal)
            g_mult[] = min(g_mult[]+0.1, 2.0)
        elseif ispressed(scene, Keyboard.minus)
            g_mult[] = max(g_mult[]-0.1, -3.0)
        end
    end
    on(scene.events.window_area) do wa
        update_cam!(state_plot, FRect(-1,-1,2,2))
    end
    display(scene);
    # push!(run_signal, true)
    nothing
end
#
