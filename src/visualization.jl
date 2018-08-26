for package in ["GLVisualize", "GeometryTypes", "GLAbstraction", "IterTools",
                "Colors","Reactive","Interact"]
    try
        sp = Symbol(package)
        @eval using $sp
    catch err
        info("Package $package not installed. Trying to...")
        Pkg.add(package)
    end
end


# This function maps the spin state to an array of colors
function color_gen(v0,basecolor)
    map(v0) do x
        if x==1
            RGB(0f0,0f0,0f0)
        elseif x==-1.
            RGB(1f0,min(1f0,Float32(basecolor)),0f0)
        end
    end
end

# Reset window
function reset_window()
try
    empty!(window)
    close(color_signal)
    close(state_map)
    close(timesignal)
    close(temperature)
catch UndefVarError
end
end

function adjust_cam!(window;eyepos_vec=Vec3f0(0,0,+1),lookat_vec=Vec3f0(0,0,0),
    up_vec=cross(lookat_vec-eyepos_vec,-Vec3f0(1,0,0))
    )

    push!(window.cameras[:perspective].eyeposition, eyepos_vec)
    push!(window.cameras[:perspective].lookat, lookat_vec)
    push!(window.cameras[:perspective].up, up_vec)
    push!(window.cameras[:perspective].fov, 90)
end
