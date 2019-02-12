# for package in ["GLVisualize", "GeometryTypes", "GLAbstraction", "IterTools",
#                 "Colors","Reactive","Interact"]
#     try
#         sp = Symbol(package)
#         @eval using $sp
#     catch err
#         info("Package $package not installed. Trying to...")
#         Pkg.add(package)
#     end
# end
import Base.Iterators: product
import LinearAlgebra: cross

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

function adjust_cam!(scene::Scene;eyepos_vec=Vec3f0(0,0,+1),lookat_vec=Vec3f0(0,0,0),
    up_vec=cross(lookat_vec-eyepos_vec,-Vec3f0(1,0,0))
    )
    cam = cameracontrols(scene)
    cam.eyeposition[] =  eyepos_vec
    cam.lookat[] = lookat_vec
    cam.upvector[] = up_vec
    cam.fov[] = 90
    update_cam!(scene,cam)
end

function get_primitives(L::Integer)
    return (Point3f0[Point3f0(2*xi/L-1.,2*yi/L-1.,0.0) for (xi,yi) in product(0:L-1,0:L-1)],
     HyperRectangle(Vec3f0(0.),Vec3f0(2/L,2/L,2/L))
     )
 end
