import Base.Iterators: product
import LinearAlgebra: cross

# This function maps the spin state to an array of colors
function color_gen(v0, color1, color2)
    map(v0) do x
        if x==1
            color1
        elseif x==-1
            color2
        end
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

# Returns a tuple of LxL positions on a square lattice, and a HyperRectangle of
# appropriate size
function get_primitives(L::Integer)
    return (Point2f0[Point2f0(2*xi/L-1.,2*yi/L-1.) for (xi,yi) in product(0:L-1,0:L-1)],
     GeometryBasics.HyperRectangle(0,0,0,2/L,2/L,0)
     )
 end
