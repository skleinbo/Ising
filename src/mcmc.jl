### Monte Carlo routines for
### the 2D Ising model

import DataStructures: Stack
import StatsBase: countmap
import StaticArrays: SVector, @SVector

"""
    randomConfiguration(L)

Returns a LxL Matrix filled randomly with ±1.
"""
function randomConfiguration(L::Int)
    return rand(Int8[-1,1], L,L)
end

"""
    frustratedConfiguration(L)

Returns a configuration where all bonds are frustrated.
"""
function frustratedConfiguration(L::Int)
    return Int8[(-1)^(i+j) for i=1:L,j=1:L]
end

"""
    H(state, J=1.,h=0.)

Energy of Ising configuration for a given coupling `J` and external field `h`.
"""
function H(state::Matrix{Int8},J=1.,h=0.)
    L = size(state,1)
    s = h*(sum(state))

    @inbounds for i in 1:L, j in 1:L
        s += -J*state[i,j]*(state[i%L+1,j]+state[i,j%L+1])
    end
    return s
end

"""
    dH(state,pos,J=1.,h=0.)

Energy difference between the current state `state`
and one spin flipped at position `(i,j)` given coupling and field strength.
"""
function dH(state::Matrix{Int8},i,j,J=1.,h=0.)
    L = size(state,1)
    @inbounds return 2*J*state[i,j]*( state[i%L+1,j]+state[i,j%L+1]+
        state[i==1 ? L : i-1,j]+state[i,j==1 ? L : j-1] ) + 2*h*state[i,j]
end

"""
    m(state)

Magnetization per spin
"""
m(state) = Float64(sum(state)/length(state))

@inline neighbors(i,j,L) = [(i%L+1,j),(i,j%L+1),(i==1 ? L : i-1,j),(i,j==1 ? L : j-1)]
@inline function neighbors(i,L)
    r = (i-1)%L
    c = div(i-1, L)
    # @info r,c
    SVector{4, Int32}(
        1 + (r+1)%L + c*L,                  # up
        1 + ((r-1)>=0 ? r-1 : L-1) + c*L,   # down
        1 + r + L*((c+1)%L),                # right
        1 + r + L*((c-1)>=0 ? c-1 : L-1)    # left
    )
end
@inline function neighbors_leftup(i,L)
    r = (i-1)%L
    c = div(i-1, L)
    # @info r,c
    SVector{2, Int32}(
        1 + ((r-1)>=0 ? r-1 : L-1) + c*L,
        1 + r + L*((c-1)>=0 ? c-1 : L-1)
    )
end


"""
    cluster_sizes(state)

A cluster is a connected region of aligned spins. Retrieve the sizes of all of them.
"""
function cluster_sizes(state)
    # Check for each element if they belong to the same cluster
    # as their left or upper neighbor. 
    cluster = zeros(Int32, size(state))
    max_cluster = 0
    pointer = 1
    new_cluster = true
    while pointer <= length(state)
        nn = neighbors_leftup(pointer, size(state, 1))
        new_cluster = true
        for n in nn
            if cluster[n]!=0 && state[pointer]==state[n]
                cluster[pointer] = cluster[n]
                new_cluster = false
            end
        end
        if new_cluster
            cluster[pointer] = max_cluster += 1
        end
        pointer += 1
    end

    return countmap(reshape(cluster, length(cluster))) |> values |> collect
end

"""
    metropolis_step!(state,beta,h)

Perform one step of the Metropolis algorithm.
`state` is mutated.
"""
function metropolis_step!(state::Matrix{Int8},beta,h)
    i = rand(1:size(state,1))
    j = rand(1:size(state,2))
    dh = dH(state,i,j,1.,h)
    if dh <= 0 || rand()<exp(-beta*dh)
        state[i,j] *= -1
    end
    return nothing
end

"""
    metropolis_sweep!(state,n,beta,h)

Perform n Metropolis steps.
"""
function metropolis_sweep!(state,n,beta,h)
    for _ in 1:n
        metropolis_step!(state,beta,h)
    end
end

"""
    init(L,beta,h,sweep)

Setup a configuration with `LxL` spins and perform an intial thermal sweep of
`sweep` timesteps.
"""
function init(L,beta,h,sweep)
    state = fill(Int8(-1), L, L)
    # Initial sweep to get into the steady state
    metropolis_sweep!(state,sweep,beta,h)
    return state
end

"""
    run_metropolis(L, beta, h;Tmax=1, sweep=0, sample_interval=1)

Sets up a random state and runs the Metropolis algorithm for a given set of parameters.
Samples in defined intervals along the Markov-Chain.

An initial thermal sweep to go to equilibrium may be specified.

Returns an array of averaged observables: [E, E^2, m, m^2, m^4] with total
energy E and magnetisation per spin m.

# Arguments
- L::Integer:            Linear system size
- beta::Float:            Inverse temperature
- h::Float:               External field
- Tmax::Integer:          Number of steps
- sweep::Integer:         Length of the initial sweep
- sample_interval::Integer: sample interval

# Example
```
julia> run_metropolis(50, 0., 0.;Tmax=50*10^3*50^2,sample_interval=10*50^2,sweep=10^3*50^2)
5-element Array{Float64,1}:
  -39.9664      # <E>
 1753.5168      # <E^2>
    0.010382079999999908 # <m>
    0.00011572633600000021 # <m^2>
    1.6924325969920386e-8 # <m^4>
```
"""
function run_metropolis(L::Int, beta,h;Tmax::Int=1,sweep::Int=0,sample_interval::Int=1)
    ## Initialise a random state
    state = init(L,beta,h,sweep)
    return _run_metropolis!(state,beta,h,Tmax=Tmax,sample_interval=sample_interval)
end

"""
    _run_metropolis!(state::Matrix{Int8},beta,h;Tmax::Int=1,sample_interval::Int=1)

Analogous to `run_metropolis` but mutates an existing state. Called by `run_metropolis`.
"""
function _run_metropolis!(state::Matrix{Int8},beta,h;Tmax::Int=1,sample_interval::Int=1)

    ## Define a matrix to record the observables.
    ## One entry for each observable, e.g E and m.
    ## Preallocating the matrix gives much better performance than
    ## constructing it on the fly.
    observables = zeros(Float64, 5)

    k = 0 #counts the number of samples

    t = 0 #simulation steps
    e = 0.
    mag = 0.

    @inbounds begin
        while(t<Tmax)
            ## Take the defined no. of steps before
            ## recording a measurement
            metropolis_sweep!(state,sample_interval,beta,h)

            ## Record observables
            e = H(state,1.,h)
            mag = m(state)
            observables[1] += e
            observables[2] += e^2
            observables[3] += mag
            observables[4] += mag^2
            observables[5] += mag^4

            ## increment counters
            k+=1
            t+=sample_interval
        end
    end
    ## Return sample means
    return observables/k
end

"""
    metropolis_timeseries(L, β, Tmax; sample_interval=L^2, sweep=1000)

Follows along a Markov chain in time.
Initialises a random configuration, sweeps it for `sweep*L^2` timsteps, and records
magnetetisation and energy every `sample_interval` until `Tmax`.

Returns a tuple of these time series.
"""
function metropolis_timeseries(L::Int, beta, Tmax; sample_interval=L^2, sweep=1000)
    state = init(L,beta,0.,sweep*L^2)
    #time series magnetization
    tsm = Vector{Float64}(undef, div(Tmax,sample_interval))
    #time series energy
    tse = Vector{Float64}(undef, div(Tmax,sample_interval))

    t=0
    k=1
    mag0 = m(state)
    while(t<Tmax)
        metropolis_sweep!(state,sample_interval,beta,0.)
        tsm[k] = m(state)
        tse[k] = H(state)
        t+=sample_interval
        k+=1
    end
    return tsm, tse
end

### Wolff algorithm ###
### --------------- ###

function bonds(i, L)
    nbs = neighbors(i, L)
    l_nb, d_nb = nbs[4], nbs[2]
    ((i,1), (d_nb, 1), (i,2), (l_nb, 2)), nbs
end

"""
    cluster!(cluster_state, state, i,j, p)

Builds a cluster around position `(i,j)` with acceptance probability `p`.

`cluster_state` is a boolean matrix of `size(state)×2`. It is set to `true` whenever
a site belongs to the cluster.
"""
function cluster!(cluster_state, state, i, p)
    nflipped = 0
    @inbounds begin
        L = size(state, 1)
        to_visit = Stack{Int}(L^2)
        cur = i
        s = state[cur]
        push!(to_visit, cur)
        while !isempty(to_visit)
            cur = pop!(to_visit)
            if state[cur] == s
                state[cur] = -s
                nflipped += 1
            end
            bnds, nbs = bonds(cur, L)
            for (b,n) in zip(bnds, nbs)
                if !cluster_state[b...]
                    cluster_state[b...] = true
                    if state[n]==s && rand()<p
                        push!(to_visit, n)
                    end
                end
            end
        end
    end
    return nflipped
end

"""
    wolff_step!(state, cluster_state, beta)

Build a cluster and flip it.
"""
function wolff_step!(state, cluster_state, beta, h=0.0)
    @inbounds begin
        L = size(state,1)
        i = rand(1:L^2)
        fill!(cluster_state, false)
        nflipped = cluster!(cluster_state, state, i, 1.0-exp(-2.0*beta))
        # state[view(cluster_state, :,:,1)] *= -1
    end
    return nflipped
end

function wolff_sweep!(state, cluster, n, args...)
    nflipped = 0
    for _ in 1:n
        nflipped += wolff_step!(state, cluster, args...)
    end
    return nflipped
end

"""
    _run_wolff!(state,cluster,beta,h;Tmax=1,sample_interval=1)

Run the Wolff cluster algorithm on a given state. You need to provide a cluster matrix, e.g
`similar(state,Bool)`.

See also: [`_run_metropolis`]
"""
function _run_wolff!(state::Matrix{Int8},cluster,beta,h;Tmax::Int=1,sample_interval::Int=1)

    ## Define a matrix in which to record the observables.
    ## One row for each observable, e.g E and m.
    ## Preallocating the matrix gives much better performance than
    ## constructing it on the fly.
    observables = zeros(Float64,5)

    k = 0 #counts the number of samples

    # pick the middle lattice site for the spin correlation

    t=0 #number of simulation steps
    r=0 #repeat counter
    e=0.
    mag=0.
    @inbounds begin
        while(t<Tmax)
            ## Do the defined no. of steps before
            ## taking a measurement
            for _ in 1:sample_interval
                wolff_step!(state, cluster, beta, h)
            end
            ## Record observables
            e = H(state,1.,0.)
            mag = m(state)|>abs
            observables[1] += e
            observables[2] += e^2
            observables[3] += mag
            observables[4] += mag^2
            observables[5] += mag^4

            ## increment counters
            k+=1
            t+=sample_interval
        end
    end
    ## Return statistics about the observables
    return observables/k
end

"""
    run_wolff(L, beta,h;Tmax=1,sweep=0,sample_interval=1)

See also: [`run_metropolis`]
"""
function run_wolff(L::Int, beta,h;Tmax::Int=1,sweep::Int=0,sample_interval::Int=1)
    ## Initialise a random state
    function init(L,beta)
        state = randomConfiguration(L)
        cluster = fill(false, L^2, 2) # (up, right) neighbor bonds for each site
        # Initial sweep to get into the steady state
        for _ in 1:sweep
            wolff_step!(state, cluster, beta, h)
        end
        return state,cluster
    end

    state, cluster = init(L,beta)
    return _run_wolff!(state,cluster,beta,0.,Tmax=Tmax,sample_interval=sample_interval)
end
