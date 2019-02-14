### Monte Carlo routines for
### the 2D Ising model

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
    s = 0.
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
    sweep!(state,n,beta,h)

Perform n Metropolis steps.
"""
function sweep!(state,n,beta,h)
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
    state = randomConfiguration(L)
    # Initial sweep to get into the steady state
    sweep!(state,sweep,beta,h)
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
    observables = zeros(Float64,5)

    k = 0 #counts the number of samples

    t = 0 #simulation steps
    e = 0.
    mag = 0.

    @inbounds begin
        while(t<Tmax)
            ## Take the defined no. of steps before
            ## recording a measurement
            sweep!(state,sample_interval,beta,h)

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
magnetetisation every `sample_interval` until `Tmax`.
"""
function metropolis_timeseries(L::Int, beta, Tmax; sample_interval=L^2, sweep=1000)
    state = init(L,beta,0.,sweep*L^2)
    ts = Vector{Float64}(div(Tmax,sample_interval))
    t=0
    k=1
    mag0 = m(state)
    while(t<Tmax)
        sweep!(state,sample_interval,beta,0.)
        ts[k] = m(state)
        k+=1
        t+=sample_interval
    end
    return ts
end

### Wolff algorithm ###
### --------------- ###

"""
    cluster!(cluster_state, state, i,j, p)

Builds a cluster around position `(i,j)` with acceptance rate `p`.

`cluster_state` is a boolean matrix of the same size as the system. It is set to `true` whenever
a site belongs to the cluster.
"""
function cluster!(cluster_state, state, i,j, p)
    @inbounds begin L = size(state, 1)
        s = state[i,j]
        cluster_state[i,j] = true
        for neighbor in [(i%L+1,j),(i,j%L+1),(i==1 ? L : i-1,j),(i,j==1 ? L : j-1)]
            if state[neighbor...] == s && !cluster_state[neighbor...] && rand()<p
                cluster_state[neighbor...]=true
                cluster!(cluster_state, state, neighbor..., p)
            end
        end
    end
end

"""
    wolff_step!(state, cluster_state, beta)

Build a cluster and flip it.
"""
function wolff_step!(state, cluster_state, beta, h)
    @inbounds begin
        L = size(state,1)
        i,j = rand(1:L,2)
        fill!(cluster_state, false)
        cluster!(cluster_state, state, i,j, 1.0-exp(-2.0*beta))
        state[cluster_state] *= -1
    end
    return nothing
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
        cluster = zeros(Bool,L,L)
        # Initial sweep to get into the steady state
        for _ in 1:sweep
            wolff_step!(state, cluster, beta, h)
        end
        return state,cluster
    end

    state,cluster = init(L,beta)
    return _run_wolff!(state,cluster,beta,0.,Tmax=Tmax,sample_interval=sample_interval)
end
