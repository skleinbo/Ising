## Start with `julia --project=.. -i -- example01_sim.jl`

# Make sure necessary packages (as defined in ../Project.toml) are installed
import Pkg; Pkg.instantiate()

using Distributed
## If you want to run on multiple cores, uncomment the following line:
#addprocs()

@everywhere begin
    using DataFrames, CSV
    using DataFramesMeta
    using ProgressMeter
    import Base.Iterators: product
    include("../src/mcmc.jl")
end
# Temperature range
Ts = range(2.0, 2.5, length=100)
# System sizes
Ls = 20:10:60

it = product(Ls, Ts)
lit = length(it)

## This gives you a dataframe filled with the full time series
raw_rows = @showprogress pmap(it) do (L,T)
    begin
        τ = 10    # time between samples in MCS
        N = 50000 # No. of samples

        (L=L, T=T, m=metropolis_timeseries(L, 1/T, N*τ*L^2; sample_interval=τ*L^2)) # take samples
    end
end
d = DataFrame(reshape(raw_rows, lit));

## Compute averages
## m, |m|, m^2, m^4, e, e^2
e = @linq d |>
    select(:L,:T, e=mean.(:e), e2=mean.(map(x->x.^2, :e)),
         m=mean.(:m), m2=mean.(map(x->x.^2, :m)), m4=mean.(map(x->x.^4,:m)), absem=mean.(map(x->abs.(x),:m)));

CSV.write("mydata.csv", e);
