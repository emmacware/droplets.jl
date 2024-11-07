
using Random
using Combinatorics
using Distributions
using Plots
using Interpolations
using BenchmarkTools
using Droplets
include("DSDvis.jl")
include("testfunctions.jl")
using CSV
using DataFrames
using Profile



Ns::Int = 2^17
scale = Ns * (Ns - 1) / 2 / (Ns / 2)
FT = Float64

coagsettings = coag_settings{FT}(Ns=Ns,scale=scale,Δt=1)
runsettings=run_settings{FT}(coag_threading =Parallel())


ξ, R, X = runsettings.init_method(coagsettings)
I = collect(1:Ns)
drops = droplets_allocations(ξ, R, X, I, zeros(FT, div(Ns, 2)), zeros(FT, div(Ns, 2)))


bins,times = coag_runtime(1,drops,coagsettings,runsettings)
plot()
plot_dsd(bins)
