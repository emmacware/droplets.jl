
using Random
using Combinatorics
using Distributions
using Plots
using Interpolations
using BenchmarkTools
using Droplets
using CPUTime
include("DSDvis.jl")
include("testfunctions.jl")
using CSV
using DataFrames
using Profile



Ns::Int = 2^15
FT = Float64

coagsettings = coag_settings{FT}(Ns=Ns,Δt=1)
runsettings=run_settings{FT}(coag_threading =Serial())


drops = runsettings.init_method(coagsettings)


bins,times = coag_runtime(1,drops,coagsettings,runsettings)
plot()
plot_dsd(bins,runsettings)

