
using Random
using Combinatorics
using Distributions
using Plots
using Interpolations
using StaticArrays, BenchmarkTools
include("Settings.jl")
include("DSDvis.jl")
using Droplets
include("testfunctions.jl")
using CSV
using DataFrames
using Profile
using ProfileVega



Ns::Int = 2^17
scale = Ns * (Ns - 1) / 2 / (Ns / 2)
FT = Float64

coagsettings = coag_settings{FT}(Ns=Ns,scale=scale,Δt=10)
runsettings=run_settings{FT}(coag_threading =Serial())

ξ, R, X = runsettings.init_method(coagsettings)

bins,coal_time = coag_runtime(1,ξ,R,X,coagsettings,runsettings)
plot_dsd(bins)


plot(bins)