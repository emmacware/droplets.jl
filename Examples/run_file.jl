
using Random
using Plots
using Droplets
using CPUTime
include("testfunctions.jl")



Ns::Int = 2^14
FT = Float64

coagsettings = coag_settings{FT}(Ns=Ns,Δt=1)
runsettings=run_settings{FT}(coag_threading =Serial())


drops = runsettings.init_method(coagsettings)


bins,times = coag_runtime(1,drops,coagsettings,runsettings)
plot()
plot_dsd(bins*kg_to_g,runsettings)

