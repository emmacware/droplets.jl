using Droplets
using Test
# import Pkg; Pkg.add("Pkg")

@testset begin

    # Pkg.add(["Combinatorics", "Distributions", "Random", "JSON", "DelimitedFiles", "CPUTime", "Plots"])
    # using Random,Combinatorics,Distributions,CPUTime
    include("../Examples/testfunctions.jl")
    # using 
        
    FT = Float64
    setup = Dict()
    
    setup["N_PART"]= 2^16
    setup["VOLUME_M3"]= 1e6
    setup["DT_SEC"]= 1.
        
        # Plotting parameters
        setup["N_BINS"]= 128
        setup["R_BINS_MIN"]= 1e-5
        setup["R_BINS_MAX"]= 5e-3
        setup["NUM_CONC_PER_M3"]= 2^23
        setup["DIAM_AT_MEAN_VOL_M"]= 2*30.531e-6
        setup["ADDITIVE_KERNEL_COEFF"]= 1500
        setup["T_MAX_SEC"]= 3600
        setup["PLOT_TIME_STEP"]= 120 # seconds
    
    setup["N_BIN_EDGES"] = setup["N_BINS"] + 1

    Ns::Int = setup["N_PART"]
    FT = Float64
    num_bins::Int = setup["N_BINS"]
    radius_bins_edges = 10 .^ range(log10(setup["R_BINS_MIN"]), log10(setup["R_BINS_MAX"]), length=setup["N_BIN_EDGES"]) 
    R0 = FT(setup["DIAM_AT_MEAN_VOL_M"])/2
    T_MAX_SEC = setup["T_MAX_SEC"]
    PLOT_TIME_STEP = setup["PLOT_TIME_STEP"]
    output_steps = collect(0:PLOT_TIME_STEP:T_MAX_SEC)

    
    coagsettings = coag_settings{FT}(Ns=Ns,Δt=setup["DT_SEC"],ΔV=setup["VOLUME_M3"],
            golovin_kernel_coeff=FT(setup["ADDITIVE_KERNEL_COEFF"]), n0=FT(setup["NUM_CONC_PER_M3"]),R0=R0)
    runsettings_numberdens =run_settings{FT}(num_bins=num_bins,radius_bins_edges=radius_bins_edges,
            smooth = false,output_steps=output_steps,init_method=init_logarithmic,binning_method = number_density)
    runsettings_mass_density =run_settings{FT}(num_bins=num_bins,radius_bins_edges=radius_bins_edges,
            smooth = false,output_steps=output_steps,init_method=init_logarithmic,binning_method = mass_density_lnr)

    drops = runsettings_numberdens.init_method(coagsettings)
    num_dens_bin,timing = coag_runtime(1,drops,coagsettings,runsettings_numberdens)

    drops = runsettings_mass_density.init_method(coagsettings)
    mass_dens_bin,timing = coag_runtime(2,drops,coagsettings,runsettings_mass_density)

    # dict = Dict()
    # dict["Number Concentration (m^-3)"] = num_dens_bin
    # dict["Mass Density (g/m^3 dlnr)"] = mass_dens_bin
    # json_string = JSON.json(dict)

    # open("output.json", "w") do file
    #     write(file, json_string)
    # end
end
