using SpecialFunctions

export golovinSCE, density_lnr_golovin_analytic

""""
    golovinSCE(x,x0,n0,b,time)

Analytic solution to the Smoluchowski Collection Equation using the Golovin (additive) kernel (Golovin, 1963).

# Arguments
- `x:: volume evaulated at`
- `x0:: initial volume`
- `n0:: initial number density`
- `b:: Golovin kernel coefficient`
- `time:: time at which to evaluate the solution`

Note the use of the besselix function rather than besseli, to solve overflow issues. The exponential term is then 
modified to compensate. 

"""
function golovinSCE(x,x0,n0,b,time)
    time = time == 0.0 ? 1e-10 : time
    τ  =  1 - exp(-x0*n0*b*time)
    sqrt_τ = sqrt(τ)
    factor = (1-τ) / (x*sqrt_τ) # (1-τ) / (x/x0*sqrt_τ)
    expterm = exp(-(1+τ - 2*sqrt_τ)*x/x0)
    besselterm = besselix(1,2*sqrt_τ*x/x0)

    return factor*besselterm*expterm
end

function density_lnr_golovin_analytic(run_settings,coag_settings::coag_settings{FT}) where FT<:AbstractFloat
    volume_bins_edges = radius_to_volume.(run_settings.radius_bins_edges)
    mids = 0.5*(volume_bins_edges[1:end-1] + volume_bins_edges[2:end])
    mids_rad =  0.5*(run_settings.radius_bins_edges[1:end-1] + run_settings.radius_bins_edges[2:end])
    
    dm = diff(volume_bins_edges)
    dr = diff(run_settings.radius_bins_edges)
    tlen = length(run_settings.output_steps)
    bins = zeros(FT, run_settings.num_bins, tlen)
    for i in 1:tlen
        test = coag_settings.n0*coag_settings.ΔV .* golovinSCE.(
            mids,
            radius_to_volume.(coag_settings.R0),
            coag_settings.n0,
            coag_settings.golovin_kernel_coeff,
            run_settings.output_steps[i])
        ys = test .* mids_rad .* dm ./ dr
        bins[:,i] = ys ./(coag_settings.ΔV) #kg
        if run_settings.binning_method == mass_density_lnr
            bins[:,i] = bins[:,i] .*(mids) .* (1000.0)
        end
    end

    return bins
end



# radius_range = 10 .^(LinRange(-6,-3,100))
# vol_range = 4pi / 3 * (radius_range.^3) 
# x0 = 4pi / 3 * (30.561e-6)^3
# n0 = 2^23
# b = 1500

# solnarr = golovinSCE.(vol_range, x0, n0, b, 1.0) .* vol_range

# plot(radius_range,solnarr,xscale=:log10, label="Golovin SCE", xlabel="Radius (m)", title="Golovin SCE Solution", legend=:topright)
