#---------------------------------------------------------
#Visualization functions for droplet size distribution
#---------------------------------------------------------


#functions:

#binning_dsd(droplets,bins,t;smooth = true,scope_init = 1),             returns radius_bin_mids, bin_water_volume
#binning_dsd(Xunsorted,ξunsorted,bins,t;smooth = true,scope_init = 1),  returns radius_bin_mids, bin_water_volume
#smoothbins!(bins,scope_init=1),                                        returns smoothed_bin_values
#error_measure(y,ytrue,x),                                              returns integrated_error
#gnormal(R,X,ξ,Ns),                                                     returns g(lnR) guassian KDE **i dont recommend
#plot_grid_with_droplets(grid, droplets),                               returns plot of grid with droplets



#####


#binning_dsd function for superdroplets, superdroplet struct input
#becuase of what it was used for, time is an input so that it would not smooth at t=0 
function binning_dsd(droplets,bins,t;smooth = true,scope_init = 1)
    droplets = sort(droplets, by = droplet -> droplet.R)
    bin_edges = bins
    # mids = 0.5*(bin_edges[1:end-1] + bin_edges[2:end])*1e6 #um
    mids = bin_edges[1:end-1]*1e6 #um
    water = zeros(length(mids))
    droplet_idx = 1
    for j in 1:length(mids)

        lo = 4*pi .*((bin_edges[j]) .^ 3) ./3 
        hi = 4*pi .*((bin_edges[j+1]) .^ 3) ./3
        leng = log(bin_edges[j+1])-log(bin_edges[j])
        for i in droplet_idx:length(droplets)

            if droplets[i].X < hi
                water[j]=water[j]+droplets[i].X*droplets[i].ξ/leng
                droplet_idx += 1
            else
                break
            end

        end

    end
    scope = scope_init
    if t != 0 && smooth == true
        new_water = copy(water)
        for j in 1:2
            # scope = scope_init
            for i in scope+1:length(water)-scope
                new_water[i] = mean(water[i-scope:i+scope])
            end
            scope = 1
            for i in scope+1:length(water)-scope
                water[i] = mean(new_water[i-scope:i+scope])
            end
        end
    end
    return mids, water
end


#binning_dsd function for superdroplets, vector input
#becuase of what it was used for, time is an input so that it would not smooth at t=0 
function binning_dsd(Xunsorted,ξunsorted,bins,t;smooth = true,scope_init = 1)
    bin_edges = bins
    i = sortperm(Xunsorted)
    X = Xunsorted[i]
    ξ = ξunsorted[i]

    # mids = 0.5*(bin_edges[1:end-1] + bin_edges[2:end])*1e6 #um
    mids = bin_edges[1:end-1]*1e6 #um

    water = zeros(length(mids))
    droplet_idx = 1
    for j in 1:length(mids)

        lo = 4*pi .*((bin_edges[j]) .^ 3) ./3 
        hi = 4*pi .*((bin_edges[j+1]) .^ 3) ./3
        leng = log(bin_edges[j+1])-log(bin_edges[j])
        for i in droplet_idx:length(X)

            if X[i] < hi
                water[j]=water[j]+X[i]*ξ[i]/leng
                droplet_idx += 1
            else
                break
            end

        end

    end

    scope = scope_init

    if t != 0 && smooth == true
        new_water = copy(water)
        for j in 1:2
            # scope = scope_init
            for i in scope+1:length(water)-scope
                new_water[i] = mean(water[i-scope:i+scope])
            end
            scope = 1
            for i in scope+1:length(water)-scope
                water[i] = mean(new_water[i-scope:i+scope])
            end
        end


    end

    return mids, water

end

# Smoothing function for bins, moving average 
function smoothbins!(bins,scope_init=1)
    scope = scope_init
    new_bins = copy(bins)
    for j in 1:2
        # scope = scope_init
        for i in scope+1:length(bins)-scope
            new_bins[i] = mean(bins[i-scope:i+scope])
        end
        scope = 1
        for i in scope+1:length(bins)-scope
            bins[i] = mean(new_bins[i-scope:i+scope])
        end
    end
    return bins
end

#integrated error function, as in PySDM
function error_measure(y,ytrue,x)
    errors = ytrue .- y
    errors1 = errors[1:end-1] .+ errors[2:end]
    dx = diff(x)
    errors1 = errors1 .* dx ./2
    err = sum(abs.(errors1))
    return err
end

# Gaussian kernel density estimator function
# Useless honestly just bin it, smoothing function is there if you want it
function gnormal(R,X,ξ,Ns)
    σ₀ = 0.62
    σ = σ₀*Ns^(-1/5)
    W(Y) = 1/(2*π).^0.5/σ * exp.(-Y.^2 ./(2*σ^2))
    Y = [log(R[i]) .- log.(R) for i = 1:lastindex(R)]
    g = 1/ΔV * sum(ξ.*X.*ρl.*1e3.*W.(Y)) # 1e3 converts kg to g
    return g
end


function plot_grid_with_droplets(grid, droplets)


    Nx, Ny = size(grid)

    bigdrops = [drop for drop in droplets if drop.R > 1e-7]

    if isempty(bigdrops)
        plot(background=:dodgerblue)
    else
        scatter([droppie.loc[1] for droppie in bigdrops], [droppie.loc[2] for droppie in bigdrops], 
            markersize= [(7e4*droppie.R) for droppie in bigdrops],#0.4, #alpha= [1*10*droppie.R^2*droppie.ξ for droppie in bigdrops],
            label=false,background=:dodgerblue,c=cgrad(:thermal, rev = true),clims=(5, 90),
            zcolor=[droppie.R*1e6 for droppie in bigdrops],
            colorbar_title="Radius(μm)")
    end

    
    xlabel!("X")
    ylabel!("Y")
    title!("Domain")
end
