using PyCall

function analytic_soln(time,radius_bins_edges,n0,X0,ΔV, golovin_kernel_coeff)
    si = pyimport("PySDM.physics").si
    rho_w = pyimport("PySDM.physics.constants_defaults").rho_w # kg/m3
    Golovin = pyimport("PySDM.dynamics.collisions.collision_kernels").Golovin
    volume_bins_edges = radius_to_volume.(radius_bins_edges)#4/3*pi*(radius_bins_edges).^3
    mids = 0.5*(volume_bins_edges[1:end-1] + volume_bins_edges[2:end])
    mids_rad =  0.5*(radius_bins_edges[1:end-1] + radius_bins_edges[2:end])


    dm = diff(volume_bins_edges) #m^3
    dr = diff(radius_bins_edges)
    collision_kernel=Golovin(b=golovin_kernel_coeff/ si.s)

    gsol = zeros(Float64,length(mids))

    for t in time
        test = n0*ΔV .*collision_kernel.analytic_solution(x=mids,t=t,x_0=X0,N_0=n0)
        ys = test .* mids_rad .* dm ./ dr
        y_true = ys .*(mids)*(rho_w)/(ΔV) #kg
        gsol = hcat(gsol,y_true)
    end

    return gsol[:,2:end]
end