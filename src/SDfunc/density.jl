########################################################################################
# This file contains the functions to calculate the density of liquid water
########################################################################################
export calc_ρw

#Only 1 function in here, calc_ρw

function calc_ρw(droplets,ΔgridV,Nx,Ny,grid_dict)
    for i in 1:Nx
        for j in 1:Ny
            grid = (i,j)

            drop = grid_dict[i,j]
            if isempty(drop)
                ρw[i,j] = 0
                continue
            else 
                ρw[i,j] = sum([droplet.ξ*droplet.X*ρl for droplet in grid_dict[i,j]])/ΔgridV
            end
        end
    end

    return ρw, grid_dict
end


"""
    calc_ρw(droplets, ΔV)

Calculate the density of water (ρw) based on the given droplet attributes and domain volume.

# Arguments
- `droplets`: an array of droplets with attributes `ξ` (multiplicity) and `X` (volume).
- `ΔV`: grid cell volume.

# Returns
- `ρw`: density of water

"""
function calc_ρw(droplets, ΔV)
    ρw = sum([droplet.ξ * droplet.X * ρl for droplet in droplets]) / ΔV
    return ρw
end

"""
    calc_ρw(ξ, X, ΔV)

Calculate the density of liquid water (ρw) based on the given droplet multiplicity(ξ) and volume(X),
and domain volume.

# Returns
- `ρw`: density of water

"""
function calc_ρw(ξ, X, ΔV)
    ρw = sum(ξ .* X .* ρl) / ΔV
    return ρw
end