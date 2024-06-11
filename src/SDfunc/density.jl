########################################################################################
# This file contains the functions to calculate the density of liquid water
########################################################################################
export calc_ρw

#Only 1 function in here, calc_ρw
#options for arguments:
#calc_ρw(droplets,ΔgridV,Nx,Ny,grid_dict), returns ρw, grid_dict
#       -this option goes through the whole system, using a grid_dict of superdroplet 
#       mutable structs to calculate the density of liquid water as a grid center value
#calc_ρw(droplets,ΔV), returns ρw
#       -this option is meant for a box or single grid cell: 
#       or just a relevant list of superdroplets
#       droplets should be a vector of superdroplet structs
#calc_ρw(ξ,X,ΔV), returns ρw
#       -this option is meant for a box or single grid cell:
#       takes in the ξ and X values as vectors

#Not sure how to implement this in 2D looping over the full Ns superdroplets
#without splitting it up by grid yet

########################################################################################


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


function calc_ρw(droplets,ΔV)
    ρw= sum([droplet.ξ*droplet.X*ρl for droplet in droplets])/ΔV
    return ρw
end

function calc_ρw(ξ,X,ΔV)
    ρw= sum(ξ.*X.*ρl)/ΔV
    return ρw
end