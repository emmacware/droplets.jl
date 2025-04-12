############################################################################################
# This file contains functions to update the position of the droplets in the domain 
############################################################################################

export update_position!

# based on superdroplet struct type and grid_dict organization



# Functions:
# limit(a, N_x)
# limity(a, N_y)
# limitvface(a, N_y)
    # these functions are used to limit the droplet location to the domain,
    # periodic horizontal boundaries and non-periodic vertical boundaries
# add_to_grid!(droplet, grid_index, grid_dict)
    # Adds a droplet to a grid_dict list
# move_to_grid!(droplet, old_grid_index, new_grid_index, grid_dict)
    # Moves a droplet from one grid_dict list to another
# update_position!(droplets,Nx,Ny,dt,ρu,ρv,ρ, grid_dict,gridbox)
    # Updates the position of the droplets based on the velocity field and grid_dict organization 
    # maybe has some bugs? I recently redid my indexing which may have caused some issues but it worked before  


############################################################################################

limit(a, N_x) = a > N_x ? a - N_x : a < 1 ? a + N_x : a
limity(a, N_y) = a > N_y ? N_y : a < 1 ? 1 : a
limitvface(a, N_y) = a > N_y +1 ? N_y +1 : a < 1 ? 1 : a

function add_to_grid!(droplet, grid_index, grid_dict)
    if haskey(grid_dict, grid_index)
        push!(grid_dict[grid_index], droplet)
    else
        grid_dict[grid_index] = [droplet]
    end
end

# Function to move a superdroplet to a new grid
function move_to_grid!(droplet, old_grid_index, new_grid_index, grid_dict)
    # Remove droplet from old grid
    filter!(d -> d !== droplet, grid_dict[old_grid_index])

    # Add droplet to new grid
    add_to_grid!(droplet, new_grid_index, grid_dict)
end


function update_position!(droplets,Nx,Ny,dt,ρu,ρv,ρ, grid_dict,gridbox)
    moved_droplets = Set()

    for i in 1:Nx
        for j in 1:Ny
            if isempty(grid_dict[i,j])
                continue
            else
                n = limitvface(j+1,Ny)
                s = j
                e = limit(i+1,Nx)
                w = i
                for droplet in grid_dict[i,j]
                    # print(droplet.loc[1])
                    if droplet in moved_droplets
                        continue
                    else
                        droplet.loc[1] += dt*((droplet.loc[1]-i*Δx)*ρu[w,j]+ (1-(droplet.loc[1]-i*Δx))*ρu[e,j])/(2*ρ[i,j])
                        droplet.loc[2] += dt*(((droplet.loc[1]-i*Δx)*ρv[i,n]+ (1-(droplet.loc[2]-i*Δy))*ρv[i,s])/(2*ρ[i,j])-terminal_v(droplet.R))
                        push!(moved_droplets, droplet)
                        # print(droplet.loc[1])

                        if droplet.loc[1] <= gridbox[i,j][1] || droplet.loc[1] >= gridbox[i,j][2] 
                            if droplet.loc[1] < 0
                                droplet.loc[1] = droplet.loc[1]+Nx*Δx
                            elseif droplet.loc[1] > Nx*Δx
                                droplet.loc[1] = droplet.loc[1]-Nx*Δx
                            end
                            move_to_grid!(droplet, (i,j), (ceil(droplet.loc[1]/Δx),j), grid_dict)
                        end

                        if droplet.loc[2] <= gridbox[i,j][3]
                            if droplet.loc[2] <= 0
                                droplet.loc[2] = 0 #gridbox[i,j][3]
                                move_to_grid!(droplet, (i,j), (0,0), grid_dict)
                            else
                                move_to_grid!(droplet, (i,j), (i,(ceil(droplet.loc[2]/Δy))), grid_dict)
                            end
                        elseif droplet.loc[2]  >= gridbox[i,j][4] 
                            if i == Ny
                                droplet.loc[2] = gridbox[i,Ny][4]
                            else
                            move_to_grid!(droplet, (i,j), (i,(ceil(droplet.loc[2]/Δy))), grid_dict)
                            end
                        end
                    end
                end
            end
        end
    end
    return droplets, grid_dict
        
end