

#---------------------------------------------------------
# KERNELS
#---------------------------------------------------------
export terminal_v,hydrodynamic,golovin

# terminal velocity of droplets (this should move..)
"""
    terminal_v(r::FT) where FT<:AbstractFloat

Compute the terminal velocity of a droplet of radius r, using tables from
Gunn and Kinzer, (1949), https://doi.org/10.1175/1520-0469(1949)006<0243:TTVOFF>2.0.CO;2

# Arguments
- `r::FT`: The radius of the droplet in meters.

# Returns
The terminal velocity of the droplet in meters/second.

"""
function terminal_v(r::FT)::FT where FT<:AbstractFloat  # terminal velocity 


    if 2*r*100<0.0078
        tv=1.2*10e6*(r*100)^2
        tv = tv/100
    else    
        d_table = [0.078,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,
            2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8]./10

        v_table = [18,27,72,117,162,206,247,287,327,367,403,464,517,565,609,649,690,
        727,757,782,806,826,844,860,872,883,892,898,903,907,909,912,914,916,917]

        interpo_extrapo = linear_interpolation(d_table,v_table,extrapolation_bc=Line())
        tv = interpo_extrapo(2*r*100)/100 # radius in meters to diameter in cm, then velocity cm->m
    end
    return tv
end

#not working correctly:
# #collision efficiency function
# function collision_efficiency(R1::FT,R2::FT)::FT where FT<:AbstractFloat
#     #Parameterization from Berry 1967
#     #https://doi.org/10.1175/1520-0469(1967)024<0688:CDGBC>2.0.CO;2
#     r = max(R1,R2)*1e6
#     rs = min(R1,R2)*1e6

#     p = rs/r
#     D = (-27)/(r^1.65)
#     E = (-58)/(r^1.9)
#     F = (15/r)^4 +1.13 
#     G = (16.7/r)^8 +1 +0.004*r

#     Y = 1+p+D/(p^F)+E/((1-p)^G)
#     if Y<0
#         Y=0
#     end

#     return Y
# end

#---------------------------------------------------------
# Coalescence Kernels
#---------------------------------------------------------
"""
    hydrodynamic(droplets::droplet_attributes, (j,k)::Tuple{Int,Int}, settings::coag_settings{FT})::FT where FT<:AbstractFloat

Hydronynamic kernel function for droplet coalescence, used to calculate the probabilities of two 
    droplets colliding.
    K(r,r') = π * (r+r')^2 * |v(r)-v(r')| * E(r,r')
    where E(r,r') is the collision efficiency function, (currently not implemented), and v is the terminal velocity of 
    a droplet of radius size r.

# Arguments
- `droplets::droplet_attributes`: The attributes of the droplets.
- `(j,k)::Tuple{Int,Int}`: The indices of the droplets.
- `settings::coag_settings{FT}`: The coagulation settings.

"""
@inline function hydrodynamic(droplets::droplet_attributes, (j,k)::Tuple{Int,Int}, settings::coag_settings{FT})::FT where FT<:AbstractFloat
    Rj, Rk = volume_to_radius(droplets.X[j]), volume_to_radius(droplets.X[k])
    if settings.hydrodynamic_collision_eff_func == true
        E = collision_efficiency(Rj,Rk)
    else    
        E = 1
    end
    Rsum = (Rj+Rk) # sum of radius is in meters
    vdif = abs(terminal_v(Rj)-terminal_v(Rk))
    return E *π * Rsum^2 *vdif
end

"""
    golovin(droplets, (j, k), settings)

Golovin (or additive) kernel function for droplet coalescence, used to calculate the probabilities of two droplets colliding.
    K(x,x') = b(x+x'),where b is the Golovin kernel coefficient, and x is the volume of the droplet.
Taken from Golovin (1963), this kernel has an analytic solution for the Smulochowski Coagulation Equation.

# Arguments
- `droplets`: droplet attributes
- `(j, k)`: indices of the droplets
- `settings`: coagulation settings

"""
@inline function golovin(droplets::droplet_attributes, (j,k)::Tuple{Int,Int}, settings::coag_settings{FT})::FT where FT<:AbstractFloat
    return settings.golovin_kernel_coeff *(droplets.X[j] + droplets.X[k])# Xsum
end
