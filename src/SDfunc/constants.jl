##################################################
# Constants
##################################################

export Constants

# ρl = Float64(1000.0) # Density of water kg/m3.
# Rd = Float64(287.0) #Dry Air Gas Constant J/kgK
# Rv = Float64(461.0) #Vapor Gas Constant J/kgK
# gconst = 9.8 # gravitational constant m/s2
# L = 22.6e5 # Latent Heat of Vaporization J/kg
# K = 0.026 # Thermal Conductivity W/mK
# D = 2.5e-6 #Diffusivity of water vapor in air m^2/s 
# Cp = 4.2e3 # Specific Heat of Water J/kgK


# #CHANGE THESE
# #various diffusion/viscoustiy constants
# λ = 1.5625#2.623*10e−2 #W m−1 K−1
# κ = 1.5625 #m2/s−1
# k = 1.5625
# μ = 1.5625
# Dv = 1.5625

struct Constants{FT<:AbstractFloat}
    λ::FT         # brief
    κ::FT         # description
    k::FT         # of 
    μ::FT         # each
    Dv::FT        # field   Molecular diffusion constant of ???water vapor in air m2/s???
    ρl::FT        # and
    Rd::FT        # its
    Rv::FT        # units
    gconst::FT    # gravitational constant, m/s2
    L::FT         # Latent Heat of Vaporization J/kg
    Cp::FT        # Specific Heat of Dry air at constant pressure J/kgK
end

function Constants(;λ=1.5625, 
                    κ=1.5625, 
                    k=0.024, 
                    μ=1.6e-5, 
                    Dv=2.26e-5, 
                    ρl=1000.0, 
                    Rd=287.0, 
                    Rv=461.0, 
                    gconst=9.8, 
                    L=22.6e5, 
                    Cp=4.2e3)
    return Constants{Float64}(λ, κ, k, μ, Dv, ρl, Rd, Rv, gconst, L, Cp)
end

constants = Constants()
export constants

