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
    λ::FT
    κ::FT
    k::FT
    μ::FT
    Dv::FT
    ρl::FT
    Rd::FT
    Rv::FT
    gconst::FT
    L::FT
    Cp::FT
end

function Constants()

return Constants{Float64}(1.5625,1.5625,1.5625,1.5625,1.5625,1000.0,287.0,461.0,9.8,22.6e5,4.2e3)
end



