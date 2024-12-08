##################################################
# Constants
##################################################

export Constants,constants #does this need to be capital somewhere else
export diffusion_constant_water_vapor, density_liquid_water, dry_air_Rd, water_vapor_Rv, gravity_const, latent_heat_vaporization, dry_air_Cp    

# ρl = Float64(1000.0) # Density of water kg/m3.
# Rd = Float64(287.0) #Dry Air Gas Constant J/kgK
# Rv = Float64(461.0) #Vapor Gas Constant J/kgK
# gconst = 9.8 # gravitational constant m/s2
# L = 22.6e5 # Latent Heat of Vaporization J/kg
# K = 0.026 # Thermal Conductivity W/mK
# D = 2.5e-6 #Diffusivity of water vapor in air m^2/s 
# Cp = 4.2e3 # Specific Heat of Water J/kgK

# AerosolActivationParameters(::Type{FT}) where {FT <: AbstractFloat} =

Base.@kwdef struct Constants{FT<:AbstractFloat}
    λ = FT(1.5625)        # brief
    κ = FT(1.5625)         # description
    k = FT(0.024)         # from clima K_therm
    μ = FT(1.6e-5)         # from clima ν_air
    Dv = FT(2.26e-5)        # from clima D_vapor water vapor in air m2/s???
    ρl = FT(1000.0)        # and
    Rd = FT(287.0)        # its
    Rv = FT(461.0)        # units
    gconst = FT(9.8)    # gravitational constant, m/s2
    L = FT(22.6e5)         # Latent Heat of Vaporization J/kg
    Cp = FT(4181)        # Specific Heat of Dry air at constant pressure J/kgK
end


# const λ = FT(1.5625)        # brief
# const κ = FT(1.5625)         # description
# const k = FT(0.024)         # of
# const μ = FT(1.6e-5)         # each
# const diffusion_constant_water_vapor = FT(2.26e-5)        # field   Molecular diffusion constant of ???water vapor in air m2/s???
# const density_liquid_water = FT(1000.0)        # and
# const dry_air_Rd = FT(287.0)        # its
# const water_vapor_Rv = FT(461.0)        # units
# const gravity_const = FT(9.8)    # gravitational constant, m/s2
# const latent_heat_vaporization = FT(22.6e5)         # Latent Heat of Vaporization J/kg
# const dry_air_Cp = FT(4.2e3)        # Specific Heat of Dry air at constant pressure J/kgK


constants = Constants{Float32}()

# function Constants(;λ=1.5625, 
#                     κ=1.5625, 
#                     k=0.024, 
#                     μ=1.6e-5, 
#                     Dv=2.26e-5, 
#                     ρl=1000.0, 
#                     Rd=287.0, 
#                     Rv=461.0, 
#                     gconst=9.8, 
#                     L=22.6e5, 
#                     Cp=4.2e3)
#     return Constants{Float64}(λ, κ, k, μ, Dv, ρl, Rd, Rv, gconst, L, Cp)
# end


