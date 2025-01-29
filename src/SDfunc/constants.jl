##################################################
# Constants
##################################################

export Constants,constants
export kg_to_g, m_to_μm, volume_to_radius, radius_to_volume

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

constants = Constants{Float32}()

#Conversions
const kg_to_g = 1e3
const m_to_μm = 1e6

@inline function volume_to_radius(V::AbstractFloat)
    return (3*V/(4*pi))^(1/3)
end

@inline function radius_to_volume(R::AbstractFloat)
    return 4/3*pi*R^3
end





