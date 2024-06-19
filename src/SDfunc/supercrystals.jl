#Crystals!


mutable struct Crystals{FT <: Real, T <: Union{Tuple{FT,FT}, Tuple{FT,FT,FT}}}
    a::FT # equatorial radius (m)
    c::FT # polar radius (m)
    # r::FT # wet radius (m)
    loc::T #location (m)
    ξ::Int #multiplicity
    M::FT # mass (g)
    ρi::FT # ice density (kg/m^3)
    ρd::FT # deposition density (kg/m^3)

end

function create_crystals(Ns,Nx,Ny,Δx,Δy,astart,cstart,ξstart,Mstart,ρistart,ρdstart)
    crystals = []
    for i in 1:Ns
        a = astart[i]
        c = cstart[i]
        loc = [rand(0.003:Nx*Δx), rand(0.003:Ny*Δy)]
        ξ = ξstart[i]
        M = Mstart[i]
        ρi = ρistart[i]
        ρd = ρdstart[i]
        crystal = Crystals(a,c,loc,ξ,M,ρi,ρd)
        push!(crystals, crystal)
    end
    return crystals 
end
