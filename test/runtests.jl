using Droplets
using Test

@testset "Droplets.jl" begin
    FT = Float64
    Ns = 2
    scale = Ns * (Ns - 1) / 2 / (Ns / 2)
    coagsettings = coag_settings{FT}(Ns=Ns,scale=scale,Δt=1)

    ξ = [Int(2),Int(3)]
    R = [FT(1.0),FT(2.0)]
    X = 3pi/4 .*R.^3
    I = [Int(1),Int(2)]
    drops = droplets_allocations(ξ, R, X, I, zeros(FT, div(Ns, 2)), zeros(FT, div(Ns, 2)))

    @test golovin(drops, (1,2), coagsettings) == 1500*(X[1]+X[2])

    pair_Ps!(1, (1,2), drops, golovin,coagsettings)
    @test drops.pαdt[1] == 3*1500*(X[1]+X[2])

    drops.pαdt[1]=1
    drops.ϕ[1]=0.5
    sdm_update!((1,2),1, drops)

    @test drops.ξ == [2,1]

end

@testset "Golovin Test" begin
    include("golovintest.jl")
end
