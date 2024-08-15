include("supercrystals.jl")
include("constants.jl")
using DifferentialEquations

dep_params = deposition_params{Float64}(deposition="Capacitance", ρ_dep_p="Constant", 
    constants=constants, const_dt=1, solver="Euler", mval=10.0, ro=1e-6, temppoly=268.15)

astart = [200e-6, 200e-6]
ξstart = [10.0, 10.0]
Mstart = [1e-16, 1e-16]
ρistart = [917.0, 917.0]
ρdstart = [200.0, 200.0]

crystals = create_crystals(2,1,1,1,1,astart,astart,ξstart,Mstart,ρistart,ρdstart)


sat(qvarray,P) = qvarray.*P./(0.378.*qvarray .+ 0.622)
pvap = sat(0.000214,30091.93)
pdry = 30091.93 - pvap
ρm = pdry/(constants.Rd*230.97)+pvap/(constants.Rv*230.97)
ssi = 0.05
qv = 0.000214
qvsi = qv/(1+ssi)
esati = pvap/(1+ssi)
statevar = statevar_struct(230.97,30091.93,qv,ρm,esati,ssi,qvsi)

for i in 1:length(crystals)
    crystals[i].vel = termvz_ice_Heymsfield_Westbrook_2010(crystals[i].a, crystals[i].c,  crystals[i].ρ, 
        constants, statevar)
end


function icedep(du::Vector, u::Vector, p, t)
    statevar, constants, dep_params = p
    da, dc, dρ = crystal_deposition(statevar, u, constants, dep_params)
    
    du[1] = da
    du[2] = dc
    du[3] = dρ
end

tspan = (0.0, 1.0)
initial_crystal = [crystals[1].a, crystals[1].c, crystals[1].ρ]  # Array for a, c, ρ
# prob = ODEProblem(icedep,[crystals[1].a,crystals[1].c,crystals[1].ρ],tspan,(statevar,constants,dep_params))
prob = ODEProblem(icedep,initial_crystal,tspan,(statevar,constants,dep_params))
sol = solve(prob, Euler(), dt=1e-2)


# for i in 1:length(crystals)
    # crystals[i].vel = termvz_ice_Heymsfield_Westbrook_2010(crystals[i].a, crystals[i].c,  crystals[i].ρ, 
    #     constants, statevar)
# end