include("supercrystals.jl")
include("constants.jl")
using DifferentialEquations

dep_params = deposition_params{Float64}(deposition="Capacitance", ρ_dep_p="Constant", 
    constants=constants, const_dt=1, solver="Euler", mval=10.0, ro=1e-6, temppoly=268.15)

astart = [200e-6, 200e-6] # 200 microns
ξstart = [10.0, 10.0]
Mstart = [1e-16, 1e-16]
ρistart = [constants.ρ_ice, constants.ρ_ice]
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



function icedep_ode(du::Vector{Float64}, u::Vector{Float64}, p, t)
    statevar, constants, dep_params = p
    # Convert the input vector back to a Supercrystal object
    crystal = crystallize_vector(u) 
    # Compute the derivatives using the crystal_deposition function
    da, dc, dρ = crystal_deposition(statevar, crystal, constants, dep_params)
    # Update the derivatives in the du vector
    du[1] = da
    du[2] = dc
    du[8] = dρ  # Assuming `ρ` is in the 8th position in the vector
end

tspan = (0.0, 1.0)
initial_crystal = vectorize_crystal(crystals[1])
# prob = ODEProblem(icedep,[crystals[1].a,crystals[1].c,crystals[1].ρ],tspan,(statevar,constants,dep_params))
prob = ODEProblem(icedep_ode,initial_crystal,tspan,(statevar,constants,dep_params))
sol = solve(prob, Euler(), dt=1e-4, verbose=true, progress=true)


# for i in 1:length(crystals)
    # crystals[i].vel = termvz_ice_Heymsfield_Westbrook_2010(crystals[i].a, crystals[i].c,  crystals[i].ρ, 
    #     constants, statevar)
# end