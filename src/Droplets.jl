module Droplets
    # under construction!
    include(joinpath("SDfunc", "constants.jl"))
    include(joinpath("SDfunc", "coalescence.jl"))
    include(joinpath("SDfunc", "condensation.jl"))
    include(joinpath("SDfunc", "setup.jl"))
    include(joinpath("SDfunc", "updateposition.jl"))

    gconst = 9.81  # gravitational constant
    Rd = 287       # gas constant for dry air
    Cp = 1005      # specific heat capacity of dry air
    #export init_ξ_const

end
