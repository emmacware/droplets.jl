using testset



@testset begin
    include("mvp.jl")
end

@testset begin
    include("golovintest.jl")
end