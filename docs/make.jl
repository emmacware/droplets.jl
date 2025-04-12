push!(LOAD_PATH, "../src/")
using Documenter, Droplets

makedocs(
    sitename = "Droplets.jl",
    # modules = [Droplets],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
    ],
    clean = true,
)

deploydocs(repo = "github.com/emmacware/droplets.jl.git", target = "build")
