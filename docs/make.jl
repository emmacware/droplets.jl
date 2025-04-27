push!(LOAD_PATH, "../")
using Documenter, Droplets


Coalescence = [
        "SDM" => "coalescence.md",
        "Kernels" => "kernels.md",
]
Condensation = [
        "Condensation" => "condensation.md",
]

makedocs(
    sitename = "Droplets.jl",
    modules = [Droplets],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Coalescence" => Coalescence,
        "Condensation" => Condensation,
    ],
    clean = true,
)

deploydocs(repo = "github.com/emmacware/Droplets.jl.git",branch = "gh-pages", target = "build",forcepush=true)
