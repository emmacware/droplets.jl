push!(LOAD_PATH, "../")
using Documenter, Droplets

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

Coalescence = [
        "SDM" => "coalescence.md",
        "Kernels" => "kernels.md",
]
Condensation = [
        "Condensation" => "condensation.md",
]

deploydocs(repo = "github.com/emmacware/Droplets.jl.git",branch = "gh-pages", target = "build",forcepush=true)
