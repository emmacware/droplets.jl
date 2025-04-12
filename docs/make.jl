push!(LOAD_PATH, "../")
using Documenter, Droplets

makedocs(
    sitename = "Droplets.jl",
    modules = [Droplets],
    format = Documenter.HTML(),
    deploydocs = [
        Documenter.HTML(),
        Documenter.GitHubActions(),
    ]
    pages = [
        "Home" => "index.md",
    ],
    clean = true,
)

deploydocs(repo = "github.com/emmacware/droplets.jl.git", target = "build")
