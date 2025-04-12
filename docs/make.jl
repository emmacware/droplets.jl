push!(LOAD_PATH, "../")
using Documenter, Droplets

makedocs(
    sitename = "Droplets.jl",
    modules = [Droplets],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
    ],
    clean = true,
)


deploydocs(repo = "github.com/emmacware/droplets.jl.git",branch = "gh-pages", target = "build",forcepush=true)
