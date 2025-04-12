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


deploydocs(repo = "github.com/emmacware/Droplets.jl.git", devbranch = "dev",branch = "gh-pages", target = "build",forcepush=true)
