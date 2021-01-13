using Documenter

push!(LOAD_PATH,"../src/")
using FindPeaks1D

makedocs(
    sitename = "FindPeaks1D.jl",
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "index.md",
        "findpeaks1d.md",
    ],
)

deploydocs(
  repo = "github.com/ymtoo/FindPeaks1D.jl",
  versions = ["stable" => "v^", "v#.#", "dev" => "master"],
  branch = "gh-pages",
)