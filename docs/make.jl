using Documenter
using SpinMonteCarlo

makedocs(
         format=:html,
         sitename="SpinMonteCarlo.jl",
         modules=[SpinMonteCarlo],
         pages = [
                  "Home" => "index.md",
                  "API" => "api.md",
                 ]
        )
