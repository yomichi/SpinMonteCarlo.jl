using Documenter
using SpinMonteCarlo

makedocs(
         format=:html,
         sitename="SpinMonteCarlo.jl",
         modules=[SpinMonteCarlo],
         pages = [
                  "Home" => "index.md",
                  "Develop your program" => "develop.md",
                  "API" => "api.md",
                 ]
        )
