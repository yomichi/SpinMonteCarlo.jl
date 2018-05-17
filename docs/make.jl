using Documenter
using SpinMonteCarlo

makedocs(
         format=:html,
         sitename="SpinMonteCarlo.jl",
         modules=[SpinMonteCarlo],
         pages = [
                  "Home" => "index.md",
                  "Develop your program" => "develop.md",
                  "Library" => Any[
                    "Public" => "lib/public.md",
                    "Interanals" => "lib/internals.md",
                   ]
                 ]
        )
