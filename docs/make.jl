using Documenter
using SpinMonteCarlo

makedocs( format=:html,
          sitename="SpinMonteCarlo.jl",
          modules=[SpinMonteCarlo],
          pages = [ "Home" => "index.md",
                    "Develop your program" => "develop.md",
                    "Library" => Any[
                      "Public" => "lib/public.md",
                      "Interanals" => "lib/internals.md",
                     ]
                  ]
        )

deploydocs( repo = "github.com/yomichi/SpinMonteCarlo.jl.git",
            julia="0.6",
            target="build",
            deps = nothing,
            make = nothing,
          )
