using Documenter
using SpinMonteCarlo

makedocs( format=:html,
          sitename="SpinMonteCarlo.jl",
          modules=[SpinMonteCarlo],
          pages = [
                   "Home" => "index.md",
                   "Manual" => Any[
                                   "Run Monte Carlo" => "runmc.md",
                                   "Develop Monte Carlo" => "develop.md",
                                  ],
                   "Library" => Any[
                                    "Public" => "lib/public.md",
                                    "Internals" => "lib/internals.md",
                                   ]
                  ]
        )

deploydocs( repo = "github.com/yomichi/SpinMonteCarlo.jl.git",
            julia="0.7",
            target="build",
            deps = nothing,
            make = nothing,
          )
