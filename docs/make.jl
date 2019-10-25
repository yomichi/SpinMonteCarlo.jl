using Documenter
using SpinMonteCarlo

makedocs( format=Documenter.HTML(),
          sitename="SpinMonteCarlo.jl",
          modules=[SpinMonteCarlo],
          pages = [
                   "Home" => "index.md",
                   "Manual" => Any[
                                   "Run Monte Carlo" => "runmc.md",
                                   "Generate lattice" => "lattice.md",
                                   "Develop Monte Carlo" => "develop.md",
                                  ],
                   "Library" => Any[
                                    "Public" => "lib/public.md",
                                    "Internals" => "lib/internals.md",
                                   ]
                  ]
        )

deploydocs(
           repo = "github.com/yomichi/SpinMonteCarlo.jl.git",
           target = "build"
          )
