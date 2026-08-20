# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

This release fixes a number of bugs that silently returned wrong physics, so
results produced with v1.2.2 and earlier should be regarded with suspicion for
the affected models and update methods. Random number streams have changed as
well, so a run cannot be reproduced bit-for-bit across this boundary.

### Breaking

#### Random number generation

- The default generator is now `Xoshiro` instead of `MersenneTwister`. It draws
  about twice as fast and carries 91 bytes of state instead of 19283, which
  matters because the generator is serialized into every checkpoint. Set
  `param["RNG"] = MersenneTwister` to keep the old generator.
- Model structs carry the generator type as a type parameter, so a model can
  hold any `AbstractRNG`. Previously the field was declared concrete and
  handing a model anything else failed on field assignment.
- `runMC(param)` no longer reseeds the generator after constructing the model.
  The constructor is now the single point where the seed is applied. Previously
  the stream was rewound right after part of it had been used to lay out the
  initial spins, so streams differ from earlier versions.
- `"Seed"` must be an integer. Non-integer seeds used to be routed through
  `Base.hash`, whose values Julia does not promise to keep stable across
  releases.
- `runMC(params)` derives a per-job child seed from each element's `"ID"`, so a
  parallel sweep explores independent streams. Previously every worker given
  the same `"Seed"` produced N copies of one simulation.
- Checkpoint files written by earlier versions cannot be restored.

#### Removed and changed APIs

- `gen_snapshot!`, `gensave_snapshot!`, and `load_snapshot` are removed and
  replaced with stubs that raise an informative error. They had been broken
  since Julia 1.0 and every call path threw immediately. A visualization
  extension is planned as a replacement (issue #35).
- `extrapolate_stderror` and `extrapolate_tau` now return a **1-sigma standard
  error** as the second element. They previously returned the half-width of a
  95% confidence interval, because `LsqFit.estimate_errors` scales the standard
  error by a Student-t quantile. To recover the old value, multiply by
  `quantile(TDist(dof), 0.975)` where `dof = nfit - 2` and `nfit` is the number
  of binning levels actually fitted.
- `extrapolate_stderror` and `extrapolate_tau` throw `ArgumentError` when fewer
  than three binning levels are available, or when a level's statistic is not
  finite. Previously `ErrorException`, `DimensionMismatch`, or `DomainError`
  leaked out of LsqFit and Distributions.
- `dim` is a function owned by this package rather than an extension of
  `Distributions.dim`. Code that does `using SpinMonteCarlo, Distributions` must
  now qualify the call as `SpinMonteCarlo.dim`.
- `<<(::VectorObservable, ...)` accepts `Vector` rather than `AbstractVector`,
  matching what the `push!` methods actually implement. Passing a view used to
  dispatch on `<<` and then hit a `MethodError` inside `push!`.
- `binning` (scalar and vector) throws `ArgumentError` on empty data, or when
  the requested number of bins exceeds the number of observations, instead of
  silently returning all-`NaN` observables.
- The helicity modulus is reported per dimension as `"Helicity Modulus x"`,
  `"Helicity Modulus y"`, ... instead of a hard-coded 2D pair. One-dimensional
  lattices no longer emit `"Helicity Modulus y"`, and three-dimensional
  lattices no longer write out of bounds.
- Potts `SW_update!` and `Wolff_update!` throw `ArgumentError` for
  antiferromagnetic couplings, which the Fortuin-Kasteleyn cluster
  construction does not support. They previously sampled at effectively
  infinite temperature without complaint.
- The correctly spelled `"Periodic Boundary Condition"` key is now honoured.
  The old misspelled `"Periodic Boudary Condition"` still works but warns.

### Fixed

- **Ising and AshkinTeller `local_update!` were not ergodic.** A deterministic
  sweep order combined with always accepting `dE = 0` flips made the Markov
  chain reducible: on an L=8 chain the energy came out as -0.485 against the
  exact -0.818. Sites (and, for AshkinTeller, species) are now chosen at
  random.
- **Ising `SW_update!` and `Wolff_update!` were wrong for `J < 0`.** Only
  satisfied bonds are activated now, with `p = 1 - exp(-2|J|/T)`, and whole
  clusters flip by a random sign. The improved estimator and the energy offset
  were generalized accordingly.
- **The XY and Clock helicity modulus ignored the coupling `J`.**
- **The AshkinTeller energy had its overall sign inverted.**
- **The QuantumXXZ improved estimator was biased** for `Sign * Magnetization^2`
  and `^4` when the sign depends on the loop flips: +3.5% at `Gamma = 0.7`, and
  -0.9% on a three-site XY ring. Linear observables and the energy were already
  unbiased.
- `interpolate!` mutated the shared `Expr` in `stdbravais`, freezing lattice
  constants after the first generation and contaminating every lattice sharing
  a Bravais definition.
- `"cubic lattice"` was missing its `:parameters` key, so all 3D lattice
  generation failed with `KeyError`.
- `LatticeParameter`'s `getindex`/`setindex!` returned the `:name` field
  instead of the requested key.
- `generatelattice_std` mutated the user's `param["L"]` array.
- `BinningVectorObservable.push!` crashed on the 256th push with the default
  `minbinnum`, and its extrapolations passed malformed data to the fitter.
  Both had been broken since Julia 1.0.
- `merge!(::SimpleVectorObservable, ...)` always threw `UndefVarError`;
  `var(::JackknifeVector)` produced `NaN`; broadcast methods for
  `JackknifeVector` silently returned an empty observable; `jackknife` on a set
  of vector observables threw a `MethodError`; `reset!(::MCObservableSet)`
  recursed infinitely.
- `confidence_interval(obs, :sigma4)` threw `FieldError`, and an unparseable
  symbol threw the same error instead of a clear one.
- The example on the documentation top page had an unbalanced parenthesis, and
  the lattice lists in `README.md` and `docs/src/lattice.md` were out of sync
  with the lattices that actually exist.

### Added

- Aqua.jl quality checks in the test suite.
- Transfer-matrix-based regression tests on chains (Ising FM/AF for local, SW
  and Wolff updates; AshkinTeller; the `J` dependence of the helicity modulus;
  3D XY), so this class of physics bug is caught automatically.
- This changelog.

### Removed

- The `LsqFit` dependency. Its only use was a straight-line fit in
  `extrapolate_detail`, now done in closed form.
- The `DataStructures` dependency. Its only use was `Stack(Deque{Int}())` in the
  `Wolff_update!` implementations, which `Vector` serves identically.
- `src/observables/second_jackknife.jl` (never included since the Julia 1.0
  port, and broken beyond the load path), `tool/latticeviewer.jl` (built on
  Makie APIs that no longer exist), and `__precompile__()` (the default since
  Julia 1.0).

### Changed

- Minimum supported Julia version is 1.10.
- Dependency compat bounds modernized; the docs environment is pinned to
  Documenter 1.

## [1.2.2] and earlier

See the git history.
