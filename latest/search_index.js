var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Home-1",
    "page": "Home",
    "title": "Home",
    "category": "section",
    "text": "This package provides Markov chain Monte Carlo solvers for lattice spin systems. Several models, lattices, and algorithms are already defined. Moreover, you can define your own model, lattice, and algorithm and combine with pre-defined ones."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "julia> Pkg.add(\"SpinMonteCarlo\")"
},

{
    "location": "index.html#Simple-example-1",
    "page": "Home",
    "title": "Simple example",
    "category": "section",
    "text": "The following simple example calculates and prints temperature dependence of specific heat for Ising model on a 16times16 square lattice by Swendsen-Wang algorithm.using SpinMonteCarlo\n\nconst model = Ising\nconst lat = square_lattice\nconst L = 16\nconst update! = SW_update!\n\nconst Tc = 2.0/log1p(sqrt(2))\nconst Ts = Tc*linspace(0.85, 1.15, 31)\nconst MCS = 8192\nconst Therm = MCS >> 3\n\nfor T in Ts\n    param = Parameter(\"Model\"=>model, \"Lattice\"=>lat,\n                      \"L\"=>L, \"T\"=>T, \"J\"=>1.0,\n                      \"Update Method\"=>update!,\n                      \"MCS\"=>MCS, \"Thermalization\"=>Therm,\n                     )\n    result = runMC(param)\n    println(@sprintf(\"%f %.15f %.15f\",\n                      T, mean(result[\"Specific Heat\"]), stderror(result[\"Specific Heat\"])))\nend"
},

{
    "location": "runmc.html#",
    "page": "Run Monte Carlo",
    "title": "Run Monte Carlo",
    "category": "page",
    "text": ""
},

{
    "location": "runmc.html#Run-Monte-Carlo-1",
    "page": "Run Monte Carlo",
    "title": "Run Monte Carlo",
    "category": "section",
    "text": ""
},

{
    "location": "develop.html#",
    "page": "Develop Monte Carlo",
    "title": "Develop Monte Carlo",
    "category": "page",
    "text": ""
},

{
    "location": "develop.html#Develop-Monte-Carlo-1",
    "page": "Develop Monte Carlo",
    "title": "Develop Monte Carlo",
    "category": "section",
    "text": ""
},

{
    "location": "develop.html#Lattice-1",
    "page": "Develop Monte Carlo",
    "title": "Lattice",
    "category": "section",
    "text": "You can define your own lattice as an instance of Lattice."
},

{
    "location": "develop.html#Model-1",
    "page": "Develop Monte Carlo",
    "title": "Model",
    "category": "section",
    "text": "Model should contain following fields: lat :: Lattice and rng :: Random.MersenneTwister. Model also should have a constructor taking param :: Parameter as a argument.You should define convert_parameter for your Model. This is a helper function which takes Model and Parameter and returns arguments of update method and estimator.@gen_convert_parameter helps you to define convert_parameter. For example, if your model::A needs a scalar T = param[\"T\"] and a vector Js = param[\"J\"] with numbondtypes(model) elements,@gen_convert_parameter(A, (\"T\", 1, 1.0), (\"J\", numbondtypes, 1.0))defines your documented and type-stable  convert_parameter(model::A, param::Parameter) which returns T and Js."
},

{
    "location": "develop.html#Note:-1",
    "page": "Develop Monte Carlo",
    "title": "Note:",
    "category": "section",
    "text": "That the second element of each tuple is not a Function means that a return value is a scalar (the case of \"T\").\nThe third element of each tuple is the default value."
},

{
    "location": "develop.html#Update-method-1",
    "page": "Develop Monte Carlo",
    "title": "Update method",
    "category": "section",
    "text": "\"Update method\" is a function which (in-place) updates model::Model under some parameters such as temperature T. For example, local_update!(model::Ising, T, Js) updates a spin configuration of model by local spin flip and Metropolice-Hasting algorithm under temperature T and coupling constants Js. \"Update method\" can return some object which will be used in \"Estimator\" as extra information."
},

{
    "location": "develop.html#Example-1",
    "page": "Develop Monte Carlo",
    "title": "Example",
    "category": "section",
    "text": "Swendsen-Wang algorithm SW_update!(::Ising, ::Parameter) returns cluster information sw::SWInfo used in improved_estimator."
},

{
    "location": "develop.html#Estimator-1",
    "page": "Develop Monte Carlo",
    "title": "Estimator",
    "category": "section",
    "text": "\"Estimator\" is a function which returns observables of a configuration model as a Dict{String, Any}. Arguments of a \"Estimator\" are model::Model, parameters (return of convert_parameter), and extra information (return of \"Update method\") in order.Default \"Estimator\" is determined by model and \"Update Method\" as return of default_estimator(model, param[\"Update Method\"])."
},

{
    "location": "develop.html#Example-2",
    "page": "Develop Monte Carlo",
    "title": "Example",
    "category": "section",
    "text": "improved_estimator(model::Ising, T::Real, Js::AbstractArray, sw::SWInfo) takes a return of SW_update!(model::Ising, T::Real, Js::AbstractArray) as the last argument."
},

{
    "location": "develop.html#Postprocess-1",
    "page": "Develop Monte Carlo",
    "title": "Postprocess",
    "category": "section",
    "text": "postproc(model::Model, param::Parameter, obs::MCObservableSet) is a post process of a Monte Carlo run. Most important objective is to calculate functions of expectation values stored in obs.For example, \"Specific Heat\" C can be calculated from \"Energy^2\" langle E^2rangle big N^2, \"Energy\" leftlangle Erightrangle big N, the number of site N, and temperature T as C = left(NbigT^2right)left leftlangle E^2 rightrangle - left langle E rightrangle^2 right. This is realized asjk = jackknife(obs)\nT = param[\"T\"]\nN = numsites(model)\njk[\"Specific Heat\"] = (N/(T*T)) * (jk[\"Energy^2\"] - jk[\"Energy\"]^2)jackknife converts MCObservableSet (e.g. BinningObservableSet) to JackknifeObservableSet, which enables to calculate functions of mean values and these statistical errors.postproc returns a MCObservableSet (usually jk::JackknifeObservableSet above), which is also the return value of runMC."
},

{
    "location": "lib/public.html#",
    "page": "Public",
    "title": "Public",
    "category": "page",
    "text": ""
},

{
    "location": "lib/public.html#Public-interfaces-1",
    "page": "Public",
    "title": "Public interfaces",
    "category": "section",
    "text": "Documentation for SpinMonteCarlo.jl\'s public interface (exported)."
},

{
    "location": "lib/public.html#SpinMonteCarlo.runMC",
    "page": "Public",
    "title": "SpinMonteCarlo.runMC",
    "category": "function",
    "text": "runMC(param::Parameter)\nrunMC(params::AbstractArray{Parameter} ; parallel::Bool=false)\n\nRuns Monte Carlo simulation(s) and returns calculated observables.\n\nKeyward aruguments\n\nparallel: If true, runs simulations in parallel (uses pmap instead of map).\n\nRequired keys in param\n\n\"Model\"\n\"Update Method\"\n\nOptional keys in param\n\n\"MCS\": The number of Monte Carlo steps after thermalization\nDefault: 8192\n\"Thermalization\": The number of Monte Carlo steps for thermalization\nDefault: MCS>>3\n\"Seed\": The initial seed of the random number generator, MersenneTwister\nDefault: determined randomly (see Random.seed!)\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Driver-1",
    "page": "Public",
    "title": "Driver",
    "category": "section",
    "text": "CurrentModule = SpinMonteCarlorunMC"
},

{
    "location": "lib/public.html#SpinMonteCarlo.Ising",
    "page": "Public",
    "title": "SpinMonteCarlo.Ising",
    "category": "type",
    "text": "Ising model with energy E = -sum_ij J_ij sigma_i sigma_j, where sigma_i takes value of 1 (up spin) or -1 (down spin).\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.Potts",
    "page": "Public",
    "title": "SpinMonteCarlo.Potts",
    "category": "type",
    "text": "Q state Potts model with energy E = -sum_ij delta_sigma_i sigma_j, where sigma_i takes an integer value from 1 to Q and delta is a Kronecker\'s delta. Order parameter (total magnetization) is defined as \\begin{equation}     M = \\frac{Q-1}{Q}N1 - \\frac{1}{Q}(N-N1), \\end{equation} where N is the number of sites and N_1 is the number of sigma=1 spins.\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.Clock",
    "page": "Public",
    "title": "SpinMonteCarlo.Clock",
    "category": "type",
    "text": "Q state clock model with energy E = -sum_ij J_ij cos(theta_i - theta_j), where theta_i = 2pi sigma_iQ and sigma_i takes an integer value from 1 to Q.\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.XY",
    "page": "Public",
    "title": "SpinMonteCarlo.XY",
    "category": "type",
    "text": "XY model with energy E = -sum_ij J_ij cos(theta_i - theta_j), where theta_i = 2pi sigma_i and sigma_i in 0 1).\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.QuantumXXZ",
    "page": "Public",
    "title": "SpinMonteCarlo.QuantumXXZ",
    "category": "type",
    "text": "Spin-S XXZ model represented as the following Hamiltonian, \\begin{equation} \\mathcal{H} = \\sum{i,j} \\left[ J{ij}^z Si^z Sj^z              + \\frac{J{ij}^{xy}}{2} (Si^+ Sj^- + Si^-Sj^+) \\right]             - \\sumi \\Gammai Si^x, \\end{equation} where S^x S^y S^z are x y and z component of spin operator with length S, and S^pm equiv S^x pm iS^y are ladder operator. A state is represented by a product state (spins at tau=0) of local S^z diagonal basis and an operator string (perturbations). A local spin with length S is represented by a symmetrical summation of 2S sub spins with length 12.\n\n\n\n"
},

{
    "location": "lib/public.html#Model-1",
    "page": "Public",
    "title": "Model",
    "category": "section",
    "text": "Ising\nPotts\nClock\nXY\nQuantumXXZ"
},

{
    "location": "lib/public.html#SpinMonteCarlo.dimer_lattice",
    "page": "Public",
    "title": "SpinMonteCarlo.dimer_lattice",
    "category": "function",
    "text": "dimer_lattice()\ndimer_lattice(param::Dict)\n\nGenerates dimer lattice (indeed, this is not a \"lattice\")\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.chain_lattice",
    "page": "Public",
    "title": "SpinMonteCarlo.chain_lattice",
    "category": "function",
    "text": "chain_lattice(L::Integer)\nchain_lattice(param::Dict)\n\nGenerates chain lattice with length L = param[\"L\"]\n\nnsitetypes is 2: all 1 (2) sites do not connect to another 1 (2) site.\n\nnbondtypes is 2: 1 bond connects 2n-1 and 2n sites and 2 bond connects 2n and 2n+1 sites.\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.square_lattice",
    "page": "Public",
    "title": "SpinMonteCarlo.square_lattice",
    "category": "function",
    "text": "square_lattice(L::Integer, W::Integer=L)\nsquare_lattice(param::Dict)\n\nGenerates square lattice with size L=param[\"L\"] times W=param[\"W\"].\n\nnsitetypes is 2: all 1 (2) sites do not connect to another 1 (2) site.\n\nnbondtypes is 2: 1 bonds are parallel to x axis and 2 are parallel to y axis.\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.triangular_lattice",
    "page": "Public",
    "title": "SpinMonteCarlo.triangular_lattice",
    "category": "function",
    "text": "triangular_lattice(L::Integer, W::Integer=L)\ntriangular_lattice(param::Dict)\n\nGenerates triangular lattice with size L=param[\"L\"] times W=param[\"W\"].\n\nnsitetypes is 3: all 1, 2, and 3 sites do not connect to another 1, 2, and 3 site, respectively.\n\nnbondtypes is 3: 1, 2, and 3 bonds make an angle of 0, 60, and 120 degree with x axis, respectively.\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.cubic_lattice",
    "page": "Public",
    "title": "SpinMonteCarlo.cubic_lattice",
    "category": "function",
    "text": "cubic_lattice(L::Integer, W::Integer=L, H::Integer=W)\ncubic_lattice(param::Dict)\n\nGenerates cubic lattice with size L=param[\"L\"] times W=param[\"W\"] times H=param[\"H\"].\n\nnsitetypes is 2: all 1 (2) sites do not connect to another 1 (2) site.\n\nnbondtypes is 2: 1, 2, and 3 bonds are parallel to x, y, and z axis, respectively.\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.fully_connected_lattice",
    "page": "Public",
    "title": "SpinMonteCarlo.fully_connected_lattice",
    "category": "function",
    "text": "fully_connected_lattice(N::Integer)\nfully_connected_lattice(param::Dict)\n\nGenerates N=param[\"N\"] site fully connected lattice\n\nBoth nsitetypes and nbondtypes are 1.\n\n\n\n"
},

{
    "location": "lib/public.html#Lattice-generator-1",
    "page": "Public",
    "title": "Lattice generator",
    "category": "section",
    "text": "dimer_lattice\nchain_lattice\nsquare_lattice\ntriangular_lattice\ncubic_lattice\nfully_connected_lattice"
},

{
    "location": "lib/public.html#SpinMonteCarlo.local_update!",
    "page": "Public",
    "title": "SpinMonteCarlo.local_update!",
    "category": "function",
    "text": "local_update!(model, param)\nlocal_update!(model, T::Real, Js::AbstractArray)\n\nUpdates spin configuration by local spin flip and Metropolice algorithm  under the temperature T = param[\"T\"] and coupling constants J = param[\"J\"]\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.SW_update!",
    "page": "Public",
    "title": "SpinMonteCarlo.SW_update!",
    "category": "function",
    "text": "SW_update!(model, param::Parameter)\nSW_update!(model, T::Real, Js::AbstractArray)\n\nUpdates spin configuration by Swendsen-Wang algorithm under temperature T=param[\"T\"] and coupling constants J=param[\"J\"]\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.Wolff_update!",
    "page": "Public",
    "title": "SpinMonteCarlo.Wolff_update!",
    "category": "function",
    "text": "Wolff_update!(model, param::Parameter)\nWolff_update!(model, T::Real, Js::AbstractArray)\n\nUpdates spin configuration by Wolff algorithm under temperature T=param[\"T\"] and coupling constants J=param[\"J\"]\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.loop_update!",
    "page": "Public",
    "title": "SpinMonteCarlo.loop_update!",
    "category": "function",
    "text": "loop_update!(model, param::Parameter)\nloop_update!(model, T::Real,\n             Jz::AbstractArray,\n             Jxy::AbstractArray,\n             Gamma:AbstractArray)\n\nUpdates spin configuration by loop algorithm  under the temperature T = param[\"T\"] and coupling constants Jz, Jxy and transverse field Gamma\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Update-method-1",
    "page": "Public",
    "title": "Update method",
    "category": "section",
    "text": "An index of model parameter (e.g., Js) is corresponding to sitetype or bondtype.local_update!\nSW_update!\nWolff_update!\nloop_update!"
},

{
    "location": "lib/public.html#SpinMonteCarlo.simple_estimator",
    "page": "Public",
    "title": "SpinMonteCarlo.simple_estimator",
    "category": "function",
    "text": "simple_estimator(model::Ising, T::Real, Js::AbstractArray)\nsimple_estimator(model::Potts, T::Real, Js::AbstractArray)\n\nReturns the following observables as Dict{String, Any}\n\nObservables\n\n\"Energy\"\nEnergy per spin (site)\n\"Energy^2\"\n\"Magnetization\"\nTotal magnetization per spin (order paremeter)\n\"|Magnetization|\"\n\"Magnetization^2\"\n\"Magnetization^4\"\n\n\n\n\n\nsimple_estimator(model::Clock, T::Real, Js::AbstractArray)\nsimple_estimator(model::XY, T::Real, Js::AbstractArray)\n\nReturns the following observables as Dict{String, Any}\n\nObservables\n\n\"Energy\"\nEnergy per spin (site)\n\"Energy^2\"\n\"|Magnetization|\"\nAbsolute value of total magnetization per spin (order paremeter)\n\"|Magnetization|^2\"\n\"|Magnetization|^4\"\n\"Magnetization x\"\nx component of total magnetization per spin (order paremeter)\n\"|Magnetization x|\"\n\"Magnetization x^2\"\n\"Magnetization x^4\"\n\"Magnetization y\"\ny component of total magnetization per spin (order paremeter)\n\"|Magnetization y|\"\n\"Magnetization y^2\"\n\"Magnetization y^4\"\n\"Helicity Modulus x\"\n\"Helicity Modulus y\"\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.improved_estimator",
    "page": "Public",
    "title": "SpinMonteCarlo.improved_estimator",
    "category": "function",
    "text": "improved_estimator(model::Ising, T::Real, Js::AbstractArray, sw::SWInfo)\n\nReturns the following observables as Dict{String, Any} using cluster information sw\n\nObservables\n\n\"Energy\"\nEnergy per spin (site)\n\"Energy^2\"\n\"Magnetization\"\nTotal magnetization per spin (site)\n\"|Magnetization|\"\n\"|Magnetization|^2\"\n\"|Magnetization|^4\"\n\n\n\n\n\nimproved_estimator(model::Potts, T::Real, Js::AbstractArray, sw::SWInfo)\n\nReturns the following observables as Dict{String, Any} using cluster information sw\n\nObservables\n\n\"Energy\"\n\"Energy^2\"\n\"Magnetization\"\n\"|Magnetization|\"\n\"|Magnetization|^2\"\n\"|Magnetization|^4\"\n\n\n\n\n\nimproved_estimator(model::QuantumXXZ, T::Real, Js::AbstractArray, uf::UnionFind)\n\nReturns the following observables as Dict{String, Any} using loop information uf\n\nObservables\n\n\"Sign\"\nSign of the weight function\n\"Sign * Energy\"\nEnergy per spin (site)\n\"Sign * Energy^2\"\n\"Sign * Magnetization\"\nTotal magnetization (Sz) per spin (site)\n\"Sign * |Magnetization|\"\n\"Sign * Magnetization^2\"\n\"Sign * Magnetization^4\"\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Estimator-1",
    "page": "Public",
    "title": "Estimator",
    "category": "section",
    "text": "simple_estimator\nimproved_estimator"
},

{
    "location": "lib/public.html#SpinMonteCarlo.convert_parameter",
    "page": "Public",
    "title": "SpinMonteCarlo.convert_parameter",
    "category": "function",
    "text": "convert_parameter(model, param)\n\nGenerates arguments of updater and estimator.\n\nExample\n\njulia> model = Ising(chain_lattice(4));\n\njulia> p = convert_parameter(model, Parameter(\"J\"=>1.0))\n(1.0, [1.0, 1.0]) # T and Js\n\njulia> p = convert_parameter(model, Parameter(\"J\"=>[1.5, 0.5]))\n(1.0, [1.5, 0.5]) # J can take a vector whose size is `numbondtypes(model)`\n\njulia> model.spins\n4-element Array{Int64,1}:\n  1\n  1\n  1\n -1\n\njulia> local_update!(model, p...);\n\njulia> model.spins\n4-element Array{Int64,1}:\n  1\n  1\n -1\n -1\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.convert_parameter-Tuple{Ising,Dict{String,Any}}",
    "page": "Public",
    "title": "SpinMonteCarlo.convert_parameter",
    "category": "method",
    "text": "convert_parameter(model::Union{Ising, Potts, Clock, XY}, param::Parameter)\n\nKeynames:\n\n\"T\": a scalar (default: 1.0).\n\"J\": a vector with numbondtypes(model) elements (default: 1.0).\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.convert_parameter-Tuple{QuantumXXZ,Dict{String,Any}}",
    "page": "Public",
    "title": "SpinMonteCarlo.convert_parameter",
    "category": "method",
    "text": "convert_parameter(model::QuantumXXZ, param::Parameter)\n\nKeynames:\n\n\"T\": a scalar (default: 1.0).\n\"Jz\": a vector with numbondtypes(model) elements (default: 1.0).\n\"Jxy\": a vector with numbondtypes(model) elements (default: 1.0).\n\"Gamma\": a vector with numsitetypes(model) elements (default: 0.0).\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Utility-1",
    "page": "Public",
    "title": "Utility",
    "category": "section",
    "text": "convert_parameter\nconvert_parameter(::Ising, ::Parameter)\nconvert_parameter(::QuantumXXZ, ::Parameter)"
},

{
    "location": "lib/internals.html#",
    "page": "Internals",
    "title": "Internals",
    "category": "page",
    "text": ""
},

{
    "location": "lib/internals.html#Internal-APIs-1",
    "page": "Internals",
    "title": "Internal APIs",
    "category": "section",
    "text": "Documentation for SpinMonteCarlo.jl\'s internals (not exported)."
},

{
    "location": "lib/internals.html#SpinMonteCarlo.accumulateObservables!",
    "page": "Internals",
    "title": "SpinMonteCarlo.accumulateObservables!",
    "category": "function",
    "text": "accumulateObservables!(model, obs::MCObservableSet, localobs::Dict)\n\nAccumulates localobs into obs. For example, obs[\"Energy\"] << localobs[\"Energy\"].\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.postproc",
    "page": "Internals",
    "title": "SpinMonteCarlo.postproc",
    "category": "function",
    "text": "postproc(model::Model, param::Parameter, obs::MCObservableSet)\n\nPost process of observables. For example, Specific heat will be calculated from energy, energy^2, and temperature.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.postproc-Tuple{Ising,Dict{String,Any},Dict{String,Obs} where Obs<:MCObservable}",
    "page": "Internals",
    "title": "SpinMonteCarlo.postproc",
    "category": "method",
    "text": "postproc(model::Union{Ising, Potts}, param::Parameter, obs::MCObservableSet)\n\nObservables to be calculated\n\nIn the following, m is total magnetization per site and epsilon is total energy per site.\n\n\"Binder Ratio\"\nR = fracleft langle m^4 right rangleleft langle m^2 rightrangle^2\n\"Susceptibility\"\nchi = fracNTleft(leftlangle m^2rightrangleright)\n\"Connected Susceptibility\"\nfracNTleft(leftlangle m^2rightrangle - leftlangle m rightrangle^2right)\n\"Specific Heat\"\nfracNT^2left(leftlangle epsilon^2rightrangle - leftlangle epsilon rightrangle^2right)\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.postproc-Tuple{Clock,Dict{String,Any},Dict{String,Obs} where Obs<:MCObservable}",
    "page": "Internals",
    "title": "SpinMonteCarlo.postproc",
    "category": "method",
    "text": "postproc(model::Union{Clock, XY}, param::Parameter, obs::MCObservableSet)\n\nObservables to be calculated\n\nIn the following, m is total magnetization per site and epsilon is total energy per site.\n\n\"Binder Ratio x\"\nfracleft langle m_x^4 right rangleleft langle m_x^2 rightrangle^2\n\"Binder Ratio y\"\nfracleft langle m_y^4 right rangleleft langle m_y^2 rightrangle^2\n\"Binder Ratio\"\nfracleft langle m^4 right rangleleft langle m^2 rightrangle^2\n\"Susceptibility x\"\nfracNTleft(leftlangle m_x^2rightrangleright)\n\"Susceptibility y\"\nfracNTleft(leftlangle m_y^2rightrangleright)\n\"Susceptibility y\"\nfracNTleft(leftlangle m^2rightrangleright)\n\"Connected Susceptibility x\"\nfracNTleft(leftlangle m_x^2rightrangle - leftlangle m_x rightrangle^2right)\n\"Connected Susceptibility y\"\nfracNTleft(leftlangle m_y^2rightrangle - leftlangle m_y rightrangle^2right)\n\"Connected Susceptibility\"\nfracNTleft(leftlangle m^2rightrangle - leftlangle m rightrangle^2right)\n\"Specific Heat\"\nfracNT^2left(leftlangle epsilon^2rightrangle - leftlangle epsilon rightrangle^2right)\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.postproc-Tuple{QuantumXXZ,Dict{String,Any},Dict{String,Obs} where Obs<:MCObservable}",
    "page": "Internals",
    "title": "SpinMonteCarlo.postproc",
    "category": "method",
    "text": "postproc(model::QuantumXXZ, param::Parameter, obs::MCObservableSet)\n\nObservables to be calculated\n\nIn the following, s is sign of weight, m is total magnetization per site, and epsilon is total energy per site.\n\n\"Magnetization\"\nleftlangle m srightrangleBigleftlangle s rightrangle\n\"|Magnetization|\"\nleftlangle m srightrangleBigleftlangle s rightrangle\n\"Magnetization^2\"\nleftlangle m^2 srightrangleBigleftlangle s rightrangle\n\"Magnetization^4\"\nleftlangle m^4 srightrangleBigleftlangle s rightrangle\n\"Energy\"\nleftlangle epsilon srightrangleBigleftlangle s rightrangle\n\"Energy^2\"\nleftlangle epsilon^2 srightrangleBigleftlangle s rightrangle\n\"Binder Ratio\"\nfracleft langle m^4 right rangleleft langle m^2 rightrangle^2\n\"Susceptibility\"\nfracNTleft(leftlangle m^2rightrangleright)\n\"Connected Susceptibility\"\nfracNTleft(leftlangle m^2rightrangle - leftlangle m rightrangle^2right)\n\"Specific Heat\"\nfracNT^2left(leftlangle epsilon^2rightrangle - leftlangle epsilon rightrangle^2right)\n\n\n\n"
},

{
    "location": "lib/internals.html#Driver-1",
    "page": "Internals",
    "title": "Driver",
    "category": "section",
    "text": "CurrentModule = SpinMonteCarloaccumulateObservables!\npostproc\npostproc(::Ising,::Parameter,::MCObservableSet)\npostproc(::Clock,::Parameter,::MCObservableSet)\npostproc(::QuantumXXZ,::Parameter,::MCObservableSet)"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.Lattice",
    "page": "Internals",
    "title": "SpinMonteCarlo.Lattice",
    "category": "type",
    "text": "Lattice\n\nFields\n\ndim :: Int\ndimension of lattice\nsize :: Vector{Int}\nlength of lattice in each dimension\nnsitetypes :: Int\nthe number of sitetypes\nnbondtypes :: Int\nthe number of bondtypes\nsites :: Vector{Vector{Int}}\nlist of site indecies for each sitetype\nbonds :: Vector{Vector{Int}}\nlist of bond indecies for each sitetype\nnsites :: Int\nthe number of sites\nnbonds :: Int\nthe number of bonds\nsitetypes :: Vector{Int}\nsitetype of each site\nbondtypes :: Vector{Int}\nbondtype of each bond\ntransvector :: Matrix{Float64}\nlattice vector represented in Cartesian system\n@assert size(transvector) == (dim, dim)\nsite_coords :: Matrix{Float64}\ncoordinate of each site represented in lattice system\n@assert size(site_coords) == (dim, nsites)\nbond_dirs :: Matrix{Float64}\ndisplacement of each bond represented in lattice system\n@assert size(bond_dirs) == (dim, nbonds)\nneighborsites :: Vector{Vector{Int}}\nlist of indecies of neighbor sites of each site\nneighborbonds :: Vector{Vector{Int}}\nlist of indecies of neighbor bonds of each site\nsource :: Vector{Int}\nthe index of an end site of each bond index\ntarget :: Vector{Int}\nthe index of the other end site of each bond index\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#Base.size-Tuple{Lattice}",
    "page": "Internals",
    "title": "Base.size",
    "category": "method",
    "text": "size(lat::Lattice, [dim::Integer])\nsize(model::Model, [dim::Integer])\n\nReturns the size of lattice.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#PDMats.dim-Tuple{Lattice}",
    "page": "Internals",
    "title": "PDMats.dim",
    "category": "method",
    "text": "dim(lat::Lattice)\ndim(model::Model)\n\nReturns the dimension of lattice.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.bonddirection-Tuple{Lattice,Integer}",
    "page": "Internals",
    "title": "SpinMonteCarlo.bonddirection",
    "category": "method",
    "text": "bonddirection(lat::Lattice, bond::Integer)\nbonddirection(model::Model, bond::Integer)\n\nReturns the direction of the bond as vector in the Cartesian system\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.bondtype-Tuple{Lattice,Integer}",
    "page": "Internals",
    "title": "SpinMonteCarlo.bondtype",
    "category": "method",
    "text": "bondtype(lat::Lattice, bond::Integer)\nbondtype(model::Model, bond::Integer)\n\nReturns the type of bond\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.lattice_bonddirection-Tuple{Lattice,Integer}",
    "page": "Internals",
    "title": "SpinMonteCarlo.lattice_bonddirection",
    "category": "method",
    "text": "lattice_bonddirection(lat::Lattice, bond::Integer)\nlattice_bonddirection(model::Model, bond::Integer)\n\nReturns the direction of the bond as vector in the lattice system\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.lattice_sitecoordinate-Tuple{Lattice,Integer}",
    "page": "Internals",
    "title": "SpinMonteCarlo.lattice_sitecoordinate",
    "category": "method",
    "text": "lattice_sitecoordinate(lat::Lattice, site::Integer)\nlattice_sitecoordinate(model::Model, site::Integer)\n\nReturns the coordinate of the site in the lattice system\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.neighbors-Tuple{Lattice,Integer}",
    "page": "Internals",
    "title": "SpinMonteCarlo.neighbors",
    "category": "method",
    "text": "neighbors(lat::Lattice, site::Integer)\nneighbors(model::Model, site::Integer)\n\nReturns the neighbor sites and bonds of site.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.numbonds-Tuple{Lattice,Integer}",
    "page": "Internals",
    "title": "SpinMonteCarlo.numbonds",
    "category": "method",
    "text": "numbonds(lat::Lattice, bondtype::Integer)\nnumbonds(model::Model, bondtype::Integer)\n\nReturns the number of bondtype bonds.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.numbonds-Tuple{Lattice}",
    "page": "Internals",
    "title": "SpinMonteCarlo.numbonds",
    "category": "method",
    "text": "numbonds(lat::Lattice)\nnumbonds(model::Model)\n\nReturns the number of all bonds.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.numsites-Tuple{Lattice,Integer}",
    "page": "Internals",
    "title": "SpinMonteCarlo.numsites",
    "category": "method",
    "text": "numsites(lat::Lattice, sitetype::Integer)\nnumsites(model::Model, sitetype::Integer)\n\nReturns the number of sitetype sites.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.numsites-Tuple{Lattice}",
    "page": "Internals",
    "title": "SpinMonteCarlo.numsites",
    "category": "method",
    "text": "numsites(lat::Lattice)\nnumsites(model::Model)\n\nReturns the number of all sites.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.sitecoordinate-Tuple{Lattice,Integer}",
    "page": "Internals",
    "title": "SpinMonteCarlo.sitecoordinate",
    "category": "method",
    "text": "sitecoordinate(lat::Lattice, site::Integer)\nsitecoordinate(model::Model, site::Integer)\n\nReturns the coordinate of the site in the Cartesian system\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.sitetype-Tuple{Lattice,Integer}",
    "page": "Internals",
    "title": "SpinMonteCarlo.sitetype",
    "category": "method",
    "text": "sitetype(lat::Lattice, site::Integer)\nsitetype(model::Model, site::Integer)\n\nReturns the type of site\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.source-Tuple{Lattice,Integer}",
    "page": "Internals",
    "title": "SpinMonteCarlo.source",
    "category": "method",
    "text": "source(lat::Lattice, bond::Integer)\nsource(model::Model, bond::Integer)\n\nReturns the source site of bond.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.target-Tuple{Lattice,Integer}",
    "page": "Internals",
    "title": "SpinMonteCarlo.target",
    "category": "method",
    "text": "target(lat::Lattice, bond::Integer)\ntarget(model::Model, bond::Integer)\n\nReturns the target site of bond.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.bonds-Tuple{Lattice,Integer}",
    "page": "Internals",
    "title": "SpinMonteCarlo.bonds",
    "category": "method",
    "text": "bonds(lat::Lattice, bondtype::Integer)\nbonds(model::Model, bondtype::Integer)\n\nReturns bonds with bondtype\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.neighborbonds-Tuple{Lattice,Integer}",
    "page": "Internals",
    "title": "SpinMonteCarlo.neighborbonds",
    "category": "method",
    "text": "neighborbonds(lat::Lattice, site::Integer)\nneighborbonds(model::Model, site::Integer)\n\nReturns the neighbor bonds of site.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.neighborsites-Tuple{Lattice,Integer}",
    "page": "Internals",
    "title": "SpinMonteCarlo.neighborsites",
    "category": "method",
    "text": "neighborsites(lat::Lattice, site::Integer)\nneighborsites(model::Model, site::Integer)\n\nReturns the neighbor sites of site.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.numbondtypes-Tuple{Lattice}",
    "page": "Internals",
    "title": "SpinMonteCarlo.numbondtypes",
    "category": "method",
    "text": "numbondtypes(lat::Lattice)\nnumbondtypes(model::Model)\n\nReturns the number of bondtypes.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.numsitetypes-Tuple{Lattice}",
    "page": "Internals",
    "title": "SpinMonteCarlo.numsitetypes",
    "category": "method",
    "text": "numsitetypes(lat::Lattice)\nnumsitetypes(model::Model)\n\nReturns the number of sitetypes.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.sites-Tuple{Lattice,Integer}",
    "page": "Internals",
    "title": "SpinMonteCarlo.sites",
    "category": "method",
    "text": "sites(lat::Lattice, sitetype::Integer)\nsites(model::Model, sitetype::Integer)\n\nReturns sites with sitetype\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#Lattice-1",
    "page": "Internals",
    "title": "Lattice",
    "category": "section",
    "text": "Modules = [SpinMonteCarlo]\nPages = [\"src/lattice.jl\"]"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.LoopElementType",
    "page": "Internals",
    "title": "SpinMonteCarlo.LoopElementType",
    "category": "type",
    "text": "Enumtype including LET_*\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.LET_Cut",
    "page": "Internals",
    "title": "SpinMonteCarlo.LET_Cut",
    "category": "constant",
    "text": "Loop element depicted as\n\n|\no\n|\n\nor matrix \n\n1 1 |+>\n1 1 |->\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.LET_FMLink",
    "page": "Internals",
    "title": "SpinMonteCarlo.LET_FMLink",
    "category": "constant",
    "text": "Loop element depicted as\n\n|  |\n|--|\n|  |\n\nor matrix \n\n1 0 0 0 |++>\n0 0 0 0 |+->\n0 0 0 0 |-+>\n0 0 0 1 |-->\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.LET_AFLink",
    "page": "Internals",
    "title": "SpinMonteCarlo.LET_AFLink",
    "category": "constant",
    "text": "Loop element depicted as\n\n|  |\n|~~|\n|  |\n\nor matrix \n\n0 0 0 0 |++>\n0 1 0 0 |+->\n0 0 1 0 |-+>\n0 0 0 0 |-->\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.LET_Vertex",
    "page": "Internals",
    "title": "SpinMonteCarlo.LET_Vertex",
    "category": "constant",
    "text": "Loop element depicted as\n\n|  |\n~~~~\n~~~~\n|  |\n\nor matrix \n\n0 0 0 0 |++>\n0 1 1 0 |+->\n0 1 1 0 |-+>\n0 0 0 0 |-->\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.LET_Cross",
    "page": "Internals",
    "title": "SpinMonteCarlo.LET_Cross",
    "category": "constant",
    "text": "Loop element depicted as\n\n|  |\n \\/\n /\\\n|  |\n\nor matrix \n\n1 0 0 0 |++>\n0 0 1 0 |+->\n0 1 0 0 |-+>\n0 0 0 1 |-->\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.LocalLoopOperator",
    "page": "Internals",
    "title": "SpinMonteCarlo.LocalLoopOperator",
    "category": "type",
    "text": "(Imaginary-temporary and spatial) local operator as a perturbation with assigned loop element.\n\nFields\n\nlet_type : assigned loop element\nisdiagonal : operator is diagonal or not\nin other words, two states connecting this perturbation are equivalent to each other or not.\ntime : imaginary time (taubeta in 01)) which this perturbation acts on.\nspace : index of subspin or subbond which this perturbation acts on.\nbottom_id :: index of node of union find assigned to a loop\ntop_id :: index of node of union find assigned to the other loop\n\n\n\n"
},

{
    "location": "lib/internals.html#Model-1",
    "page": "Internals",
    "title": "Model",
    "category": "section",
    "text": "LoopElementType\nLET_Cut\nLET_FMLink\nLET_AFLink\nLET_Vertex\nLET_Cross\nLocalLoopOperator"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.default_estimator",
    "page": "Internals",
    "title": "SpinMonteCarlo.default_estimator",
    "category": "function",
    "text": "default_estimator(model, updatemethod!)\n\nDetermines estimator to be used when param[\"Estimator\"] is not set.\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.@gen_convert_parameter",
    "page": "Internals",
    "title": "SpinMonteCarlo.@gen_convert_parameter",
    "category": "macro",
    "text": "@gen_convert_parameter(model_typename, (keyname, size_fn, default)...)\n\nGenerates convert_parameter(model::model_typename, param::Parameter).\n\nExample\n\n@gen_convert_parameter(A, (\"A\", numbondtypes, 1.0), (\"B\", 1, 1))\n\ngenerates a function equivalent to the following:\n\ndoc\"\"\"\n    convert_parameter(model::A, param::Parameter)\n\n# Keynames\n\n- \"A\": a vector with `numbondtypes(model)` elements (default: 1)\n- \"B\": a scalar (default: 1.0)\n\n\"\"\"\nfunction convert_parameter(model::A, param::Parameter)\n    ## if `size_fn` is a `Function`,\n    ## result is a vector whose size is `size_fn(model)`.\n    ## `param[\"A\"]` can take a scalar or a vector.\n    a = get(param, \"A\", 1.0)\n    as = zeros(Float64, numbondtypes(model))\n    as .= a\n\n    ## otherwise,\n    ## result is a scalar.\n    b = get(param, \"B\", 1) :: Int\n\n    return as, b\nend\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.UnionFind",
    "page": "Internals",
    "title": "SpinMonteCarlo.UnionFind",
    "category": "type",
    "text": "Union-find algorithm.\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.addnode!",
    "page": "Internals",
    "title": "SpinMonteCarlo.addnode!",
    "category": "function",
    "text": "addnode!(u::UnionFind)\n\nAdds a new node into u and returns the number of nodes including the added node.\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.unify!",
    "page": "Internals",
    "title": "SpinMonteCarlo.unify!",
    "category": "function",
    "text": "unify!(u, n1, n2)\n\nConnects n1 and n2 nodes and returns the root.\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.clusterize!",
    "page": "Internals",
    "title": "SpinMonteCarlo.clusterize!",
    "category": "function",
    "text": "clusterize!(u::UnionFind)\n\nAssigns cluster ID to each node and returns the number of clusters.\n\n\n\n"
},

{
    "location": "lib/internals.html#SpinMonteCarlo.clusterid",
    "page": "Internals",
    "title": "SpinMonteCarlo.clusterid",
    "category": "function",
    "text": "clusterid(u::UnionFind, i::Integer)\n\nReturns the index of the cluster where i node belongs.\n\n\n\n"
},

{
    "location": "lib/internals.html#Utility-1",
    "page": "Internals",
    "title": "Utility",
    "category": "section",
    "text": "default_estimator\n@gen_convert_parameter\nUnionFind\naddnode!\nunify!\nclusterize!\nclusterid"
},

]}
