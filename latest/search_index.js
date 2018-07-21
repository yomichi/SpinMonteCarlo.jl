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
    "location": "develop.html#",
    "page": "Develop your program",
    "title": "Develop your program",
    "category": "page",
    "text": ""
},

{
    "location": "develop.html#Develop-your-program-1",
    "page": "Develop your program",
    "title": "Develop your program",
    "category": "section",
    "text": ""
},

{
    "location": "develop.html#Lattice-1",
    "page": "Develop your program",
    "title": "Lattice",
    "category": "section",
    "text": ""
},

{
    "location": "develop.html#Model-1",
    "page": "Develop your program",
    "title": "Model",
    "category": "section",
    "text": "Model should contain following fields: lat :: Lattice and rng :: Random.MersenneTwister. Model also should have a constructor taking param :: Parameter as a argument.You should define convert_parameter for your Model. This is a helper function which takes Model and Parameter and returns arguments of update method and estimator.@gen_convert_parameter helps you to define convert_parameter. For example, if your model::A needs a scalar T = param[\"T\"] and a vector Js = param[\"J\"] with numbondtypes(model) elements,@gen_convert_parameter(A, (\"T\", 1, 1.0), (\"J\", numbondtypes, 1.0))defines your documented and type-stable  convert_parameter(model::A, param::Parameter) which returns T and Js."
},

{
    "location": "develop.html#Note:-1",
    "page": "Develop your program",
    "title": "Note:",
    "category": "section",
    "text": "That the second element of each tuple is not a Function means that a return value is a scalar (the case of \"T\").\nThe third element of each tuple is the default value."
},

{
    "location": "develop.html#Update-method-1",
    "page": "Develop your program",
    "title": "Update method",
    "category": "section",
    "text": "\"Update method\" is a function which (in-place) updates model::Model under some parameters such as temperature T. For example, local_update!(model::Ising, T, Js) updates a spin configuration of model by local spin flip and Metropolice-Hasting algorithm under temperature T and coupling constants Js. \"Update method\" can return some object which will be used in \"Estimator\" as extra information."
},

{
    "location": "develop.html#Example-1",
    "page": "Develop your program",
    "title": "Example",
    "category": "section",
    "text": "Swendsen-Wang algorithm SW_update!(::Ising, ::Parameter) returns cluster information sw::SWInfo used in improved_estimator."
},

{
    "location": "develop.html#Estimator-1",
    "page": "Develop your program",
    "title": "Estimator",
    "category": "section",
    "text": "\"Estimator\" is a function which returns observables of a configuration model as a Dict{String, Any}. Arguments of a \"Estimator\" are model::Model, parameters (return of convert_parameter), and extra information (return of \"Update method\") in order.Default \"Estimator\" is determined by model and \"Update Method\" as return of default_estimator(model, param[\"Update Method\"])."
},

{
    "location": "develop.html#Example-2",
    "page": "Develop your program",
    "title": "Example",
    "category": "section",
    "text": "improved_estimator(model::Ising, T::Real, Js::AbstractArray, sw::SWInfo) takes a return of SW_update!(model::Ising, T::Real, Js::AbstractArray) as the last argument."
},

{
    "location": "develop.html#Postprocess-1",
    "page": "Develop your program",
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
    "text": "Documentation for SpinMonteCarlo.jl\'s public interface."
},

{
    "location": "lib/public.html#Driver-1",
    "page": "Public",
    "title": "Driver",
    "category": "section",
    "text": "CurrentModule = SpinMonteCarlorunMC"
},

{
    "location": "lib/public.html#Model-1",
    "page": "Public",
    "title": "Model",
    "category": "section",
    "text": "Ising\nPotts\nClock\nXY\nQuantumXXZ"
},

{
    "location": "lib/public.html#Lattice-1",
    "page": "Public",
    "title": "Lattice",
    "category": "section",
    "text": "Modules = [SpinMonteCarlo]\nPages = [\"src/lattice.jl\"]"
},

{
    "location": "lib/public.html#Lattice-generator-1",
    "page": "Public",
    "title": "Lattice generator",
    "category": "section",
    "text": "dimer_lattice\nchain_lattice\nsquare_lattice\ntriangular_lattice\ncubic_lattice\nfully_connected_lattice"
},

{
    "location": "lib/public.html#Update-method-1",
    "page": "Public",
    "title": "Update method",
    "category": "section",
    "text": "An index of model parameter (e.g., Js) is corresponding to sitetype or bondtype.local_update!\nSW_update!\nWolff_update!\nloop_update!"
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
    "text": "convert_parameter(model::Union{Ising, Potts, Clock, XY}, param::Parameter)\n\nKeynames:\n\n\"T\": a scalar (default: 1.0).\n\"J\": a vector with numbondtypes(model) elements (default: 1.0).\n\n\n\nconvert_parameter(model::QuantumXXZ, param::Parameter)\n\nKeynames:\n\n\"T\": a scalar (default: 1.0).\n\"Jz\": a vector with numbondtypes(model) elements (default: 1.0).\n\"Jxy\": a vector with numbondtypes(model) elements (default: 1.0).\n\"Gamma\": a vector with numsitetypes(model) elements (default: 0.0).\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.convert_parameter-Tuple{SpinMonteCarlo.Ising,Dict{String,Any}}",
    "page": "Public",
    "title": "SpinMonteCarlo.convert_parameter",
    "category": "method",
    "text": "convert_parameter(model::Union{Ising, Potts, Clock, XY}, param::Parameter)\n\nKeynames:\n\n\"T\": a scalar (default: 1.0).\n\"J\": a vector with numbondtypes(model) elements (default: 1.0).\n\n\n\n"
},

{
    "location": "lib/public.html#SpinMonteCarlo.convert_parameter-Tuple{SpinMonteCarlo.QuantumXXZ,Dict{String,Any}}",
    "page": "Public",
    "title": "SpinMonteCarlo.convert_parameter",
    "category": "method",
    "text": "convert_parameter(model::QuantumXXZ, param::Parameter)\n\nKeynames:\n\n\"T\": a scalar (default: 1.0).\n\"Jz\": a vector with numbondtypes(model) elements (default: 1.0).\n\"Jxy\": a vector with numbondtypes(model) elements (default: 1.0).\n\"Gamma\": a vector with numsitetypes(model) elements (default: 0.0).\n\n\n\n"
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
    "text": "Documentation for SpinMonteCarlo.jl\'s internals."
},

{
    "location": "lib/internals.html#Driver-1",
    "page": "Internals",
    "title": "Driver",
    "category": "section",
    "text": "CurrentModule = SpinMonteCarloaccumulateObservables!\npostproc"
},

{
    "location": "lib/internals.html#Model-1",
    "page": "Internals",
    "title": "Model",
    "category": "section",
    "text": "LoopElementType\nLET_Cut\nLET_FMLink\nLET_AFLink\nLET_Vertex\nLET_Cross\nLocalLoopOperator"
},

{
    "location": "lib/internals.html#Utility-1",
    "page": "Internals",
    "title": "Utility",
    "category": "section",
    "text": "default_estimator\n@gen_convert_parameter\nUnionFind\naddnode!\nunify!\nclusterize!\nclusterid"
},

]}
