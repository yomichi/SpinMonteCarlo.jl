using Serialization
using Distributed
import Distributed.pmap

export runMC

@doc """
    runMC(param::Parameter)
    runMC(params::AbstractArray{Parameter}
          ;
          parallel::Bool=false,
          autoID::Bool=true)

Runs Monte Carlo simulation(s) and returns calculated observables.

# Restart

If a checkpoint file named `"\$(param["Checkpoint Filename Prefix"])_\$(param["ID"]).dat"` exists and
`param["Checkpoint Interval"] > 0.0`, `runMC` loads this file and restarts the pending simulation.
NOTE: Restart will fail if the version or the system image of julia change (see the doc of `Serialization.serialize` ).

# Keyward arguments
- `autoID`: If true, `"ID"`s will be set (overwritten) as `params[i]["ID"] = i`.
  Because a child seed is derived from `("Seed", "ID")`, this is what makes a set of
  simulations sharing one `"Seed"` use distinct random number streams.
- `parallel`: If true, runs simulations in parallel (uses `pmap` instead of `map`).

# Required keys in `param`
- "Model"
- "Lattice"
- "Update Method"

# Optional keys in `param`
- "MCS": The number of Monte Carlo steps after thermalization
    - Default: `8192`
- "Thermalization": The number of Monte Carlo steps for thermalization
    - Default: `MCS>>3`
- "Binning Size": The size of binning
    - Default: `0`
- "Number of Bins": The number of bins
    - Default: `0`
    - If both "Binning Size" and "Number of Bins" are not given, "Binning Size" is set to `floor(sqrt(MCS))`.
- "Seed": The initial seed of the random number generator.
    - Default: determined randomly (the generator is seeded from system entropy)
    - When "ID" is also given, a child seed is derived from the pair ("Seed", "ID"),
      so simulations sharing one seed still use distinct random number streams.
      The derivation is arithmetic and reproduces across julia versions; it requires
      "Seed" to be an integer.
- "RNG": The *type* (not an instance) of the random number generator.
    - Default: `Random.Xoshiro`
    - Only the constructor the other keys call for is needed: `T(seed)` when
      "Seed" is given, `T()` when it is not. `Random.RandomDevice` and
      `Random.TaskLocalRNG` provide only the latter, so they work only without "Seed".
    - `Random.MersenneTwister` selects the generator used before v1.3.
- "Checkpoint Filename Prefix": See the "Restart" section.
    - Default: `"cp"`
- "ID": Job ID. Used both for the checkpoint filename (see the "Restart" section)
  and, together with "Seed", for deriving a per-job random number stream.
    - Default: `0`. Note that no child seed is derived when the key itself is absent,
      so a single `runMC(param)` without "ID" depends on "Seed" alone.
- "Checkpoint Interval": Time interval between writing checkpoint file in seconds.
    - Default: `0.0`, this means that NO checkpoint file will be loaded and saved.
- "Snapshot Interval": Number of measurement MCS between snapshot writes. This is
    measured in MCS, not seconds. A value of `0` disables snapshots.
    - Default: `0`
    - When restarting with snapshots enabled before and after the checkpoint, changing
      this value or "Thermalization" is rejected to preserve the snapshot schedule, as
      is shrinking "MCS" below the progress the checkpoint holds. Parameters driving
      the chain itself ("Update Method", "T", "J", ...) are not checked: restarting
      with a different temperature is allowed and produces snapshots that do not
      correspond to a single simulation.
- "Snapshot Filename Prefix": Prefix for snapshot files. The filename is
    `"<prefix>_<ID>.txt"`.
    - Default: `"snapshot"`
"""
function runMC(params::AbstractArray{T}; parallel::Bool=false,
               autoID::Bool=true) where {T<:Dict}
    map_fn = ifelse(parallel, pmap, map)
    return map_fn(enumerate(params)) do (id, p)
        if autoID
            p["ID"] = id
        end
        return runMC(p)
    end
end

function runMC(param::Parameter)
    model = param["Model"](param)
    return runMC(model, param)
end

function runMC(model, param::Parameter)
    verbose = get(param, "Verbose", false)::Bool
    if verbose
        println("Start: ", param)
    end
    cp_filename = @sprintf("%s_%d.dat",
                           get(param, "Checkpoint Filename Prefix", "cp")::String,
                           get(param, "ID", 0)::Int)
    cp_interval = get(param, "Checkpoint Interval", 0.0)::Float64
    snapshot_interval = get(param, "Snapshot Interval", 0)::Int
    snapshot_filename = @sprintf("%s_%d.txt",
                                 get(param, "Snapshot Filename Prefix", "snapshot")::String,
                                 get(param, "ID", 0)::Int)
    tm = time()

    MCS = get(param, "MCS", 8192)::Int
    Therm = get(param, "Thermalization", MCS >> 3)::Int

    mcs = 0
    MCS += Therm
    # obs = BinningObservableSet()
    obs = SimpleObservableSet()
    makeMCObservable!(obs, "Time per MCS")
    makeMCObservable!(obs, "MCS per Second")

    restarted = false
    saved_snapshot_state = nothing
    if cp_interval > 0.0 && ispath(cp_filename)
        open(cp_filename) do io
            model = deserialize(io)
            obs = deserialize(io)
            mcs = deserialize(io)
            return saved_snapshot_state = deserialize(io)
        end
        restarted = true
    end

    snapshot_interval < 0 &&
        throw(ArgumentError("Snapshot Interval must be nonnegative; got " *
                            "$snapshot_interval."))
    if snapshot_interval > 0
        snapshot(model)
    end

    nwritten::Int = 0
    if restarted
        if snapshot_interval > 0 && mcs > MCS
            throw(ArgumentError("cannot shrink \"MCS\" while snapshots are enabled: the " *
                                "checkpoint has already reached step $mcs, past the $MCS " *
                                "steps (\"Thermalization\" $Therm + \"MCS\" " *
                                "$(MCS - Therm)) this run would take."))
        end
        if !isnothing(saved_snapshot_state) && snapshot_interval > 0
            saved_therm = saved_snapshot_state.therm
            saved_interval = saved_snapshot_state.interval
            if saved_therm != Therm
                throw(ArgumentError("Thermalization changed from saved value " *
                                    "$saved_therm to current value $Therm."))
            end
            if saved_interval != snapshot_interval
                throw(ArgumentError("Snapshot Interval changed from saved value " *
                                    "$saved_interval to current value " *
                                    "$snapshot_interval."))
            end
            nwritten = saved_snapshot_state.nwritten
            truncate_snapshots!(snapshot_filename, nwritten)
        elseif isnothing(saved_snapshot_state) && snapshot_interval > 0
            @warn "Snapshots were enabled while restarting; starting a new snapshot file."
            open(snapshot_filename, "w") do io
            end
        end
    elseif snapshot_interval > 0
        open(snapshot_filename, "w") do io
        end
    end

    update! = param["Update Method"]::Function
    if haskey(param, "Estimator")
        estimator = param["Estimator"]::Function
    else
        estimator = default_estimator(model, update!)
    end
    pp = get(param, "Post Process", postproc)::Function
    p = convert_parameter(model, param)

    while mcs < MCS
        if mcs < Therm
            update!(model, p...)
        else
            t = @elapsed begin
                st = update!(model, p...)
                localobs = estimator(model, p..., st)
            end
            obs["Time per MCS"] << t
            obs["MCS per Second"] << (1.0 / t)
            accumulateObservables!(model, obs, localobs)
        end
        mcs += 1
        if snapshot_interval > 0 && mcs > Therm &&
           (mcs - Therm) % snapshot_interval == 0
            open(snapshot_filename, "a") do io
                return save_snapshot(io, model)
            end
            nwritten += 1
        end
        if cp_interval > 0.0 && time() - tm > cp_interval
            open(cp_filename, "w") do io
                serialize(io, model)
                serialize(io, obs)
                serialize(io, mcs)
                snapshot_state = snapshot_interval > 0 ?
                                 (therm=Therm, interval=snapshot_interval,
                                  nwritten=nwritten) : nothing
                return serialize(io, snapshot_state)
            end
            tm += cp_interval
        end
    end

    if cp_interval > 0.0
        open(cp_filename, "w") do io
            serialize(io, model)
            serialize(io, obs)
            serialize(io, mcs)
            snapshot_state = snapshot_interval > 0 ?
                             (therm=Therm, interval=snapshot_interval,
                              nwritten=nwritten) : nothing
            return serialize(io, snapshot_state)
        end
    end

    binsize = get(param, "Binning Size", 0)::Int
    numbins = get(param, "Number of Bins", 0)::Int
    binned = binning(obs; binsize=binsize, numbins=numbins)

    jk = pp(model, param, binned)

    if verbose
        println("Finish: ", param)
    end
    return jk
end

@doc """
    accumulateObservables!(model, obs::MCObservableSet, localobs::Dict)

Accumulates `localobs` into `obs`. For example, `obs["Energy"] << localobs["Energy"]`.
"""
function accumulateObservables!(::Model, obs::MCObservableSet, localobs::Measurement)
    if length(obs) < 3
        @inbounds for key in keys(localobs)
            makeMCObservable!(obs, key)
            obs[key] << localobs[key]
        end
    else
        @inbounds for key in keys(localobs)
            obs[key] << localobs[key]
        end
    end
    return obs
end
