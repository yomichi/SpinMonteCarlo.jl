export snapshot, save_snapshot, load_snapshots
export gen_snapshot!, gensave_snapshot!, load_snapshot

@doc doc"""
    snapshot(model::Model)

Return a vector containing the model's current spin configuration. The returned vector
does not share memory with the model.

The generic implementation requires `model.spins` to be an
`AbstractMatrix{<:Number}`. It follows Julia's column-major order. In particular,
Ashkin-Teller spins are ordered as `σ₁ τ₁ σ₂ τ₂ …`; the other classical
models have one value per site. For the XY model, each value is the internal
`σ ∈ [0, 1)` representation, with the physical angle given by `θ = 2πσ`.
"""
function snapshot(model::Model)
    return copy(vec(model.spins))
end

function snapshot(model::QuantumXXZ)
    return throw(ArgumentError("QuantumXXZ snapshots are not supported because both spins " *
                               "and ops are required to represent its state."))
end

@doc doc"""
    save_snapshot(io::IO, model::Model; sep=" ")
    save_snapshot(filename::AbstractString, model::Model; sep=" ", append=false)

Write the model's current spin configuration as one line of separator-delimited text.
The filename form truncates an existing file unless `append=true`.
"""
function save_snapshot(io::IO, model::Model; sep::AbstractString=" ")
    return println(io, join(snapshot(model), sep))
end

function save_snapshot(filename::AbstractString, model::Model;
                       sep::AbstractString=" ", append::Bool=false)
    mode = append ? "a" : "w"
    return open(filename, mode) do io
        return save_snapshot(io, model; sep=sep)
    end
end

@doc doc"""
    load_snapshots(source; sep=nothing)
    load_snapshots(T::Type{<:Number}, source; sep=nothing)

Load separator-delimited snapshot configurations from an `IO` or filename. Each data
line becomes one column of the returned `snapshot_length × nconfigs` matrix. Lines
whose first non-whitespace character is `#` and lines that split into no tokens are
ignored. The default element type is `Float64`.
"""
load_snapshots(source; sep=nothing) = load_snapshots(Float64, source; sep=sep)

function load_snapshots(T::Type{<:Number}, io::IO; sep=nothing)
    configs = Vector{T}[]
    expected = nothing
    for (line_number, line) in enumerate(eachline(io))
        startswith(lstrip(line), '#') && continue
        tokens = isnothing(sep) ? split(line) : split(line, sep; keepempty=false)
        isempty(tokens) && continue
        if !isnothing(expected) && length(tokens) != expected
            throw(ArgumentError("line $line_number has $(length(tokens)) elements; " *
                                "expected $expected elements"))
        end
        isnothing(expected) && (expected = length(tokens))
        push!(configs, parse.(T, tokens))
    end
    isempty(configs) && return Matrix{T}(undef, 0, 0)
    return reduce(hcat, configs)
end

function load_snapshots(T::Type{<:Number}, filename::AbstractString; sep=nothing)
    return open(filename, "r") do io
        return load_snapshots(T, io; sep=sep)
    end
end

function truncate_snapshots!(filename::AbstractString, nlines::Integer)
    if !isfile(filename)
        nlines > 0 &&
            @warn "Snapshot file $filename does not exist; expected $nlines lines."
        return nothing
    end

    open(filename, "r+") do io
        for _ in 1:nlines
            if eof(io)
                @warn "Snapshot file $filename has fewer than $nlines lines."
                return nothing
            end
            line = readuntil(io, '\n'; keep=true)
            if isempty(line) || line[end] != '\n'
                @warn "Snapshot file $filename has fewer than $nlines lines."
                return nothing
            end
        end
        return truncate(io, position(io))
    end
    return nothing
end

"""gen_snapshot! was removed in v1.3; use the snapshot APIs instead."""
function gen_snapshot!(args...; kwargs...)
    return error("gen_snapshot! was removed in v1.3. Set param[\"Snapshot Interval\"] " *
                 "to write spin configurations from runMC, or call " *
                 "save_snapshot(io, model) directly.")
end

"""gensave_snapshot! was removed in v1.3; use the snapshot APIs instead."""
function gensave_snapshot!(args...; kwargs...)
    return error("gensave_snapshot! was removed in v1.3. Set " *
                 "param[\"Snapshot Interval\"] to write spin configurations from " *
                 "runMC, or call save_snapshot(io, model) directly.")
end

"""load_snapshot was removed in v1.3; use `load_snapshots` instead."""
function load_snapshot(args...; kwargs...)
    return error("load_snapshot was removed in v1.3. Use load_snapshots to read " *
                 "snapshot configurations.")
end
