export gen_snapshot!, gensave_snapshot!, load_snapshot

"""gen_snapshot! was removed in v1.3 because it had been broken since Julia 1.0."""
function gen_snapshot!(args...; kwargs...)
    return error("gen_snapshot! was removed in v1.3 because it had been broken since Julia 1.0. A visualization extension is planned as a replacement.")
end

"""gensave_snapshot! was removed in v1.3 because it had been broken since Julia 1.0."""
function gensave_snapshot!(args...; kwargs...)
    return error("gensave_snapshot! was removed in v1.3 because it had been broken since Julia 1.0. A visualization extension is planned as a replacement.")
end

"""load_snapshot was removed in v1.3 because it had been broken since Julia 1.0."""
function load_snapshot(args...; kwargs...)
    return error("load_snapshot was removed in v1.3 because it had been broken since Julia 1.0. A visualization extension is planned as a replacement.")
end
