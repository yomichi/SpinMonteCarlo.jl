## Behaviour tests for the snapshot feature.
## Every test file is `include`d into the same module, so the helpers below are
## prefixed with `snapshot_` to keep them out of the way of the other files.

"Runs `f(dir)` with the working directory moved into a fresh temporary directory."
function snapshot_in_tempdir(f)
    return mktempdir() do dir
        return cd(() -> f(dir), dir)
    end
end

"Returns the rendered message of the exception thrown by `f` (`\"\"` if it returns)."
function snapshot_errmsg(f)
    try
        f()
    catch err
        return sprint(showerror, err)
    end
    return ""
end

"Number of lines in `filename`, counting a file that does not exist as empty."
snapshot_countlines(filename) = isfile(filename) ? countlines(filename) : 0

"A minimal Ising parameter set (8 sites); the given pairs are merged in."
function snapshot_param(pairs::Pair...)
    param = Parameter("Model" => Ising,
                      "Lattice" => "chain lattice",
                      "L" => 8,
                      "J" => 1.0,
                      "T" => 1.5,
                      "Seed" => SEED,
                      "Update Method" => local_update!)
    for (key, value) in pairs
        param[key] = value
    end
    return param
end

snapshot_chain(L) = generatelattice(Parameter("Lattice" => "chain lattice", "L" => L))

@testset "snapshot(model)" begin
    L = 8
    Q = 5
    lat = snapshot_chain(L)
    N = numsites(lat)

    ## name, constructor, expected length, expected element type, value range
    models = (("Ising", () -> Ising(lat, SEED), N, Int, x -> x == 1 || x == -1),
              ("Potts", () -> Potts(lat, Q, SEED), N, Int, x -> 1 <= x <= Q),
              ("Clock", () -> Clock(lat, Q, SEED), N, Int, x -> 1 <= x <= Q),
              ("XY", () -> XY(lat, SEED), N, Float64, x -> 0.0 <= x < 1.0),
              ("AshkinTeller", () -> AshkinTeller(lat, SEED), 2N, Int,
               x -> x == 1 || x == -1))

    @testset "$name" for (name, make, len, T, inrange) in models
        model = make()
        conf = snapshot(model)
        @test conf isa Vector
        @test length(conf) == len
        @test eltype(conf) === T
        @test all(inrange, conf)
        @test conf == vec(model.spins)
    end

    @testset "AshkinTeller interleaves sigma and tau" begin
        model = AshkinTeller(snapshot_chain(4), SEED)
        model.spins .= [1 1 -1 -1
                        1 -1 1 -1]
        @test snapshot(model) == [1, 1, 1, -1, -1, 1, -1, -1]

        model.spins .= [-1 -1 -1 -1
                        1 1 1 1]
        @test snapshot(model) == [-1, 1, -1, 1, -1, 1, -1, 1]
    end

    @testset "$name shares no memory with the model" for (name, make, len, T,
                                                          inrange) in models

        model = make()
        conf = snapshot(model)
        expected = copy(conf)
        spins = copy(model.spins)

        ## writing into the returned vector must not reach the model
        conf[1] -= 100
        @test model.spins == spins

        ## mutating the model must not reach the already returned vector
        expected[1] -= 100
        model.spins .-= 100
        @test conf == expected
    end

    @testset "QuantumXXZ is rejected" begin
        model = QuantumXXZ(snapshot_chain(4), 1 // 2, SEED)
        @test_throws ArgumentError snapshot(model)
        @test occursin("ops", snapshot_errmsg(() -> snapshot(model)))
    end
end

@testset "save_snapshot" begin
    lat = snapshot_chain(6)
    model = Ising(lat, SEED)
    ## spelled out instead of calling `snapshot` so that these tests keep
    ## reporting on `save_snapshot` alone
    conf = vec(model.spins)

    @testset "one call writes one line" begin
        io = IOBuffer()
        save_snapshot(io, model)
        text = String(take!(io))
        @test count(==('\n'), text) == 1
        @test endswith(text, "\n")
        @test chomp(text) == join(conf, " ")
        @test parse.(Int, split(chomp(text))) == conf

        io = IOBuffer()
        for _ in 1:4
            save_snapshot(io, model)
        end
        text = String(take!(io))
        @test count(==('\n'), text) == 4
        @test all(==(join(conf, " ")), split(chomp(text), "\n"))
    end

    @testset "sep chooses the separator" begin
        io = IOBuffer()
        save_snapshot(io, model; sep=",")
        text = String(take!(io))
        @test count(==('\n'), text) == 1
        @test !occursin(" ", text)
        @test chomp(text) == join(conf, ",")
        @test parse.(Int, split(chomp(text), ",")) == conf
    end

    @testset "the filename form truncates unless append is set" begin
        snapshot_in_tempdir() do dir
            filename = "conf.txt"
            write(filename, "stale 1 2 3 4 5\nstale 1 2 3 4 5\n")

            save_snapshot(filename, model)
            @test countlines(filename) == 1
            @test !occursin("stale", read(filename, String))
            @test readlines(filename) == [join(conf, " ")]

            save_snapshot(filename, model)
            @test countlines(filename) == 1

            for _ in 1:3
                save_snapshot(filename, model; append=true)
            end
            @test countlines(filename) == 4
            @test all(==(join(conf, " ")), readlines(filename))

            save_snapshot(filename, model; append=true, sep=",")
            @test countlines(filename) == 5
            @test readlines(filename)[end] == join(conf, ",")
        end
    end
end

@testset "load_snapshots" begin
    @testset "each line of the file is one column" begin
        m = load_snapshots(IOBuffer("1 2 3 4\n5 6 7 8\n9 10 11 12\n"))
        @test m isa Matrix{Float64}
        @test size(m) == (4, 3)
        @test m[:, 1] == [1.0, 2.0, 3.0, 4.0]
        @test m[:, 3] == [9.0, 10.0, 11.0, 12.0]
        @test m == Float64[1 5 9
                           2 6 10
                           3 7 11
                           4 8 12]
    end

    @testset "the element type can be chosen" begin
        m = load_snapshots(Int, IOBuffer("1 -1 1\n-1 -1 1\n"))
        @test m isa Matrix{Int}
        @test m == [1 -1
                    -1 -1
                    1 1]
    end

    @testset "blank lines and comments are skipped" begin
        text = "# a header\n1 2 3\n\n   # an indented comment\n4 5 6\n\n"
        m = load_snapshots(IOBuffer(text))
        @test size(m) == (3, 2)
        @test m == Float64[1 4
                           2 5
                           3 6]
    end

    @testset "runs of whitespace are a single separator" begin
        m = load_snapshots(Int, IOBuffer("  11   12\t\t13  \n"))
        @test size(m) == (3, 1)
        @test m[:, 1] == [11, 12, 13]
    end

    @testset "an explicit sep drops empty fields" begin
        m = load_snapshots(Int, IOBuffer("1,2,,3,\n,4,5,6\n"); sep=",")
        @test size(m) == (3, 2)
        @test m == [1 4
                    2 5
                    3 6]
    end

    @testset "a source without data lines gives a 0x0 matrix" begin
        m = load_snapshots(IOBuffer(""))
        @test m isa Matrix{Float64}
        @test size(m) == (0, 0)

        m = load_snapshots(Int, IOBuffer("# nothing but a comment\n\n"))
        @test m isa Matrix{Int}
        @test size(m) == (0, 0)
    end

    @testset "a row of the wrong length is an error" begin
        text = repeat("11 12 13 14\n", 6) * "101 202\n"
        @test_throws ArgumentError load_snapshots(IOBuffer(text))
        msg = snapshot_errmsg(() -> load_snapshots(IOBuffer(text)))
        @test occursin(r"\b7\b", msg)   # the offending line number
        @test occursin(r"\b2\b", msg)   # the number of elements it has
    end

    @testset "a truncated last line is an error" begin
        text = repeat("11 12 13 14\n", 2) * "101 202"   # no trailing newline
        @test_throws ArgumentError load_snapshots(IOBuffer(text))
        msg = snapshot_errmsg(() -> load_snapshots(IOBuffer(text)))
        @test occursin(r"\b3\b", msg)
        @test occursin(r"\b2\b", msg)
    end

    @testset "round trip through save_snapshot" begin
        snapshot_in_tempdir() do dir
            model = Ising(snapshot_chain(8), SEED)
            Js = fill(1.0, SpinMonteCarlo.numbondtypes(model))
            filename = "ising.txt"
            confs = Vector{Int}[]
            for k in 1:3
                local_update!(model, 1.5, Js)
                push!(confs, snapshot(model))
                save_snapshot(filename, model; append=(k > 1))
            end

            m = load_snapshots(Int, filename)
            @test size(m) == (8, 3)
            for k in 1:3
                @test m[:, k] == confs[k]
            end
        end
    end

    @testset "XY values survive the round trip exactly" begin
        snapshot_in_tempdir() do dir
            model = XY(snapshot_chain(8), SEED)
            model.spins[1, :] .= [0.1, 0.2, 1 / 3, 2 / 7,
                                  0.0, nextfloat(0.0), prevfloat(1.0),
                                  0.6180339887498949]
            conf = snapshot(model)
            filename = "xy.txt"
            save_snapshot(filename, model)

            m = load_snapshots(filename)
            @test size(m) == (8, 1)
            @test isequal(m[:, 1], conf)
        end
    end
end

@testset "runMC writes snapshots" begin
    @testset "the feature is off by default" begin
        snapshot_in_tempdir() do dir
            runMC(snapshot_param("MCS" => 4, "Thermalization" => 2))
            @test !isfile("snapshot_0.txt")
            @test isempty(filter(f -> endswith(f, ".txt"), readdir()))

            runMC(snapshot_param("MCS" => 4, "Thermalization" => 2,
                                 "Snapshot Interval" => 0))
            @test !isfile("snapshot_0.txt")
            @test isempty(filter(f -> endswith(f, ".txt"), readdir()))
        end
    end

    @testset "the filename is made of the prefix and the ID" begin
        snapshot_in_tempdir() do dir
            runMC(snapshot_param("MCS" => 4, "Thermalization" => 2,
                                 "Snapshot Interval" => 2))
            @test isfile("snapshot_0.txt")
            @test countlines("snapshot_0.txt") == 2

            runMC(snapshot_param("MCS" => 4, "Thermalization" => 2,
                                 "Snapshot Interval" => 2,
                                 "Snapshot Filename Prefix" => "conf",
                                 "ID" => 7))
            @test isfile("conf_7.txt")
            @test countlines("conf_7.txt") == 2
        end
    end

    @testset "MCS=$mcs Therm=$therm interval=$interval" for (mcs, therm, interval,
                                                             expected) in
                                                            ((8, 4, 2, 4),
                                                             (9, 4, 2, 4),
                                                             (8, 0, 2, 4),
                                                             (8, 4, 1, 8),
                                                             (8, 4, 8, 1),
                                                             (4, 4, 8, 0))
        snapshot_in_tempdir() do dir
            p = snapshot_param("MCS" => mcs, "Thermalization" => therm,
                               "Snapshot Interval" => interval)
            @test_logs min_level = Test.Logging.Warn runMC(p)
            @test snapshot_countlines("snapshot_0.txt") == expected
            if expected > 0
                m = load_snapshots(Int, "snapshot_0.txt")
                @test size(m) == (8, expected)
                @test all(x -> x == 1 || x == -1, m)
            end
        end
    end

    @testset "only measurement steps are counted" begin
        snapshot_in_tempdir() do dir
            ## thermalization is longer than one interval: counting the
            ## thermalization steps as well would give div(11, 2) == 5 lines.
            p = snapshot_param("MCS" => 6, "Thermalization" => 5,
                               "Snapshot Interval" => 2)
            model = Ising(p)
            runMC(model, p)

            m = load_snapshots(Int, "snapshot_0.txt")
            @test size(m) == (8, 3)
            @test all(x -> x == 1 || x == -1, m)
            ## MCS is a multiple of the interval, so the last line is the
            ## configuration the run ended with.
            @test m[:, end] == snapshot(model)
        end
    end

    @testset "a fresh run empties an existing file" begin
        snapshot_in_tempdir() do dir
            write("snapshot_0.txt", "999 999 999 999 999 999 999 999\n" ^ 3)
            p = snapshot_param("MCS" => 4, "Thermalization" => 2,
                               "Snapshot Interval" => 2)
            runMC(p)
            @test countlines("snapshot_0.txt") == 2
            @test !occursin("999", read("snapshot_0.txt", String))
        end
    end

    @testset "MCS = 0 writes no line" begin
        snapshot_in_tempdir() do dir
            write("snapshot_0.txt", "999 999 999 999 999 999 999 999\n")
            p = snapshot_param("MCS" => 0, "Thermalization" => 4,
                               "Snapshot Interval" => 2)
            ## runMC itself cannot finish without a single measurement (binning
            ## an empty observable fails), which is none of this feature's
            ## business; all that is checked here is that no line was written.
            try
                runMC(p)
            catch
            end
            @test snapshot_countlines("snapshot_0.txt") == 0
            @test !occursin("999",
                            isfile("snapshot_0.txt") ?
                            read("snapshot_0.txt", String) : "")
        end
    end

    @testset "a negative interval is an error" begin
        snapshot_in_tempdir() do dir
            p = snapshot_param("MCS" => 4, "Thermalization" => 2,
                               "Snapshot Interval" => -1)
            @test_throws ArgumentError runMC(p)
        end
    end

    @testset "QuantumXXZ is rejected before the run starts" begin
        snapshot_in_tempdir() do dir
            tripwire = (model, args...) -> error("update! must not be called")
            p = Parameter("Model" => QuantumXXZ,
                          "Lattice" => "chain lattice",
                          "L" => 4,
                          "S" => 1 // 2,
                          "Jz" => 1.0,
                          "Jxy" => 1.0,
                          "T" => 1.0,
                          "Seed" => SEED,
                          "MCS" => 4,
                          "Thermalization" => 4,
                          "Update Method" => tripwire,
                          "Snapshot Interval" => 2)
            ## the tripwire turns "failed only after thermalization" into an
            ## ErrorException, which is not an ArgumentError.
            @test_throws ArgumentError runMC(p)
            @test !isfile("snapshot_0.txt")
        end
    end

    @testset "snapshots do not disturb the random number stream" begin
        snapshot_in_tempdir() do dir
            plain = runMC(snapshot_param("MCS" => 40, "Thermalization" => 10))
            snapped = runMC(snapshot_param("MCS" => 40, "Thermalization" => 10,
                                           "Snapshot Interval" => 3))
            @test countlines("snapshot_0.txt") == div(40, 3)
            @test issetequal(keys(plain), keys(snapped))
            for name in keys(plain)
                if name in ["MCS per Second", "Time per MCS"]
                    continue
                end
                @test isequal(mean(plain[name]), mean(snapped[name]))
                @test isequal(stderror(plain[name]), stderror(snapped[name]))
            end
        end
    end

    @testset "runMC over an array writes one file per ID" begin
        snapshot_in_tempdir() do dir
            params = [snapshot_param("MCS" => 6, "Thermalization" => 2,
                                     "Snapshot Interval" => 2) for _ in 1:3]
            runMC(params)
            @test !isfile("snapshot_0.txt")
            for id in 1:3
                filename = "snapshot_$(id).txt"
                @test isfile(filename)
                @test countlines(filename) == 3
                m = load_snapshots(Int, filename)
                @test size(m) == (8, 3)
                @test all(x -> x == 1 || x == -1, m)
            end
        end
    end
end

@testset "snapshot and checkpoint" begin
    @testset "restart reproduces an uninterrupted run (interval $cp)" for cp in (Inf,
                                                                                 1e-9)
        snapshot_in_tempdir() do dir
            function run_stages(subdir, mcss; extra=false)
                mkdir(subdir)
                return cd(subdir) do
                    for (i, mcs) in enumerate(mcss)
                        if extra && i == length(mcss)
                            ## the configurations a crash could have left behind
                            ## after the last checkpoint
                            open("snapshot_0.txt", "a") do io
                                println(io, "0 0 0 0 0 0 0 0")
                                return println(io, "1 1 1 1 1 1 1 1")
                            end
                        end
                        runMC(snapshot_param("MCS" => mcs, "Thermalization" => 4,
                                             "Snapshot Interval" => 2,
                                             "Checkpoint Interval" => cp))
                    end
                    return read("snapshot_0.txt")
                end
            end

            reference = run_stages("plain", (8,))
            @test count(==(UInt8('\n')), reference) == 4
            @test run_stages("restarted", (4, 8)) == reference
            @test run_stages("truncated", (4, 8); extra=true) == reference
        end
    end

    @testset "restarting with a changed schedule is rejected" begin
        @testset "Thermalization changed" begin
            snapshot_in_tempdir() do dir
                p = snapshot_param("MCS" => 8, "Thermalization" => 4,
                                   "Snapshot Interval" => 2,
                                   "Checkpoint Interval" => Inf)
                runMC(p)
                before = read("snapshot_0.txt")

                p["Thermalization"] = 6
                @test_throws ArgumentError runMC(p)
                msg = snapshot_errmsg(() -> runMC(p))
                @test occursin(r"\b4\b", msg)   # the saved value
                @test occursin(r"\b6\b", msg)   # the current one
                ## the point of the check is that the file is left alone
                @test read("snapshot_0.txt") == before
            end
        end

        @testset "Snapshot Interval changed" begin
            snapshot_in_tempdir() do dir
                p = snapshot_param("MCS" => 20, "Thermalization" => 4,
                                   "Snapshot Interval" => 2,
                                   "Checkpoint Interval" => Inf)
                runMC(p)

                p["Snapshot Interval"] = 5
                @test_throws ArgumentError runMC(p)
                msg = snapshot_errmsg(() -> runMC(p))
                @test occursin(r"\b2\b", msg)
                @test occursin(r"\b5\b", msg)
            end
        end

        @testset "an omitted Thermalization moves with MCS" begin
            snapshot_in_tempdir() do dir
                ## "Thermalization" defaults to MCS>>3, so it silently changes
                ## from 1 to 2 here.
                p = snapshot_param("MCS" => 8, "Snapshot Interval" => 2,
                                   "Checkpoint Interval" => Inf)
                runMC(p)
                @test countlines("snapshot_0.txt") == 4

                p["MCS"] = 16
                @test_throws ArgumentError runMC(p)
                msg = snapshot_errmsg(() -> runMC(p))
                @test occursin(r"\b1\b", msg)
                @test occursin(r"\b2\b", msg)
            end
        end

        @testset "MCS may grow but not shrink" begin
            snapshot_in_tempdir() do dir
                p = snapshot_param("MCS" => 10, "Thermalization" => 4,
                                   "Snapshot Interval" => 5,
                                   "Checkpoint Interval" => Inf)
                runMC(p)
                @test countlines("snapshot_0.txt") == 2

                p["MCS"] = 5
                @test_throws ArgumentError runMC(p)

                p["MCS"] = 20
                runMC(p)
                @test countlines("snapshot_0.txt") == 4
            end
        end
    end

    @testset "enabling snapshots on restart starts from an empty file" begin
        snapshot_in_tempdir() do dir
            p = snapshot_param("MCS" => 4, "Thermalization" => 4,
                               "Checkpoint Interval" => Inf)
            runMC(p)
            @test !isfile("snapshot_0.txt")

            p["MCS"] = 8
            p["Snapshot Interval"] = 2
            @test_logs (:warn,) match_mode = :any runMC(p)
            ## only the four measurement steps taken after the restart
            @test countlines("snapshot_0.txt") == 2
            m = load_snapshots(Int, "snapshot_0.txt")
            @test size(m) == (8, 2)
            @test all(x -> x == 1 || x == -1, m)
        end
    end

    @testset "disabling snapshots on restart leaves the file alone" begin
        snapshot_in_tempdir() do dir
            p = snapshot_param("MCS" => 4, "Thermalization" => 4,
                               "Snapshot Interval" => 2,
                               "Checkpoint Interval" => Inf)
            runMC(p)
            before = read("snapshot_0.txt")
            @test count(==(UInt8('\n')), before) == 2

            p["MCS"] = 8
            p["Snapshot Interval"] = 0
            runMC(p)
            @test isfile("snapshot_0.txt")
            @test read("snapshot_0.txt") == before
        end
    end

    @testset "re-enabling snapshots starts from an empty file again" begin
        snapshot_in_tempdir() do dir
            p = snapshot_param("MCS" => 4, "Thermalization" => 4,
                               "Snapshot Interval" => 2,
                               "Checkpoint Interval" => Inf)
            runMC(p)
            @test countlines("snapshot_0.txt") == 2

            p["MCS"] = 8
            p["Snapshot Interval"] = 0
            runMC(p)
            @test countlines("snapshot_0.txt") == 2

            p["MCS"] = 12
            p["Snapshot Interval"] = 2
            @test_logs (:warn,) match_mode = :any runMC(p)
            ## the measurements of the disabled stage have no lines, so the
            ## file must not be continued: only steps 9..12 are written.
            @test countlines("snapshot_0.txt") == 2
        end
    end

    @testset "a deleted snapshot file only warns" begin
        snapshot_in_tempdir() do dir
            p = snapshot_param("MCS" => 4, "Thermalization" => 4,
                               "Snapshot Interval" => 2,
                               "Checkpoint Interval" => Inf)
            runMC(p)
            rm("snapshot_0.txt")

            p["MCS"] = 8
            @test_logs (:warn,) match_mode = :any runMC(p)
            ## the lost lines are not restored, the new ones are written
            @test countlines("snapshot_0.txt") == 2
        end
    end

    @testset "a shortened snapshot file only warns" begin
        snapshot_in_tempdir() do dir
            p = snapshot_param("MCS" => 8, "Thermalization" => 4,
                               "Snapshot Interval" => 2,
                               "Checkpoint Interval" => Inf)
            runMC(p)
            @test countlines("snapshot_0.txt") == 4
            first_line = readlines("snapshot_0.txt")[1]
            write("snapshot_0.txt", first_line * "\n")

            p["MCS"] = 16
            @test_logs (:warn,) match_mode = :any runMC(p)
            ## the surviving line stays, the four new ones are appended
            @test countlines("snapshot_0.txt") == 5
            @test readlines("snapshot_0.txt")[1] == first_line
        end
    end

    @testset "truncate_snapshots!" begin
        function writelines(filename, n)
            open(filename, "w") do io
                for i in 1:n
                    println(io, "$i $i $i")
                end
            end
            return filename
        end

        @testset "extra lines are dropped" begin
            snapshot_in_tempdir() do dir
                filename = writelines("s.txt", 5)
                SpinMonteCarlo.truncate_snapshots!(filename, 3)
                @test readlines(filename) == ["1 1 1", "2 2 2", "3 3 3"]
                @test endswith(read(filename, String), "\n")
            end
        end

        @testset "an exact match is left alone" begin
            snapshot_in_tempdir() do dir
                filename = writelines("s.txt", 3)
                before = read(filename)
                SpinMonteCarlo.truncate_snapshots!(filename, 3)
                @test read(filename) == before
            end
        end

        @testset "too few lines warns and keeps what is there" begin
            snapshot_in_tempdir() do dir
                fn = writelines("s.txt", 2)
                before = read(fn)
                @test_logs((:warn,), match_mode = :any,
                           SpinMonteCarlo.truncate_snapshots!(fn, 5))
                @test read(fn) == before
            end
        end

        @testset "a missing file warns unless nothing is asked for" begin
            snapshot_in_tempdir() do dir
                @test_logs((:warn,), match_mode = :any,
                           SpinMonteCarlo.truncate_snapshots!("none.txt", 2))
                @test_logs(min_level = Test.Logging.Warn,
                           SpinMonteCarlo.truncate_snapshots!("none.txt", 0))
            end
        end

        @testset "truncating to zero empties the file" begin
            snapshot_in_tempdir() do dir
                filename = writelines("s.txt", 4)
                SpinMonteCarlo.truncate_snapshots!(filename, 0)
                @test isfile(filename)
                @test filesize(filename) == 0
            end
        end

        @testset "a partial last line does not survive" begin
            snapshot_in_tempdir() do dir
                filename = "s.txt"
                write(filename, "1 1 1\n2 2 2\n3 3")
                SpinMonteCarlo.truncate_snapshots!(filename, 2)
                @test read(filename, String) == "1 1 1\n2 2 2\n"
            end
        end

        @testset "appending after truncation continues the file" begin
            snapshot_in_tempdir() do dir
                model = Ising(snapshot_chain(3), SEED)
                conf = snapshot(model)
                filename = writelines("s.txt", 5)
                SpinMonteCarlo.truncate_snapshots!(filename, 2)
                save_snapshot(filename, model; append=true)

                lines = readlines(filename)
                @test length(lines) == 3
                @test lines[1:2] == ["1 1 1", "2 2 2"]
                @test lines[3] == join(conf, " ")

                m = load_snapshots(Int, filename)
                @test size(m) == (3, 3)
                @test m[:, 3] == conf
            end
        end
    end

    @testset "a checkpoint of the old format cannot be restored" begin
        snapshot_in_tempdir() do dir
            open("cp_0.dat", "w") do io
                return SpinMonteCarlo.serialize(io, Ising(snapshot_chain(8), SEED))
            end
            p = snapshot_param("MCS" => 4, "Thermalization" => 2,
                               "Snapshot Interval" => 2,
                               "Checkpoint Interval" => Inf)
            @test_throws Exception runMC(p)
        end
    end
end

@testset "removed snapshot API" begin
    @testset "$name" for (name, f) in (("gen_snapshot!", () -> gen_snapshot!()),
                                       ("gensave_snapshot!", () -> gensave_snapshot!()),
                                       ("load_snapshot", () -> load_snapshot()))
        @test_throws ErrorException f()
        msg = snapshot_errmsg(f)
        @test occursin("removed in v1.3", msg)
        ## the word boundaries keep the name of the removed function itself
        ## ("gensave_snapshot!") from passing as a pointer to the new API
        @test occursin("Snapshot Interval", msg) ||
              occursin(r"\bsave_snapshot\b", msg) ||
              occursin(r"\bload_snapshots\b", msg)
    end

    ## the stubs swallow any argument list
    @test_throws ErrorException gen_snapshot!(nothing, 1.0; MCS=1)
    @test_throws ErrorException gensave_snapshot!(IOBuffer(), nothing, 1.0)
    @test_throws ErrorException load_snapshot(IOBuffer("1 2 3"))

    ## the singular form points at the plural one
    @test occursin(r"\bload_snapshots\b", snapshot_errmsg(() -> load_snapshot()))
end
