import ExaPowerIO
import PowerIO

_case_path(name) = joinpath(@__DIR__, "..", "data", name)

_matches(lhs::Real, rhs::Real) = isapprox(lhs, rhs; rtol = 1e-9, atol = 1e-10, nans = true)
_matches(lhs, rhs) = lhs == rhs

function _test_matching_rows(actual, expected)
    @test length(actual) == length(expected)
    for (a, e) in zip(actual, expected), field in propertynames(e)
        @test _matches(getproperty(a, field), getproperty(e, field))
    end
end

function _test_matches_exapowerio(actual, expected)
    @test actual.baseMVA == [expected.baseMVA]
    for table in (:bus, :gen, :branch, :arc, :storage)
        _test_matching_rows(getproperty(actual, table), getproperty(expected, table))
    end
end

# Backend-free: everything here runs on the serial path, so `runtests` calls it from
# one slice rather than letting it repeat in every CI job.
function powerio_parser_tests()
    @testset "PowerIO parser backend" begin
        matpower_cases = [
            "pglib_opf_case3_lmbd.m",
            "pglib_opf_case5_pjm.m",
            "pglib_opf_case14_ieee.m",
            "pglib_opf_case3_lmbd_mod.m",
            "pglib_opf_case5_pjm_mod.m",
        ]

        for case in matpower_cases
            path = _case_path(case)
            expected = ExaPowerIO.parse_matpower(Float64, path)
            actual = ExaModelsPower.parse_ac_power_data(path, Float64)
            _test_matches_exapowerio(actual, expected)
        end

        # Every model here is built from the parse, and a single abstract field in it
        # leaves `ExaModel(CORE, data)` unresolved under `--trim=safe`. The row types
        # must also be isbits, or the GPU backends refuse the arrays built from them.
        #
        # Broken against PowerIO 0.9.0: its table accessors return the JSON3 value
        # union, so containers built from them infer as `Any`, and `T` arrives as a
        # keyword-typed `Type{<:Real}` that nothing downstream can specialize on. The
        # rows ARE concrete at run time, which the isbits checks below assert; only
        # inference cannot prove it. Flips to a failure when PowerIO pins them.
        let path = _case_path("pglib_opf_case14_ieee.m")
            @test_broken isconcretetype(
                Base.infer_return_type(ExaModelsPower.parse_ac_power_data, (String, Type{Float64})),
            )
            p = ExaModelsPower.parse_ac_power_data(path, Float64)
            for table in (:bus, :gen, :branch, :arc)
                @test isbitstype(eltype(getproperty(p, table)))
            end
        end

        # The multiperiod busarray is built immutably from the LoadSeries; its element
        # type must stay concrete/isbits so the GPU kernels accept it (mirrors the
        # SCOPF boundary test, and surfaces a PowerIO regression without a GPU runner).
        let mp_net = ExaModelsPower.parse_mp_network(_case_path("pglib_opf_case3_lmbd.m"))
            mp_series = PowerIO.LoadSeries(mp_net, [1.0, 0.95]; T = Float64)
            mp_data = ExaModelsPower.parse_mp_power_data(mp_net, mp_series, 2, 0.1, Float64)
            @test isconcretetype(eltype(mp_data.busarray))
            @test isbitstype(eltype(mp_data.busarray))

            # A horizon longer than the series is a clear error rather than a
            # BoundsError from inside the comprehension.
            @test_throws DimensionMismatch ExaModelsPower.parse_mp_power_data(
                mp_net, mp_series, 3, 0.1, Float64)
        end

        # A load curve scales the LOADS. A fixed bus shunt is fixed: case14 carries
        # bs = 19.0 at bus 9, and it is the same in every period.
        let net = ExaModelsPower.parse_mp_network(_case_path("pglib_opf_case14_ieee.m"))
            series = PowerIO.LoadSeries(net, [1.0, 0.5]; T = Float64)
            d = ExaModelsPower.parse_mp_power_data(net, series, 2, 0.1, Float64)
            shunt = [b.b.bs for b in d.busarray[:, 1]]
            @test any(!iszero, shunt)
            @test shunt == [b.b.bs for b in d.busarray[:, 2]]
            @test [b.b.pd for b in d.busarray[:, 2]] ≈ 0.5 .* [b.b.pd for b in d.busarray[:, 1]]
        end

        mktempdir() do dir
            matpower = _case_path("pglib_opf_case3_lmbd.m")
            psse_path = joinpath(dir, "case3.raw")
            powerworld_path = joinpath(dir, "case3.aux")
            powermodels_path = joinpath(dir, "case3.json")

            write(psse_path, first(PowerIO.convert_file(matpower, "psse")))
            write(powerworld_path, first(PowerIO.convert_file(matpower, "powerworld")))
            write(powermodels_path, first(PowerIO.convert_file(matpower, "powermodels-json")))

            for path in (psse_path, powerworld_path)
                data = ExaModelsPower.parse_ac_power_data(path, Float64)
                @test length(data.bus) == 3
                @test length(data.branch) == 3
                @test length(data.gen) == 3
            end

            pm_data = ExaModelsPower.parse_ac_power_data(
                powermodels_path,
                Float64;
                from = "powermodels",
            )
            @test length(pm_data.bus) == 3
            @test length(pm_data.branch) == 3
            @test length(pm_data.gen) == 3

            @test !isnothing(ac_opf_model(psse_path; form = Polar())[1])
            @test !isnothing(ac_opf_model(powerworld_path; form = Rect())[1])
            @test !isnothing(dcopf_model(psse_path)[1])
            @test !isnothing(
                mpopf_model(
                    powermodels_path,
                    [1.0, 0.95];
                    from = "powermodels",
                )[1],
            )

            # Passing only one of pd/qd is a clear error, not a MethodError from a
            # `nothing` reaching PowerIO.LoadSeries. Both are read from file, or both
            # are supplied as matrices; the files here are never read (the error is
            # raised before load construction).
            @test_throws ArgumentError mpopf_model(
                matpower, "unused.Pd", "unused.Qd"; pd = zeros(3, 2))
            @test_throws ArgumentError mpopf_model(
                matpower, "unused.Pd", "unused.Qd"; qd = zeros(3, 2))
        end
    end
end
