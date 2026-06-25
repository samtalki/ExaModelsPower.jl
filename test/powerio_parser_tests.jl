import ExaPowerIO
import PowerIO

function _case_path(name)
    return joinpath(@__DIR__, "..", "data", name)
end

function _matches(lhs, rhs)
    if lhs isa Real && rhs isa Real
        return isapprox(lhs, rhs; rtol = 1e-9, atol = 1e-10, nans = true)
    end
    return lhs == rhs
end

function _test_matching_fields(lhs, rhs, fields)
    for field in fields
        @test _matches(getproperty(lhs, field), getproperty(rhs, field))
    end
end

function _test_matching_rows(actual, expected, fields)
    @test length(actual) == length(expected)
    for i in 1:min(length(actual), length(expected))
        _test_matching_fields(actual[i], expected[i], fields)
    end
end

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

        @test actual.baseMVA == [expected.baseMVA]
        _test_matching_rows(actual.bus, expected.bus,
            (:i, :bus_i, :type, :pd, :qd, :gs, :bs, :area, :vm, :va,
             :baseKV, :zone, :vmax, :vmin))
        _test_matching_rows(actual.gen, expected.gen,
            (:i, :bus, :pg, :qg, :qmax, :qmin, :vg, :mbase, :status,
             :pmax, :pmin, :model_poly, :startup, :shutdown, :n, :c))
        _test_matching_rows(actual.branch, expected.branch,
            (:i, :f_bus, :t_bus, :br_r, :br_x, :b_fr, :b_to, :g_fr,
             :g_to, :rate_a, :rate_b, :rate_c, :tap, :shift, :status,
             :angmin, :angmax, :f_idx, :t_idx, :c1, :c2, :c3, :c4,
             :c5, :c6, :c7, :c8))
        _test_matching_rows(actual.arc, expected.arc, (:i, :bus, :rate_a))
        _test_matching_rows(actual.storage, expected.storage,
            (:i, :storage_bus, :Pexts, :Qexts, :energy, :energy_rating,
             :charge_rating, :discharge_rating, :charge_efficiency,
             :discharge_efficiency, :thermal_rating, :qmin, :qmax, :Zr,
             :Zim, :p_loss, :q_loss, :status))
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

        pm_data = ExaModelsPower.parse_ac_power_data(powermodels_path, Float64; from = "powermodels")
        @test length(pm_data.bus) == 3
        @test length(pm_data.branch) == 3
        @test length(pm_data.gen) == 3

        @test !isnothing(ac_opf_model(psse_path; form = :polar)[1])
        @test !isnothing(ac_opf_model(powerworld_path; form = :rect)[1])
        @test !isnothing(dcopf_model(psse_path)[1])
        @test !isnothing(mpopf_model(powermodels_path, [1.0, 0.95]; from = "powermodels")[1])
    end
end
