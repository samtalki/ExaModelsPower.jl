import JSON
import PowerIO

# On a CUDA backend, MadNLP's defaults (SparseKKTSystem + MUMPS) are CPU-only and cannot
# assemble the KKT matrix from device arrays, so solve the condensed system with cuDSS on
# the GPU. On the CPU, pin the linear solver to Umfpack: MadNLP 0.10 changed the sparse
# default to MUMPS, which drives these OPFs into restoration, whereas Umfpack (the default
# through MadNLP 0.8, and what main uses) converges cleanly.
function exasolve(model, backend; kwargs...)
    opts = backend isa CUDABackend ?
        (; kkt_system = MadNLP.SparseCondensedKKTSystem, linear_solver = MadNLPGPU.CUDSSSolver) :
        (; linear_solver = MadNLP.UmfpackSolver)
    return madnlp(model; opts..., kwargs...)
end

function test_case3(result, result_pm, result_nlp_pm, pg, qg, p, q)
    test_static_case(result, result_pm, result_nlp_pm, pg, qg)

    #Branches are encoded differently in solutions, so matches are hard coded
    vars_dict =  Dict("p" => p, "q" => q)
    for st_var in ["p", "q"]
        var = vars_dict[st_var]
        for i in 1:length(result_pm["solution"]["branch"])
            @test isapprox(Array(solution(result, var))[i], result_pm["solution"]["branch"][string(i)][string(st_var, "f")], atol = result.options.tol*100)
        end
    end
end

function test_static_case(result, result_pm, result_nlp_pm, pg, qg)
    @test result.status == result_nlp_pm.status
    @test isapprox(result.objective, result_pm["objective"], rtol = result.options.tol*100)
    for i in 1:length(result_pm["solution"]["gen"])
        @test isapprox(Array(solution(result, pg))[i], result_pm["solution"]["gen"][string(i)]["pg"], atol = result.options.tol*1000)
        @test isapprox(Array(solution(result, qg))[i], result_pm["solution"]["gen"][string(i)]["qg"], atol = result.options.tol*1000)
    end
end

function test_polar_voltage(result, result_pm, va, vm)
    @info "ppolar"
    @info result_pm
    for i in 1:length(result_pm["solution"]["bus"])
        @test isapprox(Array(solution(result, va))[i], result_pm["solution"]["bus"][string(i)]["va"], atol = result.options.tol*100)
        @test isapprox(Array(solution(result, vm))[i], result_pm["solution"]["bus"][string(i)]["vm"], rtol = result.options.tol*100)
    end
end

function test_rect_voltage(result, result_pm, vr, vim)
    @info "rrect"
    @info result_pm
    for i in 1:length(result_pm["solution"]["bus"])
        @test isapprox(Array(solution(result, vr))[i], result_pm["solution"]["bus"][string(i)]["vr"], rtol = result.options.tol*100)
        @test isapprox(Array(solution(result, vim))[i], result_pm["solution"]["bus"][string(i)]["vi"], atol = result.options.tol*100)
    end
end

function test_case5(result, result_pm, result_nlp_pm, pg, qg, p, q)
    test_static_case(result, result_pm, result_nlp_pm, pg, qg)
end

function test_case14(result, result_pm, result_nlp_pm, pg, qg, p, q)
    test_static_case(result, result_pm, result_nlp_pm, pg, qg)
end

function test_float32(m, m64, result, backend)
    x1 = result.solution
    tol = 2.71828^(log(result.options.tol)/2)
    x2 = x1 .* (1 .+ 0.01 .* (2 .* rand(size(x1)) .- 1))
    x3 = x1 .* (1 .+ 0.01 .* (2 .* rand(size(x1)) .- 1))
    for x in [x1, x2, x3]
        @test isapprox(NLPModelsJuMP.obj(m, x), NLPModelsJuMP.obj(m64, x), rtol = tol)
        @test isapprox(NLPModelsJuMP.cons(m, x), NLPModelsJuMP.cons(m64, x), rtol = tol)
        if backend != CUDABackend()
            @test isapprox(NLPModelsJuMP.grad(m, x), NLPModelsJuMP.grad(m64, x), rtol = tol)
            @test isapprox(NLPModelsJuMP.jac(m, x), NLPModelsJuMP.jac(m64, x), rtol = tol)
        end
    end
end

function test_mp_case(result, true_sol)
    @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
    @test isapprox(result.objective, true_sol, rtol = result.options.tol*100)
end

function sc_tests(filename, backend, T)
    uc_filename = filename*"_solution.json"
    filename = filename*".json"
    model, cons, vars, lengths, sc_data_array = ExaModelsPower.goc3_model(filename, uc_filename; backend=backend, T=T)
    test_goc3_parser_boundary(sc_data_array, lengths)
    backend === nothing && test_goc3_reactive_capability_rows(filename, uc_filename)
    result = exasolve(model, backend; max_iter=1, tol=1e-2)
end

function test_goc3_parser_boundary(sc_data, lengths)
    L_J_ln = lengths[2]
    L_J_ac = lengths[3]
    L_N_p = lengths[12]

    @test all(r -> r.j == r.j_ac && r.j_ac == r.j_ln, sc_data.aclbrancharray)
    @test all(r -> r.j == r.j_ac && r.j_ac == r.j_xf + L_J_ln, sc_data.acxbrancharray)
    @test all(r -> r.j == r.j_ac && r.j_ac == r.j_xf + L_J_ln, sc_data.fpdarray)
    @test all(r -> r.j == r.j_ac && r.j_ac == r.j_xf + L_J_ln, sc_data.fwrarray)
    @test all(r -> r.j == r.j_ac && r.j_ac == r.j_xf + L_J_ln, sc_data.vpdarray)
    @test all(r -> r.j == r.j_ac && r.j_ac == r.j_xf + L_J_ln, sc_data.vwrarray)
    @test all(r -> r.j == r.j_ac && r.j_ac == r.j_ln, sc_data.jtk_ln_flattened)
    @test all(r -> r.j == r.j_ac && r.j_ac == r.j_xf + L_J_ln, sc_data.jtk_xf_flattened)
    @test all(r -> r.j == r.j_dc + L_J_ac, sc_data.dclinearray)
    @test all(r -> r.j == r.j_dc + L_J_ac, sc_data.jtk_dc_flattened)

    @test all(r -> r.n == r.n_p, sc_data.preservearray)
    @test all(r -> r.n == r.n_p, sc_data.preservesetarray_pr)
    @test all(r -> r.n == r.n_p, sc_data.preservesetarray_cs)
    @test all(r -> r.n == r.n_q + L_N_p, sc_data.qreservearray)
    @test all(r -> r.n == r.n_q + L_N_p, sc_data.qreservesetarray_pr)
    @test all(r -> r.n == r.n_q + L_N_p, sc_data.qreservesetarray_cs)

    if !isempty(sc_data.aclbrancharray)
        @test propertynames(first(sc_data.aclbrancharray))[1:3] == (:j, :j_ac, :j_ln)
    end
    if !isempty(sc_data.acxbrancharray)
        @test propertynames(first(sc_data.acxbrancharray))[1:3] == (:j, :j_ac, :j_xf)
    end
    if !isempty(sc_data.jtk_xf_flattened)
        @test propertynames(first(sc_data.jtk_xf_flattened))[1:5] == (:flat_jtk_xf, :ctg, :j, :j_ac, :j_xf)
    end
end

function set_q_bound_cap!(dev, beta_ub, beta_lb, q_0_ub, q_0_lb)
    dev["q_bound_cap"] = 1
    dev["q_linear_cap"] = 0
    dev["beta_ub"] = beta_ub
    dev["beta_lb"] = beta_lb
    dev["q_0_ub"] = q_0_ub
    dev["q_0_lb"] = q_0_lb
    return dev
end

function set_q_linear_cap!(dev, beta, q_0)
    dev["q_bound_cap"] = 0
    dev["q_linear_cap"] = 1
    dev["beta"] = beta
    dev["q_0"] = q_0
    return dev
end

function test_goc3_reactive_capability_rows(filename, uc_filename)
    data_json = JSON.parsefile(filename)
    devices = data_json["network"]["simple_dispatchable_device"]
    producers = [dev for dev in devices if dev["device_type"] == "producer"]
    consumers = [dev for dev in devices if dev["device_type"] == "consumer"]

    set_q_bound_cap!(producers[1], 0.31, -0.17, 0.41, -0.29)
    set_q_linear_cap!(producers[2], 0.13, 0.07)
    set_q_bound_cap!(consumers[1], 0.23, -0.19, 0.37, -0.31)
    set_q_linear_cap!(consumers[2], 0.11, 0.05)

    data = PowerIO.parse_goc3_json(data_json)
    uc_data = JSON.parsefile(uc_filename)
    sc_data, lengths, _ = ExaModelsPower.parse_sc_data(data, uc_data, data_json)
    L_T = lengths[11]

    @test length(sc_data.prarray_pqbounds) == L_T
    @test length(sc_data.prarray_pqe) == L_T
    @test length(sc_data.csarray_pqbounds) == L_T
    @test length(sc_data.csarray_pqe) == L_T
    @test isconcretetype(eltype(sc_data.prarray_pqbounds))
    @test isconcretetype(eltype(sc_data.prarray_pqe))
    @test isconcretetype(eltype(sc_data.csarray_pqbounds))
    @test isconcretetype(eltype(sc_data.csarray_pqe))
    @test propertynames(first(sc_data.prarray_pqbounds))[1:3] == (:j, :jprcs, :j_pr)
    @test propertynames(first(sc_data.prarray_pqe))[1:3] == (:j, :jprcs, :j_pr)
    @test propertynames(first(sc_data.csarray_pqbounds))[1:3] == (:j, :jprcs, :j_cs)
    @test propertynames(first(sc_data.csarray_pqe))[1:3] == (:j, :jprcs, :j_cs)
end

function test_dcopf_case(result, result_pm, pg, pf)
    @test result.status == MadNLP.SOLVE_SUCCEEDED || result.status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
    @test isapprox(result.objective, result_pm["objective"], rtol = result.options.tol*100)
    for i in 1:length(result_pm["solution"]["gen"])
        @test isapprox(Array(solution(result, pg))[i], result_pm["solution"]["gen"][string(i)]["pg"], atol = result.options.tol*1000)
    end
    for i in 1:length(result_pm["solution"]["branch"])
        @test isapprox(Array(solution(result, pf))[i], result_pm["solution"]["branch"][string(i)]["pf"], atol = result.options.tol*1000)
    end
end
