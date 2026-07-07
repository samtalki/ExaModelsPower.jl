# PowerIO (src/goc3.jl) returns the general GOC3 topology and time-series rows,
# keyed by uid and per-class index. The stacked global variable numbering used
# below (j, j_pr, j_cs, j_prcs, j_sh: offsets into one variable vector) is
# specific to this optimization model, so it is defined and threaded on here.
is_pr(uid::Int, L_J_pr::Int, L_J_cs::Int, producers_first::Bool)::Bool =
    producers_first ? uid < L_J_pr : uid >= L_J_cs
function is_pr(uid_str::String, L_J_pr::Int, L_J_cs::Int, producers_first::Bool)::Bool
    is_pr(parse(Int, match(r"\d+", uid_str).match), L_J_pr, L_J_cs, producers_first)
end

get_j_prcs(uid_str::String, L_J_pr::Int, L_J_cs::Int, producers_first::Bool) =
    parse(Int, match(r"\d+", uid_str).match) + 1

function get_j_pr(uid_str::String, L_J_pr::Int, L_J_cs::Int, producers_first::Bool)
    offset = Int64(producers_first ? 0 : (-L_J_cs))
    parse(Int, match(r"\d+", uid_str).match) + 1 + offset
end

function get_j_cs(uid_str::String, L_J_pr::Int, L_J_cs::Int, producers_first::Bool)
    offset = Int64(producers_first ? (-L_J_pr) : 0)
    parse(Int, match(r"\d+", uid_str).match) + 1 + offset
end

_uidnum(s) = parse(Int, match(r"\d+", s).match)

const goc3_bus_id = PowerIO.goc3_bus_id
const get_as = PowerIO.goc3_interval_bounds

# Shutdown power capability p_sdpc[j_prcs, t, t_prime]. Model-specific: it is
# indexed in this model's stacked producer/consumer variable space.
function goc3_shutdown_power_cap(data, lengths, producers_first)
    L_J_cs = lengths.L_J_cs
    L_J_pr = lengths.L_J_pr
    L_J_cspr = lengths.L_J_cspr
    L_T = lengths.L_T
    periods = data.periods
    dt = Float64.(data.dt)
    # Interval end/start times precomputed once from the cumulative durations,
    # replacing the per-cell get_as(dt, ·) slice-and-sum.
    cdt = cumsum(dt)
    a_end = cdt
    a_start = cdt .- dt

    p_sdpc = zeros(L_J_cspr, L_T, L_T)
    for (key, val) in data.sdd_lookup
        jprcs = get_j_prcs(val["uid"], L_J_pr, L_J_cs, producers_first)
        p0 = val["initial_status"]["p"]
        ramp = val["p_shutdown_ramp_ub"]
        p_lb = data.sdd_ts_lookup[key]["p_lb"]
        for t in periods
            for t_prime in periods
                if t_prime == 1 && t >= t_prime
                    p_sdpc[jprcs, t, t_prime] = p0 - ramp*(a_end[t] - a_start[t_prime])
                elseif t >= t_prime
                    p_sdpc[jprcs, t, t_prime] = p_lb[t_prime-1] - ramp*(a_end[t] - a_start[t_prime])
                end
            end
        end
    end
    return p_sdpc
end

function parse_sc_data(data, uc_data, data_json)
    producers_first = data.sdd_lookup[minimum(keys(data.sdd_lookup))]["device_type"] == "producer"
    sc_data, lengths, cost_vector_pr, cost_vector_cs = PowerIO.goc3_static_data(data)

    (L_J_xf, L_J_ln, L_J_ac, L_J_dc, L_J_br, L_J_cs,
    L_J_pr, L_J_cspr, L_J_sh, I, L_T, L_N_p, L_N_q) = lengths

    # Thread this model's stacked global variable indices onto PowerIO's general
    # rows, reproducing the numbering the model layer expects.
    _j_pr(uid) = (j = _uidnum(uid) + L_J_br + 1, j_prcs = get_j_prcs(uid, L_J_pr, L_J_cs, producers_first), j_pr = get_j_pr(uid, L_J_pr, L_J_cs, producers_first))
    _j_cs(uid) = (j = _uidnum(uid) + L_J_br + 1, j_prcs = get_j_prcs(uid, L_J_pr, L_J_cs, producers_first), j_cs = get_j_cs(uid, L_J_pr, L_J_cs, producers_first))
    sc_data = (
        bus = sc_data.bus,
        shunt = [(j = _uidnum(s.uid) + L_J_br + L_J_cspr + 1, j_sh = _uidnum(s.uid) + 1, s...) for s in sc_data.shunt],
        acl_branch = [(j = b.j_ln, j_ac = b.j_ln, b...) for b in sc_data.acl_branch],
        acx_branch = [(j = b.j_xf + L_J_ln, j_ac = b.j_xf + L_J_ln, b...) for b in sc_data.acx_branch],
        vpd = [(j = b.j_xf + L_J_ln, j_ac = b.j_xf + L_J_ln, b...) for b in sc_data.vpd],
        fpd = [(j = b.j_xf + L_J_ln, j_ac = b.j_xf + L_J_ln, b...) for b in sc_data.fpd],
        vwr = [(j = b.j_xf + L_J_ln, j_ac = b.j_xf + L_J_ln, b...) for b in sc_data.vwr],
        fwr = [(j = b.j_xf + L_J_ln, j_ac = b.j_xf + L_J_ln, b...) for b in sc_data.fwr],
        dc_branch = [(j = b.j_dc + L_J_ac, b...) for b in sc_data.dc_branch],
        prod = [(j = _uidnum(p.uid) + L_J_br + 1, j_pr = get_j_pr(p.uid, L_J_pr, L_J_cs, producers_first), j_prcs = get_j_prcs(p.uid, L_J_pr, L_J_cs, producers_first), p...) for p in sc_data.prod],
        cons = [(j = _uidnum(p.uid) + L_J_br + 1, j_cs = get_j_cs(p.uid, L_J_pr, L_J_cs, producers_first), j_prcs = get_j_prcs(p.uid, L_J_pr, L_J_cs, producers_first), p...) for p in sc_data.cons],
        active_reserve = [(n = r.n_p, r...) for r in sc_data.active_reserve],
        reactive_reserve = [(n = r.n_q + L_N_p, r...) for r in sc_data.reactive_reserve],
        active_reserve_set_pr = [(i = r.i, j = _uidnum(r.uid) + L_J_br + 1, n = r.n, n_p = r.n_p, j_pr = get_j_pr(r.uid, L_J_pr, L_J_cs, producers_first), j_prcs = get_j_prcs(r.uid, L_J_pr, L_J_cs, producers_first)) for r in sc_data.active_reserve_set_pr],
        active_reserve_set_cs = [(i = r.i, j = _uidnum(r.uid) + L_J_br + 1, n = r.n, n_p = r.n_p, j_cs = get_j_cs(r.uid, L_J_pr, L_J_cs, producers_first), j_prcs = get_j_prcs(r.uid, L_J_pr, L_J_cs, producers_first)) for r in sc_data.active_reserve_set_cs],
        reactive_reserve_set_pr = [(i = r.i, j = _uidnum(r.uid) + L_J_br + 1, n = r.n, n_q = r.n_q, j_pr = get_j_pr(r.uid, L_J_pr, L_J_cs, producers_first), j_prcs = get_j_prcs(r.uid, L_J_pr, L_J_cs, producers_first)) for r in sc_data.reactive_reserve_set_pr],
        reactive_reserve_set_cs = [(i = r.i, j = _uidnum(r.uid) + L_J_br + 1, n = r.n, n_q = r.n_q, j_cs = get_j_cs(r.uid, L_J_pr, L_J_cs, producers_first), j_prcs = get_j_prcs(r.uid, L_J_pr, L_J_cs, producers_first)) for r in sc_data.reactive_reserve_set_cs],
    )

    periods = data.periods
    dt = Float64.(data.dt)
    K = length(data_json["reliability"]["contingency"])
    PowerIO.goc3_add_status_flags!(uc_data["time_series_output"]["ac_line"], data.ac_line_lookup)
    PowerIO.goc3_add_status_flags!(uc_data["time_series_output"]["two_winding_transformer"], data.twt_lookup)
    PowerIO.goc3_add_status_flags!(uc_data["time_series_output"]["simple_dispatchable_device"], data.sdd_lookup)

    T_supc_pr = [
            (j = parse(Int, match(r"\d+", val["uid"]).match) + L_J_br + 1,
            j_pr=get_j_pr(val["uid"], L_J_pr, L_J_cs, producers_first),
            j_prcs=get_j_prcs(val["uid"], L_J_pr, L_J_cs, producers_first),
            t = t,
            t_prime = t_prime, 
            p_supc = data.sdd_ts_lookup[key]["p_lb"][t_prime] - val["p_startup_ramp_ub"]*(get_as(dt, t_prime)[3] - get_as(dt, t)[3]),
            u_su = uc["su_status"][t_prime]
        )
        for (key, val) in data.sdd_lookup
        if is_pr(val["uid"], L_J_pr, L_J_cs, producers_first)
        for t in periods
        for t_prime in periods
        if t_prime > t && data.sdd_ts_lookup[key]["p_lb"][t_prime] - val["p_startup_ramp_ub"]*(get_as(dt, t_prime)[3] - get_as(dt, t)[3]) > 0
        for uc in uc_data["time_series_output"]["simple_dispatchable_device"]
        if val["uid"] == uc["uid"]
        ]
    
    #This sum corresponds to constraint 69 (summing p_supc*u_su)
    sum_T_supc_pr = zeros(L_J_pr, length(periods))
    #This sum corresponds to constraint 112/113 (summing u_su)
    sum2_T_supc_pr = zeros(L_J_pr, length(periods))

    for b in T_supc_pr
        sum_T_supc_pr[b.j_pr, b.t] += b.p_supc*b.u_su
        sum2_T_supc_pr[b.j_pr, b.t] += b.u_su
    end

    T_supc_cs = [
            (j = parse(Int, match(r"\d+", val["uid"]).match) + L_J_br + 1,
            j_cs=get_j_cs(val["uid"], L_J_pr, L_J_cs, producers_first),
            j_prcs=get_j_prcs(val["uid"], L_J_pr, L_J_cs, producers_first),
            t = t,
            t_prime = t_prime, 
            p_supc = data.sdd_ts_lookup[key]["p_lb"][t_prime] - val["p_startup_ramp_ub"]*(get_as(dt, t_prime)[3] - get_as(dt, t)[3]),
            u_su = uc["su_status"][t_prime]
        )
        for (key, val) in data.sdd_lookup
        if !is_pr(val["uid"], L_J_pr, L_J_cs, producers_first)
        for t in periods
        for t_prime in periods
        if t_prime > t && data.sdd_ts_lookup[key]["p_lb"][t_prime] - val["p_startup_ramp_ub"]*(get_as(dt, t_prime)[3] - get_as(dt, t)[3]) > 0
        for uc in uc_data["time_series_output"]["simple_dispatchable_device"]
        if val["uid"] == uc["uid"]
        ]
    #This sum corresponds to constraint 69 (p_supc*u_su)
    sum_T_supc_cs = zeros(L_J_cs, L_T)
    #This sum corresponds to constraints 122-126 (u_su)
    sum2_T_supc_cs = zeros(L_J_cs, L_T)
    for b in T_supc_cs
        sum_T_supc_cs[b.j_cs, b.t] += b.p_supc*b.u_su
        sum2_T_supc_cs[b.j_cs, b.t] += b.u_su
    end


    # p_sdpc (shutdown power capability), indexed [j_prcs, t, t_prime]
    p_sdpc = goc3_shutdown_power_cap(data, lengths, producers_first)

    T_sdpc_pr = [
            (j = parse(Int, match(r"\d+", val["uid"]).match) + L_J_br + 1,
            j_pr=get_j_pr(val["uid"], L_J_pr, L_J_cs, producers_first),
            j_prcs=get_j_prcs(val["uid"], L_J_pr, L_J_cs, producers_first),
            t = t,
            t_prime = t_prime,
            p_sdpc = p_sdpc[parse(Int, match(r"\d+", val["uid"]).match)+1, t, t_prime],
            u_sd = uc["sd_status"][t_prime]
            )
        for (key, val) in data.sdd_lookup
        if is_pr(val["uid"], L_J_pr, L_J_cs, producers_first)
        for t in periods
        for t_prime in periods
        if t_prime <= t && p_sdpc[parse(Int, match(r"\d+", val["uid"]).match)+1, t, t_prime] > 0
        for uc in uc_data["time_series_output"]["simple_dispatchable_device"]
        if val["uid"] == uc["uid"]
        ]

    #This sum corresponds to constraint 70 (summing p_sdpc*u_sd)
    sum_T_sdpc_pr = zeros(L_J_pr, L_T)
    #This sum corresponds to constraints 112 and 113 (summing u_sd)
    sum2_T_sdpc_pr = zeros(L_J_pr, L_T)
    for b in T_sdpc_pr
        sum_T_sdpc_pr[b.j_pr, b.t] += b.p_sdpc*b.u_sd
        sum2_T_sdpc_pr[b.j_pr, b.t] += b.u_sd
    end

    T_sdpc_cs = [
            (j = parse(Int, match(r"\d+", val["uid"]).match) + L_J_br + 1,
            j_cs=get_j_cs(val["uid"], L_J_pr, L_J_cs, producers_first),
            j_prcs=get_j_prcs(val["uid"], L_J_pr, L_J_cs, producers_first),
            t = t,
            t_prime = t_prime,
            p_sdpc = p_sdpc[parse(Int, match(r"\d+", val["uid"]).match)+1, t, t_prime],
            u_sd = uc["sd_status"][t_prime]
            )
        for (key, val) in data.sdd_lookup
        if !is_pr(val["uid"], L_J_pr, L_J_cs, producers_first)
        for t in periods
        for t_prime in periods
        if t_prime <= t && p_sdpc[parse(Int, match(r"\d+", val["uid"]).match)+1, t, t_prime] > 0
        for uc in uc_data["time_series_output"]["simple_dispatchable_device"]
        if val["uid"] == uc["uid"]
        ]
    #This sum corresponds to constraint 70 (p_sdpc*u_sd)
    sum_T_sdpc_cs = zeros(L_J_cs, L_T)
    #This sum corresponds to constraint 122-126 (u_sd)
    sum2_T_sdpc_cs = zeros(L_J_cs, L_T)
    for b in T_sdpc_cs
        sum_T_sdpc_cs[b.j_cs, b.t] += b.p_sdpc*b.u_sd
        sum2_T_sdpc_cs[b.j_cs, b.t] += b.u_sd
    end
    

    # Multi-interval energy requirement windows and their per-period membership.
    ew = PowerIO.goc3_energy_windows(data)
    W_en_max_pr = [(w_en_max_pr_ind = r.w_en_max_pr_ind, _j_pr(r.uid)..., a_en_max_start = r.a_en_max_start, a_en_max_end = r.a_en_max_end, e_max = r.e_max) for r in ew.W_en_max_pr]
    W_en_max_cs = [(w_en_max_cs_ind = r.w_en_max_cs_ind, _j_cs(r.uid)..., a_en_max_start = r.a_en_max_start, a_en_max_end = r.a_en_max_end, e_max = r.e_max) for r in ew.W_en_max_cs]
    W_en_min_pr = [(w_en_min_pr_ind = r.w_en_min_pr_ind, _j_pr(r.uid)..., a_en_min_start = r.a_en_min_start, a_en_min_end = r.a_en_min_end, e_min = r.e_min) for r in ew.W_en_min_pr]
    W_en_min_cs = [(w_en_min_cs_ind = r.w_en_min_cs_ind, _j_cs(r.uid)..., a_en_min_start = r.a_en_min_start, a_en_min_end = r.a_en_min_end, e_min = r.e_min) for r in ew.W_en_min_cs]
    T_w_en_max_pr = [(w_en_max_pr_ind = r.w_en_max_pr_ind, _j_pr(r.uid)..., t = r.t, dt = r.dt) for r in ew.T_w_en_max_pr]
    T_w_en_max_cs = [(w_en_max_cs_ind = r.w_en_max_cs_ind, _j_cs(r.uid)..., t = r.t, dt = r.dt) for r in ew.T_w_en_max_cs]
    T_w_en_min_pr = [(w_en_min_pr_ind = r.w_en_min_pr_ind, _j_pr(r.uid)..., t = r.t, dt = r.dt) for r in ew.T_w_en_min_pr]
    T_w_en_min_cs = [(w_en_min_cs_ind = r.w_en_min_cs_ind, _j_cs(r.uid)..., t = r.t, dt = r.dt) for r in ew.T_w_en_min_cs]
    L_W_en_max_pr = length(W_en_max_pr)
    L_W_en_max_cs = length(W_en_max_cs)
    L_W_en_min_pr = length(W_en_min_pr)
    L_W_en_min_cs = length(W_en_min_cs)

    pb_pr, pb_cs = PowerIO.goc3_price_blocks(cost_vector_pr, cost_vector_cs)
    p_jtm_flattened_pr = [(flat_k = r.flat_k, _j_pr(r.uid)..., t = r.t, m = r.m, c_en = r.c_en, p_max = r.p_max) for r in pb_pr]
    p_jtm_flattened_cs = [(flat_k = r.flat_k, _j_cs(r.uid)..., t = r.t, m = r.m, c_en = r.c_en, p_max = r.p_max) for r in pb_cs]

    # Post-contingency surviving AC branches. PowerIO enumerates the pure survivor
    # rows per contingency; here we attach the UC on_status and expand over
    # periods, preserving the (ctg, t, component, branch) flat ordering.
    ctg_ac_survivors = PowerIO.goc3_ac_contingency_survivors(data, lengths)
    uc_ln_lookup = Dict(uc["uid"] => uc for uc in uc_data["time_series_output"]["ac_line"])
    uc_xf_lookup = Dict(uc["uid"] => uc for uc in uc_data["time_series_output"]["two_winding_transformer"])

    jtk_ln_flattened = Vector{@NamedTuple{flat_jtk_ln::Int, ctg::Int, j::Int, j_ac::Int, j_ln::Int, to_bus::Int, fr_bus::Int, b_sr::Float64, s_max_ctg::Float64, u_on::Int, t::Int, dt::Float64}}()
    flat_jtk_ln = 1
    for rows in ctg_ac_survivors.ln
        for t in periods
            for r in rows
                haskey(uc_ln_lookup, r.uid) || continue
                uc = uc_ln_lookup[r.uid]
                push!(jtk_ln_flattened, (flat_jtk_ln=flat_jtk_ln, ctg=r.ctg, j=r.j_ln, j_ac=r.j_ln, j_ln=r.j_ln,
                to_bus=r.to_bus, fr_bus=r.fr_bus, b_sr=r.b_sr, s_max_ctg=r.s_max_ctg, u_on=uc["on_status"][t], t=t, dt=dt[t]))
                flat_jtk_ln += 1
            end
        end
    end

    jtk_xf_flattened = Vector{@NamedTuple{flat_jtk_xf::Int, ctg::Int, j::Int, j_ac::Int, j_xf::Int, to_bus::Int, fr_bus::Int, b_sr::Float64, s_max_ctg::Float64, u_on::Int, t::Int, dt::Float64}}()
    flat_jtk_xf = 1
    for rows in ctg_ac_survivors.xf
        for t in periods
            for r in rows
                haskey(uc_xf_lookup, r.uid) || continue
                uc = uc_xf_lookup[r.uid]
                push!(jtk_xf_flattened, (flat_jtk_xf=flat_jtk_xf, ctg=r.ctg, j=r.j_xf + L_J_ln, j_ac=r.j_xf + L_J_ln, j_xf=r.j_xf,
                to_bus=r.to_bus, fr_bus=r.fr_bus, b_sr=r.b_sr, s_max_ctg=r.s_max_ctg, u_on=uc["on_status"][t], t=t, dt=dt[t]))
                flat_jtk_xf += 1
            end
        end
    end

    jtk_dc_flattened = [(flat_jtk_dc = r.flat_jtk_dc, ctg = r.ctg, j = r.j_dc + L_J_ac, j_dc = r.j_dc, to_bus = r.to_bus, fr_bus = r.fr_bus, t = r.t, dt = r.dt) for r in PowerIO.goc3_dc_contingency_flows(data)]

    empty_vpd = Vector{NamedTuple{(:j, :j_ac, :j_xf, :phi_min, :phi_max, :t), Tuple{Int64, Int64, Int64, Float64, Float64, Int64}}}()
    empty_fpd = Vector{NamedTuple{(:j, :j_ac, :j_xf, :phi_o, :t), Tuple{Int64, Int64, Int64, Float64, Int64}}}()
    empty_vwr = Vector{NamedTuple{(:j, :j_ac, :j_xf, :tau_min, :tau_max, :t), Tuple{Int64, Int64, Int64, Float64, Float64, Int64}}}()
    empty_fwr = Vector{NamedTuple{(:j, :j_ac, :j_xf, :tau_o, :t), Tuple{Int64, Int64, Int64, Float64, Int64}}}()

    sc_time_data = (
        ;
        periods = periods,

        c_p = [data.violation_cost["p_bus_vio_cost"]],
        c_q = [data.violation_cost["q_bus_vio_cost"]],
        c_s = [data.violation_cost["s_vio_cost"]],
        c_e = [data.violation_cost["e_vio_cost"]],
        busarray = [(;b..., t=t, dt=dt[t]) for b in sc_data.bus, t in periods],
        k_busarray = [(;b..., t=t, k=k) for b in sc_data.bus, t in periods, k in 1:length(data_json["reliability"]["contingency"])],
        shuntarray = [
            (;s..., t=t, u_sh = uc["step"][t])
            for s in sc_data.shunt, t in periods
            for uc in uc_data["time_series_output"]["shunt"]
            if s.uid == uc["uid"]
                ],
        k_shuntarray = [(;b..., t=t, k=k) for b in sc_data.shunt, t in periods, k in 1:length(data_json["reliability"]["contingency"])],
       

        preservearray = [(;n=r.n, n_p=r.n_p, uid=r.uid, c_rgu=r.c_rgu, c_rgd=r.c_rgd, c_scr=r.c_scr, c_nsc=r.c_nsc, c_rru=r.c_rru, c_rrd=r.c_rrd,
        σ_rgu=r.σ_rgu, σ_rgd=r.σ_rgd, σ_scr=r.σ_scr, σ_nsc=r.σ_nsc, p_rru_min=r.p_rru_min[t], p_rrd_min=r.p_rrd_min[t],
        t=t, dt=dt[t]) for r in sc_data.active_reserve, t in periods],

        qreservearray = [(;n=q.n, n_q=q.n_q, uid=q.uid, c_qru=q.c_qru, c_qrd=q.c_qrd, q_qru_min=q.q_qru_min[t], q_qrd_min=q.q_qrd_min[t], t=t, dt = dt[t])
        for q in sc_data.reactive_reserve, t in periods],

        k_prarray = [(;j_prcs=b.j_prcs, j_pr=b.j_pr, bus = b.bus, uid=b.uid, t=t, k=k) for b in sc_data.prod, t in periods, k in 1:length(data_json["reliability"]["contingency"])],
        prarray = [(;j=p.j, j_prcs=p.j_prcs, j_pr=p.j_pr, bus=p.bus, uid=p.uid, c_on=p.c_on, c_sd=p.c_sd, c_su = p.c_su, p_ru=p.p_ru, p_rd=p.p_rd,  
        p_ru_su=p.p_ru_su, p_rd_sd=p.p_rd_sd, c_rgu=p.c_rgu[t], c_rgd=p.c_rgd[t], c_scr=p.c_scr[t], c_nsc=p.c_nsc[t], c_rru_on=p.c_rru_on[t],
        c_rru_off=p.c_rru_off[t], c_rrd_on=p.c_rrd_on[t], c_rrd_off=p.c_rrd_off[t], c_qru=p.c_qru[t], c_qrd=p.c_qrd[t], p_rgu_max=p.p_rgu_max,
        p_rgd_max=p.p_rgd_max, p_scr_max=p.p_scr_max, p_nsc_max=p.p_nsc_max, p_rru_on_max=p.p_rru_on_max, p_rru_off_max=p.p_rru_off_max, 
        p_rrd_on_max=p.p_rrd_on_max, p_rrd_off_max=p.p_rrd_off_max, p_0=p.p_0, q_0=p.q_0, p_max=p.p_max[t], p_min=p.p_min[t], q_max=p.q_max[t], q_min=p.q_min[t], #sus = p.sus,
        u_on = uc["on_status"][t], u_su = uc["su_status"][t], u_sd = uc["sd_status"][t], t=t,
        sum_T_supc_pr_jt = sum_T_supc_pr[p.j_pr, t], sum_T_sdpc_pr_jt = sum_T_sdpc_pr[p.j_pr, t], sum2_T_supc_pr_jt=sum2_T_supc_pr[p.j_pr, t], sum2_T_sdpc_pr_jt=sum2_T_sdpc_pr[p.j_pr, t], dt = dt[t])
        for p in sc_data.prod, t in periods
        for uc in uc_data["time_series_output"]["simple_dispatchable_device"]
        if p.uid == uc["uid"]],

        prarray_pqbounds = isempty(val for val in values(data.sdd_lookup) if is_pr(val["uid"], L_J_pr, L_J_cs, producers_first) && val["q_bound_cap"]==1) ? 
                            empty_data = Vector{NamedTuple{(:j, :jprcs, :j_pr, :u_on, :sum2_T_supc_pr_jt, :sum2_T_sdpc_pr_jt, :beta_max, :beta_min, :q_max_p0, :q_min_p0, :t), Tuple{Int64, Int64, Int64, Int64, Int64, Int64, Float64, Float64, Float64, Float64, Int64}}}() : [
                            (;j = parse(Int, match(r"\d+", val["uid"]).match) + L_J_br + 1,
                            j_prcs=get_j_prcs(val["uid"], L_J_pr, L_J_cs, producers_first),
                            j_pr=get_j_pr(val["uid"], L_J_pr, L_J_cs, producers_first),
                            u_on = uc["on_status"][t],
                            sum2_T_supc_pr_jt=sum2_T_supc_pr[p.j_pr, t], 
                            sum2_T_sdpc_pr_jt=sum2_T_sdpc_pr[p.j_pr, t], 
                            beta_max = val["beta_ub"],
                            beta_min = val["beta_lb"],
                            q_max_p0 = val["q_0_ub"],
                            q_min_p0 = val["q_0_lb"],
                            t=t
        )
        for val in values(data.sdd_lookup), t in periods
        if is_pr(val["uid"], L_J_pr, L_J_cs, producers_first) && val["q_bound_cap"]==1
        for uc in uc_data["time_series_output"]["simple_dispatchable_device"]
        if p.uid == uc["uid"]],

        prarray_pqe = isempty(val for val in values(data.sdd_lookup) if is_pr(val["uid"], L_J_pr, L_J_cs, producers_first) && val["q_bound_cap"]==1) ? 
                            empty_data = Vector{NamedTuple{(:j, :jprcs, :j_pr, :u_on, :sum2_T_supc_pr_jt, :sum2_T_sdpc_pr_jt, :beta, :q_p0, :t), Tuple{Int64, Int64, Int64, Int64, Int64, Int64, Float64, Float64, Int64}}}() : [
                            (;j = parse(Int, match(r"\d+", val["uid"]).match) + L_J_br + 1,
                            j_prcs=get_j_prcs(val["uid"], L_J_pr, L_J_cs, producers_first),
                            j_pr=get_j_pr(val["uid"], L_J_pr, L_J_cs, producers_first),
                            u_on = uc["on_status"][t],
                            sum2_T_supc_pr_jt=sum2_T_supc_pr[p.j_pr, t], 
                            sum2_T_sdpc_pr_jt=sum2_T_sdpc_pr[p.j_pr, t], 
                            beta = val["beta"],
                            q_p0 = val["q_0"],
                            t=t
        )
        for val in values(data.sdd_lookup), t in periods
        if is_pr(val["uid"], L_J_pr, L_J_cs, producers_first) && val["q_linear_cap"]==1
        for uc in uc_data["time_series_output"]["simple_dispatchable_device"]
        if p.uid == uc["uid"]],

        k_csarray = [(;j_prcs=b.j_prcs, j_cs=b.j_cs, bus = b.bus, uid=b.uid, t=t, k=k) for b in sc_data.cons, t in periods, k in 1:length(data_json["reliability"]["contingency"])],

        csarray = [(;j=p.j, j_prcs=p.j_prcs, j_cs=p.j_cs, bus=p.bus, uid=p.uid, c_on=p.c_on, c_sd=p.c_sd, c_su=p.c_su, p_ru=p.p_ru, p_rd=p.p_rd,  
        p_ru_su=p.p_ru_su, p_rd_sd=p.p_rd_sd, c_rgu=p.c_rgu[t], c_rgd=p.c_rgd[t], c_scr=p.c_scr[t], c_nsc=p.c_nsc[t], c_rru_on=p.c_rru_on[t],
        c_rru_off=p.c_rru_off[t], c_rrd_on=p.c_rrd_on[t], c_rrd_off=p.c_rrd_off[t], c_qru=p.c_qru[t], c_qrd=p.c_qrd[t], p_rgu_max=p.p_rgu_max,
        p_rgd_max=p.p_rgd_max, p_scr_max=p.p_scr_max, p_nsc_max=p.p_nsc_max, p_rru_on_max=p.p_rru_on_max, p_rru_off_max=p.p_rru_off_max, 
        p_rrd_on_max=p.p_rrd_on_max, p_rrd_off_max=p.p_rrd_off_max, p_0=p.p_0, q_0=p.q_0, p_max=p.p_max[t], p_min=p.p_min[t], q_max=p.q_max[t], q_min=p.q_min[t], #sus=p.sus, 
        u_on = uc["on_status"][t], u_su = uc["su_status"][t], u_sd = uc["sd_status"][t], t=t,
        sum_T_supc_cs_jt = sum_T_supc_cs[p.j_cs, t], sum_T_sdpc_cs_jt = sum_T_sdpc_cs[p.j_cs, t], sum2_T_supc_cs_jt = sum2_T_supc_cs[p.j_cs, t], sum2_T_sdpc_cs_jt = sum2_T_sdpc_cs[p.j_cs, t], dt = dt[t])
        for p in sc_data.cons, t in periods
        for uc in uc_data["time_series_output"]["simple_dispatchable_device"]
        if p.uid == uc["uid"]],

        csarray_pqbounds = isempty(val for val in values(data.sdd_lookup) if !is_pr(val["uid"], L_J_pr, L_J_cs, producers_first) && val["q_bound_cap"]==1) ? 
                            empty_data = Vector{NamedTuple{(:j, :jprcs, :j_cs, :u_on, :sum2_T_supc_cs_jt, :sum2_T_sdpc_cs_jt, :beta_max, :beta_min, :q_max_p0, :q_min_p0, :t), Tuple{Int64, Int64, Int64, Int64, Int64, Int64, Float64, Float64, Float64, Float64, Int64}}}() : [
                            (;j = parse(Int, match(r"\d+", val["uid"]).match) + L_J_br + 1,
                            j_cs=get_j_cs(val["uid"], L_J_pr, L_J_cs, producers_first),
                            j_prcs=get_j_prcs(val["uid"], L_J_pr, L_J_cs, producers_first),
                            u_on = uc["on_status"][t],
                            sum2_T_supc_cs_jt=sum2_T_supc_cs[p.j_cs, t], 
                            sum2_T_sdpc_cs_jt=sum2_T_sdpc_cs[p.j_cs, t], 
                            beta_max = val["beta_ub"],
                            beta_min = val["beta_lb"],
                            q_max_p0 = val["q_0_ub"],
                            q_min_p0 = val["q_0_lb"],
                            t=t
        )
        for val in values(data.sdd_lookup), t in periods
        if !is_pr(val["uid"], L_J_pr, L_J_cs, producers_first) && val["q_bound_cap"]==1
        for uc in uc_data["time_series_output"]["simple_dispatchable_device"]
        if p.uid == uc["uid"]],

        csarray_pqe = isempty(val for val in values(data.sdd_lookup) if !is_pr(val["uid"], L_J_pr, L_J_cs, producers_first) && val["q_bound_cap"]==1) ? 
                            empty_data = Vector{NamedTuple{(:j, :jprcs, :j_cs, :u_on, :sum2_T_supc_cs_jt, :sum2_T_sdpc_cs_jt, :beta, :q_p0, :t), Tuple{Int64, Int64, Int64, Int64, Int64, Int64, Float64, Float64, Int64}}}() : [
                            (;j = parse(Int, match(r"\d+", val["uid"]).match) + L_J_br + 1,
                            j_cs=get_j_cs(val["uid"], L_J_pr, L_J_cs, producers_first),
                            j_prcs=get_j_prcs(val["uid"], L_J_pr, L_J_cs, producers_first),
                            u_on = uc["on_status"][t],
                            sum2_T_supc_cs_jt=sum2_T_supc_cs[p.j_cs, t], 
                            sum2_T_sdpc_cs_jt=sum2_T_sdpc_cs[p.j_cs, t], 
                            beta = val["beta"],
                            q_p0 = val["q_0"],
                            t=t
        )
        for val in values(data.sdd_lookup), t in periods
        if !is_pr(val["uid"], L_J_pr, L_J_cs, producers_first) && val["q_linear_cap"]==1
        for uc in uc_data["time_series_output"]["simple_dispatchable_device"]
        if p.uid == uc["uid"]],

        acxbrancharray = [
            (;j=b.j, j_ac=b.j_ac, j_xf=b.j_xf, uid=b.uid, to_bus=b.to_bus, fr_bus=b.fr_bus, c_su=b.c_su, c_sd=b.c_sd, s_max=b.s_max, g_sr=b.g_sr, b_sr=b.b_sr, b_ch=b.b_ch,
            g_fr=b.g_fr, g_to=b.g_to, b_fr=b.b_fr, b_to=b.b_to, u_on=uc["on_status"][t], u_su=uc["su_status"][t], u_sd=uc["sd_status"][t], t=t, dt = dt[t])
            for b in sc_data.acx_branch, t in periods
            for uc in uc_data["time_series_output"]["two_winding_transformer"]
            if b.uid == uc["uid"]],

        aclbrancharray = [
            (;j=b.j, j_ac=b.j_ac, j_ln=b.j_ln, uid=b.uid, to_bus=b.to_bus, fr_bus=b.fr_bus, c_su=b.c_su, c_sd=b.c_sd, s_max=b.s_max, g_sr=b.g_sr, b_sr=b.b_sr, b_ch=b.b_ch,
            g_fr=b.g_fr, g_to=b.g_to, b_fr=b.b_fr, b_to=b.b_to, u_on=uc["on_status"][t], u_su=uc["su_status"][t], u_sd=uc["sd_status"][t], t=t, dt = dt[t])
            for b in sc_data.acl_branch, t in periods
            for uc in uc_data["time_series_output"]["ac_line"]
            if b.uid == uc["uid"]],

        fpdarray = isempty(sc_data.fpd) ? empty_data = empty_fpd : [(;b..., t=t) for b in sc_data.fpd, t in periods],
        fwrarray = isempty(sc_data.fwr) ? empty_data = empty_fwr : [(;b..., t=t) for b in sc_data.fwr, t in periods],
        vpdarray = isempty(sc_data.vpd) ? empty_data = empty_vpd : [(;b..., t=t) for b in sc_data.vpd, t in periods],
        vwrarray = isempty(sc_data.vwr) ? empty_data = empty_vwr : [(;b..., t=t) for b in sc_data.vwr, t in periods],
        dclinearray = [(;b..., t=t) for b in sc_data.dc_branch, t in periods],

        p_jt_fr_dc_max = [dc.pdc_max for dc in sc_data.dc_branch, t in periods],
        p_jt_to_dc_max = [dc.pdc_max for dc in sc_data.dc_branch, t in periods],
        q_jt_fr_dc_lvar = [dc.qdc_fr_min for dc in sc_data.dc_branch, t in periods],
        q_jt_fr_dc_uvar = [dc.qdc_fr_max for dc in sc_data.dc_branch, t in periods],
        q_jt_to_dc_lvar = [dc.qdc_to_min for dc in sc_data.dc_branch, t in periods],
        q_jt_to_dc_uvar = [dc.qdc_to_max for dc in sc_data.dc_branch, t in periods],

        preservesetarray_pr = [(;b..., t=t) for b in sc_data.active_reserve_set_pr, t in periods],
        preservesetarray_cs = [(;b..., t=t) for b in sc_data.active_reserve_set_cs, t in periods],
        qreservesetarray_pr = [(;b..., t=t) for b in sc_data.reactive_reserve_set_pr, t in periods],
        qreservesetarray_cs = [(;b..., t=t) for b in sc_data.reactive_reserve_set_cs, t in periods],

        v_lvar = repeat([b.v_min for b in sc_data.bus], 1, L_T),
        v_uvar = repeat([b.v_max for b in sc_data.bus], 1, L_T),

        W_en_max_pr=W_en_max_pr,
        W_en_max_cs=W_en_max_cs,
        T_w_en_max_pr=T_w_en_max_pr,
        T_w_en_max_cs=T_w_en_max_cs,
        W_en_min_pr=W_en_min_pr,
        W_en_min_cs=W_en_min_cs,
        T_w_en_min_pr=T_w_en_min_pr,
        T_w_en_min_cs=T_w_en_min_cs,

        tk_index = [(t=t, k=k) for t in periods, k in 1:length(data_json["reliability"]["contingency"])],

        p_jtm_flattened_pr=p_jtm_flattened_pr,
        p_jtm_flattened_cs=p_jtm_flattened_cs,
        p_jtm_pr_uvar = [b.p_max for b in p_jtm_flattened_pr],
        p_jtm_cs_uvar = [b.p_max for b in p_jtm_flattened_cs],
        jtk_ln_flattened=jtk_ln_flattened,
        jtk_xf_flattened=jtk_xf_flattened,
        jtk_dc_flattened=jtk_dc_flattened,

    )

    lengths = (L_J_xf, L_J_ln, L_J_ac, L_J_dc, L_J_br, L_J_cs,
    L_J_pr, L_J_cspr, L_J_sh, I, L_T, L_N_p, L_N_q, L_W_en_min_pr, L_W_en_min_cs, L_W_en_max_pr, L_W_en_max_cs, K)

    return sc_time_data, lengths, producers_first
end

function save_go3_solution(uc_filename, solution_name, result, vars, lengths, producers_first)
    uc_data = JSON.parsefile(uc_filename)
    (L_J_xf, L_J_ln, L_J_ac, L_J_dc, L_J_br, L_J_cs,
    L_J_pr, L_J_cspr, L_J_sh, I, L_T, L_N_p, L_N_q, L_W_en_min_pr, L_W_en_min_cs, L_W_en_max_pr, L_W_en_max_cs, K) = lengths
    #Update simple dispatchable devices
    for line in uc_data["time_series_output"]["simple_dispatchable_device"]
        raw_uid = parse(Int, match(r"\d+", line["uid"]).match)
        solution_index = raw_uid + 1 #This corresponds to j_prcs
        if !is_pr(raw_uid, L_J_pr, L_J_cs, producers_first)
            #This section corresponds to consuming devices
            if producers_first
                solution_index -= L_J_pr
            end
            line["p_syn_res"] = Array(solution(result, vars.p_jt_scr_cs))[solution_index,:]
            line["p_ramp_res_up_online"] = Array(solution(result, vars.p_jt_rru_on_cs))[solution_index,:]
            line["p_nsyn_res"] = zeros(L_T)
            line["p_reg_res_up"] = Array(solution(result, vars.p_jt_rgu_cs))[solution_index,:]
            line["p_ramp_res_down_online"] = Array(solution(result, vars.p_jt_rrd_on_cs))[solution_index,:]
            line["p_on"] = Array(solution(result, vars.p_jt_on_cs))[solution_index,:]
            line["q"] = Array(solution(result, vars.q_jt_cs))[solution_index,:]
            line["p_reg_res_down"] = Array(solution(result, vars.p_jt_rgd_cs))[solution_index,:]
            line["p_ramp_res_up_offline"] = zeros(L_T)
            line["q_res_down"] = Array(solution(result, vars.q_jt_qrd_cs))[solution_index,:]
            line["q_res_up"] = Array(solution(result, vars.q_jt_qru_cs))[solution_index,:]
            line["p_ramp_res_down_offline"] = Array(solution(result, vars.p_jt_rrd_off_cs))[solution_index,:]
        else
            #Producing devices
            if !producers_first
                solution_index -= L_J_cs
            end
            line["p_syn_res"] = Array(solution(result, vars.p_jt_scr_pr))[solution_index,:]
            line["p_ramp_res_up_online"] = Array(solution(result, vars.p_jt_rru_on_pr))[solution_index,:]
            line["p_nsyn_res"] = Array(solution(result, vars.p_jt_nsc_pr))[solution_index,:]
            line["p_reg_res_up"] = Array(solution(result, vars.p_jt_rgu_pr))[solution_index,:]
            line["p_ramp_res_down_online"] = Array(solution(result, vars.p_jt_rrd_on_pr))[solution_index,:]
            line["p_on"] = Array(solution(result, vars.p_jt_on_pr))[solution_index,:]
            line["q"] = Array(solution(result, vars.q_jt_pr))[solution_index,:]
            line["p_reg_res_down"] = Array(solution(result, vars.p_jt_rgd_pr))[solution_index,:]
            line["p_ramp_res_up_offline"] = Array(solution(result, vars.p_jt_rru_off_pr))[solution_index,:]
            line["q_res_down"] = Array(solution(result, vars.q_jt_qrd_pr))[solution_index,:]
            line["q_res_up"] = Array(solution(result, vars.q_jt_qru_pr))[solution_index,:]
            line["p_ramp_res_down_offline"] = zeros(L_T)
        end
    end
    #Update two winding transformers
    for line in uc_data["time_series_output"]["two_winding_transformer"]
        solution_index = parse(Int, match(r"\d+", line["uid"]).match) + 1 #This corresponds to j_xf
        line["tm"] = Array(solution(result, vars.τ_jt_xf))[solution_index,:]
        line["ta"] = Array(solution(result, vars.φ_jt_xf))[solution_index,:] 
    end
    #ac line flows can be inferred from voltage and angle levels

    #Update DC lines
    for line in uc_data["time_series_output"]["dc_line"]
        solution_index = parse(Int, match(r"\d+", line["uid"]).match) + 1 #This corresponds to j_dc
        line["qdc_fr"] = Array(solution(result, vars.q_jt_fr_dc))[solution_index,:]
        line["pdc_fr"] = Array(solution(result, vars.p_jt_fr_dc))[solution_index,:] #pdc_to is implied from energy conservation
        line["qdc_to"] = Array(solution(result, vars.q_jt_to_dc))[solution_index,:]
    end

    #Update buses
    for line in uc_data["time_series_output"]["bus"]
        solution_index = parse(Int, match(r"\d+", line["uid"]).match) + 1 #This corresponds to i
        line["vm"] = Array(solution(result, vars.v_it))[solution_index,:]
        line["va"] = Array(solution(result, vars.θ_it))[solution_index,:]
    end

    open(solution_name, "w") do io
        JSON.print(io, uc_data, 4)
    end
end
