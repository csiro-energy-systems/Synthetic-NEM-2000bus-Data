""
function solve_mn_tnep(file, model_type::Type, optimizer; kwargs...)
    return _PM.solve_model(file, model_type, optimizer, build_mn_tnep; multinetwork=true, ref_extensions=[_PM.ref_add_on_off_va_bounds!,_PM.ref_add_ne_branch!], kwargs...)
end

"the general form of the tnep optimization model"
function build_mn_tnep(pm::_PM.AbstractPowerModel)
    variable_gen_power_real_slack(pm, nw=1)
    variable_gen_power_real_slack_abs(pm, nw=1)
    # variable_gen_power_imaginary_slack(pm, nw=1)
    # variable_gen_power_imaginary_slack_abs(pm, nw=1)
    for i in _PM.ids(pm, :gen, nw=1)
        constraint_gen_slack_abs(pm, i, nw=1)
    end

    _PM.variable_ne_branch_indicator(pm, nw=1)

    for (n, network) in _PM.nws(pm)
        _PM.variable_bus_voltage(pm, nw=n)
        _PM.variable_gen_power(pm, nw=n)
        _PM.variable_branch_power(pm, nw=n)
        _PM.variable_dcline_power(pm, nw=n)

        # _PM.variable_ne_branch_indicator(pm, nw=n)
        _PM.variable_ne_branch_power(pm, nw=n)
        _PM.variable_ne_branch_voltage(pm, nw=n)
        _PM.constraint_model_voltage(pm, nw=n)
        _PM.constraint_ne_model_voltage(pm, nw=n)

        for i in _PM.ids(pm, :ref_buses, nw=n)
            _PM.constraint_theta_ref(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :bus, nw=n)
            # _PM.constraint_ne_power_balance(pm, i, nw=n)
            constraint_ne_power_balance_gen_slack(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :branch, nw=n)
            _PM.constraint_ohms_yt_from(pm, i, nw=n)
            _PM.constraint_ohms_yt_to(pm, i, nw=n)
            _PM.constraint_voltage_angle_difference(pm, i, nw=n)
            _PM.constraint_thermal_limit_from(pm, i, nw=n)
            _PM.constraint_thermal_limit_to(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :ne_branch, nw=n)
            constraint_ne_ohms_yt_from(pm, i, nw=n)
            constraint_ne_ohms_yt_to(pm, i, nw=n)
            constraint_ne_voltage_angle_difference(pm, i, nw=n)
            constraint_ne_thermal_limit_from(pm, i, nw=n)
            constraint_ne_thermal_limit_to(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :dcline, nw=n)
            _PM.constraint_dcline_power_losses(pm, i, nw=n)
        end
    end

    objective_min_tnep_cost(pm)

end


"variable: `0 <= gen_ne[g] <= 1` for `g` in `gen`s"
function variable_ne_gen_indicator(pm::_PM.AbstractPowerModel; nw::Int=nw_id_default, relax::Bool=false, report::Bool=true)
    if !relax
        z_gen_ne = _PM.var(pm, nw)[:gen_ne] = JuMP.@variable(pm.model,
            [i in _PM.ids(pm, nw, :gen)], base_name="$(nw)_gen_ne",
            binary = true,
            start = _PM.comp_start_value(_PM.ref(pm, nw, :gen, i), "gen_tnep_start", 1.0)
        )
    else
        z_gen_ne = _PM.var(pm, nw)[:gen_ne] = JuMP.@variable(pm.model,
            [i in _PM.ids(pm, nw, :gen)], base_name="$(nw)_gen_ne",
            lower_bound = 0,
            upper_bound = 1,
            start = _PM.comp_start_value(_PM.ref(pm, nw, :gen, i), "gen_tnep_start", 1.0)
        )
    end

    report && _PM.sol_component_value(pm, nw, :gen, :gen_status, _PM.ids(pm, nw, :gen), z_gen_ne)
end


function objective_min_tnep_cost(pm::_PM.AbstractPowerModel)
    gen_cost = Dict()
    for (n, nw_ref) in _PM.nws(pm)
        for (i,gen) in nw_ref[:gen]
            pg = _PM.var(pm, n, :pg, i)

            if length(gen["cost"]) == 1
                gen_cost[(n,i)] = gen["cost"][1]
            elseif length(gen["cost"]) == 2
                gen_cost[(n,i)] = gen["cost"][1]*pg + gen["cost"][2]
            elseif length(gen["cost"]) == 3
                gen_cost[(n,i)] = gen["cost"][1]*pg^2 + gen["cost"][2]*pg + gen["cost"][3]
            else
                gen_cost[(n,i)] = 0.0
            end
        end
    end
    
    return JuMP.@objective(pm.model, Min,
        sum(branch["construction_cost"]*_PM.var(pm, 1, :branch_ne, i) for (i,branch) in _PM.nws(pm)[1][:ne_branch])
        +
        sum(gen["construction_cost"]*_PM.var(pm, 1, :pgslackabs, i) for (i,gen) in _PM.nws(pm)[1][:gen])
        +
        sum(
            sum( 
                gen_cost[(n,i)] 
                for (i,gen) in nw_ref[:gen] 
            )
            for (n, nw_ref) in _PM.nws(pm)
        )
    )
end

""
function constraint_ne_ohms_yt_from(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = _PM.ref(pm, nw, :ne_branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PM.calc_branch_y(branch)
    tr, ti = _PM.calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]

    vad_min = _PM.ref(pm, nw, :off_angmin)
    vad_max = _PM.ref(pm, nw, :off_angmax)

    constraint_ne_ohms_yt_from(pm, nw, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
end


function constraint_ne_ohms_yt_from(pm::_PM.AbstractDCPModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, vad_min, vad_max)
    p_fr  = _PM.var(pm, n, :p_ne, f_idx)
    va_fr = _PM.var(pm, n,   :va, f_bus)
    va_to = _PM.var(pm, n,   :va, t_bus)
    z = _PM.var(pm, 1, :branch_ne, i)

    JuMP.@constraint(pm.model, p_fr <= -b*(va_fr - va_to + vad_max*(1-z)) )
    JuMP.@constraint(pm.model, p_fr >= -b*(va_fr - va_to + vad_min*(1-z)) )
end

""
function constraint_ne_ohms_yt_to(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = _PM.ref(pm, nw, :ne_branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PM.calc_branch_y(branch)
    tr, ti = _PM.calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    vad_min = _PM.ref(pm, nw, :off_angmin)
    vad_max = _PM.ref(pm, nw, :off_angmax)

    constraint_ne_ohms_yt_to(pm, nw, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
end

"nothing to do, this model is symetric"
function constraint_ne_ohms_yt_to(pm::_PM.AbstractAPLossLessModels, n::Int, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, vad_min, vad_max)
end

""
function constraint_ne_voltage_angle_difference(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = _PM.ref(pm, nw, :ne_branch, i)
    f_idx = (i, branch["f_bus"], branch["t_bus"])

    vad_min = _PM.ref(pm, nw, :off_angmin)
    vad_max = _PM.ref(pm, nw, :off_angmax)

    constraint_ne_voltage_angle_difference(pm, nw, f_idx, branch["angmin"], branch["angmax"], vad_min, vad_max)
end

"`angmin*branch_ne[i] + vad_min*(1-branch_ne[i]) <= t[f_bus] - t[t_bus] <= angmax*branch_ne[i] + vad_max*(1-branch_ne[i])`"
function constraint_ne_voltage_angle_difference(pm::_PM.AbstractDCPModel, n::Int, f_idx, angmin, angmax, vad_min, vad_max)
    i, f_bus, t_bus = f_idx

    va_fr = _PM.var(pm, n, :va, f_bus)
    va_to = _PM.var(pm, n, :va, t_bus)
    z = _PM.var(pm, 1, :branch_ne, i)

    JuMP.@constraint(pm.model, va_fr - va_to <= angmax*z + vad_max*(1-z))
    JuMP.@constraint(pm.model, va_fr - va_to >= angmin*z + vad_min*(1-z))
end


""
function constraint_ne_thermal_limit_from(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = _PM.ref(pm, nw, :ne_branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    if !haskey(branch, "rate_a")
        Memento.error(_LOGGER, "constraint_thermal_limit_from_ne requires a rate_a value on all branches, calc_thermal_limits! can be used to generate reasonable values")
    end

    constraint_ne_thermal_limit_from(pm, nw, i, f_idx, branch["rate_a"])
end

""
function constraint_ne_thermal_limit_from(pm::_PM.AbstractActivePowerModel, n::Int, i, f_idx, rate_a)
    p_fr = _PM.var(pm, n, :p_ne, f_idx)
    z = _PM.var(pm, 1, :branch_ne, i)

    JuMP.@constraint(pm.model, p_fr <=  rate_a*z)
    JuMP.@constraint(pm.model, p_fr >= -rate_a*z)
end

# "`p_ne[f_idx]^2 + q_ne[f_idx]^2 <= (rate_a * branch_ne[i])^2`"
# function constraint_ne_thermal_limit_from(pm::_PM.AbstractPowerModel, n::Int, i, f_idx, rate_a)
#     p_fr = _PM.var(pm, n, :p_ne, f_idx)
#     q_fr = _PM.var(pm, n, :q_ne, f_idx)
#     z = _PM.var(pm, 1, :branch_ne, i)
#
#     JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= rate_a^2*z^2)
# end

""
function constraint_ne_thermal_limit_to(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = _PM.ref(pm, nw, :ne_branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    if !haskey(branch, "rate_a")
        Memento.error(_LOGGER, "constraint_thermal_limit_to_ne requires a rate_a value on all branches, calc_thermal_limits! can be used to generate reasonable values")
    end

    constraint_ne_thermal_limit_to(pm, nw, i, t_idx, branch["rate_a"])
end

""
function constraint_ne_thermal_limit_to(pm::_PM.AbstractActivePowerModel, n::Int, i, t_idx, rate_a)
    p_to = _PM.var(pm, n, :p_ne, t_idx)
    z = _PM.var(pm, 1, :branch_ne, i)

    JuMP.@constraint(pm.model, p_to <=  rate_a*z)
    JuMP.@constraint(pm.model, p_to >= -rate_a*z)
end

# "`p_ne[t_idx]^2 + q_ne[t_idx]^2 <= (rate_a * branch_ne[i])^2`"
# function constraint_ne_thermal_limit_to(pm::_PM.AbstractPowerModel, n::Int, i, t_idx, rate_a)
#     p_to = _PM.var(pm, n, :p_ne, t_idx)
#     q_to = _PM.var(pm, n, :q_ne, t_idx)
#     z = _PM.var(pm, 1, :branch_ne, i)
#
#     JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= rate_a^2*z^2)
# end


"variable: `pgslack[i]` for `i` in `gen`s"
function variable_gen_power_real_slack(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=false, report::Bool=true)
    pgslack = _PM.var(pm, nw)[:pgslack] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :gen)], base_name="$(nw)_pgslack",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :gen, i), "pgslack_start")
    )
    report && _PM.sol_component_value(pm, nw, :gen, :pgslack, _PM.ids(pm, nw, :gen), pgslack)
end


function variable_gen_power_real_slack_abs(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=false, report::Bool=true)
    pgslackabs = _PM.var(pm, nw)[:pgslackabs] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :gen)], base_name="$(nw)_pgslackabs",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :gen, i), "pgslackabs_start", 0.0)
    )
    report && _PM.sol_component_value(pm, nw, :gen, :pgslackabs, _PM.ids(pm, nw, :gen), pgslackabs)
end


"variable: `qgslack[i]` for `i` in `gen`s"
function variable_gen_power_imaginary_slack(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=false, report::Bool=true)
    qgslack = _PM.var(pm, nw)[:qgslack] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :gen)], base_name="$(nw)_qgslack",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :gen, i), "qgslack_start")
    )
    report && _PM.sol_component_value(pm, nw, :gen, :qgslack, _PM.ids(pm, nw, :gen), qgslack)
end

function variable_gen_power_imaginary_slack_abs(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=false, report::Bool=true)
    qgslackabs = _PM.var(pm, nw)[:qgslackabs] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :gen)], base_name="$(nw)_qgslackabs",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :gen, i), "qgslackabs_start", 0.0)
    )
    report && _PM.sol_component_value(pm, nw, :gen, :qgslackabs, _PM.ids(pm, nw, :gen), qgslackabs)
end


""
function constraint_gen_slack_abs(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    constraint_gen_slack_abs(pm, nw, i)
end

""
function constraint_gen_slack_abs(pm::_PM.AbstractPowerModel, n::Int, i::Int)
    pgslack = _PM.var(pm, n, :pgslack, i)
    pgslackabs = _PM.var(pm, n, :pgslackabs, i)
    # qgslack = _PM.var(pm, n, :qgslack, i)
    # qgslackabs = _PM.var(pm, n, :qgslackabs, i)
    JuMP.@constraint(pm.model, pgslack <= pgslackabs)
    JuMP.@constraint(pm.model, -pgslack <= pgslackabs)
    # JuMP.@constraint(pm.model, qgslack <= qgslackabs)
    # JuMP.@constraint(pm.model, -qgslack <= qgslackabs)
end


""
function constraint_ne_power_balance_gen_slack(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = _PM.ref(pm, nw, :bus, i)
    bus_arcs = _PM.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = _PM.ref(pm, nw, :bus_arcs_dc, i)
    bus_arcs_ne = _PM.ref(pm, nw, :ne_bus_arcs, i)
    bus_arcs_sw = _PM.ref(pm, nw, :bus_arcs_sw, i)
    bus_gens = _PM.ref(pm, nw, :bus_gens, i)
    bus_loads = _PM.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PM.ref(pm, nw, :bus_shunts, i)
    bus_storage = _PM.ref(pm, nw, :bus_storage, i)

    bus_pd = Dict(k => _PM.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => _PM.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_ne_power_balance_gen_slack(pm, nw, i, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_arcs_ne, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
end


"""
sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - sum(pd[d] for d in bus_loads) - sum(gs[s] for s in bus_shunts)*v^2
sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - sum(qd[d] for d in bus_loads) + sum(bs[s] for s in bus_shunts)*v^2
"""
function constraint_ne_power_balance_gen_slack(pm::_PM.AbstractPowerModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_arcs_ne, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vm   = 1.0 #_PM.var(pm, n, :vm, i)
    p    = get(_PM.var(pm, n),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    # q    = get(_PM.var(pm, n),    :q, Dict()); _PM._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(_PM.var(pm, n),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    # qg   = get(_PM.var(pm, n),   :qg, Dict()); _PM._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PM.var(pm, n),   :ps, Dict()); _PM._check_var_keys(ps, bus_storage, "active power", "storage")
    # qs   = get(_PM.var(pm, n),   :qs, Dict()); _PM._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(_PM.var(pm, n),  :psw, Dict()); _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    # qsw  = get(_PM.var(pm, n),  :qsw, Dict()); _PM._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = get(_PM.var(pm, n), :p_dc, Dict()); _PM._check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    # q_dc = get(_PM.var(pm, n), :q_dc, Dict()); _PM._check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")
    p_ne = get(_PM.var(pm, n), :p_ne, Dict()); _PM._check_var_keys(p_ne, bus_arcs_ne, "active power", "ne_branch")
    # q_ne = get(_PM.var(pm, n), :q_ne, Dict()); _PM._check_var_keys(q_ne, bus_arcs_ne, "reactive power", "ne_branch")
    
    pgslack   = get(_PM.var(pm, 1),   :pgslack, Dict()); _PM._check_var_keys(pgslack, bus_gens, "active power", "generator")
    # qgslack   = get(_PM.var(pm, 1),   :qgslack, Dict()); _PM._check_var_keys(qgslack, bus_gens, "reactive power", "generator")

    # the check "typeof(p[arc]) <: JuMP.NonlinearExpression" is required for the
    # case when p/q are nonlinear expressions instead of decision variables
    # once NLExpressions are first order in JuMP it should be possible to
    # remove this.
    nl_form = length(bus_arcs) > 0 && (typeof(p[iterate(bus_arcs)[1]]) <: JuMP.NonlinearExpression)

    # if !nl_form
        cstr_p = JuMP.@constraint(pm.model,
            sum(p[a] for a in bus_arcs)
            + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
            + sum(psw[a_sw] for a_sw in bus_arcs_sw)
            + sum(p_ne[a] for a in bus_arcs_ne)
            ==
            sum(pg[g] for g in bus_gens)
            - sum(ps[s] for s in bus_storage)
            - sum(pd for (i,pd) in bus_pd)
            - sum(gs for (i,gs) in bus_gs)*vm^2
            + sum(pgslack[g] for g in bus_gens)
        )
    # else
    #     cstr_p = JuMP.@NLconstraint(pm.model,
    #         sum(p[a] for a in bus_arcs)
    #         + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
    #         + sum(psw[a_sw] for a_sw in bus_arcs_sw)
    #         + sum(p_ne[a] for a in bus_arcs_ne)
    #         ==
    #         sum(pg[g] for g in bus_gens)
    #         - sum(ps[s] for s in bus_storage)
    #         - sum(pd for (i,pd) in bus_pd)
    #         - sum(gs for (i,gs) in bus_gs)*vm^2
    #         + sum(pgslack[g] for g in bus_gens)
    #     )
    # end

    # if !nl_form
    #     cstr_q = JuMP.@constraint(pm.model,
    #         sum(q[a] for a in bus_arcs)
    #         + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
    #         + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
    #         + sum(q_ne[a] for a in bus_arcs_ne)
    #         ==
    #         sum(qg[g] for g in bus_gens)
    #         - sum(qs[s] for s in bus_storage)
    #         - sum(qd for (i,qd) in bus_qd)
    #         + sum(bs for (i,bs) in bus_bs)*vm^2
    #         + sum(qgslack[g] for g in bus_gens)
    #     )
    # else
    #     cstr_q = JuMP.@NLconstraint(pm.model,
    #         sum(q[a] for a in bus_arcs)
    #         + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
    #         + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
    #         + sum(q_ne[a] for a in bus_arcs_ne)
    #         ==
    #         sum(qg[g] for g in bus_gens)
    #         - sum(qs[s] for s in bus_storage)
    #         - sum(qd for (i,qd) in bus_qd)
    #         + sum(bs for (i,bs) in bus_bs)*vm^2
    #         + sum(qgslack[g] for g in bus_gens)
    #     )
    # end

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = NaN
        # sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end
