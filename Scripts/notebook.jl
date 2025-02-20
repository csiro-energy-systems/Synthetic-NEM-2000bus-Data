### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ ec6d8bbd-057d-4fff-8555-225a5af5140f
begin
	using Pkg
	Pkg.add("PowerModels")
	Pkg.add("PowerModelsACDC")
	Pkg.add("JuMP")
	Pkg.add("CSV")
	Pkg.add("DataFrames")
	Pkg.add("Ipopt")
	Pkg.add("Gurobi")
	Pkg.add("Plots")
end

# ╔═╡ 2e749793-92b6-48be-bcf0-260699abb915
begin
	using PowerModels
	using PowerModelsACDC
	using JuMP
	# using InfrastructureModels
	# using CbaOPF
	using Ipopt
	using Plots
	# using PlotlyJS
	using DataFrames
	using CSV
	using Gurobi
	
	# Create short hands for the most important ones
	const _PM = PowerModels
	const _PMACDC = PowerModelsACDC
	# const _IM = InfrastructureModels

	### Assign solvers
	ac_solver =  JuMP.optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1000, "print_level" => 0)
	dc_solver =  JuMP.optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0, "method" => 2) #  https://www.gurobi.com/documentation/current/refman/method.html#parameter:Method 
	lpac_solver =  JuMP.optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0)
	soc_solver =  JuMP.optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1000, "print_level" => 0)
	
	ENV["GUROBI_HOME"]="/Library/gurobi902/mac64"
	ENV["GRB_LICENSE_FILE"]="/Users/hei06j/gurobi/gurobi_11.lic"
end

# ╔═╡ aea5ce82-c53c-4830-ab03-72407778ac92
begin
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
	
	    if _IM.report_duals(pm)
	        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
	        sol(pm, n, :bus, i)[:lam_kcl_i] = NaN
	        # sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
	    end
	end

end

# ╔═╡ 3888baba-53e4-4c87-b1b4-332063eaf7b4
begin
	function add_area_dict!(data_nem)
		    data_nem["areas"] = Dict{String, Any}()
		    data_nem["areas"]["1"] = "NSW"
		    data_nem["areas"]["2"] = "VIC"
		    data_nem["areas"]["3"] = "QLD"
		    data_nem["areas"]["4"] = "SA"
		    data_nem["areas"]["5"] = "TAS"
		end
	
	function fix_hvdc_data_issues!(data; no_bass = false, no_terra = false, no_murray = false)
		
		    if no_bass == false
		        # BASS LINK 
		        data["branch"]["543"]["br_status"] = 0
		        delete!(data["bus"], "2113")
		        delete!(data["gen"], "264")
		        delete!(data["gen"], "265")
		    end
		
		    if no_terra == false
		        # TERRANORALINK
		        data["branch"]["1568"]["br_status"] = 0
		        data["branch"]["1569"]["br_status"] = 0
		        data["branch"]["1570"]["br_status"] = 0
		    end
		
		    if no_murray ==false
		        # MURRAYLINK
		        data["branch"]["1949"]["br_status"] = 0
		        data["gen"]["222"]["gen_status"] = 0
		    end
		
		    return data
		end
	
	function add_demand_data!(data)
		    for (l, load) in data["load"]
		        # Superior bound on voluntary load reduction (not consumed power) as a fraction of the total reference demand (0 ≤ pred_rel_max ≤ 1)
		        load["pred_rel_max"] = 0
		
		        # Compensation for consuming less (i.e. voluntary demand reduction) (€/MWh)
		        load["cost_red"] = 1000
		
		        # Compensation for load curtailment (i.e. involuntary demand reduction) (€/MWh)
		        load["cost_curt"] = 1e4 
		
		        # Whether load is flexible (boolean)
		        load["flex"] = 1
		
		        # Power factor angle θ, giving the reactive power as Q = P ⨉ tan(θ)
		        load["pf_angle"] = atan(load["qd"]/load["pd"])
		
		        # Rescale cost and power input values to the p.u. values used internally in the model
		        rescale_cost = x -> x*data["baseMVA"]
		        _PM._apply_func!(load, "cost_red", rescale_cost)
		        _PM._apply_func!(load, "cost_curt", rescale_cost)
		    end
		    return data
		end
	
	function merge_all_parallel_lines(data)
		    ft_buses = Dict{String, Any}()
		    for (b, branch) in data["branch"]
		        ftbus = join([branch["f_bus"], "-", branch["t_bus"]])
		        tfbus = join([branch["t_bus"], "-", branch["f_bus"]])
		        if !haskey(ft_buses, ftbus) && !haskey(ft_buses, tfbus)
		            push!(ft_buses, ftbus => b)
		        elseif haskey(ft_buses, ftbus)
		            branch_id = ft_buses[ftbus]
		            merge_branches!(data, branch_id, branch)
		            delete!(data["branch"], b)
		        elseif haskey(ft_buses, tfbus)
		            branch_id = ft_buses[tfbus]
		            merge_branches!(data, branch_id, branch)
		            delete!(data["branch"], b)
		        end
		    end
		    
		    return data
		end
	
	function merge_branches!(data, b_idx, parallel_br)
		    data["branch"][b_idx]["rate_a"] = data["branch"][b_idx]["rate_a"] + parallel_br["rate_a"]
		    data["branch"][b_idx]["rate_b"] = data["branch"][b_idx]["rate_b"] + parallel_br["rate_b"]
		    data["branch"][b_idx]["rate_c"] = data["branch"][b_idx]["rate_c"] + parallel_br["rate_c"]
		    if haskey(data["branch"][b_idx], "c_rating") && haskey(parallel_br, "c_rating")
		        data["branch"][b_idx]["c_rating"] = data["branch"][b_idx]["c_rating"] + parallel_br["c_rating"]
		    end
		    data["branch"][b_idx]["br_r"] = (data["branch"][b_idx]["br_r"] * parallel_br["br_r"]) / (data["branch"][b_idx]["br_r"] + parallel_br["br_r"])
		    data["branch"][b_idx]["br_x"] = (data["branch"][b_idx]["br_x"] * parallel_br["br_x"]) / (data["branch"][b_idx]["br_x"] + parallel_br["br_x"])
		    data["branch"][b_idx]["b_fr"] = data["branch"][b_idx]["b_fr"] + parallel_br["b_fr"]
		    data["branch"][b_idx]["b_to"] = data["branch"][b_idx]["b_to"] + parallel_br["b_to"]
		    data["branch"][b_idx]["g_fr"] = data["branch"][b_idx]["g_fr"] + parallel_br["g_fr"]
		    data["branch"][b_idx]["g_to"] = data["branch"][b_idx]["g_to"] + parallel_br["g_to"]
		end
	
	function get_demand_data(scenario, year; verbose = true)
		    data_folder = joinpath("data", "2022 Final ISP Model", scenario, "Traces", "demand")
		    all_demand_files = readdir(data_folder)
		
		    demand = Dict{String, Any}()
		    demand_series = Dict{String, Any}()
		
		    for file in all_demand_files
		        region = file[1:(collect(findfirst("_", file)).-1)[1]]
		        verbose ? print("Reading demand trace for ", region, "\n") : nothing
		        demand_region = _DF.DataFrame(CSV.File(joinpath(data_folder, file)))
		        demand[region] = demand_region[demand_region[!, :Year] .== year, :]
		    end
		
		    demand_states = aggregate_demand(demand)
		    
		
		    for (d, demand) in demand_states
		        demand_series[d] = zeros(1, 48 * _DF.nrow(demand))
		        for row in 1:_DF.nrow(demand)
		            demand_series[d][(((row-1) * 48) + 1) : (row * 48)] = collect(values(demand[row, 4:end]))
		        end
		    end
		
		    return demand_series
		end
	
	function prepare_data(data_path; merge_parallel_lines = true, hvdc = true, no_bass = false,  no_terra = false,  no_murray = false)
		    # Get grid data from the NEM 2000 bus model m-file 
		    data = _PM.parse_file(data_path)
		    # Assign buses to states
		    add_area_dict!(data)
		
		    if hvdc
		        # Process data to fit into PMACDC model
		        _PMACDC.process_additional_data!(data)
		        # Delete DC lines which have been modelled as AC lines
		        fix_hvdc_data_issues!(data, no_bass = no_bass, no_terra = no_terra, no_murray = no_murray)
		    end
		
		    # Extend data model with flexible demand, to be able to do demand shedding if required
		    add_demand_data!(data)
		
		    # Aggregate demand data per state to modulate with hourly traces
		    # aggregate_demand_data!(data)
		
		    # fix data issues, e.g. putting generation cost in € / pu:
		    # fix_data!(data; hvdc = hvdc)
		
		    if merge_parallel_lines
		        merge_all_parallel_lines(data)
		    end
		
		    return data
		end
	
	function prepare_hourly_opf_data!(hourly_data, grid_data, demand_series, pv_series, wind_series, time_stamp)
		    hour = time_stamp
		    hourly_data["pst"] = Dict{String, Any}()
		
		    for (l, load) in hourly_data["load"]
		        load_bus = load["load_bus"]
		        area_code = hourly_data["bus"]["$load_bus"]["area"]
		        area = hourly_data["areas"]["$area_code"]
		        demand_ratio = demand_series[!,area][hour]
		        if load["pd"] >= 0 
		            load["pd"] = grid_data["load"][l]["pd"] * demand_ratio
		            load["qd"] = grid_data["load"][l]["qd"] * demand_ratio
		        end
		    end
		
		    for (g, gen) in hourly_data["gen"]
		        trace = 1
		        gen_bus = gen["gen_bus"]
		        area_code = hourly_data["bus"]["$gen_bus"]["area"]
		        area = hourly_data["areas"]["$area_code"]
		        if gen["type"] == "Wind"
		            trace = wind_series[!,area][hour]
		        elseif gen["type"] == "Solar"
		            if !isempty(pv_series[!,area])
		                trace = pv_series[!,area][hour]
		            end
		        end
		        gen["pmax"] = grid_data["gen"][g]["pmax"] * trace
		    end
		    
		    return hourly_data
		end
	
	function prepare_mn_opf_data!(opf_mn_data, grid_data, demand_series, pv_series, wind_series, hours)
		    for hour in hours
		        opf_mn_data["nw"]["$hour"]["per_unit"] = true
		
		        for (l, load) in opf_mn_data["nw"]["$hour"]["load"]
		            load_bus = load["load_bus"]
		            area_code = opf_mn_data["nw"]["$hour"]["bus"]["$load_bus"]["area"]
		            area = opf_mn_data["nw"]["$hour"]["areas"]["$area_code"]
		            demand_ratio = demand_series[!,area][hour]
		            if load["pd"] >= 0 
		                load["pd"] = grid_data["load"][l]["pd"] * demand_ratio * 0.9
		                load["qd"] = grid_data["load"][l]["qd"] * demand_ratio * 0.9
		            end
		        end
		
		        for (g, gen) in opf_mn_data["nw"]["$hour"]["gen"]
		            trace = 1
		            gen_bus = gen["gen_bus"]
		            area_code = opf_mn_data["nw"]["$hour"]["bus"]["$gen_bus"]["area"]
		            area = opf_mn_data["nw"]["$hour"]["areas"]["$area_code"]
		            if gen["type"] == "Wind"
		                trace = wind_series[!,area][hour] * 1.1
		            elseif gen["type"] == "Solar"
		                if !isempty(pv_series[!,area])
		                    trace = pv_series[!,area][hour] * 1.1
		                end
		            end
		            gen["pmax"] = grid_data["gen"][g]["pmax"] * trace
		        end
		    end
		    return opf_mn_data
		end
	
		
		# Run hourly OPF calcuations
	function run_mn_opf(opf_data, hours, demand_series, pv_series, wind_series; formulation="DC", verbose = true)
		    hourly_data = deepcopy(opf_data)
		    # Create dictionaries for inspection of results
		    pf = Dict{String, Any}(["$hour" => zeros(1, maximum(parse.(Int, collect(keys(hourly_data["branch"]))))) for hour in hours])
		    pf_mw = Dict{String, Any}(["$hour" => zeros(length(opf_data["branch"])) for hour in hours])
		    pfdc = Dict{String, Any}(["$hour" => zeros(length(opf_data["branchdc"])) for hour in hours])
		    pcurt = Dict{String, Any}(["$hour" => zeros(1, maximum(parse.(Int, collect(keys(hourly_data["branch"]))))) for hour in hours])
		    pd = Dict{String, Any}(["$hour" => zeros(length(opf_data["load"])) for hour in hours])
		    pflex = Dict{String, Any}(["$hour" => zeros(length(opf_data["load"])) for hour in hours])
		    pgmax = Dict{String, Any}(["$hour" => zeros(length(opf_data["load"])) for hour in hours])
		    pg = Dict{String, Any}(["$hour" => zeros(length(opf_data["load"])) for hour in hours])
		    pf_tot = Dict{String, Any}(["$hour" => zeros(length(opf_data["branch"])) for hour in hours])
		    pc_tot = Dict{String, Any}(["$hour" => zeros(length(opf_data["convdc"])) for hour in hours])
		    branch_duals = Dict{String, Any}(["$hour" => zeros(length(opf_data["branch"])) for hour in hours])
		    bus_duals = Dict{String, Any}(["$hour" => zeros(length(opf_data["bus"])) for hour in hours])
		    bus_ids = []
		    branch_ids = []
		
		    for hour in hours
		        # Write hourly pv, wind and demand traces into opf data
		        prepare_hourly_opf_data!(hourly_data, opf_data, demand_series, pv_series, wind_series, hour)
		
		        # Solve OPF
		        if formulation == "AC"
		            opf_result = CbaOPF.solve_cbaopf(hourly_data, _PM.ACPPowerModel, ac_solver, setting = s)
		        elseif formulation == "LPAC"
		            opf_result = CbaOPF.solve_cbaopf(hourly_data, _PM.LPACCPowerModel, lpac_solver, setting = s)
		        elseif formulation == "SOC"
		            opf_result = CbaOPF.solve_cbaopf(hourly_data, _PM.SOCWRPowerModel, soc_solver, setting = s)
		        elseif formulation == "DC"
		            opf_result = CbaOPF.solve_cbaopf(hourly_data, _PM.DCPPowerModel, dc_solver, setting = s)
		        end
		        # calculate and print some more information based on results
		        if haskey(opf_result["solution"], "load")
		            pflex["$hour"] = [opf_result["solution"]["load"]["$l"]["pflex"] for l in sort(parse.(Int, collect(keys(opf_result["solution"]["load"]))))]
		            pd["$hour"] = [hourly_data["load"]["$l"]["pd"] for l in sort(parse.(Int, collect(keys(opf_result["solution"]["load"]))))]
		            pgmax["$hour"] = [hourly_data["gen"]["$g"]["pmax"] for g in sort(parse.(Int, collect(keys(opf_result["solution"]["gen"]))))]
		            for l in sort(collect(parse.(Int, keys(opf_result["solution"]["load"]))))
		                pcurt["$hour"][1, l] = opf_result["solution"]["load"]["$l"]["pcurt"] 
		            end
		            pcurt_tot = sum([opf_result["solution"]["load"]["$l"]["pcurt"] for l in sort(parse.(Int, collect(keys(opf_result["solution"]["load"]))))]) * hourly_data["baseMVA"]
		            for b in sort(collect(parse.(Int, keys(opf_result["solution"]["branch"]))))
		                pf["$hour"][1, b] = opf_result["solution"]["branch"]["$b"]["pf"] ./ hourly_data["branch"]["$b"]["rate_a"]
		            end
		            pf_mw["$hour"] = [opf_result["solution"]["branch"]["$b"]["pf"] for b in sort(collect(parse.(Int, keys(opf_result["solution"]["branch"]))))]
		            pf_tot["$hour"] = [opf_result["solution"]["branch"]["$b"]["pf"] + opf_result["solution"]["branch"]["$b"]["pt"] for b in sort(collect(parse.(Int, keys(opf_result["solution"]["branch"]))))]
		            pc_tot["$hour"] = [opf_result["solution"]["convdc"]["$c"]["pgrid"] + opf_result["solution"]["convdc"]["$c"]["pdc"] for c in sort(collect(parse.(Int, keys(opf_result["solution"]["convdc"]))))]
		            pfdc["$hour"] = [opf_result["solution"]["branchdc"]["$b"]["pf"] for b in sort(parse.(Int, collect(keys(opf_result["solution"]["branchdc"]))))] ./ [hourly_data["branchdc"]["$b"]["rateA"] for b in sort(parse.(Int, collect(keys(opf_result["solution"]["branchdc"]))))]
		            pg["$hour"] = [opf_result["solution"]["gen"]["$g"]["pg"] for g in sort(parse.(Int, collect(keys(opf_result["solution"]["gen"]))))]
		            branch_duals["$hour"] = [opf_result["solution"]["branch"]["$b"]["mu_sm_to"] for b in sort(collect(parse.(Int, keys(opf_result["solution"]["branch"]))))]
		            bus_duals["$hour"] = [opf_result["solution"]["bus"]["$b"]["lam_kcl_r"] for b in sort(collect(parse.(Int, keys(opf_result["solution"]["bus"]))))]
		        else
		            pcurt_tot = 0
		        end
		        if verbose
		            # calculate some charactersitics for result inspecttion
		            pd_max = sum([load["pd"] for (l, load) in opf_data["load"]]) * hourly_data["baseMVA"]
		            pg_max = sum([gen["pmax"] for (g, gen) in opf_data["gen"]]) * hourly_data["baseMVA"]
		            pdh_max = sum([load["pd"] for (l, load) in hourly_data["load"]]) * hourly_data["baseMVA"]
		            pgh_max = sum([gen["pmax"] for (g, gen) in hourly_data["gen"]]) * hourly_data["baseMVA"]
		            # print charactersiticss
		            print("Grid Data: Total demand = ", pd_max, " MW, Total generation = ", pg_max, " MW","\n")
		            print("Hour ", hour,": Total demand = ", pdh_max, " MW, Total generation = ", pgh_max, " MW","\n")
		            # Write out general information
		            print("Hour: ", hour, " -> ", opf_result["termination_status"], " in ", opf_result["solve_time"], " seconds.", "\n")
		            print("Total curtailed load = ", pcurt_tot, " MW, ", pcurt_tot / pdh_max * 100,"%", "\n")
		        end
		
		        bus_ids = sort(collect(parse.(Int, keys(opf_result["solution"]["bus"]))))
		        branch_ids = sort(collect(parse.(Int, keys(opf_result["solution"]["branch"]))))
		    end
		
		    return pf, pf_mw, pfdc, pcurt, pd, pflex, pgmax, pg, pf_tot, pc_tot, bus_duals, branch_duals, bus_ids, branch_ids
		end
end

# ╔═╡ 3d3aa08d-fd2c-4837-8530-e7d4f35ad692
begin
	# You can choose select certain hours or a full year for the analysis: full year is 17520 entries
	# selected_hours = Dict{String, Any}("hour_range" => start hour:end hour)
	selected_hours = Dict{String, Any}("hour_range" => 1:48)
	timeseries_folder = "/Users/hei06j/Documents/Repositories/Remote/Synthetic-NEM-2000bus-Data/Data/Timeseries"
	demand_series = CSV.read(timeseries_folder*"/timeseries_demand.csv", DataFrame)
	wind_series = CSV.read(timeseries_folder*"/timeseries_wind.csv", DataFrame)
	pv_series = CSV.read(timeseries_folder*"/timeseries_pv.csv", DataFrame)
	
	# Select OPF method opf ∈ {"AC", "DC", "LPAC", "SOC"}
	opf = "DC"
		
	data_file_tnep = "/Users/hei06j/Documents/Repositories/Remote/Synthetic-NEM-2000bus-Data/snem2000_tnep.m"
	
	parse_file(data_file_tnep)
	# tnep_data = prepare_data(data_file_tnep; merge_parallel_lines = false, hvdc = false)
	
	# for (i, gen) in tnep_data["gen"]
	#     if gen["type"] == "Fossil" || gen["type"] == "Hydro"
	#         gen["construction_cost"] = 1E8 * tnep_data["baseMVA"]
	#     elseif gen["fuel"] == "Solar"
	#         gen["construction_cost"] = 4E6 * tnep_data["baseMVA"]  # https://www.mcgqs.com.au/media/australian-solar-farms/
	#     elseif gen["fuel"] == "Wind"
	#         gen["construction_cost"] = 8E6 * tnep_data["baseMVA"]  
	#     end
	# end

end

# ╔═╡ Cell order:
# ╟─ec6d8bbd-057d-4fff-8555-225a5af5140f
# ╟─2e749793-92b6-48be-bcf0-260699abb915
# ╟─aea5ce82-c53c-4830-ab03-72407778ac92
# ╟─3888baba-53e4-4c87-b1b4-332063eaf7b4
# ╠═3d3aa08d-fd2c-4837-8530-e7d4f35ad692
