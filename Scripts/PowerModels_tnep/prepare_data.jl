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



# function aggregate_demand_data!(grid_data)

#     grid_data["aggregated_data"] = Dict{String, Any}([area => Dict{String, Any}("demand" => 0) for (key, area) in grid_data["areas"]])

#     for (l, load) in grid_data["load"]
#         load_bus = load["load_bus"]
#         area_code = grid_data["bus"]["$load_bus"]["area"]
#         area = grid_data["areas"]["$area_code"]
#         grid_data["aggregated_data"][area]["demand"] = grid_data["aggregated_data"][area]["demand"] + load["pd"] * grid_data["baseMVA"]
#     end

#     return grid_data
# end



# function fix_data!(data; hvdc = true)
#     # Find isolated buses and put their demand zero
#     bus_arcs = Dict{String, Any}([b => [] for (b, bus) in data["bus"]])
#     for (b, branch) in data["branch"]
#         fbus = branch["f_bus"]
#         tbus = branch["t_bus"]
#         if haskey(bus_arcs, "$fbus")
#             push!(bus_arcs["$fbus"], parse(Int, b))
#         end
#         if haskey(bus_arcs, "$tbus")
#             push!(bus_arcs["$tbus"], parse(Int, b))
#         end
#         branch["tap"] = 1.0
#         branch["shift"] = 0.0
#     end

#     if hvdc == true
#         for (c, conv) in data["convdc"]
#             cbus = conv["busac_i"]
#             push!(bus_arcs["$cbus"], parse(Int, c))
#         end
#     end

#     for (l, load) in data["load"]
#         load_bus = load["load_bus"]
#         if isempty(bus_arcs["$load_bus"])
#             load["pd"] = 0.0
#             load["qd"] = 0.0
#         end
#     end

#     # generator data comming from matlab model seems two orders of magnitude too small
#     for (g, gen) in data["gen"] 
#         gen["cost"] = gen["cost"] .* data["baseMVA"]
#     end

#     return data
# end



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
