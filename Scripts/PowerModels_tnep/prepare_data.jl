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
