# Define all packages that are needed
# using Pkg
# Pkg.activate("./")
# using ISPhvdc
using PowerModels
using PowerModelsACDC
using JuMP
using InfrastructureModels
using CbaOPF
using Ipopt
using Plots
# using PlotlyJS
using DataFrames
using CSV
using Gurobi

# Create short hands for the most important ones
const _PM = PowerModels
const _PMACDC = PowerModelsACDC
const _IM = InfrastructureModels

###########################################################
################## INPUT SECTION ##########################
include("tnep_PM_functions.jl")
include("prepare_data.jl")

# You can choose select certain hours or a full year for the analysis: full year is 17520 entries
# selected_hours = Dict{String, Any}("hour_range" => start hour:end hour)
selected_hours = Dict{String, Any}("hour_range" => 1:48)
timeseries_folder = "Data/Timeseries"
demand_series = CSV.read(timeseries_folder*"/timeseries_demand.csv", DataFrame)
wind_series = CSV.read(timeseries_folder*"/timeseries_wind.csv", DataFrame)
pv_series = CSV.read(timeseries_folder*"/timeseries_pv.csv", DataFrame)

# Select OPF method opf ∈ {"AC", "DC", "LPAC", "SOC"}
opf = "DC"

### Assign solvers
ac_solver =  JuMP.optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1000, "print_level" => 0)
dc_solver =  JuMP.optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0, "method" => 2) #  https://www.gurobi.com/documentation/current/refman/method.html#parameter:Method 
lpac_solver =  JuMP.optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0)
soc_solver =  JuMP.optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1000, "print_level" => 0)

ENV["GUROBI_HOME"]="/Library/gurobi902/mac64"
ENV["GRB_LICENSE_FILE"]="/Users/hei06j/gurobi/gurobi.lic"
############ END INPUT SECTION ##############################
#############################################################


#############################################################
#################### START METHODOLOGY ######################
# ### Optimisation settings for CbaOPF.jl
# s = Dict("output" => Dict("branch_flows" => true, "duals" => true), "conv_losses_mp" => true)

# # Test case data
# data_path_hvdc = "snem2000_ACDC.m"

# function prepare_data(data_path; merge_parallel_lines = true, hvdc = true, no_bass = false,  no_terra = false,  no_murray = false)
#     # Get grid data from the NEM 2000 bus model m-file 
#     data = _PM.parse_file(data_path)
#     # Assign buses to states
#     add_area_dict!(data)

#     if hvdc
#         # Process data to fit into PMACDC model
#         _PMACDC.process_additional_data!(data)
#         # Delete DC lines which have been modelled as AC lines
#         fix_hvdc_data_issues!(data, no_bass = no_bass, no_terra = no_terra, no_murray = no_murray)
#     end

#     # Extend data model with flexible demand, to be able to do demand shedding if required
#     add_demand_data!(data)

#     # Aggregate demand data per state to modulate with hourly traces
#     # aggregate_demand_data!(data)

#     # fix data issues, e.g. putting generation cost in € / pu:
#     # fix_data!(data; hvdc = hvdc)

#     if merge_parallel_lines
#         merge_all_parallel_lines(data)
#     end

#     return data
# end


# opf_data = prepare_data(data_path_hvdc; merge_parallel_lines=true)

# ##
# # Select hours
# # hours = select_hours(year, selection = selected_hours)
# hours = selected_hours["hour_range"]

# # Run hourly OPF calcuations
# function run_mn_opf(opf_data, hours, demand_series, pv_series, wind_series; formulation="DC", verbose = true)
#     hourly_data = deepcopy(opf_data)
#     # Create dictionaries for inspection of results
#     pf = Dict{String, Any}(["$hour" => zeros(1, maximum(parse.(Int, collect(keys(hourly_data["branch"]))))) for hour in hours])
#     pf_mw = Dict{String, Any}(["$hour" => zeros(length(opf_data["branch"])) for hour in hours])
#     pfdc = Dict{String, Any}(["$hour" => zeros(length(opf_data["branchdc"])) for hour in hours])
#     pcurt = Dict{String, Any}(["$hour" => zeros(1, maximum(parse.(Int, collect(keys(hourly_data["branch"]))))) for hour in hours])
#     pd = Dict{String, Any}(["$hour" => zeros(length(opf_data["load"])) for hour in hours])
#     pflex = Dict{String, Any}(["$hour" => zeros(length(opf_data["load"])) for hour in hours])
#     pgmax = Dict{String, Any}(["$hour" => zeros(length(opf_data["load"])) for hour in hours])
#     pg = Dict{String, Any}(["$hour" => zeros(length(opf_data["load"])) for hour in hours])
#     pf_tot = Dict{String, Any}(["$hour" => zeros(length(opf_data["branch"])) for hour in hours])
#     pc_tot = Dict{String, Any}(["$hour" => zeros(length(opf_data["convdc"])) for hour in hours])
#     branch_duals = Dict{String, Any}(["$hour" => zeros(length(opf_data["branch"])) for hour in hours])
#     bus_duals = Dict{String, Any}(["$hour" => zeros(length(opf_data["bus"])) for hour in hours])
#     bus_ids = []
#     branch_ids = []

#     for hour in hours
#         # Write hourly pv, wind and demand traces into opf data
#         prepare_hourly_opf_data!(hourly_data, opf_data, demand_series, pv_series, wind_series, hour)

#         # Solve OPF
#         if formulation == "AC"
#             opf_result = CbaOPF.solve_cbaopf(hourly_data, _PM.ACPPowerModel, ac_solver, setting = s)
#         elseif formulation == "LPAC"
#             opf_result = CbaOPF.solve_cbaopf(hourly_data, _PM.LPACCPowerModel, lpac_solver, setting = s)
#         elseif formulation == "SOC"
#             opf_result = CbaOPF.solve_cbaopf(hourly_data, _PM.SOCWRPowerModel, soc_solver, setting = s)
#         elseif formulation == "DC"
#             opf_result = CbaOPF.solve_cbaopf(hourly_data, _PM.DCPPowerModel, dc_solver, setting = s)
#         end
#         # calculate and print some more information based on results
#         if haskey(opf_result["solution"], "load")
#             pflex["$hour"] = [opf_result["solution"]["load"]["$l"]["pflex"] for l in sort(parse.(Int, collect(keys(opf_result["solution"]["load"]))))]
#             pd["$hour"] = [hourly_data["load"]["$l"]["pd"] for l in sort(parse.(Int, collect(keys(opf_result["solution"]["load"]))))]
#             pgmax["$hour"] = [hourly_data["gen"]["$g"]["pmax"] for g in sort(parse.(Int, collect(keys(opf_result["solution"]["gen"]))))]
#             for l in sort(collect(parse.(Int, keys(opf_result["solution"]["load"]))))
#                 pcurt["$hour"][1, l] = opf_result["solution"]["load"]["$l"]["pcurt"] 
#             end
#             pcurt_tot = sum([opf_result["solution"]["load"]["$l"]["pcurt"] for l in sort(parse.(Int, collect(keys(opf_result["solution"]["load"]))))]) * hourly_data["baseMVA"]
#             for b in sort(collect(parse.(Int, keys(opf_result["solution"]["branch"]))))
#                 pf["$hour"][1, b] = opf_result["solution"]["branch"]["$b"]["pf"] ./ hourly_data["branch"]["$b"]["rate_a"]
#             end
#             pf_mw["$hour"] = [opf_result["solution"]["branch"]["$b"]["pf"] for b in sort(collect(parse.(Int, keys(opf_result["solution"]["branch"]))))]
#             pf_tot["$hour"] = [opf_result["solution"]["branch"]["$b"]["pf"] + opf_result["solution"]["branch"]["$b"]["pt"] for b in sort(collect(parse.(Int, keys(opf_result["solution"]["branch"]))))]
#             pc_tot["$hour"] = [opf_result["solution"]["convdc"]["$c"]["pgrid"] + opf_result["solution"]["convdc"]["$c"]["pdc"] for c in sort(collect(parse.(Int, keys(opf_result["solution"]["convdc"]))))]
#             pfdc["$hour"] = [opf_result["solution"]["branchdc"]["$b"]["pf"] for b in sort(parse.(Int, collect(keys(opf_result["solution"]["branchdc"]))))] ./ [hourly_data["branchdc"]["$b"]["rateA"] for b in sort(parse.(Int, collect(keys(opf_result["solution"]["branchdc"]))))]
#             pg["$hour"] = [opf_result["solution"]["gen"]["$g"]["pg"] for g in sort(parse.(Int, collect(keys(opf_result["solution"]["gen"]))))]
#             branch_duals["$hour"] = [opf_result["solution"]["branch"]["$b"]["mu_sm_to"] for b in sort(collect(parse.(Int, keys(opf_result["solution"]["branch"]))))]
#             bus_duals["$hour"] = [opf_result["solution"]["bus"]["$b"]["lam_kcl_r"] for b in sort(collect(parse.(Int, keys(opf_result["solution"]["bus"]))))]
#         else
#             pcurt_tot = 0
#         end
#         if verbose
#             # calculate some charactersitics for result inspecttion
#             pd_max = sum([load["pd"] for (l, load) in opf_data["load"]]) * hourly_data["baseMVA"]
#             pg_max = sum([gen["pmax"] for (g, gen) in opf_data["gen"]]) * hourly_data["baseMVA"]
#             pdh_max = sum([load["pd"] for (l, load) in hourly_data["load"]]) * hourly_data["baseMVA"]
#             pgh_max = sum([gen["pmax"] for (g, gen) in hourly_data["gen"]]) * hourly_data["baseMVA"]
#             # print charactersiticss
#             print("Grid Data: Total demand = ", pd_max, " MW, Total generation = ", pg_max, " MW","\n")
#             print("Hour ", hour,": Total demand = ", pdh_max, " MW, Total generation = ", pgh_max, " MW","\n")
#             # Write out general information
#             print("Hour: ", hour, " -> ", opf_result["termination_status"], " in ", opf_result["solve_time"], " seconds.", "\n")
#             print("Total curtailed load = ", pcurt_tot, " MW, ", pcurt_tot / pdh_max * 100,"%", "\n")
#         end

#         bus_ids = sort(collect(parse.(Int, keys(opf_result["solution"]["bus"]))))
#         branch_ids = sort(collect(parse.(Int, keys(opf_result["solution"]["branch"]))))
#     end

#     return pf, pf_mw, pfdc, pcurt, pd, pflex, pgmax, pg, pf_tot, pc_tot, bus_duals, branch_duals, bus_ids, branch_ids
# end

# pf, pf_mw, pfdc, pcurt, pd, pflex, pgmax, pg, pf_tot, pc_tot, bus_duals, branch_duals, bus_ids, branch_ids = run_mn_opf(opf_data, hours, demand_series, pv_series, wind_series; formulation="DC", verbose = false);


# ## to find out which bus is a candidate for branch connection, determine which bus has the most significant pos/neg dual variation
# sensitive_buses = []
# for (i, bus_id) in enumerate(bus_ids)
#     bus_daily_duals = [duals[i] for (h, duals) in bus_duals]
#     bus_max_abs = maximum(abs.(bus_daily_duals))
#     if !isempty(bus_daily_duals[bus_daily_duals.>0])
#         bus_max_pos_normalised = maximum(bus_daily_duals[bus_daily_duals.>0]) / bus_max_abs
#     else 
#         bus_max_pos_normalised = 0
#     end
#     if !isempty(bus_daily_duals[bus_daily_duals.<0])
#         bus_min_neg_normalised = minimum(bus_daily_duals[bus_daily_duals.<0]) / bus_max_abs
#     else
#         bus_min_neg_normalised = 0
#     end
#     bus_mean_variation = abs((bus_max_pos_normalised + bus_min_neg_normalised) / 2)

#     if abs(bus_max_pos_normalised) .> 0.07 && abs(bus_min_neg_normalised) .> 0.07
#         if opf_data["bus"]["$bus_id"]["base_kv"] > 132
#             # @show (bus_id, opf_data["bus"]["$bus_id"]["area"], bus_max_pos_normalised*bus_max_abs, bus_min_neg_normalised*bus_max_abs)
#             push!(sensitive_buses, (bus_id, opf_data["bus"]["$bus_id"]["base_kv"], opf_data["bus"]["$bus_id"]["area"], bus_max_pos_normalised*bus_max_abs, bus_min_neg_normalised*bus_max_abs))
#         end
#     end
# end
# @show sensitive_buses

# total_variation_sorted =  [x[3]/maximum([abs(x[3]),abs(x[4])]) + x[4]/maximum([abs(x[3]),abs(x[4])]) for x in sensitive_buses]
# sorted_buses = sort(sensitive_buses, by=x->x[3]/maximum([abs(x[3]),abs(x[4])]) + x[4]/maximum([abs(x[3]),abs(x[4])]))

##
""" 
Transmission line expansion cost: OHL: 1 m$/km for double circuit OHL 400 kV; lifetime: 60 years -> 1E6/60 per km per year -> 1E6/60/365 per km per day
"""

data_file_tnep = "snem2000_tnep.m"
tnep_data = prepare_data(data_file_tnep; merge_parallel_lines = false, hvdc = false)

for (i, gen) in tnep_data["gen"]
    if gen["type"] == "Fossil" || gen["type"] == "Hydro"
        gen["construction_cost"] = 1E8 * tnep_data["baseMVA"]
    elseif gen["fuel"] == "Solar"
        gen["construction_cost"] = 4E6 * tnep_data["baseMVA"]  # https://www.mcgqs.com.au/media/australian-solar-farms/
    elseif gen["fuel"] == "Wind"
        gen["construction_cost"] = 8E6 * tnep_data["baseMVA"]  # https://aemo.com.au/-/media/files/stakeholder_consultation/consultations/nem-consultations/2022/2023-inputs-assumptions-and-scenarios-consultation/supporting-materials-for-2023/aurecon-2022-cost-and-technical-parameter-review.pdf
    end
end

##
hours = 96

if hours > 1
    tnep_mn_data = _PM.replicate(tnep_data, hours)
    prepare_mn_opf_data!(tnep_mn_data, tnep_data, demand_series, pv_series, wind_series, collect(1:hours))
    pm = _PM.instantiate_model(tnep_mn_data, _PM.DCPPowerModel, build_mn_tnep; ref_extensions=[_PM.ref_add_on_off_va_bounds!, _PM.ref_add_ne_branch!]);
    result = _PM.optimize_model!(pm, relax_integrality=false, optimizer=dc_solver, solution_processors=[])
else
    ### test a single period for validation purposes
    i = hours
    tnep_i_data = _PM.replicate(tnep_data, 1)
    tnep_i_data["nw"]["1"] = deepcopy(tnep_mn_data["nw"]["$i"])
    tnep_i_data["nw"]["1"]["per_unit"] = true
    pm = _PM.instantiate_model(tnep_i_data, _PM.DCPPowerModel, build_mn_tnep; ref_extensions=[_PM.ref_add_on_off_va_bounds!, _PM.ref_add_ne_branch!]);
    result = _PM.optimize_model!(pm, relax_integrality=false, optimizer=dc_solver, solution_processors=[])
    Plots.plot([gen["pgslack"] for (i, gen) in result["solution"]["nw"]["1"]["gen"]])
    Plots.plot([load["pdslack"] for (i, load) in result["solution"]["nw"]["1"]["load"]])
end

## Plots
mkdir("./Scripts/PowerModels_TransmissionNetworkExpansionPlanning/Figures")

buses_NSW_xy = [(bus["x"], bus["y"]) for (i, bus) in opf_data["bus"] if bus["area"]==1 && haskey(bus, "x")]
buses_Vic_xy = [(bus["x"], bus["y"]) for (i, bus) in opf_data["bus"] if bus["area"]==2 && haskey(bus, "x")]
buses_QLD_xy = [(bus["x"], bus["y"]) for (i, bus) in opf_data["bus"] if bus["area"]==3 && haskey(bus, "x")]
buses_SA_xy  = [(bus["x"], bus["y"]) for (i, bus) in opf_data["bus"] if bus["area"]==4 && haskey(bus, "x")]
buses_TAS_xy = [(bus["x"], bus["y"]) for (i, bus) in opf_data["bus"] if bus["area"]==5 && haskey(bus, "x")]

sorted_bus_ids = [758, 93, 135, 10004, 807, 10005, 119, 1566, 844, 1557, 1569, 1570, 1742, 1745, 1764, 1523]
bus_iaxy = [(i, opf_data["bus"]["$i"]["area"], opf_data["bus"]["$i"]["base_kv"], opf_data["bus"]["$i"]["x"], opf_data["bus"]["$i"]["y"]) for i in sorted_bus_ids]
bus_ixy = [(i, opf_data["bus"]["$i"]["x"], opf_data["bus"]["$i"]["y"]) for i in sorted_bus_ids]
bus_xy = [(opf_data["bus"]["$i"]["x"], opf_data["bus"]["$i"]["y"]) for i in sorted_bus_ids]

# plotlyjs()
map_buses = Plots.scatter(buses_NSW_xy, color=1, label="NSW")
Plots.scatter!(buses_Vic_xy, color=2, label="VIC")
Plots.scatter!(buses_QLD_xy, color=3, label="QLD")
Plots.scatter!(buses_SA_xy, color=4, label="SA")
Plots.scatter!(buses_TAS_xy, color=5, label="TAS")
Plots.scatter!(bus_xy, color=:black, label="Candidate buses", xticks=false, yticks=false)
# Plots.savefig(map_buses, "Figures/map_buses.pdf")

for (i, ne_branch) in tnep_data["ne_branch"]
    ft_bus_x = [tnep_data["bus"]["$(ne_branch["f_bus"])"]["x"], tnep_data["bus"]["$(ne_branch["t_bus"])"]["x"]]
    ft_bus_y = [tnep_data["bus"]["$(ne_branch["f_bus"])"]["y"], tnep_data["bus"]["$(ne_branch["t_bus"])"]["y"]]
    if i == "1"
        Plots.plot!(ft_bus_x, ft_bus_y, color=:black, linewidth=2, label="Candidate lines")
    else
        Plots.plot!(ft_bus_x, ft_bus_y, color=:black, linewidth=2, label=false)
    end
end
map_buses
# Plots.savefig(map_buses, "Figures/map_buses_branches.pdf")


# plotlyjs()
map_buses_built_lines = Plots.scatter(buses_NSW_xy, color=1, label="NSW")
Plots.scatter!(buses_Vic_xy, color=2, label="VIC")
Plots.scatter!(buses_QLD_xy, color=3, label="QLD")
Plots.scatter!(buses_SA_xy, color=4, label="SA")
Plots.scatter!(buses_TAS_xy, color=5, label="TAS")
Plots.scatter!(bus_xy, color=:black, label="Candidate buses", xticks=false, yticks=false)
for (i, ne_branch) in tnep_data["ne_branch"]
    if result["solution"]["nw"]["1"]["ne_branch"]["$i"]["built"] == 1
        ft_bus_x = [tnep_data["bus"]["$(ne_branch["f_bus"])"]["x"], tnep_data["bus"]["$(ne_branch["t_bus"])"]["x"]]
        ft_bus_y = [tnep_data["bus"]["$(ne_branch["f_bus"])"]["y"], tnep_data["bus"]["$(ne_branch["t_bus"])"]["y"]]
        if i == "1"
            Plots.plot!(ft_bus_x, ft_bus_y, color=:black, linewidth=2, label="Built lines")
        else
            Plots.plot!(ft_bus_x, ft_bus_y, color=:black, linewidth=2, label=false)
        end
    end
end
map_buses_built_lines
# Plots.savefig(map_buses_built_lines, "Figures/map_buses_branches_built.pdf")


##
pv_average = Plots.plot(pv_series[!,"NSW"][1:48], label="NSW", linewidth=2)
Plots.plot!(pv_series[!,"VIC"][1:48], label="VIC", linewidth=2)
Plots.plot!(pv_series[!,"QLD"][1:48], label="QLD", linewidth=2)
Plots.plot!(pv_series[!,"SA"][1:48], label="SA", linewidth=2)
Plots.plot!(pv_series[!,"TAS"][1:48], label="TAS", linewidth=2)
xlabel!("Average day")
ylabel!("Ratio")
title!("PV Average Ratio per State")
# Plots.savefig(pv_average, "Figures/pv_series_average.pdf")


wind_average = Plots.plot(wind_series[!,"NSW"][1:48], label="NSW", linewidth=2)
Plots.plot!(wind_series[!,"VIC"][1:48], label="VIC", linewidth=2)
Plots.plot!(wind_series[!,"QLD"][1:48], label="QLD", linewidth=2)
Plots.plot!(wind_series[!,"SA"][1:48], label="SA", linewidth=2)
Plots.plot!(wind_series[!,"TAS"][1:48], label="TAS", linewidth=2, legend=:top)
xlabel!("Average day")
ylabel!("Ratio")
title!("Wind Average Ratio per State")
# Plots.savefig(wind_average, "Figures/wind_series_average.pdf")

demand_average = Plots.plot(demand_series[!,"NSW"][1:48], label="NSW", linewidth=2)
Plots.plot!(demand_series[!,"VIC"][1:48], label="VIC", linewidth=2)
Plots.plot!(demand_series[!,"QLD"][1:48], label="QLD", linewidth=2)
Plots.plot!(demand_series[!,"SA"][1:48], label="SA", linewidth=2)
Plots.plot!(demand_series[!,"TAS"][1:48], label="TAS", linewidth=2)
xlabel!("Average day")
ylabel!("Ratio")
title!("Demand Average Ratio per State")
# Plots.savefig(demand_average, "Figures/demand_series_average.pdf")

