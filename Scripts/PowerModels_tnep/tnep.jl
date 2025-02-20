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

ENV["GUROBI_HOME"]="/Library/gurobi1103/macos_universal2"
ENV["GRB_LICENSE_FILE"]="/Users/hei06j/gurobi/gurobi_11.lic"
############ END INPUT SECTION ##############################
#############################################################


#############################################################
#################### START METHODOLOGY ######################
### Optimisation settings for CbaOPF.jl
s = Dict("output" => Dict("branch_flows" => true, "duals" => true), "conv_losses_mp" => true)

# Test case data
data_path_hvdc = "snem2000_ACDC.m"

opf_data = prepare_data(data_path_hvdc; merge_parallel_lines=true)

##
# Select hours
# hours = select_hours(year, selection = selected_hours)
hours = selected_hours["hour_range"]

pf, pf_mw, pfdc, pcurt, pd, pflex, pgmax, pg, pf_tot, pc_tot, bus_duals, branch_duals, bus_ids, branch_ids = run_mn_opf(opf_data, hours, demand_series, pv_series, wind_series; formulation="DC", verbose = false);


## to find out which bus is a candidate for branch connection, determine which bus has the most significant pos/neg dual variation
sensitive_buses = []
for (i, bus_id) in enumerate(bus_ids)
    bus_daily_duals = [duals[i] for (h, duals) in bus_duals]
    bus_max_abs = maximum(abs.(bus_daily_duals))
    if !isempty(bus_daily_duals[bus_daily_duals.>0])
        bus_max_pos_normalised = maximum(bus_daily_duals[bus_daily_duals.>0]) / bus_max_abs
    else 
        bus_max_pos_normalised = 0
    end
    if !isempty(bus_daily_duals[bus_daily_duals.<0])
        bus_min_neg_normalised = minimum(bus_daily_duals[bus_daily_duals.<0]) / bus_max_abs
    else
        bus_min_neg_normalised = 0
    end
    bus_mean_variation = abs((bus_max_pos_normalised + bus_min_neg_normalised) / 2)

    if abs(bus_max_pos_normalised) .> 0.02 && abs(bus_min_neg_normalised) .> 0.02
        if opf_data["bus"]["$bus_id"]["base_kv"] > 132
            # @show (bus_id, opf_data["bus"]["$bus_id"]["area"], bus_max_pos_normalised*bus_max_abs, bus_min_neg_normalised*bus_max_abs)
            push!(sensitive_buses, (bus_id, opf_data["bus"]["$bus_id"]["base_kv"], opf_data["bus"]["$bus_id"]["area"], bus_max_pos_normalised*bus_max_abs, bus_min_neg_normalised*bus_max_abs))
        end
    end
end
@show sensitive_buses

total_variation_sorted =  [x[3]/maximum([abs(x[3]),abs(x[4])]) + x[4]/maximum([abs(x[3]),abs(x[4])]) for x in sensitive_buses]
sorted_buses = sort(sensitive_buses, by=x->x[3]/maximum([abs(x[3]),abs(x[4])]) + x[4]/maximum([abs(x[3]),abs(x[4])]))

##
# """ 
# Transmission line expansion cost: OHL: 1 m$/km for double circuit OHL 400 kV; lifetime: 60 years -> 1E6/60 per km per year -> 1E6/60/365 per km per day
# """

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
hours = 144

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
mkpath("./Scripts/PowerModels_tnep/Figures")

buses_NSW_xy = [(bus["x"], bus["y"]) for (i, bus) in tnep_data["bus"] if bus["area"]==1 && haskey(bus, "x")]
buses_Vic_xy = [(bus["x"], bus["y"]) for (i, bus) in tnep_data["bus"] if bus["area"]==2 && haskey(bus, "x")]
buses_QLD_xy = [(bus["x"], bus["y"]) for (i, bus) in tnep_data["bus"] if bus["area"]==3 && haskey(bus, "x")]
buses_SA_xy  = [(bus["x"], bus["y"]) for (i, bus) in tnep_data["bus"] if bus["area"]==4 && haskey(bus, "x")]
buses_TAS_xy = [(bus["x"], bus["y"]) for (i, bus) in tnep_data["bus"] if bus["area"]==5 && haskey(bus, "x")]


gen_Fossil_xy = [(tnep_data["bus"]["$(gen["gen_bus"])"]["x"], tnep_data["bus"]["$(gen["gen_bus"])"]["y"]) for (i,gen) in tnep_data["gen"] if gen["type"]=="Fossil" && haskey(tnep_data["bus"]["$(gen["gen_bus"])"], "x")]
gen_Wind_xy = [(tnep_data["bus"]["$(gen["gen_bus"])"]["x"], tnep_data["bus"]["$(gen["gen_bus"])"]["y"]) for (i,gen) in tnep_data["gen"] if gen["type"]=="Wind" && haskey(tnep_data["bus"]["$(gen["gen_bus"])"], "x")]
gen_Solar_xy = [(tnep_data["bus"]["$(gen["gen_bus"])"]["x"], tnep_data["bus"]["$(gen["gen_bus"])"]["y"]) for (i,gen) in tnep_data["gen"] if gen["type"]=="Solar" && haskey(tnep_data["bus"]["$(gen["gen_bus"])"], "x")]
gen_Hydro_xy = [(tnep_data["bus"]["$(gen["gen_bus"])"]["x"], tnep_data["bus"]["$(gen["gen_bus"])"]["y"]) for (i,gen) in tnep_data["gen"] if gen["type"]=="Hydro" && haskey(tnep_data["bus"]["$(gen["gen_bus"])"], "x")]

plotlyjs()
map_buses = Plots.scatter(buses_NSW_xy, color=:white, label="NSW")
Plots.scatter!(buses_Vic_xy, color=:white, label="VIC")
Plots.scatter!(buses_QLD_xy, color=:white, label="QLD")
Plots.scatter!(buses_SA_xy, color=:white, label="SA")
Plots.scatter!(buses_TAS_xy, color=:white, label="TAS")
# Plots.scatter!(gen_Fossil_xy, color=:black, label="Fossil")
Plots.scatter!(gen_Wind_xy, color=:green, label="Wind")
Plots.scatter!(gen_Solar_xy, color=:red, label="Solar")
Plots.scatter!(gen_Hydro_xy, color=:blue, label="Hydro")

sorted_bus_ids = [758, 93, 135, 10004, 807, 10005, 119, 1566, 844, 1557, 1569, 1570, 1742, 1745, 1764, 1523]
bus_iaxy = [(i, tnep_data["bus"]["$i"]["area"], tnep_data["bus"]["$i"]["base_kv"], tnep_data["bus"]["$i"]["x"], tnep_data["bus"]["$i"]["y"]) for i in sorted_bus_ids]
bus_ixy = [(i, tnep_data["bus"]["$i"]["x"], tnep_data["bus"]["$i"]["y"]) for i in sorted_bus_ids]
bus_xy = [(tnep_data["bus"]["$i"]["x"], tnep_data["bus"]["$i"]["y"]) for i in sorted_bus_ids]

# plotlyjs()
map_buses = Plots.scatter(buses_NSW_xy, color=1, label="NSW")
Plots.scatter!(buses_Vic_xy, color=2, label="VIC")
Plots.scatter!(buses_QLD_xy, color=3, label="QLD")
Plots.scatter!(buses_SA_xy, color=4, label="SA")
Plots.scatter!(buses_TAS_xy, color=5, label="TAS")
Plots.scatter!(bus_xy, color=:black, label="Candidate buses", xticks=false, yticks=false)
Plots.savefig(map_buses, "./Scripts/PowerModels_tnep/Figures/map_buses.pdf")

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
Plots.savefig(map_buses, "./Scripts/PowerModels_tnep/Figures/map_buses_branches.pdf")


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
Plots.savefig(map_buses_built_lines, "./Scripts/PowerModels_tnep/Figures/map_buses_branches_built.pdf")


[ne_branch["built"] for (n, ne_branch) in result["solution"]["nw"]["1"]["ne_branch"]]
ne1 = [nw["ne_branch"]["1"]["pf"] for (n, nw) in result["solution"]["nw"]] .* tnep_data["baseMVA"]
ne2 = [nw["ne_branch"]["2"]["pf"] for (n, nw) in result["solution"]["nw"]] .* tnep_data["baseMVA"]
ne3 = [nw["ne_branch"]["3"]["pf"] for (n, nw) in result["solution"]["nw"]] .* tnep_data["baseMVA"]
ne4 = [nw["ne_branch"]["4"]["pf"] for (n, nw) in result["solution"]["nw"]] .* tnep_data["baseMVA"]
ne5 = [nw["ne_branch"]["5"]["pf"] for (n, nw) in result["solution"]["nw"]] .* tnep_data["baseMVA"]
ne6 = [nw["ne_branch"]["6"]["pf"] for (n, nw) in result["solution"]["nw"]] .* tnep_data["baseMVA"]
ne7 = [nw["ne_branch"]["7"]["pf"] for (n, nw) in result["solution"]["nw"]] .* tnep_data["baseMVA"]
ne8 = [nw["ne_branch"]["8"]["pf"] for (n, nw) in result["solution"]["nw"]] .* tnep_data["baseMVA"]
ne9 = [nw["ne_branch"]["9"]["pf"] for (n, nw) in result["solution"]["nw"]] .* tnep_data["baseMVA"]
ne_plot = Plots.scatter(ne1, label="NSW-NSW 1")
Plots.scatter!(ne2, label="NSW-NSW 2")
Plots.scatter!(ne3, label="QLD-QLD 1", color=8)
Plots.scatter!(ne4, label="QLD-QLD 2")
Plots.scatter!(ne5, label="NSW-VIC 1")
Plots.scatter!(ne6, label="NSW-VIC 2")
Plots.scatter!(ne7, label="NSW-QLD")
Plots.scatter!(ne8, label="NSW-SA", color=3)
Plots.scatter!(ne9, label="VIC-SA")
xlabel!("Time (Hour)")
ylabel!("Power (MVA)")
title!("Power Transfer Across Candidate Lines")
Plots.savefig(ne_plot, "./Scripts/PowerModels_tnep/Figures/Network_expansion.pdf")


##
pv_average = Plots.plot(pv_series[!,"NSW"][1:48], label="NSW", linewidth=2)
Plots.plot!(pv_series[!,"VIC"][1:48], label="VIC", linewidth=2)
Plots.plot!(pv_series[!,"QLD"][1:48], label="QLD", linewidth=2)
Plots.plot!(pv_series[!,"SA"][1:48], label="SA", linewidth=2)
Plots.plot!(pv_series[!,"TAS"][1:48], label="TAS", linewidth=2)
xlabel!("Average day")
ylabel!("Ratio")
title!("PV Average Ratio per State")
# Plots.savefig(pv_average, "./Scripts/PowerModels_tnep/Figures/pv_series_average.pdf")


wind_average = Plots.plot(wind_series[!,"NSW"][1:48], label="NSW", linewidth=2)
Plots.plot!(wind_series[!,"VIC"][1:48], label="VIC", linewidth=2)
Plots.plot!(wind_series[!,"QLD"][1:48], label="QLD", linewidth=2)
Plots.plot!(wind_series[!,"SA"][1:48], label="SA", linewidth=2)
Plots.plot!(wind_series[!,"TAS"][1:48], label="TAS", linewidth=2, legend=:top)
xlabel!("Average day")
ylabel!("Ratio")
title!("Wind Average Ratio per State")
# Plots.savefig(wind_average, "./Scripts/PowerModels_tnep/Figures/wind_series_average.pdf")

demand_average = Plots.plot(demand_series[!,"NSW"][1:48], label="NSW", linewidth=2)
Plots.plot!(demand_series[!,"VIC"][1:48], label="VIC", linewidth=2)
Plots.plot!(demand_series[!,"QLD"][1:48], label="QLD", linewidth=2)
Plots.plot!(demand_series[!,"SA"][1:48], label="SA", linewidth=2)
Plots.plot!(demand_series[!,"TAS"][1:48], label="TAS", linewidth=2)
xlabel!("Average day")
ylabel!("Ratio")
title!("Demand Average Ratio per State")
# Plots.savefig(demand_average, "./Scripts/PowerModels_tnep/Figures/demand_series_average.pdf")

