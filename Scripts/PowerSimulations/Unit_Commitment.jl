# using Pkg
using PowerSystems
using PowerSimulations
using HydroPowerSimulations
using PowerSystemCaseBuilder
using PowerNetworkMatrices
using PowerGraphics
using HiGHS
import DataStructures: SortedDict
using DataFrames
using TimeSeries
# using Ipopt
# using SCIP
# using Juniper
# using PowerAnalytics
# using InfrastructureSystems
# using Dates

const PSY = PowerSystems
const PSI = PowerSimulations
const PSB = PowerSystemCaseBuilder
const PNM = PowerNetworkMatrices
# const PA = PowerAnalytics
# const IS = InfrastructureSystems


## build the system and add timeseries data
include("./Descripter.jl")
# include("./Timeseries_data.jl")
include("./Timeseries_data_week.jl")

file_path  = "snem2000.m"
sys_da = build_snem2000_bus_matpower_DA(file_path)
# PSY.transform_single_time_series!(sys_da, 24, Hour(1))
# PSY.transform_single_time_series!(sys_da, 48, Minute(30))

##
PTDF_matrix = PNM.PTDF(sys_da) # linear_solver="KLU"
template_uc = ProblemTemplate()
set_network_model!(template_uc, PSI.NetworkModel(PTDFPowerModel, PTDF_matrix=PTDF_matrix, duals=[CopperPlateBalanceConstraint], use_slacks=true))
set_device_model!(template_uc, Line, StaticBranch) # StaticBranchUnbounded
set_device_model!(template_uc, TapTransformer, StaticBranch)
set_device_model!(template_uc, Transformer2W, StaticBranch)
set_device_model!(template_uc, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_uc, ThermalStandard, ThermalStandardUnitCommitment)
set_device_model!(template_uc, HydroDispatch, HydroDispatchRunOfRiver)
set_device_model!(template_uc, PowerLoad, StaticPowerLoad)
for (k, v) in template_uc.branches
    v.duals = [NetworkFlowConstraint]
end

##
# solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
# solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>nl_solver)
# solver = optimizer_with_attributes(SCIP.Optimizer)
solver = optimizer_with_attributes(HiGHS.Optimizer)#, "mip_rel_gap" => 0.5)
problem = DecisionModel(template_uc, sys_da; optimizer = solver, horizon = Hour(24))
build!(problem, output_dir = mktempdir())
solve!(problem)

##
res = OptimizationProblemResults(problem)

plot_folder = "./plots"
mkpath(plot_folder)
plt_demand = plot_demand(res)
Plots.savefig(plt_demand, plot_folder*"/demand.pdf")

plt_fuel = plot_fuel(res)
Plots.savefig(plt_fuel, plot_folder*"/fuel.pdf")

duals = Dict([k => read_dual(res, k) for k in list_dual_keys(res) if PSI.get_entry_type(k) == NetworkFlowConstraint])
λ = read_dual(res, "CopperPlateBalanceConstraint__System")[:, 2]
# flow_duals = outerjoin(values(duals)..., on = :DateTime)
flow_duals = values(duals).dict[ PSI.ConstraintKey{NetworkFlowConstraint, Line}("") ]
flow_duals1 = values(duals).dict[ PSI.ConstraintKey{NetworkFlowConstraint, TapTransformer}("") ]
flow_duals2 = values(duals).dict[ PSI.ConstraintKey{NetworkFlowConstraint, Transformer2W}("") ]
flow_duals = innerjoin(flow_duals, flow_duals1, on = :DateTime)
flow_duals = innerjoin(flow_duals, flow_duals2, on = :DateTime)

μ = Matrix(flow_duals[:, PTDF_matrix.axes[2]])
LMP = flow_duals[:, [:DateTime]]
for bus in get_components(Bus, sys_da)
    bus_index = findall(x->x==get_number(bus), PTDF_matrix.axes[1])
    LMP[:, get_name(bus)] = λ .+ μ * PTDF_matrix[:, bus_index[1]]
end
@show LMP

##
timestamps = get_realized_timestamps(res)
variable_thermal_power = read_variable(res, "ActivePowerVariable__ThermalStandard")
variable_hydro_power = read_variable(res, "ActivePowerVariable__HydroDispatch")
variable_renewable_power = read_variable(res, "ActivePowerVariable__RenewableDispatch")
variable_line_flows = read_variable(res, "FlowActivePowerVariable__Line")
variable_thermal_on = read_variable(res, "OnVariable__ThermalStandard")
variable_thermal_start = read_variable(res, "StartVariable__ThermalStandard")
variable_thermal_stop = read_variable(res, "StopVariable__ThermalStandard")
parameter_load_power = read_parameter(res, "ActivePowerTimeSeriesParameter__PowerLoad")
parameter_renewable_power = read_parameter(res, "ActivePowerTimeSeriesParameter__RenewableDispatch")
prices = read_dual(res, "CopperPlateBalanceConstraint__System")

##
# plotlyjs() # loads the plotlyjs backend
plot_dataframe(LMP, timestamps)
plt_thermal_gens = plot_dataframe(variable_thermal_power, timestamps; stack = true)
plt_renewables = plot_dataframe(variable_renewable_power, timestamps; stack = true)
plt_hydro_power = plot_dataframe(variable_hydro_power, timestamps; stack = true)
plt_line_flows = plot_dataframe(variable_line_flows, timestamps; stack = false, legend=false)
plt_thermal_gens_on = plot_dataframe(variable_thermal_on, timestamps; stack = true)
plt_thermal_gens_start = plot_dataframe(variable_thermal_start, timestamps; stack = true)
plt_thermal_gens_stop = plot_dataframe(variable_thermal_stop, timestamps; stack = true)


Plots.savefig(plt_thermal_gens, plot_folder*"/plt_thermal_gens.pdf")
Plots.savefig(plt_renewables, plot_folder*"/plt_renewables.pdf")
Plots.savefig(plt_hydro_power, plot_folder*"/plt_hydro_power.pdf")
Plots.savefig(plt_line_flows, plot_folder*"/plt_line_flows.pdf")
Plots.savefig(plt_thermal_gens_on, plot_folder*"/plt_thermal_gens_on.pdf")
Plots.savefig(plt_thermal_gens_start, plot_folder*"/plt_thermal_gens_start.pdf")
Plots.savefig(plt_thermal_gens_stop, plot_folder*"/plt_thermal_gens_stop.pdf")
