using Pkg
Pkg.activate("./")
using PowerSystems
using PowerSimulations
using PowerSimulationsDynamics
using PowerSystemCaseBuilder
using PowerNetworkMatrices
using PowerFlows
using PowerGraphics
using PowerAnalytics
using HiGHS
using Ipopt
using Sundials
import DataStructures: SortedDict
using Dates
using DataFrames
using TimeSeries
using OrdinaryDiffEq
using Plots
using XLSX


const PSY = PowerSystems
# const PSI = PowerSimulations
const PSD = PowerSimulationsDynamics
const PSB = PowerSystemCaseBuilder
# const PNM = PowerNetworkMatrices
# const PA = PowerAnalytics


function get_dataframe(data_folder, filename)
    file = data_folder * "/" * filename
    sheetnames = XLSX.sheetnames(XLSX.readxlsx(file))
    df = DataFrames.DataFrame(XLSX.readtable(file, sheetnames[1]))
    colnames_01 = [df[1,k] for k=1:length(df[1,:])]
    colnames = [split(name,"\r")[1] for name in colnames_01]
    colnames = [split(name,"_x000D")[1] for name in colnames]
    colnames = [split(name,"\n")[1] for name in colnames]
    DataFrames.rename!(df, colnames)

    return df[3:end, :]
end

data_folder = joinpath(pwd(), "Data/original")
data_dict = Dict{String,Any}()
for filename in ["Gen", "Gen_Vref", "Gen_Exciters", "Gen_Governors_HYGOV", "Gen_Governors_TGOV", "Gen_Stabilizers"]
    data_dict[filename] = get_dataframe(data_folder, filename*".xlsx")
end

##
# include("./Descripter_1803.jl")
# # include("./Timeseries_data.jl")
# # # sys_rt = build_snem1803_bus_matpower_RT()
# sys_da = build_snem1803_bus_matpower_DA()

# PSY.transform_single_time_series!(sys_rt, 12, Minute(5))
# PSY.transform_single_time_series!(sys_da, 24, Hour(1))


# file_path  = "/Users/hei06j/Documents/repositories/remote/NEM2000synthetic/data/matpower/snem1803.m"
# file_path  = joinpath(pwd(), "snem197.m") # not stable
file_path  = joinpath(pwd(), "snemNSW.m") # stable - converging
# file_path  = joinpath(pwd(), "snemVIC.m") # not stable
# file_path  = joinpath(pwd(), "snemQLD.m") # stable - not converging
# file_path  = joinpath(pwd(), "snemSA.m")  # not stable

# sys_kwargs = PSB.filter_kwargs(; kwargs...)
# file_path = PSB.get_raw_data(; kwargs...)
# data_dir = dirname(dirname(file_path))
pm_data = PSY.PowerModelsData(file_path)

for (i, branch) in pm_data.data["branch"]
    if branch["br_x"] == 0
        branch["br_x"] = 1E-3
    end
end


hydro_gens = [name[7:end] for name in data_dict["Gen_Governors_HYGOV"][!,"Machine Controls"]]
thermal_gens = [name[7:end] for name in data_dict["Gen_Governors_TGOV"][!,"Machine Controls"]]

for (i, gen) in pm_data.data["gen"]
    if gen["name"] ∈ hydro_gens
        gen["fuel"] = "Hydro"
        gen["type"] = "HY"
    elseif gen["name"] ∈ thermal_gens
        gen["fuel"] = "Gas"
        gen["type"] = "CC"
    else
        # gen["fuel"] = "Solar"
        # gen["type"] = "PV"
        gen["fuel"] = "Gas"
        gen["type"] = "CC"
    end
end
sys = PSY.System(pm_data)

for l in get_components(PSY.StandardLoad, sys)
    transform_load_to_constant_impedance(l)
end

### Read and assing dynamic parameters 
parse_float(x) = parse.(Float64, x)

include("dynamic_gen_components.jl")
include("dynamic_inverter_components.jl")


gen_list = [gen.name for gen in collect(get_components(Generator, sys))]

for (i, gen) in pm_data.data["gen"]
    static_gen = get_component(Generator, sys, gen["name"])

    Gen_id = findall(x->x==gen["name"], data_dict["Gen"][!,"Machine"])
    if ~isempty(Gen_id)
        Gen_id = Gen_id[1]
        Gen_data = data_dict["Gen"][Gen_id,:]

        Exciter_id = findall(x->x=="IEEET1_"*gen["name"], data_dict["Gen_Exciters"][!,"Machine Controls"])[1]
        Exciter_data = data_dict["Gen_Exciters"][Exciter_id,:]
        # vref = data_dict["Gen_Vref"][Exciter_id, "A"]
        vref = string(pm_data.data["bus"]["$(gen["gen_bus"])"]["vm"])
        
        PSS_id = findall(x->x=="PSS2B_"*gen["name"], data_dict["Gen_Stabilizers"][!,"Machine Controls"])[1]
        PSS_data = data_dict["Gen_Stabilizers"][PSS_id,:]

        ### if the generator is a synchronous generator
        HYGOV_id = findall(x->x=="HYGOV_"*gen["name"], data_dict["Gen_Governors_HYGOV"][!,"Machine Controls"])
        TGOV1_id = findall(x->x=="TGOV1_"*gen["name"], data_dict["Gen_Governors_TGOV"][!,"Machine Controls"])

        if ~isempty(HYGOV_id)
            HYGOV_id = HYGOV_id[1]
            HYGOV_data = data_dict["Gen_Governors_HYGOV"][HYGOV_id,:]
            # dyn_gen = dyn_gen_hygov(static_gen, Gen_data, Exciter_data, vref, HYGOV_data, PSS_data)
            # dyn_gen = dyn_gen_simple(static_gen)
            # dyn_gen = dyn_gen_simple_v1(static_gen, Gen_data)
            dyn_gen = dyn_gen_simple_hygov(static_gen, Gen_data, HYGOV_data, Exciter_data, vref, PSS_data)

        elseif ~isempty(TGOV1_id)
            # @show gen["name"]
            TGOV1_id = TGOV1_id[1]
            TGOV1_data = data_dict["Gen_Governors_TGOV"][TGOV1_id,:]
            dyn_gen = dyn_gen_tgov1(static_gen, Gen_data, Exciter_data, vref, TGOV1_data, PSS_data)
            # dyn_gen = dyn_gen_simple(static_gen)
            # dyn_gen = dyn_gen_simple_v1(static_gen, Gen_data)
            # dyn_gen = dyn_gen_simple_tgov1(static_gen, Gen_data, TGOV1_data, Exciter_data, vref, PSS_data)

        end

        add_component!(sys, dyn_gen, static_gen)

    else
        # @show gen["gen_bus"]
        
        bus_name = pm_data.data["bus"]["$(gen["gen_bus"])"]["name"]
        bus = sys.data.components.data[ACBus][bus_name]
        sta_src = statis_source(gen, bus)

        remove_component!(sys, static_gen)
        add_component!(sys, sta_src)

        # dyn_gen = dyn_gen_simple(static_gen)
        # add_component!(sys, dyn_gen, static_gen)
        
        # gen["pmax"]*pm_data.data["baseMVA"]
        # dyn_inverter = dyn_inv(static_gen, bus.base_voltage, gen["pmax"]*pm_data.data["baseMVA"]/bus.base_voltage; Vdc=2*bus.base_voltage)
        # add_component!(sys, dyn_inverter, static_gen)
    end
    
end

time_span = (0.0, 1.0)
sim = PSD.Simulation(ResidualModel, sys, pwd(), time_span)
res_small_signal_before_fault = small_signal_analysis(sim)
plot(res_small_signal_before_fault.eigenvalues, seriestype=:scatter, xlim=(-1,1))

# """ 
# A good sanity check it running a power flow on the system to make sure all the components are properly scaled
# and that the system is properly balanced. We can use PowerSystems to perform this check. 
# We can get the results back and perform a sanity check.
# """

# res_pf = solve_powerflow(PowerFlows.ACPowerFlow(), sys)
# res_pf["flow_results"]
# res_pf["bus_results"]
# plot(res_pf["bus_results"][!,"Vm"], seriestype=:scatter)

##

time_span = (0.0, 20.0)

case_gen = collect(PSY.get_components(PSY.DynamicGenerator, sys))[1]
gen_trip = GeneratorTrip(1.5, case_gen)
[g for (g,gen) in sys.data.components.data[HydroDispatch]]

sim_rodas = PSD.Simulation(MassMatrixModel, sys, pwd(), time_span, gen_trip)
# sim_rodas = PSD.Simulation(MassMatrixModel, sys, pwd(), time_span)
PSD.execute!(sim_rodas, Rodas4(), saveat = 0.00001, atol = 1e-10, rtol = 1e-10)#, initializealg = NoInit())
res_rodas = read_results(sim_rodas)

res_small_signal_after_fault = small_signal_analysis(sim_rodas)
plot!(res_small_signal_after_fault.eigenvalues, seriestype=:scatter, xlim=(-1,0))


plt = plot()
for bus_id in [id for (id,bus) in res_rodas.bus_lookup][1:20]
    v_rodas = get_voltage_magnitude_series(res_rodas, bus_id)
    plot!(v_rodas, label=false)
end
display(plt)
title!("Bus voltage magnitues (V pu)")
display(plt)
savefig(plt, "./scripts/PowerSimulations/dynamics_fault_bus_voltage.pdf")


plt = plot()
c = 0
for (i, line) in sys.data.components.data[Line]
    freq = get_real_current_series(res_rodas, line.name)
    plot!(freq, label=false)
end
display(plt)

plt = plot()
for (i, gen) in sys.data.components.data[DynamicGenerator{RoundRotorQuadratic, SingleMass, AVRTypeI, SteamTurbineGov1, PSS2B}]
    freq = get_frequency_series(res_rodas, gen.name)
    plot!(freq.*50, label=false)
end
display(plt)

plt = plot()
for (i, gen) in sys.data.components.data[ThermalStandard]
    freq = get_frequency_series(res_rodas, gen.name)
    plot!(freq.*50, label=false, xticks=(1:200:2001, 0:2:20))
end
title!("Thermal generators frequency (Hz)")
display(plt)
savefig(plt, "./scripts/PowerSimulations/dynamics_fault_freq_thermal.pdf")


#Obtain data for angles
t, θ_oc = PSD.get_state_series(res_rodas, ("gen_1037_4", :θ_oc))
length(t) # == length(θ_oc)
t, ω_oc = PSD.get_state_series(res_rodas, ("generator-103-1", :ω_oc))
length(t) # == length(ω_oc)

##
time_span = (0.0, 20.0)
line_tripped = [l for (l,line) in sys.data.components.data[Line]][161]
perturbation_trip = BranchTrip(1.0, Line, line_tripped)
sim_rodas = PSD.Simulation(MassMatrixModel, sys, pwd(), time_span, perturbation_trip)
# sim_rodas = PSD.Simulation(MassMatrixModel, sys, pwd(), time_span)
PSD.execute!(sim_rodas, Rodas4(), saveat = 0.01, atol = 1e-10, rtol = 1e-10, initializealg = NoInit())
res_rodas = read_results(sim_rodas)
bus_id = [id for (id,bus) in res_rodas.bus_lookup][10]
v_rodas = get_voltage_magnitude_series(res_rodas, bus_id);
plot(v_rodas)

##
time_span = (0.0, 20.0)
line_tripped = [l for (l,line) in sys.data.components.data[Line]][162]
perturbation_trip = BranchTrip(1.0, Line, line_tripped)
sim = PSD.Simulation(ResidualModel, sys, pwd(), time_span, perturbation_trip)
# sim = PSD.Simulation(ResidualModel, sys, pwd(), time_span)
PSD.execute!(sim, IDA(), dtmax = 0.01)
res_ida = read_results(sim)
bus_id = [id for (id,bus) in res_ida.bus_lookup][10]
v_ida = get_voltage_magnitude_series(res_ida, bus_id);
plot!(v_ida)



## small signal analysis
res_small_signal = small_signal_analysis(sim_rodas)
plot(res_small_signal.eigenvalues, seriestype=:scatter, xlim=(-1,1))


## frequency plots
res = res_rodas

plt = plot()
for (i, gen) in sys.data.components.data[DynamicGenerator{RoundRotorQuadratic, SingleMass, AVRTypeI, SteamTurbineGov1, PSS2B}]
    freq = get_frequency_series(res, gen.name)
    plot!(freq)
end
display(plt)

# for (i, gen) in sys.data.components.data[DynamicGenerator{SalientPoleQuadratic, SingleMass, AVRTypeI, HydroTurbineGov, PSS2B}]
#     freq = get_frequency_series(res, gen.name)
#     plot!(freq)
# end
for (i, gen) in sys.data.components.data[DynamicGenerator{SalientPoleQuadratic, SingleMass, AVRFixed, HydroTurbineGov, PSS2B}]
    freq = get_frequency_series(res, gen.name)
    plot!(freq)
end
display(plt)


##
p = plot();
for b in get_components(Bus, sys)
    voltage_series = get_voltage_magnitude_series(result, get_number(b))
    plot!(
        p,
        voltage_series;
        xlabel = "Time",
        ylabel = "Voltage Magnitude [pu]",
        label = "Bus - $(get_name(b))",
    );
end


p2 = plot();
for g in get_components(ThermalStandard, sys)
    state_series = get_state_series(result, (get_name(g), :ω))
    plot!(
        p2,
        state_series;
        xlabel = "Time",
        ylabel = "Speed [pu]",
        label = "$(get_name(g)) - ω",
    );
end


# freq_gen_5040_4 = get_frequency_series(res, "gen_5040_4")
# plot(freq_gen_5040_4)

##
# sim = PSD.Simulation(
#     ResidualModel, #Type of model used
#     threebus_sys, #system
#     pwd(), #folder to output results
#     tspan, #time span
#     Ybus_change, #Type of perturbation
# )

###
# # Will print the initial states. It also give the symbols used to describe those states.
# # show_states_initial_value(sim)

# #Will export a dictionary with the initial condition values to explore
# # x0_init = PSD.get_initial_conditions(sim)

# # # Step 4: Run the simulation of the Static Lines System
# # Run the simulation
# PSD.execute!(sim, IDA(), dtmax = 0.01)




##
time_span = (0.0, 9.0)
perturbation_trip = BranchTrip(1.0, Line, "bus_5095-bus_5111-i_line_5095_to_5111_1_cp")
sim_rodas = PSD.Simulation(MassMatrixModel, sys, pwd(), time_span, perturbation_trip)

execute!(sim_rodas, Rodas4(), saveat = 0.01, atol = 1e-10, rtol = 1e-10, initializealg = NoInit())

##

using Logging
file_dir = joinpath(pkgdir(NEM2000synthetic), "scripts", "PowerSimulations", "data")
sys1 = System(joinpath(file_dir, "14bus.raw"), joinpath(file_dir, "dyn_data.dyr"))
sim1 = PSD.Simulation(ResidualModel, sys1, file_dir, (0.0, 20.0),
    # BranchTrip(1.0, Line, "BUS 02-BUS 04-i_4");
    console_level = Logging.Info,
)
show_states_initial_value(sim1)
x0_init = PSD.get_initial_conditions(sim1)
PSD.execute!(sim1, IDA(); abstol = 1e-8)
