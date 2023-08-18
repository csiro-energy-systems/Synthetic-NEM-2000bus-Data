using PowerModelsACDC
using PowerModels
using JuMP
using Ipopt
using Plots
using DataFrames
using SCS
const _PMACDC = PowerModelsACDC

optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>10000)

##
data_file = "snem2000.m"
snem2000 = parse_file(data_file)
[i for (i,bus) in snem2000["bus"] if bus["name"]=="bus_1304"]
[i for (i,branch) in snem2000["branch"] if branch["name"]=="line_2202_to_5001_1_cp"]
[i for (i,gen) in snem2000["gen"] if gen["name"]=="gen_1218_1"]


## 
data_file = "snem2000_acdc.m"
result_nlp_opf = _PMACDC.run_acdcopf(data_file, ACPPowerModel, optimizer)
result_soc_opf = _PMACDC.run_acdcopf(data_file, SOCWRPowerModel, optimizer)
result_dcp_opf = _PMACDC.run_acdcopf(data_file, DCPPowerModel, optimizer)
obj_snem2000_acdc = [result_nlp_opf["objective"] result_soc_opf["objective"] result_dcp_opf["objective"]]
time_snem2000_acdc = [result_nlp_opf["solve_time"] result_soc_opf["solve_time"] result_dcp_opf["solve_time"]]

##
pg_nlp = Dict(gen["pg"] => id for (id, gen) in result_nlp_opf["solution"]["gen"])
qg_nlp = Dict(gen["qg"] => id for (id, gen) in result_nlp_opf["solution"]["gen"])
pg_soc = Dict(gen["pg"] => id for (id, gen) in result_soc_opf["solution"]["gen"])
qg_soc = Dict(gen["qg"] => id for (id, gen) in result_soc_opf["solution"]["gen"])
pg_dc = Dict(gen["pg"] => id for (id, gen) in result_dcp_opf["solution"]["gen"])
qg_dc = Dict(gen["qg"] => id for (id, gen) in result_dcp_opf["solution"]["gen"])

pg_plot = plot([p for (p, i) in pg_nlp], [p for (p, i) in pg_nlp], seriestype=:scatter, linestyle=:dash, label="ACP", xlabel="ACP", ylabel="SOC/DC", title="Comparison of Generator Active Power Setpoints (MW)", titlefotnsize=12)
plot!([p for (p, i) in pg_nlp], [p for (p, i) in pg_soc], seriestype=:scatter, label="SOC")
plot!([p for (p, i) in pg_nlp], [p for (p, i) in pg_dc], seriestype=:scatter, label="DC")

qg_plot = plot([q for (q, i) in qg_nlp], [q for (q, i) in qg_nlp], linestyle=:dash, label=false, xlabel="ACP", ylabel="SOC/DC", title="Comparison of Generator Reactive Power Setpoints (MVAR)", titlefotnsize=12)
plot!([q for (q, i) in qg_nlp], [q for (q, i) in qg_soc], seriestype=:scatter, label="SOC")
plot!([q for (q, i) in qg_nlp], [q for (q, i) in qg_dc], seriestype=:scatter, label="DC")