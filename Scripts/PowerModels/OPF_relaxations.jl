using PowerModels
using JuMP
using Ipopt
using Plots
using DataFrames
using SCS
const _PM = PowerModels

optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>10000)

## 
data_file = "snem1803.m"
result_nlp_opf = _PM.solve_ac_opf(data_file, optimizer)
result_soc_opf = _PM.solve_opf(data_file, SOCWRPowerModel, optimizer)
result_dcp_opf = _PM.solve_opf(data_file, DCPPowerModel, optimizer)
obj_snem1803 = [result_nlp_opf["objective"] result_soc_opf["objective"] result_dcp_opf["objective"]]
time_snem1803 = [result_nlp_opf["solve_time"] result_soc_opf["solve_time"] result_dcp_opf["solve_time"]]

## 
data_file = "snem197.m"
result_nlp_opf = _PM.solve_ac_opf(data_file, optimizer)
result_soc_opf = _PM.solve_opf(data_file, SOCWRPowerModel, optimizer)
result_dcp_opf = _PM.solve_opf(data_file, DCPPowerModel, optimizer)
obj_snem197 = [result_nlp_opf["objective"] result_soc_opf["objective"] result_dcp_opf["objective"]]
time_snem197 = [result_nlp_opf["solve_time"] result_soc_opf["solve_time"] result_dcp_opf["solve_time"]]


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