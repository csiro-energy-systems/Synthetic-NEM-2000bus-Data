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
data_file = "snem2000_acdc.m"
s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true)

result_nlp_opf = _PMACDC.run_acdcopf(data_file, ACPPowerModel, optimizer; setting=s)
result_soc_opf = _PMACDC.run_acdcopf(data_file, SOCWRPowerModel, optimizer; setting=s)
result_dcp_opf = _PMACDC.run_acdcopf(data_file, DCPPowerModel, optimizer; setting=s )
obj_snem2000_acdc = [result_nlp_opf["objective"] result_soc_opf["objective"] result_dcp_opf["objective"]]
time_snem2000_acdc = [result_nlp_opf["solve_time"] result_soc_opf["solve_time"] result_dcp_opf["solve_time"]]
