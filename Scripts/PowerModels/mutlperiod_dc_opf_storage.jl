using PowerModels
using JuMP
using Ipopt
using Plots
using DataFrames
using HiGHS
const _PM = PowerModels

optimizer = JuMP.optimizer_with_attributes(HiGHS.Optimizer)

## 
data_file = "snem1803_storage.m"
data = PowerModels.parse_file(data_file)
result_dcp_opf = _PM._solve_opf_strg(data, DCPPowerModel, optimizer, setting = Dict("output" => Dict("duals" => true)))
obj_snem1803 = result_dcp_opf["objective"]
time_snem1803 = result_dcp_opf["solve_time"]

data["time_series"] = Dict(
    "num_steps" => 3,
    "load" => Dict(
        "1" => Dict("pd" => [3.0, 3.5, 4.0]),
        "2" => Dict("pd" => [3.0, 3.5, 4.0]),
        "3" => Dict("pd" => [4.0, 3.5, 3.0])
    )
)
mn_data = make_multinetwork(data)
result_dc_storage_opf = _PM.solve_mn_opf_strg(mn_data, DCPPowerModel, optimizer,  setting = Dict("output" => Dict("duals" => true)))
obj_snem1803_storage = result_dc_storage_opf["objective"]
time_snem1803_storage = result_dc_storage_opf["solve_time"]

##
pg_dc = Dict(gen["pg"] => id for (id, gen) in result_dcp_opf["solution"]["gen"])
qg_dc = Dict(gen["qg"] => id for (id, gen) in result_dcp_opf["solution"]["gen"])
