using PowerModels
using JuMP
using Plots
using DataFrames
using HiGHS
const _PM = PowerModels

optimizer_highs = JuMP.optimizer_with_attributes(HiGHS.Optimizer)

## 
data_file = "snem1803_storage.m"
data = PowerModels.parse_file(data_file)
result_dcp_opf = _PM._solve_opf_strg_mi(data, DCPPowerModel, optimizer_highs, setting = Dict("output" => Dict("duals" => true)))
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
result_dc_storage_opf = _PM.solve_mn_opf_strg(mn_data, DCPPowerModel, optimizer_highs,  setting = Dict("output" => Dict("duals" => true)))
obj_snem1803_storage = result_dc_storage_opf["objective"]
time_snem1803_storage = result_dc_storage_opf["solve_time"]

# duals of kcl at buses are stored in solution -> nw is time step index -> bus index -> lam_kcl_r, 
duals_ts_1 = Dict(b=> bus["lam_kcl_r"] for (b,bus)  in result_dc_storage_opf["solution"]["nw"]["1"]["bus"])
# but they will be zeros in presence of binary variables. see https://jump.dev/JuMP.jl/stable/tutorials/linear/mip_duality/
##
pg_dc = Dict(gen["pg"] => id for (id, gen) in result_dc_storage_opf["solution"]["nw"]["1"]["gen"])
