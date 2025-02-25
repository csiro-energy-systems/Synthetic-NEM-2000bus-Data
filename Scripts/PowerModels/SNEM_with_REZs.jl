using Pkg
Pkg.activate("./Scripts/PowerModels")
using PowerModels
using JuMP
using Ipopt
# Pkg.add("PlotlyJS")
# Pkg.build("PlotlyJS")
using PlotlyJS

const _PM = PowerModels

optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>10000)

## 
data_file = "snem1803_rez_2022 ISP Step Change_2032.m"
data = _PM.parse_file(data_file)

## 


bus_coords = [(bus["x"], bus["y"]) for (i,bus) in data["bus"]]
load_coords = [(data["bus"]["$(load["load_bus"])"]["x"], data["bus"]["$(load["load_bus"])"]["y"]) for (i,load) in data["load"]]
gen_coords = [(data["bus"]["$(gen["gen_bus"])"]["x"], data["bus"]["$(gen["gen_bus"])"]["y"]) for (i,gen) in data["gen"]]
rez_coords =  [(data["bus"]["$(gen["gen_bus"])"]["x"], data["bus"]["$(gen["gen_bus"])"]["y"]) for (i,gen) in data["gen"] if !startswith(gen["name"], "gen")]

buses_without_coords =  [i for (i,bus) in data["bus"] if isnan(bus["x"])]

plot(
    [scatter(
        mode="markers",
        x=first.(bus_coords),
        y=last.(bus_coords),
        name="bus"
    ),
    scatter(
        mode="markers",
        x=first.(load_coords),
        y=last.(load_coords),
        name="load"
    ),
    # scatter(
    #     mode="markers",
    #     x=first.(gen_coords),
    #     y=last.(gen_coords),
    #     label="gen"
    # ),
    scatter(
        mode="markers",
        x=first.(rez_coords),
        y=last.(rez_coords),
        color="black",
        name="rez",
    )],
    config=PlotConfig(scrollZoom=true),  
    Layout(
        width=1000, height=800,
        margin=attr(l=20,r=20,t=20,b=20),
        paper_bgcolor="LightSteelBlue"
    )
)