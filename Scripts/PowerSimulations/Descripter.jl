function build_snem2000_bus_matpower_RT(file_path; kwargs...)
    sys_kwargs = PSB.filter_kwargs(; kwargs...)
    # file_path = PSB.get_raw_data(; kwargs...)
    data_dir = dirname(dirname(file_path))
    pm_data = PSY.PowerModelsData(file_path)
    
    for (i, branch) in pm_data.data["branch"]
        if branch["br_x"] == 0
            branch["br_x"] = 1E-3
        end
    end
    
    ### generator types should be aligned with PowerSystems.jl enumerated typs
    gen_types = Dict(
            "Coal"=>"CT", 
            "Gas"=>"CC",
            "Solar"=>"PV",
            "Wind"=>"WT",
            "Hydro"=>"HY",
            # "Storage"=>"BA",
            )

    ### PowerSystems.jl does not recognize distillate, biomass, and capbanks
    ### Battery storage should be modelled as a PowerModels.jl storage, but here we assume batteries are solar to run this case
    ### CapBank/SVC/StatCom/SynCon generators will be replaced by gas, but with zero costs.
    for (i, gen) in pm_data.data["gen"]
        if gen["fuel"] == "Water"
            gen["fuel"] = "Hydro"
        elseif gen["fuel"] ∈ ["Natural Gas" "CapBank/SVC/StatCom/SynCon"]
            gen["fuel"] = "Gas"
        elseif gen["fuel"] ∈ ["Black Coal" "Brown Coal"]
            gen["fuel"] = "Coal"
        elseif gen["fuel"] == "Battery"
            gen["fuel"] = "Solar"
        end
        gen["type"] = gen_types[gen["fuel"]]
    end
    
    sys = PSY.System(pm_data; sys_kwargs...)
    # sys = System(pm_data; sys_kwargs...)


    if get(kwargs, :add_forecasts, true)
        # total_active_power_PowerLoad = sum(abs.([load.active_power for (i, load) in sys.data.components.data[PowerLoad]]))
        for (ix, l) in enumerate(PSY.get_components(PSY.PowerLoad, sys))
            forecast_data = SortedDict{Dates.DateTime, TimeSeries.TimeArray}()
            for t in 1:2 # loop over days
                mult = 1 # abs(l.active_power)/total_active_power_PowerLoad
                ta = load_timeseries_DA[t][1]
                values(ta[1]) .*= mult
                for i in 1:length(ta) # loop over hours
                    ini_time = timestamp(ta[i]) # get the hour
                    tra = load_timeseries_RT[t] # get the time series
                    values(tra[1]) .*= mult       # multiply the time series values with mult
                    data = when(tra[1], hour, hour(ini_time[1])) # get the subset ts for that hour
                    forecast_data[ini_time[1]] = data
                end
            end
            PSY.add_time_series!(
                sys,
                l,
                PSY.Deterministic("max_active_power", forecast_data),
            )
        end
        
        for (ix, l) in enumerate(PSY.get_components(PSY.RenewableGen, sys))
            forecast_data = SortedDict{Dates.DateTime, TimeSeries.TimeArray}()
            for t in 1:2
                if l.prime_mover_type.value == 21      # Solar
                    ta = solar_timeseries_DA[t][1]
                elseif l.prime_mover_type.value == 22  # Wind
                    ta = wind_timeseries_DA[t][1]
                end
                for i in 1:length(ta)
                    ini_time = timestamp(ta[i])
                    if l.prime_mover_type.value == 21      # Solar
                        data = when(solar_timeseries_RT[t][1], hour, hour(ini_time[1]))
                    elseif l.prime_mover_type.value == 22  # Wind
                        data = when(wind_timeseries_RT[t][1], hour, hour(ini_time[1]))
                    end
                    forecast_data[ini_time[1]] = data
                end
            end
            PSY.add_time_series!(
                sys,
                l,
                PSY.Deterministic("max_active_power", forecast_data),
            )
        end

        for (ix, l) in enumerate(PSY.get_components(PSY.HydroDispatch, sys))
            forecast_data = SortedDict{Dates.DateTime, TimeSeries.TimeArray}()
            for t in 1:2
                ta = hydro_timeseries_DA[t][1]
                for i in 1:length(ta)
                    ini_time = timestamp(ta[i])
                    data = when(hydro_timeseries_RT[t][1], hour, hour(ini_time[1]))
                    forecast_data[ini_time[1]] = data
                end
            end
            PSY.add_time_series!(
                sys,
                l,
                PSY.Deterministic("max_active_power", forecast_data),
            )
        end

    end

    return sys
end



function build_snem2000_bus_matpower_DA(file_path; kwargs...)
    sys_kwargs = PSB.filter_kwargs(; kwargs...)
    # file_path = PSB.get_raw_data(; kwargs...)
    data_dir = dirname(dirname(file_path))
    pm_data = PSY.PowerModelsData(file_path)
    
    for (i, gen) in pm_data.data["gen"]
        gen["pmax"] *= 1
        gen["mbase"] *= 1
        gen["qmin"] *= 1
        gen["qmax"] *= 1
    end

    for (i, branch) in pm_data.data["branch"]
        if branch["br_x"] == 0
            branch["br_x"] = 1E-3
        end
    end
    
    ### generator types should be aligned with PowerSystems.jl enumerated typs
    gen_types = Dict(
            "Coal"=>"CT",
            "Gas"=>"CC",
            "Solar"=>"PV",
            "Wind"=>"WT",
            "Hydro"=>"HY",
            # "Storage"=>"BA",
            )

    ### PowerSystems.jl does not recognize distillate, biomass, and capbanks
    ### Battery storage should be modelled as a PowerModels.jl storage, but here we assume batteries are solar to run this case
    ### CapBank/SVC/StatCom/SynCon generators will be replaced by gas, but with zero costs.
    for (i, gen) in pm_data.data["gen"]
        if gen["fuel"] == "Water"
            gen["fuel"] = "Hydro"
        elseif gen["fuel"] ∈ ["Natural Gas" "CapBank/SVC/StatCom/SynCon"]
            gen["fuel"] = "Gas"
        elseif gen["fuel"] ∈ ["Black Coal" "Brown Coal"]
            gen["fuel"] = "Coal"
        elseif gen["fuel"] == "Battery"
            gen["fuel"] = "Solar"
        end
        gen["type"] = gen_types[gen["fuel"]]
    end

    sys = PSY.System(pm_data; sys_kwargs...)

    if get(kwargs, :add_forecasts, true)
        total_active_power_PowerLoad = sum(abs.([load.active_power for (i, load) in sys.data.components.data[PowerLoad]]))
        for (ix, l) in enumerate(PSY.get_components(PSY.PowerLoad, sys))
            forecast_data = SortedDict{Dates.DateTime, TimeSeries.TimeArray}()
            for t in 1:2 # loop over days
                mult = 1 #abs(l.active_power)/total_active_power_PowerLoad
                ini_time = timestamp(load_timeseries_DA[t][1])[1]
                forecast_data[ini_time] = load_timeseries_DA[t][1] .* mult
            end
            PSY.add_time_series!(
                sys,
                l,
                PSY.Deterministic("max_active_power", forecast_data),
            )
        end

        for (ix, l) in enumerate(PSY.get_components(PSY.RenewableGen, sys))
            forecast_data = SortedDict{Dates.DateTime, TimeSeries.TimeArray}()
            for t in 1:2
                if l.prime_mover_type.value == 21 
                    ini_time = timestamp(solar_timeseries_DA[t][1])[1]
                    forecast_data[ini_time] = solar_timeseries_DA[t][1]
                    # ta = solar_timeseries_DA[t][1]
                elseif l.prime_mover_type.value == 22  # Wind
                    ini_time = timestamp(wind_timeseries_DA[t][1])[1]
                    forecast_data[ini_time] = wind_timeseries_DA[t][1]
                    # ta = wind_timeseries_DA[t][1]
                end
            end
            PSY.add_time_series!(
                sys,
                l,
                PSY.Deterministic("max_active_power", forecast_data),
            )
        end

        for (ix, l) in enumerate(PSY.get_components(PSY.HydroDispatch, sys))
            forecast_data = SortedDict{Dates.DateTime, TimeSeries.TimeArray}()
            for t in 1:2
                ini_time = timestamp(hydro_timeseries_DA[t][1])[1]
                forecast_data[ini_time] = hydro_timeseries_DA[t][1]
            end
            PSY.add_time_series!(
                sys,
                l,
                PSY.Deterministic("max_active_power", forecast_data),
            )
        end

    end

    return sys
end



const NEM_SYSTEM_CATALOG = [
    SystemDescriptor(
        name = "snem2000_bus_matpower_DA",
        description = "matpower 2000-bus synthetic NEM system with DA time series",
        category = PSITestSystems,
        raw_data = "snem2000.m",
        build_function = build_snem2000_bus_matpower_DA,
    )
    SystemDescriptor(
        name = "snem2000_bus_matpower_RT",
        description = "matpower 2000-bus synthetic NEM system with RT time series",
        category = PSITestSystems,
        raw_data = "snem2000.m",
        build_function = build_snem2000_bus_matpower_RT,
    )
]