dyn_inv(static_gen, Vrated, Irated; Vdc=600.0) = DynamicInverter(
           name = get_name(static_gen),
           ω_ref = 1.0, # ω_ref,
           converter = AverageConverter(rated_voltage = Vrated, rated_current = Irated),
           outer_control = OuterControl(
               VirtualInertia(Ta = 2.0, kd = 400.0, kω = 20.0),
               ReactivePowerDroop(kq = 0.2, ωf = 1000.0),
           ),
           inner_control = VoltageModeControl(
               kpv = 0.59,     #Voltage controller proportional gain
               kiv = 736.0,    #Voltage controller integral gain
               kffv = 0.0,     #Binary variable enabling the voltage feed-forward in output of current controllers
               rv = 0.0,       #Virtual resistance in pu
               lv = 0.2,       #Virtual inductance in pu
               kpc = 1.27,     #Current controller proportional gain
               kic = 14.3,     #Current controller integral gain
               kffi = 0.0,     #Binary variable enabling the current feed-forward in output of current controllers
               ωad = 50.0,     #Active damping low pass filter cut-off frequency
               kad = 0.2,
           ),
           dc_source = FixedDCSource(voltage = Vdc),
           freq_estimator = KauraPLL(
               ω_lp = 500.0, #Cut-off frequency for LowPass filter of PLL filter.
               kp_pll = 0.084,  #PLL proportional gain
               ki_pll = 4.69,   #PLL integral gain
           ),
           filter = LCLFilter(lf = 0.08, rf = 0.003, cf = 0.074, lg = 0.2, rg = 0.01),
       )
