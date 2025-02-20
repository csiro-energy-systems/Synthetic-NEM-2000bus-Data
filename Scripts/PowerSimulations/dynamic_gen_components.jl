machine_AFM(data) = AndersonFouadMachine(
            R = parse_float(data["Rs"]),
            Xd = parse_float(data["Xd"]),
            Xq = parse_float(data["Xq"]),
            Xd_p = parse_float(data["Xpd"]),
            Xq_p = parse_float(data["Xpq"]),
            Xd_pp = parse_float(data["Xppd"]),
            Xq_pp = parse_float(data["Xppq"]),
            Td0_p = parse_float(data["Tpd0"]),
            Tq0_p = parse_float(data["Tpq0"]),
            Td0_pp = parse_float(data["Tppd0"]),
            Tq0_pp = parse_float(data["Tppq0"])
        )
machine_RRM(data) = RoundRotorQuadratic(
            R = parse_float(data["Rs"]),
            Td0_p = parse_float(data["Tpd0"]),
            Td0_pp = parse_float(data["Tppd0"]),
            Tq0_p = parse_float(data["Tpq0"]),
            Tq0_pp = parse_float(data["Tppq0"]),
            Xd = parse_float(data["Xd"]),
            Xq = parse_float(data["Xq"]),
            Xd_p = parse_float(data["Xpd"]),
            Xq_p = parse_float(data["Xpq"]),
            Xd_pp = parse_float(data["Xppd"]),
            Xl = parse_float(data["Xls"]),
            # Se = (0.0, 0.0)
            Se = (parse_float.(split(data["Im_Points"][2:end-1], " ")[parse(Int, data["NdSat"])-1]), parse_float.(split(data["Im_Points"][2:end-1], " ")[parse(Int, data["NdSat"])])),
            # γ_d1 = (parse_float(data["Xppd"]) - parse_float(data["Xls"])) / (parse_float(data["Xpd"]) - parse_float(data["Xls"])),
            # γ_q1 = (parse_float(data["Xppq"]) - parse_float(data["Xls"])) / (parse_float(data["Xpq"]) - parse_float(data["Xls"])),
            # γ_d2 = (parse_float(data["Xpd"]) - parse_float(data["Xppd"])) / (parse_float(data["Xpd"]) - parse_float(data["Xls"]))^2,
            # γ_q2 = (parse_float(data["Xpq"]) - parse_float(data["Xppq"])) / (parse_float(data["Xpd"]) - parse_float(data["Xls"]))^2,
            # γ_qd = (parse_float(data["Xq"]) - parse_float(data["Xls"])) / (parse_float(data["Xd"]) - parse_float(data["Xls"])),
        )
machine_SPM(data) = SalientPoleQuadratic(
            R = parse_float(data["Rs"]),
            Td0_p = parse_float(data["Tpd0"]),
            Td0_pp = parse_float(data["Tppd0"]),
            Tq0_pp = parse_float(data["Tppq0"]),
            Xd = parse_float(data["Xd"]),
            Xq = parse_float(data["Xq"]),
            Xd_p = parse_float(data["Xpd"]),
            Xd_pp = parse_float(data["Xppd"]),
            Xl = parse_float(data["Xls"]),
            # Se = (0.0, 0.0)
            Se = (parse_float.(split(data["Im_Points"][2:end-1], " ")[parse(Int, data["NdSat"])-1]), parse_float.(split(data["Im_Points"][2:end-1], " ")[parse(Int, data["NdSat"])])),
            # γ_d1 = (parse_float(data["Xppd"]) - parse_float(data["Xls"])) / (parse_float(data["Xpd"]) - parse_float(data["Xls"])),
            # γ_q1 = (parse_float(data["Xppq"]) - parse_float(data["Xls"])) / (parse_float(data["Xpq"]) - parse_float(data["Xls"])),
            # γ_d2 = (parse_float(data["Xpd"]) - parse_float(data["Xppd"])) / (parse_float(data["Xpd"]) - parse_float(data["Xls"]))^2,
        )
machine_classic() = BaseMachine(
            0.0, #R
            0.2995, #Xd_p
            0.7087, #eq_p
        )
shaft_damping(data) = SingleMass(
            H = parse_float(data["InertiaConstant"]),
            D = parse_float(data["AbsoluteDamping"])
        )
shaft_damping_simple() = SingleMass(
    3.148, #H
    2.0, #D
    )
avr_ieeet1(data, vref) = IEEET1(
            Tr = parse_float(data["Tr_param"]),
            Ka = parse_float(data["Ka_param"]),
            Ta = parse_float(data["Ta_param"]),
            Vr_lim = (parse_float(data["VRmin_param"]), parse_float(data["VRmax_param"])),
            Ke = parse_float(data["Ke_param"]),
            Te = parse_float(data["Te_param"]),
            Kf = parse_float(data["Kf_param"]),
            Tf = parse_float(data["Tf_param"]),
            switch = 1,
            E_sat = (parse_float.(split(data["Efd12_param"][2:end-1], " ")[2]), parse_float.(split(data["Efd12_param"][2:end-1], " ")[3])),
            Se = (parse_float.(split(data["SeEfd12_param"][2:end-1], " ")[2]), parse_float.(split(data["SeEfd12_param"][2:end-1], " ")[3])),
            V_ref = parse_float(vref)
            # saturation_coeffs::Tuple{Float64, Float64} # this might be "SeEfd12_param"
        )
avr_type1(data, vref) = AVRTypeI(
            Ka = parse_float(data["Ka_param"]),
            Ke = parse_float(data["Ke_param"]),
            Kf = parse_float(data["Kf_param"]),
            Ta = parse_float(data["Ta_param"]),
            Te = parse_float(data["Te_param"]),
            Tf = parse_float(data["Tf_param"]),
            Tr = 0.00000001, # parse_float(data["Tr_param"]), most of these values are 0 in the table
            Va_lim = (min=-50.0, max=50.0),
            ### Ae = s1 / ( exp( ln(s1/s2) * e1 / (e1-e2) ) )
            Ae = parse_float.(split(data["SeEfd12_param"][2:end-1], " ")[2]) / (exp( parse_float.(split(data["SeEfd12_param"][2:end-1], " ")[2]) * (log(parse_float.(split(data["SeEfd12_param"][2:end-1], " ")[2]) / parse_float.(split(data["SeEfd12_param"][2:end-1], " ")[3])) / (parse_float.(split(data["Efd12_param"][2:end-1], " ")[2]) - parse_float.(split(data["Efd12_param"][2:end-1], " ")[3])))  ) ),
            ### Be = ln(s1/s2) / (e1-e2)
            Be = log(parse_float.(split(data["SeEfd12_param"][2:end-1], " ")[2]) / parse_float.(split(data["SeEfd12_param"][2:end-1], " ")[3])) / (parse_float.(split(data["Efd12_param"][2:end-1], " ")[2]) - parse_float.(split(data["Efd12_param"][2:end-1], " ")[3])),
            # Ae = 0.0039,
            # Be = 1.555,
            V_ref = parse_float(vref)
            # Vr_lim = (parse_float(data["VRmin_param"]), parse_float(data["VRmax_param"])),
            # E_sat = (parse_float.(split(data["Efd12_param"][2:end-1], " ")[2]), parse_float.(split(data["Efd12_param"][2:end-1], " ")[3])),
            # Se = (parse_float.(split(data["SeEfd12_param"][2:end-1], " ")[2]), parse_float.(split(data["SeEfd12_param"][2:end-1], " ")[3])),
        )
# avr_sexs(data, vref) = SEXS(
#             Ta_Tb::Float64
#             Tb::Float64
#             K::Float64
#             Te::Float64
#             V_lim::MinMax
#             V_ref = parse_float(vref)
avr_none() = AVRFixed(0.0)
tg_hygov(data) = HydroTurbineGov(
            R = parse_float(data["R_param"]),
            r = parse_float(data["r_param"]),
            Tr = parse_float(data["TR_param"]),
            Tf = parse_float(data["TF_param"]),
            Tg = parse_float(data["TG_param"]),
            VELM = parse_float(data["VELM_param"]),
            gate_position_limits = (parse_float(data["GMIN_param"]), parse_float(data["GMAX_param"])),
            Tw = parse_float(data["TW_param"]),
            At = parse_float(data["At_param"]),
            D_T = parse_float(data["Dturb_param"]),
            q_nl = parse_float(data["qNL_param"]),
            P_ref = parse_float(data["Pm0_param"])
        )

tg_tgov1(data, pss_data) = SteamTurbineGov1(
            R = parse_float(data["R_param"]),
            T1 = parse_float(data["T1_param"]),
            valve_position_limits = (parse_float(data["vmin_param"]), parse_float(data["vmax_param"])), 
            T2 = parse_float(data["T2_param"]),
            T3 = parse_float(data["T3_param"]),
            D_T = parse_float(data["Dt_param"]),
            DB_h = 0,
            DB_l = 0,
            T_rate = 0,
            P_ref = parse_float(pss_data["Pe0_param"])
        )

tg_none() = TGFixed(1.0) #efficiency
pss_pss2b(data) = PSS2B(
            input_code_1 = 1,
            remote_bus_control_1 = 0,
            input_code_2 = 1,
            remote_bus_control_2 = 0,
            M_rtf = parse(Int, data["M_param"]),
            N_rtf = parse(Int, data["N_param"]),
            Tw1 = parse_float(data["Tw1_param"]),
            Tw2 = parse_float(data["Tw2_param"]), # !== 0.0 ? parse_float(data["Tw2_param"]) : 1,
            T6 = parse_float(data["T6_param"]), # !== 0.0 ? parse_float(data["T6_param"]) : 0.001,
            Tw3 = parse_float(data["Tw3_param"]),
            Tw4 = parse_float(data["Tw4_param"]),
            T7 = parse_float(data["T7_param"]), # !== 0.0 ? parse_float(data["T7_param"]) : 4,
            Ks2 = parse_float(data["Ks2_param"]),
            Ks3 = parse_float(data["Ks3_param"]),
            T8 = parse_float(data["T8_param"]), # !== 0.0 ? parse_float(data["T8_param"]) : 0.1,
            T9 = parse_float(data["T9_param"]), # !== 0.0 ? parse_float(data["T9_param"]) : 0.1,
            Ks1 = parse_float(data["Ks1_param"]),
            T1 = parse_float(data["T1_param"]),
            T2 = parse_float(data["T2_param"]),
            T3 = parse_float(data["T3_param"]),
            T4 = parse_float(data["T4_param"]),
            T10 = parse_float(data["T10_param"]),
            T11 = 0.0001, #parse_float(data["T11_param"]),
            Vs1_lim = (parse_float(data["VSI1min_param"]), parse_float(data["VSI1max_param"])),
            Vs2_lim = (parse_float(data["VSI2min_param"]), parse_float(data["VSI2max_param"])),
            Vst_lim = (parse_float(data["VSTmin_param"]), parse_float(data["VSTmax_param"]))
        )

pss_none() = PSSFixed(0.0)

#######################################################

dyn_gen_hygov(static_gen, Gen_data, Exciter_data, vref, HYGOV_data, PSS_data) = DynamicGenerator(
            name = PSY.get_name(static_gen),
            ω_ref = 1.0,
            # machine = machine_AFM(Gen_data),
            machine = machine_SPM(Gen_data),
            shaft = shaft_damping(Gen_data),
            avr = avr_type1(Exciter_data, vref),
            prime_mover = tg_hygov(HYGOV_data),
            pss = pss_pss2b(PSS_data)
        )

dyn_gen_tgov1(static_gen, Gen_data, Exciter_data, vref, TGOV1_data, PSS_data) = DynamicGenerator(
            name = PSY.get_name(static_gen),
            ω_ref = 1.0,
            # machine = machine_AFM(Gen_data),
            machine = machine_RRM(Gen_data),
            shaft = shaft_damping(Gen_data),
            avr = avr_type1(Exciter_data, vref),
            prime_mover = tg_tgov1(TGOV1_data, PSS_data),
            pss = pss_pss2b(PSS_data)
        )

dyn_gen_simple_hygov(static_gen, Gen_data, HYGOV_data, Exciter_data, vref, PSS_data) = DynamicGenerator(
            name = PSY.get_name(static_gen),
            ω_ref = 1.0,
            # machine = machine_classic(),
            # shaft = shaft_damping_simple(),
            # machine = machine_AFM(Gen_data),
            machine = machine_SPM(Gen_data),
            shaft = shaft_damping(Gen_data),
            # avr = avr_none(),
            avr = avr_type1(Exciter_data, vref),
            # prime_mover = tg_none(),
            prime_mover = tg_hygov(HYGOV_data),
            # pss = pss_none(),
            pss = pss_pss2b(PSS_data)
        )

dyn_gen_simple_tgov1(static_gen, Gen_data, TGOV1_data, Exciter_data, vref, PSS_data) = DynamicGenerator(
            name = PSY.get_name(static_gen),
            ω_ref = 1.0,
            # machine = machine_classic(),
            # machine = machine_AFM(Gen_data),
            machine = machine_RRM(Gen_data),
            shaft = shaft_damping(Gen_data),
            # avr = avr_none(),
            avr = avr_type1(Exciter_data, vref),
            prime_mover = tg_tgov1(TGOV1_data, PSS_data),
            # pss = pss_none(),
            pss = pss_pss2b(PSS_data)
        )

dyn_gen_simple_v1(static_gen, Gen_data) = DynamicGenerator(
            name = PSY.get_name(static_gen),
            ω_ref = 1.0,
            # machine = machine_classic(),
            machine = machine_AFM(Gen_data),
            shaft = shaft_damping(Gen_data),
            avr = avr_none(),
            prime_mover = tg_none(),
            pss = pss_none(),
        )

dyn_gen_simple(static_gen) = DynamicGenerator(
            name = PSY.get_name(static_gen),
            ω_ref = 1.0,
            machine = machine_classic(),
            shaft = shaft_damping_simple(),
            avr = avr_none(),
            prime_mover = tg_none(),
            pss = pss_none(),
        )

statis_source(gen, bus) = Source(
            name = String(gen["name"]),
            available = true,
            bus = bus,
            active_power = gen["pg"],
            reactive_power = gen["qg"],
            R_th = 0.,
            X_th = 0.001,
            internal_voltage = bus.magnitude,
            internal_angle = bus.angle,
            dynamic_injector = nothing
        )