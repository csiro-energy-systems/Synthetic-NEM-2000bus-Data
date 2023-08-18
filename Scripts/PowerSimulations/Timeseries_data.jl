using TimeSeries
using Dates
using Random
Random.seed!(123)
using PowerSystems
const PSY = PowerSystems

DayAhead = collect(
    DateTime("1/1/2024  0:00:00", "d/m/y  H:M:S"):Hour(1):DateTime(
        "1/1/2024  23:00:00",
        "d/m/y  H:M:S",
    ),
)
#Dispatch_11am =  collect(DateTime("1/1/2024  0:11:00", "d/m/y  H:M:S"):Minute(15):DateTime("1/1/2024  12::00", "d/m/y  H:M:S"))

solar_ts_DA = [
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0.351105684
    0.632536266
    0.99463925
    1
    0.944237283
    0.396681234
    0.366511428
    0.155125829
    0.040872694
    0
    0
    0
    0
    0
    0
]

wind_ts_DA = [
    0.985205412
    0.991791369
    0.997654144
    1
    0.998663733
    0.995497149
    0.992414567
    0.98252418
    0.957203427
    0.927650911
    0.907181989
    0.889095913
    0.848186718
    0.766813846
    0.654052531
    0.525336131
    0.396098004
    0.281771509
    0.197790004
    0.153241012
    0.131355854
    0.113688144
    0.099302656
    0.069569628
]

hydro_inflow_ts_DA = [
    0.314300
    0.386684
    0.228582
    0.226677
    0.222867
    0.129530
    0.144768
    0.365731
    0.207628
    0.622885
    0.670507
    0.676221
    0.668602
    0.407638
    0.321919
    0.369541
    0.287632
    0.449544
    0.630505
    0.731462
    0.777178
    0.712413
    0.780988
    0.190485
];

# hydro_inflow_ts_DA /= maximum(hydro_inflow_ts_DA)

loadbus2_ts_DA = [
    0.792729978
    0.723201574
    0.710952098
    0.677672816
    0.668249175
    0.67166919
    0.687608809
    0.711821241
    0.756320618
    0.7984057
    0.827836527
    0.840362459
    0.84511032
    0.834592803
    0.822949221
    0.816941743
    0.824079963
    0.905735139
    0.989967048
    1
    0.991227765
    0.960842114
    0.921465115
    0.837001437
]

loadbus3_ts_DA = [
    0.831093782
    0.689863228
    0.666058513
    0.627033103
    0.624901388
    0.62858924
    0.650734211
    0.683424321
    0.750876413
    0.828347191
    0.884248576
    0.888523615
    0.87752169
    0.847534405
    0.8227661
    0.803809323
    0.813282799
    0.907575962
    0.98679848
    1
    0.990489904
    0.952520972
    0.906611479
    0.624307054
    # 0.824307054
]

loadbus4_ts_DA = [
    0.871297342
    0.670489749
    0.642812243
    0.630092987
    0.652991383
    0.671971681
    0.716278493
    0.770885833
    0.810075243
    0.85562361
    0.892440566
    0.910660449
    0.922135467
    0.898416969
    0.879816542
    0.896390855
    0.978598576
    0.96523761
    1
    0.969626503
    0.901212601
    0.81894251
    0.771004923
    0.717847996
]

### Operating Reserve Demand Curve
ORDC_cost = [(9000.0, 0.0), (6000.0, 0.2), (500.0, 0.4), (10.0, 0.6), (0.0, 0.8)]

ORDC_cost_ts = [
    TimeSeries.TimeArray(DayAhead, repeat([ORDC_cost], 24)),
    TimeSeries.TimeArray(DayAhead + Day(1), repeat([ORDC_cost], 24)),
]

# TODO: add a sensible cost for hybrid devices
# hybrid_cost_ts = [
#     TimeSeries.TimeArray(DayAhead, repeat([25.0], 24)),
#     TimeSeries.TimeArray(DayAhead + Day(1), repeat([25.0], 24)),
# ]

# Reserve_ts = [TimeSeries.TimeArray(DayAhead, rand(24)), TimeSeries.TimeArray(DayAhead + Day(1), rand(24))]

hydro_timeseries_DA = [
    [TimeSeries.TimeArray(DayAhead, hydro_inflow_ts_DA)],
    [TimeSeries.TimeArray(DayAhead + Day(1), ones(24) * 0.1 + hydro_inflow_ts_DA)],
];

storage_target = zeros(24)
storage_target[end] = 0.1
storage_target_DA = [
   [TimeSeries.TimeArray(DayAhead, storage_target)],
   [TimeSeries.TimeArray(DayAhead + Day(1), storage_target)],
];

hydro_budget_DA = [
    [TimeSeries.TimeArray(DayAhead, hydro_inflow_ts_DA * 0.8)],
    [TimeSeries.TimeArray(DayAhead + Day(1), hydro_inflow_ts_DA * 0.8)],
];

RealTime = collect(
    DateTime("1/1/2024 0:00:00", "d/m/y H:M:S"):Minute(5):DateTime(
        "1/1/2024 23:55:00",
        "d/m/y H:M:S",
    ),
)

hydro_timeseries_RT = [
    [TimeSeries.TimeArray(RealTime, repeat(hydro_inflow_ts_DA, inner = 12))],
    [TimeSeries.TimeArray(RealTime + Day(1), ones(288) * 0.1 + repeat(hydro_inflow_ts_DA, inner = 12))],
];

storage_target_RT = [
    [TimeSeries.TimeArray(RealTime, repeat(storage_target, inner = 12))],
    [TimeSeries.TimeArray(RealTime + Day(1), repeat(storage_target, inner = 12))],
];

hydro_budget_RT = [
    [TimeSeries.TimeArray(RealTime, repeat(hydro_inflow_ts_DA  * 0.8, inner = 12))],
    [TimeSeries.TimeArray(RealTime + Day(1), repeat(hydro_inflow_ts_DA  * 0.8, inner = 12))],
];

hybrid_cost_ts_RT = [
    [TimeSeries.TimeArray(RealTime, repeat([25.0], 288))],
    [TimeSeries.TimeArray(RealTime + Day(1), ones(288) * 0.1 + repeat([25.0], 288))],
];


load_timeseries_DA = [
    [TimeSeries.TimeArray(DayAhead, loadbus3_ts_DA)],
    [TimeSeries.TimeArray(DayAhead + Day(1), rand(24) * 0.1 + loadbus3_ts_DA)],
];

load_timeseries_RT = [
    [TimeSeries.TimeArray(RealTime, repeat(loadbus4_ts_DA, inner = 12))],
    [TimeSeries.TimeArray(RealTime + Day(1), rand(288) * 0.01 + repeat(loadbus4_ts_DA, inner = 12))],
]


solar_timeseries_DA = [
    [TimeSeries.TimeArray(DayAhead, solar_ts_DA)],
    [TimeSeries.TimeArray(DayAhead + Day(1), (1 .- (rand(24) * 0.1)) .* solar_ts_DA)],
];

solar_timeseries_RT = [
    [TimeSeries.TimeArray(RealTime, repeat(solar_ts_DA, inner = 12))],
    [TimeSeries.TimeArray(RealTime + Day(1), (1 .- (rand(288) * 0.1)) .* repeat(solar_ts_DA, inner = 12))],
]

wind_timeseries_DA = [
    [TimeSeries.TimeArray(DayAhead, wind_ts_DA)],
    [TimeSeries.TimeArray(DayAhead + Day(1), rand(24) * 0.1 + wind_ts_DA)],
];

wind_timeseries_RT = [
    [TimeSeries.TimeArray(RealTime, repeat(wind_ts_DA, inner = 12))],
    [TimeSeries.TimeArray(RealTime + Day(1), rand(288) * 0.1 + repeat(wind_ts_DA, inner = 12))],
]

Iload_timeseries_DA = [
    [TimeSeries.TimeArray(DayAhead, loadbus4_ts_DA)],
    [TimeSeries.TimeArray(DayAhead + Day(1), loadbus4_ts_DA + 0.1 * rand(24))],
]

Iload_timeseries_RT = [
    [TimeSeries.TimeArray(RealTime, repeat(loadbus4_ts_DA, inner = 12))],
    [TimeSeries.TimeArray(RealTime + Day(1), rand(288) * 0.1 + repeat(loadbus4_ts_DA, inner = 12))],
]



