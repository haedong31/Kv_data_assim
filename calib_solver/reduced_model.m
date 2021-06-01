function [ykto, ykslow1, ykslow2, ykss, yksum] = reduced_model(p, hold_volt, volt, time_space, Ek)
    % default values
    kto0 = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
        0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
    kslow10 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
    kslow20 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 5334, 4912, 0.05766];
    kss0 = [22.5, 40.0, 7.7, 0.0862, 1235.5, 13.17, 0.0428];

    param_kto = zeros(17,1);
    param_kslow1 = zeros(13,1);
    param_kslow2 = zeros(11,1);
    param_kss = zeros(7,1);

    % fixed parameters for Ikto: p4, p7, p8, p9, p11, p12, fkto
    fixed_kto_idx = [4, 7, 8, 9, 11, 12, 15];
    param_kto(fixed_kto_idx) = kto0(fixed_kto_idx);
    
    % tuning parameters for Ikto
    tune_kto_idx = setdiff(1:17, fixed_kto_idx);
    param_kto(tune_kto_idx) = p(1:10);

    % fixed parameters for Ikslow1: p4, p5, p6, p7, p10, fkslow1
    fixed_kslow1_idx = [4, 5, 6, 7, 10, 11];
    param_kslow1(fixed_kslow1_idx) = kslow10(fixed_kslow1_idx);

    % tuning parameters for Ikslow1
    tune_kslow1_idx = setdiff(1:13, fixed_kslow1_idx);
    param_kslow1(tune_kslow1_idx) = p(11:17);

    % fixed parameters for Ikslow2
    fixed_kslow2_idx = [4, 5, 6, 7, 10];
    param_kslow2(fixed_kslow2_idx) = kslow20(fixed_kslow2_idx);

    % tuning parameters for Ikslow2
    tune_kslow2_idx = setdiff(1:11, fixed_kslow2_idx);
    param_kslow2(tune_kslow2_idx) = p([11, 12, 13, 14, 18, 16]);
    
    % fixed parameters for Ikss
    fixed_kss_idx = [3, 4, 5];
    param_kss(fixed_kss_idx) = kss0(fixed_kss_idx);

    % tuning parameters for Ikss
    tune_kss_idx = setdiff(1:7, fixed_kss_idx);
    param_kss(tune_kss_idx) = p([11, 13, 19, 20]);

    % generate K+ currents
    ykto = ikto(param_kto, hold_volt, volt, time_space, Ek);
    ykslow1 = ikslow1(param_kslow1, hold_volt, volt, time_space, Ek);
    ykslow2 = ikslow2(param_kslow2, hold_volt, volt, time_space, Ek);
    ykss = ikss(param_kss, hold_volt, volt, time_space, Ek);
    yksum = ykto + ykslow1 + ykslow2 + ykss;
end
