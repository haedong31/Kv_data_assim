function [ykto, ykslow, ykss, yksum] = kcurrent_model2(p, hold_volt, volt, time_space, Ek)
    % default values
    kto0 = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
        0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
    kslow0 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
    kss0 = [22.5, 40.0, 7.7, 0.0862, 1235.5, 13.17, 0.0428];

    param_kto = zeros(17,1);
    param_kslow = zeros(13,1);
    param_kss = zeros(7,1);

    % fixed parameters for ikto: p4, p7, p8, p9, p11, p12, fkto
    fixed_kto_idx = [4, 7, 8, 9, 11, 12, 15];
    param_kto(fixed_kto_idx) = kto0(fixed_kto_idx);
    
    % tuning parameters for ikto
    tune_kto_idx = setdiff(1:17, fixed_kto_idx);
    param_kto(tune_kto_idx) = p(1:10);

    % some parameters of ikslow are shared with ikss
    % ikslow: p1, p2, p3, p4, p5, p6, p7, p8
    % ikss: p1, p3, p4 (* for p1, p2, p3 of ikss)

    % fixed parameters for ikslow: p4, p6, p7, p10, fkslow1
    fixed_kslow_idx = [4, 6, 7, 10, 11];
    param_kslow(fixed_kslow_idx) = kslow0(fixed_kslow_idx);

    % tuning parameters for ikslow
    tune_kslow_idx = setdiff(1:13, fixed_kslow_idx);
    param_kslow(tune_kslow_idx) = p(11:18);
    
    % shared parameters for ikss
    shared_ikss_idx = 1:3;
    param_kss(shared_ikss_idx) = param_kslow([1, 3, 4]);

    % fixed parameters for ikss
    fixed_kss_idx = [4, 5];
    param_kss(fixed_kss_idx) = kss0(fixed_kss_idx);

    % tuning parameters for ikss; p6, GKss
    tune_kss_idx = setdiff(1:7, [shared_ikss_idx, fixed_kss_idx]);
    param_kss(tune_kss_idx) = p([19, 20]);

    % generate K+ currents
    ykto = ikto(param_kto, hold_volt, volt, time_space, Ek);
    ykslow = ikslow1(param_kslow, hold_volt, volt, time_space, Ek);
    ykss = ikss(param_kss, hold_volt, volt, time_space, Ek);
    yksum = ykto + ykslow + ykss;
end
