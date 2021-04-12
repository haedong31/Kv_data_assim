function [ykto, ykslow1, ykslow2, ykss, yksum] = full_model(p, hold_volt, volt, time_space, Ek)
    % default values
    kto0 = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
        0.000152, 0.067083, 0.00095, 0.051335, 0.2087704319, 0.14067, 0.387];
    kslow10 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
    kslow20 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 5334, 4912, 0.05766];
    kss0 = [22.5, 40.0, 7.7, 0.0862, 1235.5, 13.17, 0.0428];

    % assign parameters
    param_kto = zeros(17,1);
    param_kslow1 = zeros(13,1);
    param_kslow2 = zeros(11,1);
    param_kss = zeros(7,1);

    % for Iktof
    param_kto(1:14) = p(1:14);
    param_kto(15) = kto0(15);
    param_kto(16:17) = p(15:16);

    % for Ikslow1
    param_kslow1(1:10) = p(17:26);
    param_kslow1(11) = kslow10(11);
    param_kslow1(12:13) = p(27:28);

    % for Ikslow2
    param_kslow2(1:8) = p(17:24);
    param_kslow2(9:10) = p(29:30);
    param_kslow2(11) = p(27);

    % for Ikss
    param_kss(1) = p(17);
    param_kss(2) = p(19);
    param_kss(3) = p(20);
    param_kss(4:7) = p(31:34);

    % generate K+ currents
    ykto = Ikto(param_kto, hold_volt, volt, time_space, Ek);
    ykslow1 = Ikslow1(param_kslow1, hold_volt, volt, time_space, Ek);
    ykslow2 = Ikslow2(param_kslow2, hold_volt, volt, time_space, Ek);
    ykss = Ikss(param_kss, hold_volt, volt, time_space, Ek);
    yksum = ykto + ykslow1 + ykslow2 + ykss;
end
