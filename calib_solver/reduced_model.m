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

    % fixed parameters for Ikto
    % p4, p7, p8, p9, p11, p12, fkto


    % fixed parameters for Ikslow1

    % fixed parameters for Ikslow2

    % fixed parameters for Ikss

    % tuning parameters for Ikto

    % tuning parameters for Ikslow1

    % tuning parameters for Ikslow2

    % tuning parameters for Ikss

end