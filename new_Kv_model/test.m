clc
close all
clear variables

% time space
hold_pt = 100;
end_pt = 5000;
time_space = cell(1, 3);

hold_t = 0:hold_pt;
pulse_t = (hold_pt + 1):end_pt;
pulse_t_adj = pulse_t - pulse_t(1);
t = [hold_t, pulse_t];

time_space{1} = t;
time_space{2} = hold_t;
time_space{3} = pulse_t_adj;

% voltages
v = -50:10:50;

% Ikslow1
p1 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 803.0, 18.0, 0.9214774521, 0.05766, 0.07496];
y = Ikslow1(p1, -70, v(1), time_space, -91.1);

figure(1)
plot(t, y)
hold on
for i=2:length(v)
   y = Ikslow1(p1, -70, v(i), time_space, -91.1);
   plot(t, y)
end
hold off

% Ikslow2
p2 = [22.5, 45.2, 40.0, 7.7, 5.7, 6.1, 0.0629, 2.058, 5334, 4912, 0.05766];
y = Ikslow2(p2, -70, v(1), time_space, -91.1);

figure(2)
plot(t, y)
hold on
for i=2:length(v)
   y = Ikslow2(p2, -70, v(i), time_space, -91.1);
   plot(t, y)
end
hold off

% Ikss
p3 = [22.5, 40.0, 7.7, 0.0862, 1235.5, 13.17, 0.0428];
y = Ikss(p3, -70, v(1), time_space, -91.1);

figure(3)
plot(t, y)
hold on
for i=2:length(v)
   y = Ikss(p3, -70, v(i), time_space, -91.1);
   plot(t, y)
end
hold off

% Ikto
p0 = [33, 15.5, 20, 16, 8, 7, 0.03577, 0.06237, 0.04516, 4, ...
    0.0989, 0.0019, 0.067083, 0.051335, 0.08, 0.5, 0.2087704319, 0.14067, 0.387];

y =Ikto(p0, -70, v(1), time_space, -91.1);
figure(4)
plot(t, y)
hold on
for i=2:length(v)
    y = Ikto(p0, -70, v(i), time_space, -91.1);
    plot(t,y)
end
hold off
