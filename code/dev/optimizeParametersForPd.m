
% create data
T = 60;
dt = 0.001;
time = (dt : dt : T)';
t_min = 10;
t_max = 120;

step_frequency = 2;
oscillation_amplitude = 0;
%     oscillation_amplitude = 0.25;

t_acc_start = 5;
t_acc_duration = 2;

speed_base = 1;
speed_delta = 1;

g_p = 0.15;
g_v = 0.25;

vel_base = speed_base + speed_delta * sinusSigmoid(time, t_acc_start, t_acc_start + t_acc_duration);
vel_oscillation = oscillation_amplitude * sin(2 * pi * time * step_frequency);
heelstrike_oscillation = sin(2 * pi * time * step_frequency);
[~, heelstrike_indices] = findpeaks(heelstrike_oscillation);

vel = vel_base + vel_oscillation;
pos = cumtrapz(time, vel);
v_init = vel(1);
time_indices_to_analyze = (time > t_min & time < t_max);
    
[pos_paced_sms, pace_sms] = selfPacer_sms(pos, time, v_init, g_p, g_v, heelstrike_indices);
[pos_paced_pidAvg, pace_pidAvg] = selfPacer_pidWithAveraging(pos, time, 0.1, 0.25, v_init, heelstrike_indices);

indices_to_compare = (time >= t_acc_start);

options = optimoptions('fminunc', 'MaxIterations', 200, 'OptimalityTolerance', 1e-4);
x0 = [0.2, 0.4];
fun = @(x)(sum((selfPacer_pidWithAveraging(pos, time, x(1), x(2), v_init, heelstrike_indices) - pos_paced_sms).^2));

% x_solution = fminunc(fun,x0)

% x_solution = [0.215711430241940   0.446652431903695]; % for g_p, g_v = 0.1, 0.25, speed_delta = 1, t_acc_duration = 0.5
% x_solution = [0.215380849943874   0.446725643127156]; % for g_p, g_v = 0.1, 0.25, speed_delta = 1, t_acc_duration = 2
% x_solution = [0.213735518485160   0.447292360822209]; % for g_p, g_v = 0.1, 0.25, speed_delta = 1, t_acc_duration = 5
% x_solution = [0.205617543670900   0.448987670029429]; % for g_p, g_v = 0.1, 0.25, speed_delta = 1, t_acc_duration = 15

% x_solution = [0.449757541893018   0.490173726158898]; % for g_p, g_v = 0.2, 0.3, speed_delta = 2, t_acc_duration = 0.5
% x_solution = [0.448410398034812   0.490363366842182]; % for g_p, g_v = 0.2, 0.3, speed_delta = 2, t_acc_duration = 2
% x_solution = [0.441404379661514   0.492067201558008]; % for g_p, g_v = 0.2, 0.3, speed_delta = 2, t_acc_duration = 5
% x_solution = [0.408390574276085   0.498082800447377]; % for g_p, g_v = 0.2, 0.3, speed_delta = 2, t_acc_duration = 15

% trying to find g_p, g_v values with corresponding k_p, k_v values similar to what we usually use, k_p = 0.3, k_v = 

x_solution = [0.2154, 0.4467];

[pos_paced_pidAvg, pace_pidAvg] = selfPacer_pidWithAveraging(pos, time, x_solution(1), x_solution(2), v_init, heelstrike_indices);
error = fun(x_solution);

figure; hold on
plot(time, pos_paced_pidAvg, 'linewidth', 2)
plot(time, pos_paced_sms, 'linewidth', 2)

