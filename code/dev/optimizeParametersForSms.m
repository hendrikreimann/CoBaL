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

% inverse problem to check algorithm, the solution to this should be g_p = 0.1, g_v = 0.25
alpha = 0.2154;
beta = 0.4467;

% find SMS gains close to what we usually use
alpha = 0.3;
beta = 0.54;


vel_base = speed_base + speed_delta * sinusSigmoid(time, t_acc_start, t_acc_start + t_acc_duration);
vel_oscillation = oscillation_amplitude * sin(2 * pi * time * step_frequency);
heelstrike_oscillation = sin(2 * pi * time * step_frequency);
[~, heelstrike_indices] = findpeaks(heelstrike_oscillation);

vel = vel_base + vel_oscillation;
pos = cumtrapz(time, vel);
v_init = vel(1);
time_indices_to_analyze = (time > t_min & time < t_max);
    
[pos_paced_pidAvg, pace_pidAvg] = selfPacer_pidWithAveraging(pos, time, alpha, beta, v_init, heelstrike_indices);
[pos_paced_sms, pace_sms] = selfPacer_sms(pos, time, v_init, 0.1, 0.25, heelstrike_indices);

indices_to_compare = (time >= t_acc_start);

options = optimoptions('fminunc', 'MaxIterations', 200, 'OptimalityTolerance', 1e-4);
x0 = [0.1, 0.25];
fun = @(x)(sum((selfPacer_sms(pos, time, v_init, x(1), x(2), heelstrike_indices) - pos_paced_pidAvg).^2));

x_solution = fminunc(fun,x0)

% trying to find g_p, g_v values with corresponding k_p, k_v values similar to what we usually use, alpha = 0.3, beta = 0.54
% x_solution = [0.1370, 0.3071];

[pos_paced_sms, pace_sms] = selfPacer_sms(pos, time, v_init, x_solution(1), x_solution(2), heelstrike_indices);
error = fun(x_solution);

figure; hold on
plot(time, pos_paced_pidAvg, 'linewidth', 2)
plot(time, pos_paced_sms, 'linewidth', 2)

