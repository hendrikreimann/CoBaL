% test the PID self pacer
dt = 0.01;
time = 0 : dt : 20;

omega = 2;
offset = 1;
vel = sin(time * 2 * pi * omega) + offset;
pos = cumtrapz(time, vel);

belt_velocity_timeseries = pidSelfPacer(pos, time, 1, 0, 0.5, vel(1));
belt_pos_timeseries = cumtrapz(time, belt_velocity_timeseries);
pos_beltspace = pos - belt_pos_timeseries;

% figure; hold on; title('world space')
% plot(time, pos)

figure; hold on; title('belt space')
plot(time, pos_beltspace)

figure; hold on; title('velocity')
plot(time, belt_velocity_timeseries)
