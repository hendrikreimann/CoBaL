% load data
load dev_data.mat

% pace - linear
[pos_paced_linear, pace_linear] = selfPacer_linear(pos, time);

% pace - pid
[pos_paced_pid, pace_pid] = selfPacer_pid(pos, time, 0.5, 0, 0.5, v_init);

% pace - ideal (splined)
[pos_paced_splined, pace_splined] = selfPacer_spline(pos, time, heelstrike_indices);

% pace - pid with position averaged across previous step
[pos_paced_pidAvg, pace_pidAvg] = selfPacer_pidWithAveraging(pos, time, 0.5, 0, 0.5, v_init, heelstrike_indices);

% pace - minimum jerk
[pos_paced_mmj, pace_mmj] = selfPacer_mmj(pos, time, v_init, heelstrike_indices);


% calculate derivatives
pace_accel_linear = deriveByTime(pace_linear, time);
pace_accel_pid = deriveByTime(pace_pid, time);
pace_accel_splined = deriveByTime(pace_splined, time);
pace_accel_pidAvg = deriveByTime(pace_pidAvg, time);
pace_accel_mmj = deriveByTime(pace_mmj, time);

pace_jerk_linear = deriveByTime(pace_accel_linear, time);
pace_jerk_pid = deriveByTime(pace_accel_pid, time);
pace_jerk_splined = deriveByTime(pace_accel_splined, time);
pace_jerk_pidAvg = deriveByTime(pace_accel_pidAvg, time);
pace_jerk_mmj = deriveByTime(pace_accel_mmj, time);

% calculate mean jerks
mean_accel_linear = rms(pace_accel_linear);
mean_accel_pid = rms(pace_accel_pid);
mean_accel_splined = rms(pace_accel_splined);
mean_accel_pidAvg = rms(pace_accel_pidAvg);
mean_accel_mmj = rms(pace_accel_mmj);
mean_jerk_linear = rms(pace_jerk_linear);
mean_jerk_pid = rms(pace_jerk_pid);
mean_jerk_splined = rms(pace_jerk_splined);
mean_jerk_pidAvg = rms(pace_jerk_pidAvg);
mean_jerk_mmj = rms(pace_jerk_mmj);

% plot
figure; hold on; title('positions - world space')
plot(time, pos_paced_linear, 'displayname', 'linear');
plot_pid = plot(time, pos_paced_pid, 'displayname', 'pid', 'linewidth', 2);
plot_splined = plot(time, pos_paced_splined, 'displayname', 'splined', 'linewidth', 2);
plot_pidAvg = plot(time, pos_paced_pidAvg, 'displayname', 'pid with averaging', 'linewidth', 2);
plot_mmj = plot(time, pos_paced_mmj, 'displayname', 'minimum jerk', 'linewidth', 2);
legend('show')

figure; hold on; title('pace velocity')
plot(time, pace_linear, 'displayname', 'linear');
plot(time, pace_pid, 'displayname', 'pid', 'linewidth', 2);
plot(time, pace_splined, 'displayname', 'splined', 'linewidth', 2);
plot(time, pace_pidAvg, 'displayname', 'pid with averaging', 'linewidth', 2);
plot(time, pace_mmj, 'displayname', 'minimum jerk', 'linewidth', 2);
legend('show')

figure; hold on; title('pace acceleration')
plot(time, pace_accel_linear, 'displayname', 'linear');
plot(time, pace_accel_pid, 'displayname', 'pid', 'linewidth', 2);
plot(time, pace_accel_splined, 'displayname', 'splined', 'linewidth', 2);
plot(time, pace_accel_pidAvg, 'displayname', 'pid with averaging', 'linewidth', 2);
plot(time, pace_accel_mmj, 'displayname', 'minimum jerk', 'linewidth', 2);
legend('show')

figure; hold on; title('pace jerk')
plot(time, pace_jerk_linear, 'displayname', 'linear');
plot(time, pace_jerk_pid, 'displayname', 'pid', 'linewidth', 2);
plot(time, pace_jerk_splined, 'displayname', 'splined', 'linewidth', 2);
plot(time, pace_jerk_pidAvg, 'displayname', 'pid with averaging', 'linewidth', 2);
plot(time, pace_jerk_mmj, 'displayname', 'minimum jerk', 'linewidth', 2);
legend('show')

figure; 
subplot(2, 1, 1); hold on; title('mean acceleration')
bar(1, mean_accel_pid, 'facecolor', get(plot_pid, 'color'))
bar(2, mean_accel_splined, 'facecolor', get(plot_splined, 'color'))
bar(3, mean_accel_pidAvg, 'facecolor', get(plot_pidAvg, 'color'))
bar(4, mean_accel_mmj, 'facecolor', get(plot_mmj, 'color'))
set(gca, 'xtick', []);
subplot(2, 1, 2); hold on; title('mean jerk')
bar(1, mean_jerk_pid, 'facecolor', get(plot_pid, 'color'))
bar(2, mean_jerk_splined, 'facecolor', get(plot_splined, 'color'))
bar(3, mean_jerk_pidAvg, 'facecolor', get(plot_pidAvg, 'color'))
bar(4, mean_jerk_mmj, 'facecolor', get(plot_mmj, 'color'))
set(gca, 'xtick', [1 2 3], 'xticklabels', {'pid', 'splined', 'pid with averaging', 'minimum jerk'});

return
% plot
figure; hold on
plot(time, pos)
plot(time(heelstrike_indices), pos(heelstrike_indices), 'v');