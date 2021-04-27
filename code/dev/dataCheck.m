% extract pelvis position in anterior-posterior direction from inverse kinematics data
labels = strsplit(data_ik.header{9}, '\t');
pelvis_ap_data_indicator = strcmp(labels(2:end), 'pelvis_tx');
pelvis_pos_belt = data_ik.trajectories(:, pelvis_ap_data_indicator);
pelvis_time = data_ik.time;

% extract pelvis position in anterior-posterior direction from marker data
LPSI_y_indicator = strcmp(data_mocap.marker_labels, 'LPSI_y');
marker_time = data_mocap.time_mocap;
LPSI_pos_world = data_mocap.marker_trajectories(:, LPSI_y_indicator);

% use actual treadmill velocity
treadmill_vel_world = data_plc.belt_speed_left_trajectory;
treadmill_time = data_plc.time;
treadmill_pos_world = cumtrapz(treadmill_time, treadmill_vel_world); % this is the translation between world and treadmill
belt_to_world_treadmilltime = -treadmill_pos_world;
belt_to_world_kinematicstime = spline(treadmill_time, belt_to_world_treadmilltime, pelvis_time);

% calculate position of pelvis in world space
pelvis_pos_world = pelvis_pos_belt + belt_to_world_kinematicstime;
pelvis_vel_world = deriveByTime(pelvis_pos_world, pelvis_time);

% find linear fit - best constant treadmill velocity
p_fit = polyfit(pelvis_time, pelvis_pos_belt, 1);
linear_fit_pos_belt = p_fit(1) * pelvis_time + p_fit(2);
linear_fit_time = pelvis_time;
linear_fit_pos_world = linear_fit_pos_belt + belt_to_world_kinematicstime;
linear_fit_vel_world = p_fit(1) * ones(size(pelvis_time));

% simulate PID controller
pid_controlled_treadmill_vel = pidSelfPacer(pelvis_pos_belt, pelvis_time, 0.5, 0, 3.5, treadmill_vel_world(1));
pid_controlled_treadmill_pos = cumtrapz(pelvis_time, pid_controlled_treadmill_vel);
pid_controlled_pelvis_pos_world = pelvis_pos_belt - pid_controlled_treadmill_pos';

% filter belt space position
filter_order = 6;
cutoff_frequency = 1;
sampling_rate = 200;
[b_filter, a_filter] = butter(filter_order, cutoff_frequency/(sampling_rate/2));
filtered_pelvis_pos_belt = filtfilt(b_filter, a_filter, pelvis_pos_belt);
filtered_pelvis_pos_world = filtered_pelvis_pos_belt + belt_to_world_kinematicstime;
filtered_pelvis_vel_world = deriveByTime(filtered_pelvis_pos_belt, pelvis_time);

% filter out step-frequency oscillations by hand to get ideal data
heelstrike_left_times = data_events.event_data{strcmp(data_events.event_labels, 'left_touchdown')};
heelstrike_left_indices = findClosestIndex(heelstrike_left_times, pelvis_time);
heelstrike_right_times = data_events.event_data{strcmp(data_events.event_labels, 'right_touchdown')};
heelstrike_right_indices = findClosestIndex(heelstrike_right_times, pelvis_time);
heelstrike_all_indices = sort([heelstrike_left_indices; heelstrike_right_indices]);
pelvis_pos_world_splined = spline(heelstrike_all_indices, pelvis_pos_world(heelstrike_all_indices), 1 : length(pelvis_time));
pelvis_pos_belt_splined = spline(heelstrike_all_indices, pelvis_pos_belt(heelstrike_all_indices), 1 : length(pelvis_time));
pelvis_pos_belt_splined_world = pelvis_pos_belt_splined - pid_controlled_treadmill_pos;
pelvis_vel_belt_splined = deriveByTime(pelvis_pos_belt_splined, pelvis_time);


% save dev data
dev_data.pos = pelvis_pos_belt;
dev_data.time = pelvis_time;
dev_data.heelstrike_indices = heelstrike_all_indices;
dev_data.v_init = treadmill_vel_world(1);
% save('dev_data', '-struct', 'dev_data');


% plot in belt space
figure; hold on; title('position - belt space')
plot(pelvis_time, pelvis_pos_belt)
plot(treadmill_time, treadmill_pos_world, 'displayname', 'actual self-pacing')
plot(pelvis_time(heelstrike_all_indices), pelvis_pos_belt(heelstrike_all_indices), 'v', 'displayname', 'heelstrikes')
plot(pelvis_time, pelvis_pos_belt_splined, 'displayname', 'pelvis position, splined', 'linewidth', 2)
plot(linear_fit_time, linear_fit_pos_belt, 'displayname', 'optimal constant speed')
plot(pelvis_time, pid_controlled_treadmill_pos, 'displayname', 'pid')
legend('show')

% plot in belt space minus constant
figure; hold on; title('position - belt space minus optimal constant speed')
plot(pelvis_time, pelvis_pos_belt - linear_fit_pos_belt, 'displayname', 'actual pelvis position')
plot(pelvis_time(heelstrike_all_indices), pelvis_pos_belt(heelstrike_all_indices) - linear_fit_pos_belt(heelstrike_all_indices), 'v', 'displayname', 'heelstrikes')
plot(pelvis_time, pelvis_pos_belt_splined - linear_fit_pos_belt', 'displayname', 'pelvis position, splined', 'linewidth', 2)
plot(pelvis_time, pid_controlled_treadmill_pos - linear_fit_pos_belt', 'displayname', 'pid')
legend('show')
return

% plot in world space
f_pos_world = figure; a_pos_world = axes; hold on; title('position - world space')
plot(pelvis_time, pelvis_pos_world, 'displayname', 'pelvis position with actual self-pacing')
plot(pelvis_time(heelstrike_all_indices), pelvis_pos_world(heelstrike_all_indices), 'v', 'displayname', 'heelstrikes')
plot(pelvis_time, pelvis_pos_world_splined, 'displayname', 'pelvis position, splined')
plot(pelvis_time, pelvis_pos_belt_splined_world, 'displayname', 'pelvis position belt, splined, transformed')
plot(pelvis_time, linear_fit_pos_world, 'displayname', 'pelvis position with constant speed')
plot(pelvis_time, pid_controlled_pelvis_pos_world, 'displayname', 'pid')
plot(pelvis_time, filtered_pelvis_pos_world, 'displayname', 'filtered')
legend('show')

% plot velocity
f_vel = figure; a_vel = axes; hold on; title('velocities')
plot(treadmill_time, treadmill_vel_world, 'displayname', 'actual treadmill speed')
plot(linear_fit_time, linear_fit_vel_world, 'displayname', 'optimal constant speed')
plot(pelvis_time, pid_controlled_treadmill_vel, 'displayname', 'pid')
plot(pelvis_time, filtered_pelvis_vel_world, 'displayname', 'filtered')
plot(pelvis_time, pelvis_vel_belt_splined, 'displayname', 'splined')
legend('show')

% layout
set(f_pos_world, 'units', 'normalized', 'position', [0 0.5 0.8 0.45])
set(f_vel, 'units', 'normalized', 'position', [0 0 0.8 0.45])
linkaxes([a_pos_world, a_vel], 'x')



