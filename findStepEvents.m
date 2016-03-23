% findStepEvents

trial_number = 3;

%% prepare

% load data
load subjectInfo.mat;
model_file_name = makeFileName(date, subject_id, 'model');
load(model_file_name);
marker_trajectories_file_name = makeFileName(date, subject_id, 'walking', trial_number, 'markerTrajectories');
load(marker_trajectories_file_name);
marker_trajectories_file_name = makeFileName(date, subject_id, 'walking', trial_number, 'forcePlateData');
load(marker_trajectories_file_name);

% extract data
left_heel_marker = 34;
left_toes_marker = 35;
right_heel_marker = 42;
right_toes_marker = 43;
left_heel_marker_indices = reshape([(left_heel_marker - 1) * 3 + 1; (left_heel_marker - 1) * 3 + 2; (left_heel_marker - 1) * 3 + 3], 1, length(left_heel_marker)*3);
left_toes_marker_indices = reshape([(left_toes_marker - 1) * 3 + 1; (left_toes_marker - 1) * 3 + 2; (left_toes_marker - 1) * 3 + 3], 1, length(left_toes_marker)*3);
right_heel_marker_indices = reshape([(right_heel_marker - 1) * 3 + 1; (right_heel_marker - 1) * 3 + 2; (right_heel_marker - 1) * 3 + 3], 1, length(right_heel_marker)*3);
right_toes_marker_indices = reshape([(right_toes_marker - 1) * 3 + 1; (right_toes_marker - 1) * 3 + 2; (right_toes_marker - 1) * 3 + 3], 1, length(right_toes_marker)*3);
left_heel_marker_z_trajectory = marker_trajectories(:, left_heel_marker_indices(3));
left_toes_marker_z_trajectory = marker_trajectories(:, left_toes_marker_indices(3));
right_heel_marker_z_trajectory = marker_trajectories(:, right_heel_marker_indices(3));
right_toes_marker_z_trajectory = marker_trajectories(:, right_toes_marker_indices(3));

% calculate derivatives
filter_order = 2;
cutoff_frequency = 10; % cutoff frequency, in Hz

[b, a] = butter(filter_order, cutoff_frequency/(samplingRate/2));	% set filter parameters for butterworth filter: 2=order of filter;
left_toes_marker_z_vel_trajectory = centdiff(filtfilt(b, a, left_toes_marker_z_trajectory), 1/samplingRate);
right_toes_marker_z_vel_trajectory = centdiff(filtfilt(b, a, right_toes_marker_z_trajectory), 1/samplingRate);

%% find events

% find touch down indices as negative peaks of the heel marker z-position
peak_width_threshold = samplingRate * 0.25;
peak_prominence_threshold = 0.5;
[~, left_heel_peak_locations, left_heel_peak_widths] = findpeaks(-left_heel_marker_z_trajectory);
left_touchdown_indices_mocap = left_heel_peak_locations(left_heel_peak_widths > peak_width_threshold);
[~, right_heel_peak_locations, right_heel_peak_widths] = findpeaks(-right_heel_marker_z_trajectory);
right_touchdown_indices_mocap = right_heel_peak_locations(right_heel_peak_widths > peak_width_threshold);


% for touchdown, find the first significant toes z-velocity peak after each heelstrike
[~, left_toes_vel_peak_locations, left_toes_vel_peak_widths, left_toes_vel_peak_prominence] = findpeaks(left_toes_marker_z_vel_trajectory);
left_toes_vel_peak_indices = left_toes_vel_peak_locations(left_toes_vel_peak_prominence > peak_prominence_threshold);
left_pushoff_indices_mocap = [];
for i_touchdown = 1 : length(left_touchdown_indices_mocap)
    left_toes_vel_peak_indices_after_touchdown = left_toes_vel_peak_indices(left_toes_vel_peak_indices > left_touchdown_indices_mocap(i_touchdown));
    if ~isempty(left_toes_vel_peak_indices_after_touchdown)
        left_pushoff_index = left_toes_vel_peak_indices_after_touchdown(1);
        left_pushoff_indices_mocap = [left_pushoff_indices_mocap left_pushoff_index];
    end
end

[~, right_toes_vel_peak_locations, right_toes_vel_peak_widths, right_toes_vel_peak_prominence] = findpeaks(right_toes_marker_z_vel_trajectory);
right_toes_vel_peak_indices = right_toes_vel_peak_locations(right_toes_vel_peak_prominence > peak_prominence_threshold);
right_pushoff_indices_mocap = [];
for i_touchdown = 1 : length(right_touchdown_indices_mocap)
    right_toes_vel_peak_indices_after_touchdown = right_toes_vel_peak_indices(right_toes_vel_peak_indices > right_touchdown_indices_mocap(i_touchdown));
    if ~isempty(right_toes_vel_peak_indices_after_touchdown)
        right_pushoff_index = right_toes_vel_peak_indices_after_touchdown(1);
        right_pushoff_indices_mocap = [right_pushoff_indices_mocap right_pushoff_index];
    end
end

% transform to force plate time
left_pushoff_indices_force_plate = zeros(size(left_pushoff_indices_mocap));
for i_index = 1 : length(left_pushoff_indices_mocap)
    [~, index_force_plate] = min(abs(time_force_plate - time_mocap(left_pushoff_indices_mocap(i_index))));
    left_pushoff_indices_force_plate(i_index) = index_force_plate;
end
left_touchdown_indices_force_plate = zeros(size(left_touchdown_indices_mocap));
for i_index = 1 : length(left_touchdown_indices_mocap)
    [~, index_force_plate] = min(abs(time_force_plate - time_mocap(left_touchdown_indices_mocap(i_index))));
    left_touchdown_indices_force_plate(i_index) = index_force_plate;
end
right_pushoff_indices_force_plate = zeros(size(right_pushoff_indices_mocap));
for i_index = 1 : length(right_pushoff_indices_mocap)
    [~, index_force_plate] = min(abs(time_force_plate - time_mocap(right_pushoff_indices_mocap(i_index))));
    right_pushoff_indices_force_plate(i_index) = index_force_plate;
end
right_touchdown_indices_force_plate = zeros(size(right_touchdown_indices_mocap));
for i_index = 1 : length(right_touchdown_indices_mocap)
    [~, index_force_plate] = min(abs(time_force_plate - time_mocap(right_touchdown_indices_mocap(i_index))));
    right_touchdown_indices_force_plate(i_index) = index_force_plate;
end



%% form contact indicators
number_of_time_steps = length(time_mocap);

left_contact_indicators_mocap = formContactIndicatorTrajectory(left_pushoff_indices_mocap, left_touchdown_indices_mocap, number_of_time_steps);
right_contact_indicators_mocap = formContactIndicatorTrajectory(right_pushoff_indices_mocap, right_touchdown_indices_mocap, number_of_time_steps);
left_contact_indicators_force_plate = formContactIndicatorTrajectory(left_pushoff_indices_force_plate, left_touchdown_indices_force_plate, number_of_time_steps);
right_contact_indicators_force_plate = formContactIndicatorTrajectory(right_pushoff_indices_force_plate, right_touchdown_indices_force_plate, number_of_time_steps);

% form contact trajectories
left_heel_contact_trajectories = left_heel_marker_z_trajectory; left_heel_contact_trajectories(~left_contact_indicators_mocap, :) = NaN;
left_toes_contact_trajectories = left_toes_marker_z_trajectory; left_toes_contact_trajectories(~left_contact_indicators_mocap, :) = NaN;
left_toes_velocity_contact_trajectories = left_toes_marker_z_vel_trajectory; left_toes_velocity_contact_trajectories(~left_contact_indicators_mocap, :) = NaN;

right_heel_contact_trajectories = right_heel_marker_z_trajectory; right_heel_contact_trajectories(~right_contact_indicators_mocap, :) = NaN;
right_toes_contact_trajectories = right_toes_marker_z_trajectory; right_toes_contact_trajectories(~right_contact_indicators_mocap, :) = NaN;
right_toes_velocity_contact_trajectories = right_toes_marker_z_vel_trajectory; right_toes_velocity_contact_trajectories(~right_contact_indicators_mocap, :) = NaN;












step_events_file_name = makeFileName(date, subject_id, 'walking', trial_number, 'stepEvents');
save ...
  ( ...
    step_events_file_name, ...
    'left_touchdown_indices_mocap', ...
    'right_touchdown_indices_mocap', ...
    'left_pushoff_indices_mocap', ...
    'right_pushoff_indices_mocap', ...
    'left_contact_indicators_mocap', ...
    'right_contact_indicators_mocap', ...
    'left_touchdown_indices_force_plate', ...
    'right_touchdown_indices_force_plate', ...
    'left_pushoff_indices_force_plate', ...
    'right_pushoff_indices_force_plate', ...
    'left_contact_indicators_force_plate', ...
    'right_contact_indicators_force_plate' ...
  );









%% visualize
force_scaler = 2e-4;
vel_scaler = 5e-2;
acc_scaler = 1e-2;



% left events
figure; axes_left = axes; hold on
plot(time_force_plate, fzl_trajectory*force_scaler)
plot(time_mocap(left_touchdown_indices_mocap), left_toes_marker_z_trajectory(left_touchdown_indices_mocap)*0, 'o', 'linewidth', 2);
plot(time_mocap(left_pushoff_indices_mocap), left_toes_marker_z_trajectory(left_pushoff_indices_mocap)*0, 'o', 'linewidth', 2);
plot(time_force_plate(left_pushoff_indices_force_plate), zeros(size(left_pushoff_indices_force_plate)), 'x', 'linewidth', 2)
plot(time_force_plate(left_touchdown_indices_force_plate), zeros(size(left_touchdown_indices_force_plate)), 'x', 'linewidth', 2)

% right events
figure; axes_right = axes; hold on
plot(time_force_plate, fzr_trajectory*force_scaler)
plot(time_mocap(right_touchdown_indices_mocap), right_toes_marker_z_trajectory(right_touchdown_indices_mocap)*0, 'o', 'linewidth', 2);
plot(time_mocap(right_pushoff_indices_mocap), right_toes_marker_z_trajectory(right_pushoff_indices_mocap)*0, 'o', 'linewidth', 2);
plot(time_force_plate(right_pushoff_indices_force_plate), zeros(size(right_pushoff_indices_force_plate)), 'x', 'linewidth', 2)
plot(time_force_plate(right_touchdown_indices_force_plate), zeros(size(right_touchdown_indices_force_plate)), 'x', 'linewidth', 2)

linkaxes([axes_left, axes_right], 'x')
distFig('rows', 2)

return
% right events
figure; axes; hold on
plot(time_force_plate, fzr_trajectory*force_scaler)
plot(time_mocap, right_heel_marker_z_trajectory, 'linewidth', 1);
% plot(time_mocap, right_heel_contact_trajectories, 'linewidth', 2);
plot(time_mocap, right_toes_marker_z_trajectory, 'linewidth', 1);
% plot(time_mocap, right_toes_contact_trajectories, 'linewidth', 2);
% plot(time_mocap, right_toes_marker_z_vel_trajectory * vel_scaler, 'linewidth', 1);
% plot(right_touchdown_indices, right_toes_marker_z_vel_trajectory(right_touchdown_indices) * vel_scaler, 'o', 'linewidth', 2);
% plot(right_toes_vel_peak_indices, right_toes_marker_z_vel_trajectory(right_toes_vel_peak_indices) * vel_scaler, 'o', 'linewidth', 2);
% plot(time_mocap(right_touchdown_indices), right_toes_marker_z_trajectory(right_touchdown_indices)*0, 'o', 'linewidth', 2);
% plot(time_mocap(right_pushoff_indices), right_toes_marker_z_trajectory(right_pushoff_indices)*0, 'o', 'linewidth', 2);



% plot(right_toes_vel_peak_locations, right_toes_vel_peak_prominence, 'o', 'linewidth', 2);

return


% left heel position
figure; axes; hold on
plot(time_force_plate, fzl_trajectory*force_scaler)
plot(time_mocap, left_heel_marker_z_trajectory, 'linewidth', 1);
plot(time_mocap, left_heel_marker_z_vel_trajectory * vel_scaler, 'linewidth', 2);
plot(time_mocap, left_heel_marker_z_acc_trajectory * acc_scaler, 'linewidth', 2);
plot(time_mocap(left_touchdown_indices), left_heel_marker_z_trajectory(left_touchdown_indices), 'o');
legend('force', 'pos', 'vel', 'acc')

% left toes position
figure; axes; hold on
plot(time_force_plate, fzl_trajectory*force_scaler)
plot(time_mocap, left_toes_marker_z_trajectory, 'linewidth', 1);
plot(time_mocap, left_toes_marker_z_vel_trajectory * vel_scaler, 'linewidth', 2);
plot(time_mocap(left_touchdown_indices), left_toes_marker_z_vel_trajectory(left_touchdown_indices) * vel_scaler, 'o', 'linewidth', 2);
% plot(time_mocap, left_toes_marker_z_acc_trajectory * acc_scaler, 'linewidth', 2);
% legend('force', 'pos', 'vel', 'acc')
return



% left heel position
figure; axes; hold on
% plot(time_force_plate, fzl_trajectory*force_scaler)
plot(time_mocap, left_heel_marker_z_trajectory);
% plot(time_mocap(left_heel_peak_locations), left_heel_marker_z_trajectory(left_heel_peak_locations), 'o');
plot(time_mocap(left_heel_touchdown_indices), left_heel_marker_z_trajectory(left_heel_touchdown_indices), 'o');

% left toes position
figure; axes; hold on
% plot(time_force_plate, fzl_trajectory*force_scaler)
plot(time_mocap, left_toes_marker_z_trajectory);
% plot(time_mocap(left_toes_peak_locations), left_toes_marker_z_trajectory(left_toes_peak_locations), 'o');
plot(time_mocap(left_toes_pushoff_indices), left_toes_marker_z_trajectory(left_toes_pushoff_indices), 'o');

% right heel position
figure; axes; hold on
% plot(time_force_plate, fzl_trajectory*force_scaler)
plot(time_mocap, right_heel_marker_z_trajectory);
% plot(time_mocap(right_heel_peak_locations), right_heel_marker_z_trajectory(right_heel_peak_locations), 'o');
plot(time_mocap(right_heel_touchdown_indices), right_heel_marker_z_trajectory(right_heel_touchdown_indices), 'o');

% right toes position
figure; axes; hold on
% plot(time_force_plate, fzl_trajectory*force_scaler)
plot(time_mocap, right_toes_marker_z_trajectory);
% plot(time_mocap(right_toes_peak_locations), right_toes_marker_z_trajectory(right_toes_peak_locations), 'o');
plot(time_mocap(right_toes_pushoff_indices), right_toes_marker_z_trajectory(right_toes_pushoff_indices), 'o');

return

% left foot velocity
figure; axes; hold on
plot(time_mocap, left_heel_marker_z_trajectory);
plot(time_mocap, left_heel_marker_z_velocity_trajectory);
plot(time_mocap(left_heel_marker_upward_indices), left_heel_marker_z_velocity_trajectory(left_heel_marker_upward_indices), 'o');
plot(time_mocap(left_heel_marker_downward_indices), left_heel_marker_z_velocity_trajectory(left_heel_marker_downward_indices), 'o');

% plot(time_mocap, left_heel_vel_above_pos_threshold);
% plot(time_mocap, -left_heel_vel_below_neg_threshold);
plot(time_mocap, left_heel_marker_z_velocity_pos_threshold*ones(size(time_mocap)))
plot(time_mocap, left_heel_marker_z_velocity_neg_threshold*ones(size(time_mocap)))










return
