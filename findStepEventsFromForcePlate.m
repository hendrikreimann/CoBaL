% find step events from force plate

trials_to_process = 6 : 25;
trials_to_process = 3;

for i_trial = trials_to_process

    % load data
    load subjectInfo.mat;
    model_file_name = makeFileName(date, subject_id, 'model');
    load(model_file_name);
    marker_trajectories_file_name = makeFileName(date, subject_id, 'walking', i_trial, 'markerTrajectories');
    load(marker_trajectories_file_name);
    marker_trajectories_file_name = makeFileName(date, subject_id, 'walking', i_trial, 'forcePlateData');
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

    % find events
    cop_threshold = 0.1;
    threshold_side = sign(abs(copxl_trajectory) - cop_threshold);
    left_touchdown_indices_force_plate = find(diff(sign(abs(copxl_trajectory) - cop_threshold)) > 0) + 1;
    left_pushoff_indices_force_plate = find(diff(sign(abs(copxl_trajectory) - cop_threshold)) < 0);
    right_touchdown_indices_force_plate = find(diff(sign(abs(copxr_trajectory) - cop_threshold)) > 0) + 1;
    right_pushoff_indices_force_plate = find(diff(sign(abs(copxr_trajectory) - cop_threshold)) < 0);

    % transform to mocap time
    left_pushoff_indices_mocap = zeros(size(left_pushoff_indices_force_plate));
    for i_index = 1 : length(left_pushoff_indices_force_plate)
        [~, index_mocap] = min(abs(time_mocap - time_force_plate(left_pushoff_indices_force_plate(i_index))));
        left_pushoff_indices_mocap(i_index) = index_mocap;
    end
    left_touchdown_indices_mocap = zeros(size(left_touchdown_indices_force_plate));
    for i_index = 1 : length(left_touchdown_indices_force_plate)
        [~, index_mocap] = min(abs(time_mocap - time_force_plate(left_touchdown_indices_force_plate(i_index))));
        left_touchdown_indices_mocap(i_index) = index_mocap;
    end
    right_pushoff_indices_mocap = zeros(size(right_pushoff_indices_force_plate));
    for i_index = 1 : length(right_pushoff_indices_force_plate)
        [~, index_mocap] = min(abs(time_mocap - time_force_plate(right_pushoff_indices_force_plate(i_index))));
        right_pushoff_indices_mocap(i_index) = index_mocap;
    end
    right_touchdown_indices_mocap = zeros(size(right_touchdown_indices_force_plate));
    for i_index = 1 : length(right_touchdown_indices_force_plate)
        [~, index_mocap] = min(abs(time_mocap - time_force_plate(right_touchdown_indices_force_plate(i_index))));
        right_touchdown_indices_mocap(i_index) = index_mocap;
    end

    %% form contact indicators
    number_of_time_steps = length(time_mocap);

    left_contact_indicators_mocap = formContactIndicatorTrajectory(left_pushoff_indices_mocap, left_touchdown_indices_mocap, number_of_time_steps);
    right_contact_indicators_mocap = formContactIndicatorTrajectory(right_pushoff_indices_mocap, right_touchdown_indices_mocap, number_of_time_steps);
    left_contact_indicators_force_plate = formContactIndicatorTrajectory(left_pushoff_indices_force_plate, left_touchdown_indices_force_plate, number_of_time_steps);
    right_contact_indicators_force_plate = formContactIndicatorTrajectory(right_pushoff_indices_force_plate, right_touchdown_indices_force_plate, number_of_time_steps);

    step_events_file_name = makeFileName(date, subject_id, 'walking', i_trial, 'stepEvents');
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
    
    disp(['Trial ' num2str(i_trial) ' completed']);
end



% left events
force_scaler = 2e-4;
figure; axes_left = axes; hold on
plot(time_force_plate, fzl_trajectory*force_scaler)
plot(time_mocap(left_touchdown_indices_mocap), left_toes_marker_z_trajectory(left_touchdown_indices_mocap)*0, 'o', 'linewidth', 2);
plot(time_mocap(left_pushoff_indices_mocap), left_toes_marker_z_trajectory(left_pushoff_indices_mocap)*0, 'o', 'linewidth', 2);
plot(time_mocap, left_heel_marker_z_trajectory);

% right events
figure; axes_right = axes; hold on
plot(time_force_plate, fzr_trajectory*force_scaler)
plot(time_mocap(right_touchdown_indices_mocap), right_toes_marker_z_trajectory(right_touchdown_indices_mocap)*0, 'o', 'linewidth', 2);
plot(time_mocap(right_pushoff_indices_mocap), right_toes_marker_z_trajectory(right_pushoff_indices_mocap)*0, 'o', 'linewidth', 2);
plot(time_mocap, right_heel_marker_z_trajectory);

linkaxes([axes_left, axes_right], 'x')
distFig('rows', 2)
return



figure; axes; hold on
plot(time_force_plate, copxl_trajectory);
plot(time_force_plate(left_touchdown_indices_force_plate), copxl_trajectory(left_touchdown_indices_force_plate), 'o');
plot(time_force_plate(left_pushoff_indices_force_plate), copxl_trajectory(left_pushoff_indices_force_plate), 'o');

figure; axes; hold on
plot(time_force_plate, copxr_trajectory);
plot(time_force_plate(right_touchdown_indices_force_plate), copxr_trajectory(right_touchdown_indices_force_plate), 'o');
plot(time_force_plate(right_pushoff_indices_force_plate), copxr_trajectory(right_pushoff_indices_force_plate), 'o');

% figure; axes; hold on
% plot(time_force_plate, fzl_trajectory)
% plot(time_force_plate(left_touchdown_indices_force_plate), fzl_trajectory(left_touchdown_indices_force_plate), 'o');
% plot(time_force_plate(left_pushoff_indices_force_plate), fzl_trajectory(left_pushoff_indices_force_plate), 'o');

% plot(time_force_plate, copxr_trajectory);
