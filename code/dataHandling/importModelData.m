% EarlyLeft: 1 in early swing, 0 in late swing
% EarlyRight: 1 in early swing, 0 in late swing
% Contact: stance/swing

% load data
loaded_data = load(['out' filesep 'out_normal']); % TODO: generalize this
model_results_raw = loaded_data.out_normal;

if ~directoryExists('processed')
    mkdir('processed')
end
if ~directoryExists('analysis')
    mkdir('analysis')
end

subject_settings = loadSettingsFromFile('subject');
collection_date = subject_settings.get('collection_date');
subject_id = subject_settings.get('subject_id');
trial_type = 'walking';
trial_number = 1;

% determine heel-strikes
left_heelstrike_time_indices = diff(model_results_raw.Contact(:, 1)) == 1;
right_heelstrike_time_indices = diff(model_results_raw.Contact(:, 2)) == 1;
left_heelstrike_times = model_results_raw.tout(left_heelstrike_time_indices);
right_heelstrike_times = model_results_raw.tout(right_heelstrike_time_indices);

left_pushoff_time_indices = diff(model_results_raw.Contact(:, 1)) == -1;
right_pushoff_time_indices = diff(model_results_raw.Contact(:, 2)) == -1;
left_pushoff_times = model_results_raw.tout(left_pushoff_time_indices);
right_pushoff_times = model_results_raw.tout(right_pushoff_time_indices);

variables_to_save.event_labels = {'left_pushoff';'left_touchdown';'right_pushoff';'right_touchdown';};
variables_to_save.event_data = {left_pushoff_times; left_heelstrike_times; right_pushoff_times; right_heelstrike_times};

step_events_file_name = ['analysis' filesep makeFileName(collection_date, subject_id, trial_type, trial_number, 'events.mat')];
saveDataToFile(step_events_file_name, variables_to_save);

% resample
delta_t = 0.005;
time = 0 : delta_t : 100;
model_results = struct;
model_results.time = time';
model_results.trunk_pos = spline(model_results_raw.tout', model_results_raw.xyzHAT', time)';
model_results.joint_angles_left_leg = spline(model_results_raw.tout', model_results_raw.Llegq', time)';
model_results.joint_angles_right_leg = spline(model_results_raw.tout', model_results_raw.Rlegq', time)';
model_results.left_ankle = spline(model_results_raw.tout', model_results_raw.AnklePosL', time)';
model_results.right_ankle = spline(model_results_raw.tout', model_results_raw.AnklePosR', time)';

% markers
ml_index = 3;
ap_index = 1;
vert_index = 2;
marker_trajectories = ...
  [ ...
    model_results.left_ankle(:, ml_index), model_results.left_ankle(:, ap_index), model_results.left_ankle(:, vert_index), ...
    model_results.right_ankle(:, ml_index), model_results.right_ankle(:, ap_index), model_results.right_ankle(:, vert_index) ...
  ];
marker_labels = {'LANK_x', 'LANK_y', 'LANK_z', 'RANK_x', 'RANK_y', 'RANK_z'};

number_of_marker_trajectories = size(marker_trajectories, 2);
marker_directions = cell(2, number_of_marker_trajectories);
[marker_directions{1, 1 : 3 : number_of_marker_trajectories}] = deal('right');
[marker_directions{2, 1 : 3 : number_of_marker_trajectories}] = deal('left');
[marker_directions{1, 2 : 3 : number_of_marker_trajectories}] = deal('forward');
[marker_directions{2, 2 : 3 : number_of_marker_trajectories}] = deal('backward');
[marker_directions{1, 3 : 3 : number_of_marker_trajectories}] = deal('up');
[marker_directions{2, 3 : 3 : number_of_marker_trajectories}] = deal('down');

save_folder = 'processed';
time_marker = time;
sampling_rate_marker = delta_t^(-1);
save_file_name = makeFileName(collection_date, subject_id, trial_type, trial_number, 'markerTrajectories.mat');
save ...
    ( ...
    [save_folder filesep save_file_name], ...
    'marker_trajectories', ...
    'time_marker', ...
    'sampling_rate_marker', ...
    'marker_labels', ...
    'marker_directions' ...
    );
addAvailableData('marker_trajectories', 'time_marker', 'sampling_rate_marker', '_marker_labels', '_marker_directions', save_folder, save_file_name);



% joint angles
hip_flexion_sign = -1;
hip_flexion_offset = 180;
hip_abduction_sign = -1;
hip_abduction_offset = 0;
knee_flexion_sign = -1;
knee_flexion_offset = 180;
ankle_flexion_sign = -1;
ankle_flexion_offset = 90;
left_hip_flexion_trajectory = hip_flexion_sign * rad2deg(model_results.joint_angles_left_leg(:, 1)) + hip_flexion_offset;
left_hip_abduction_trajectory = hip_abduction_sign * rad2deg(model_results.joint_angles_left_leg(:, 2)) + hip_abduction_offset;
left_knee_flexion_trajectory = knee_flexion_sign * rad2deg(model_results.joint_angles_left_leg(:, 3)) + knee_flexion_offset;
left_ankle_flexion_trajectory = ankle_flexion_sign * rad2deg(model_results.joint_angles_left_leg(:, 4)) + ankle_flexion_offset;
right_hip_flexion_trajectory = hip_flexion_sign * rad2deg(model_results.joint_angles_right_leg(:, 1)) + hip_flexion_offset;
right_hip_abduction_trajectory = hip_abduction_sign * rad2deg(model_results.joint_angles_right_leg(:, 2)) + hip_abduction_offset;
right_knee_flexion_trajectory = knee_flexion_sign * rad2deg(model_results.joint_angles_right_leg(:, 3)) + knee_flexion_offset;
right_ankle_flexion_trajectory = ankle_flexion_sign * rad2deg(model_results.joint_angles_right_leg(:, 4)) + ankle_flexion_offset;

joint_angle_trajectories = ...
  [ ...
    left_hip_flexion_trajectory, ...
    left_hip_abduction_trajectory, ...
    left_knee_flexion_trajectory, ...
    left_ankle_flexion_trajectory, ...
    right_hip_flexion_trajectory, ...
    right_hip_abduction_trajectory, ...
    right_knee_flexion_trajectory, ...
    right_ankle_flexion_trajectory ...
  ];
joint_labels = ...
  { ...
    'left_hip_flexion', 'left_hip_adduction', 'left_knee_flexion', 'left_ankle_flexion', ...
    'right_hip_flexion', 'right_hip_adduction', 'right_knee_flexion', 'right_ankle_flexion', ...
  };
joint_directions = ...
  {
    'flexion', 'adduction', 'flexion', 'dorsiflexion', 'flexion', 'adduction', 'flexion', 'dorsiflexion'; ...
    'extension', 'abduction', 'extension', 'plantarflexion', 'extension', 'abduction', 'extension', 'plantarflexion';
  };

% save joint angles
variables_to_save = struct;
variables_to_save.joint_angle_trajectories = joint_angle_trajectories;
variables_to_save.time = time;
variables_to_save.sampling_rate = delta_t^(-1);
variables_to_save.joint_labels = joint_labels;
variables_to_save.joint_directions = joint_directions;

save_folder = 'processed';
save_file_name = makeFileName(collection_date, subject_id, trial_type, trial_number, 'kinematicTrajectories.mat');
save([save_folder filesep save_file_name], '-struct', 'variables_to_save');
addAvailableData ...
  ( ...
    'joint_angle_trajectories', ...
    'time', ...
    'sampling_rate', ...
    '_joint_labels', ...
    '_joint_directions', ...
    save_folder, ...
    save_file_name ...
  );











