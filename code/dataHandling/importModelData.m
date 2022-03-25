% EarlyLeft: 1 in early swing, 0 in late swing
% EarlyRight: 1 in early swing, 0 in late swing
% Contact: stance/swing

% prep
if ~directoryExists('processed')
    mkdir('processed')
end
if ~directoryExists('analysis')
    mkdir('analysis')
end

study_settings = loadSettingsFromFile('study');
subject_settings = loadSettingsFromFile('subject');
collection_date = subject_settings.get('collection_date');
subject_id = subject_settings.get('subject_id');
trial_type = 'walking';
data_to_import = study_settings.get('data_to_import', 1);

% load data
source_dir = 'out';
filename_label_trunk = subject_settings.get('filename_label_trunk');
data_dir_mat = dir([source_dir filesep '*mat']);
[file_name_list{1:length(data_dir_mat)}] = deal(data_dir_mat.name);

number_of_files = length(file_name_list);
for i_file = 1 : number_of_files
    % import this file
    data_file_name = file_name_list{i_file};
    loaded_data = load(['out' filesep data_file_name]);
    file_name_split = strsplit(data_file_name, '.');
    trial_number = str2double(file_name_split{1}((length(filename_label_trunk)+1):end)); % this assumes the naming scheme of "walkN", where N is the trial number
    model_results_raw = loaded_data.(file_name_split{1});

    % prepare resampling
    delta_t = 0.005;
    time = 0 : delta_t : 100;
    model_results = struct;
    model_results.time = time';
    
    if any(strcmp(data_to_import, 'gait_events'))
        % determine heel-strikes
        left_heelstrike_time_indices = diff(model_results_raw.Contact(:, 1)) == 1;
        right_heelstrike_time_indices = diff(model_results_raw.Contact(:, 2)) == 1;
        left_heelstrike_times = model_results_raw.tout(left_heelstrike_time_indices);
        right_heelstrike_times = model_results_raw.tout(right_heelstrike_time_indices);

        left_pushoff_time_indices = diff(model_results_raw.Contact(:, 1)) == -1;
        right_pushoff_time_indices = diff(model_results_raw.Contact(:, 2)) == -1;
        left_pushoff_times = model_results_raw.tout(left_pushoff_time_indices);
        right_pushoff_times = model_results_raw.tout(right_pushoff_time_indices);

        left_midswing_time_indices = diff(model_results_raw.EarlyLeft) == -1;
        right_midswing_time_indices = diff(model_results_raw.EarlyRight) == -1;
        left_midswing_times = model_results_raw.tout(left_midswing_time_indices);
        right_midswing_times = model_results_raw.tout(right_midswing_time_indices);

        variables_to_save.event_labels = {'left_pushoff';'left_midswing';'left_touchdown';'right_pushoff';'right_midswing';'right_touchdown';};
        variables_to_save.event_data = {left_pushoff_times; left_midswing_times; left_heelstrike_times; right_pushoff_times; right_midswing_times; right_heelstrike_times};

        step_events_file_name = ['analysis' filesep makeFileName(collection_date, subject_id, trial_type, trial_number, 'events.mat')];
        saveDataToFile(step_events_file_name, variables_to_save);
    end
    
    if any(strcmp(data_to_import, 'muscle_forces'))
        % assemble trajectories
        muscle_forces_left = spline(model_results_raw.tout', model_results_raw.LForce', time)';
        muscle_forces_right = spline(model_results_raw.tout', model_results_raw.RForce', time)';
        muscle_force_trajectories = [muscle_forces_left muscle_forces_right];
        number_of_muscles_per_leg = size(muscle_forces_left, 2);
        
        % assemble labels
        muscle_force_labels_left = cell(1, number_of_muscles_per_leg);
        muscle_force_labels_right = cell(1, number_of_muscles_per_leg);
        muscle_labels = subject_settings.get('muscle_labels');
        for i_muscle = 1 : number_of_muscles_per_leg
            muscle_force_labels_left{i_muscle} = [muscle_labels{i_muscle} '_left'];
            muscle_force_labels_right{i_muscle} = [muscle_labels{i_muscle} '_right'];
        end
        muscle_force_labels = [muscle_force_labels_left muscle_force_labels_right];
        
        % assemble directions
        muscle_force_directions_left = repmat({'+'; '-'}, 1, number_of_muscles_per_leg);
        muscle_force_directions_right = repmat({'+'; '-'}, 1, number_of_muscles_per_leg);
        muscle_force_directions = [muscle_force_directions_left muscle_force_directions_right];
        
        % save
        variables_to_save = struct;
        variables_to_save.muscle_force_trajectories = muscle_force_trajectories;
        variables_to_save.time = time;
        variables_to_save.sampling_rate = delta_t^(-1);
        variables_to_save.muscle_force_labels = muscle_force_labels;
        variables_to_save.muscle_force_directions = muscle_force_directions;

        save_folder = 'processed';
        save_file_name = makeFileName(collection_date, subject_id, trial_type, trial_number, 'muscleForceTrajectories.mat');
        save([save_folder filesep save_file_name], '-struct', 'variables_to_save');
        addAvailableData ...
          ( ...
            'muscle_force_trajectories', ...
            'time', ...
            'sampling_rate', ...
            '_muscle_force_labels', ...
            '_muscle_force_directions', ...
            save_folder, ...
            save_file_name ...
          );
    end
    
    if any(strcmp(data_to_import, 'muscle_activations'))
        % assemble trajectories
        muscle_activations_left = spline(model_results_raw.tout', model_results_raw.LForce', time)';
        muscle_activations_right = spline(model_results_raw.tout', model_results_raw.RForce', time)';
        muscle_activation_trajectories = [muscle_activations_left muscle_activations_right];
        number_of_muscles_per_leg = size(muscle_activations_left, 2);
        
        % assemble labels
        muscle_activation_labels_left = cell(1, number_of_muscles_per_leg);
        muscle_activation_labels_right = cell(1, number_of_muscles_per_leg);
        muscle_labels = subject_settings.get('muscle_labels');
        for i_muscle = 1 : number_of_muscles_per_leg
            muscle_activation_labels_left{i_muscle} = [muscle_labels{i_muscle} '_left'];
            muscle_activation_labels_right{i_muscle} = [muscle_labels{i_muscle} '_right'];
        end
        muscle_activation_labels = [muscle_activation_labels_left muscle_activation_labels_right];
        
        % assemble directions
        muscle_activation_directions_left = repmat({'+'; '-'}, 1, number_of_muscles_per_leg);
        muscle_activation_directions_right = repmat({'+'; '-'}, 1, number_of_muscles_per_leg);
        muscle_activation_directions = [muscle_activation_directions_left muscle_activation_directions_right];
        
        % save
        variables_to_save = struct;
        variables_to_save.muscle_activation_trajectories = muscle_activation_trajectories;
        variables_to_save.time = time;
        variables_to_save.sampling_rate = delta_t^(-1);
        variables_to_save.muscle_activation_labels = muscle_activation_labels;
        variables_to_save.muscle_activation_directions = muscle_activation_directions;

        save_folder = 'processed';
        save_file_name = makeFileName(collection_date, subject_id, trial_type, trial_number, 'muscleActivationTrajectories.mat');
        save([save_folder filesep save_file_name], '-struct', 'variables_to_save');
        addAvailableData ...
          ( ...
            'muscle_activation_trajectories', ...
            'time', ...
            'sampling_rate', ...
            '_muscle_activation_labels', ...
            '_muscle_activation_directions', ...
            save_folder, ...
            save_file_name ...
          );
    end
    
    if any(strcmp(data_to_import, 'joint_angles'))
        % get raw data
        model_results.joint_angles_left = spline(model_results_raw.tout', model_results_raw.Llegq', time)';
        model_results.joint_angles_right = spline(model_results_raw.tout', model_results_raw.Rlegq', time)';
        
        % process
        joint_signs = subject_settings.get('joint_signs');
        joint_angle_offsets = subject_settings.get('joint_angle_offsets');
        joint_labels = subject_settings.get('joint_labels');
        joint_directions_pos = subject_settings.get('joint_directions_pos');
        joint_directions_neg = subject_settings.get('joint_directions_neg');
        number_of_joints_per_leg = size(model_results.joint_angles_left, 2);
        
        % transform and assemble trajectories
        joint_angle_trajectories_left = zeros(size(model_results.joint_angles_left));
        joint_angle_trajectories_right = zeros(size(model_results.joint_angles_left));
        for i_joint = 1 : number_of_joints_per_leg
            %
            this_sign = joint_signs(i_joint);
            this_offset = joint_angle_offsets(i_joint);
            this_trajectory_left = model_results.joint_angles_left(:, i_joint);
            joint_angle_trajectories_left(:, i_joint) = this_sign * rad2deg(this_trajectory_left) + this_offset;
            this_trajectory_right = model_results.joint_angles_right(:, i_joint);
            joint_angle_trajectories_right(:, i_joint) = this_sign * rad2deg(this_trajectory_right) + this_offset;
        end
        joint_angle_trajectories = [joint_angle_trajectories_left joint_angle_trajectories_right];
        
        % assemble labels
        joint_angle_labels_left = cell(1, number_of_joints_per_leg);
        joint_angle_labels_right = cell(1, number_of_joints_per_leg);
        for i_joint = 1 : number_of_joints_per_leg
            joint_angle_labels_left{i_joint} = [joint_labels{i_joint} '_left'];
            joint_angle_labels_right{i_joint} = [joint_labels{i_joint} '_right'];
        end
        joint_angle_labels = [joint_angle_labels_left joint_angle_labels_right];
        
        % assemble directions
        joint_angle_directions_left = cell(2, number_of_joints_per_leg);
        joint_angle_directions_right = cell(2, number_of_joints_per_leg);
        for i_joint = 1 : number_of_joints_per_leg
            joint_angle_directions_left(:, i_joint) = {joint_directions_pos{i_joint}; joint_directions_neg{i_joint}};
            joint_angle_directions_right(:, i_joint) = {[joint_directions_pos{i_joint}]; [joint_directions_neg{i_joint}]};
        end
        joint_directions = [joint_angle_directions_left joint_angle_directions_right];

        % save joint angles
        variables_to_save = struct;
        variables_to_save.joint_angle_trajectories = joint_angle_trajectories;
        variables_to_save.time = time;
        variables_to_save.sampling_rate = delta_t^(-1);
        variables_to_save.joint_labels = joint_angle_labels;
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
    end

    if any(strcmp(data_to_import, 'joint_torques'))
        % get raw data
        model_results.joint_torques_left = spline(model_results_raw.tout', model_results_raw.Ltorque', time)';
        model_results.joint_torques_right = spline(model_results_raw.tout', model_results_raw.Rtorque', time)';
        
        % process
        joint_signs = subject_settings.get('joint_signs');
        joint_labels = subject_settings.get('joint_labels');
        joint_directions_pos = subject_settings.get('joint_directions_pos');
        joint_directions_neg = subject_settings.get('joint_directions_neg');
        number_of_joints_per_leg = size(model_results.joint_torques_left, 2);
        
        % transform and assemble trajectories
        joint_torque_trajectories_left = zeros(size(model_results.joint_torques_left));
        joint_torque_trajectories_right = zeros(size(model_results.joint_torques_left));
        for i_joint = 1 : number_of_joints_per_leg
            %
            this_sign = joint_signs(i_joint);
            this_trajectory_left = model_results.joint_torques_left(:, i_joint);
            joint_torque_trajectories_left(:, i_joint) = this_sign * this_trajectory_left;
            this_trajectory_right = model_results.joint_torques_right(:, i_joint);
            joint_torque_trajectories_right(:, i_joint) = this_sign * this_trajectory_right;
        end
        joint_torque_trajectories = [joint_torque_trajectories_left joint_torque_trajectories_right];
        
        % assemble labels
        joint_torque_labels_left = cell(1, number_of_joints_per_leg);
        joint_torque_labels_right = cell(1, number_of_joints_per_leg);
        for i_joint = 1 : number_of_joints_per_leg
            joint_torque_labels_left{i_joint} = [joint_labels{i_joint} '_left'];
            joint_torque_labels_right{i_joint} = [joint_labels{i_joint} '_right'];
        end
        joint_torque_labels = [joint_torque_labels_left joint_torque_labels_right];
        
        % assemble directions
        joint_torque_directions_left = cell(2, number_of_joints_per_leg);
        joint_torque_directions_right = cell(2, number_of_joints_per_leg);
        for i_joint = 1 : number_of_joints_per_leg
            joint_torque_directions_left(:, i_joint) = {joint_directions_pos{i_joint}; joint_directions_neg{i_joint}};
            joint_torque_directions_right(:, i_joint) = {[joint_directions_pos{i_joint}]; [joint_directions_neg{i_joint}]};
        end
        joint_directions = [joint_torque_directions_left joint_torque_directions_right];

        % save joint torques
        variables_to_save = struct;
        variables_to_save.joint_torque_trajectories = joint_torque_trajectories;
        variables_to_save.time = time;
        variables_to_save.sampling_rate = delta_t^(-1);
        variables_to_save.joint_labels = joint_torque_labels;
        variables_to_save.joint_directions = joint_directions;

        save_folder = 'processed';
        save_file_name = makeFileName(collection_date, subject_id, trial_type, trial_number, 'jointTorques.mat');
        save([save_folder filesep save_file_name], '-struct', 'variables_to_save');
        addAvailableData ...
          ( ...
            'joint_torque_trajectories', ...
            'time', ...
            'sampling_rate', ...
            '_joint_labels', ...
            '_joint_directions', ...
            save_folder, ...
            save_file_name ...
          );
    end
    
    if any(strcmp(data_to_import, 'marker'))
        model_results.left_ankle = spline(model_results_raw.tout', model_results_raw.AnklePosL', time)';
        model_results.right_ankle = spline(model_results_raw.tout', model_results_raw.AnklePosR', time)';
        
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
    end
    
    
    





end













