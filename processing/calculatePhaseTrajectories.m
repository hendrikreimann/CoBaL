function calculatePhaseTrajectories(varargin)
    % parse arguments
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    
    load('subjectInfo.mat', 'date', 'subject_id');
     % load settings
    study_settings_file = '';
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    subject_settings = SettingsCustodian('subjectSettings.txt');
   
    arm_angle_method = study_settings.get('arm_angle');
    
    % calculate variables that depend upon the step events to be identified correctly
    stretch_variables = study_settings.get('stretch_variables');
    variables_to_save = struct;
    variables_to_prune_for = {}; 
   
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            save_folder = 'processed';
            save_file_name = makeFileName(date, subject_id, condition_list{i_condition}, i_trial, 'phaseTrajectories.mat');
                        
            % marker data
            [marker_trajectories, time_mocap, sampling_rate_mocap, marker_labels, marker_directions] = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'marker_trajectories');
            
            if strcmp('elbow', arm_angle_method)
                proximal_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'LELB', 'trajectories');
            elseif strcmp('shoulder', arm_angle_method)
                proximal_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'LSHO', 'trajectories');
            end
            
            LWRA_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'LWRA', 'trajectories');
            LWRB_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'LWRB', 'trajectories');
            left_wrist_center_trajectory = (LWRA_trajectory + LWRB_trajectory) * 0.5;
            
            % calculate vectors
            
            left_arm_vector_trajectory = proximal_trajectory - left_wrist_center_trajectory;

            % calculate angles
            left_arm_angle_ap = rad2deg(atan2(-left_arm_vector_trajectory(:, 2), left_arm_vector_trajectory(:, 3)));

            % find negative peaks
            [~, left_arm_peak_locations] = findpeaks(-left_arm_angle_ap, 'MinPeakProminence', subject_settings.get('left_armswing_peak_prominence_threshold'), 'MinPeakDistance', subject_settings.get('left_armswing_peak_distance_threshold') * sampling_rate_mocap);

            % normalize
            [larm_angle_ap_normalized, larm_angle_ap_dot_normalized] = normalizePeriodicVariable(left_arm_angle_ap, time_mocap, left_arm_peak_locations);

            % calculate phase
            left_arm_phase = atan2(larm_angle_ap_dot_normalized, -larm_angle_ap_normalized);

            RANK_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'RANK', 'trajectories');
            RPSI_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'RPSI', 'trajectories');
            RASI_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'RASI', 'trajectories');
            
            % calculate vectors
            right_pelvis_center_trajectory = (RPSI_trajectory + RASI_trajectory) * 0.5;
            right_leg_vector_trajectory = right_pelvis_center_trajectory - RANK_trajectory;

            % calculate angles
            right_leg_angle = rad2deg(atan2(-right_leg_vector_trajectory(:, 2), right_leg_vector_trajectory(:, 3)));

            % find negative peaks
            [~, right_leg_peak_locations] = findpeaks(-right_leg_angle, 'MinPeakProminence', subject_settings.get('right_legswing_peak_prominence_threshold'), 'MinPeakDistance', subject_settings.get('right_legswing_peak_distance_threshold') * sampling_rate_mocap);

            % normalize
            [rleg_angle_ap_normalized, rleg_angle_ap_dot_normalized] = normalizePeriodicVariable(right_leg_angle, time_mocap, right_leg_peak_locations);

            % calculate phase
            right_leg_phase = atan2(rleg_angle_ap_dot_normalized, -rleg_angle_ap_normalized);
                     
            if strcmp('elbow', arm_angle_method)
                proximal_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'RELB', 'trajectories');
            elseif strcmp('shoulder', arm_angle_method)
                proximal_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'RSHO', 'trajectories');
            end
            
            RWRA_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'RWRA', 'trajectories');
            RWRB_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'RWRB', 'trajectories');
            right_wrist_center_trajectory = (RWRA_trajectory + RWRB_trajectory) * 0.5;
            % calculate vectors
            
            right_arm_vector_trajectory = proximal_trajectory - right_wrist_center_trajectory;

            % calculate angles
            right_arm_angle_ap = rad2deg(atan2(-right_arm_vector_trajectory(:, 2), right_arm_vector_trajectory(:, 3)));

            % find negative peaks
            [~, right_arm_peak_locations] = ...
                findpeaks ...
                  ( ...
                    -right_arm_angle_ap, ...
                    'MinPeakProminence', subject_settings.get('right_armswing_peak_prominence_threshold'), ...
                    'MinPeakDistance', subject_settings.get('right_armswing_peak_distance_threshold') * sampling_rate_mocap...
                  );

            % normalize
            [rarm_angle_ap_normalized, rarm_angle_ap_dot_normalized] = normalizePeriodicVariable(right_arm_angle_ap, time_mocap, right_arm_peak_locations);

            % calculate phase
            right_arm_phase = atan2(rarm_angle_ap_dot_normalized, -rarm_angle_ap_normalized);

            LANK_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'LANK', 'trajectories');
            LPSI_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'LPSI', 'trajectories');
            LASI_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'LASI', 'trajectories');

            % calculate vectors
            left_pelvis_center_trajectory = (LPSI_trajectory + LASI_trajectory) * 0.5;
            left_leg_vector_trajectory = left_pelvis_center_trajectory - LANK_trajectory;

            % calculate angles
            left_leg_angle_ap = rad2deg(atan2(-left_leg_vector_trajectory(:, 2), left_leg_vector_trajectory(:, 3)));

            % find negative peaks
            [~, left_leg_peak_locations] = ...
                findpeaks ...
                  ( ...
                    -left_leg_angle_ap, ...
                    'MinPeakProminence', subject_settings.get('left_legswing_peak_prominence_threshold'), ...
                    'MinPeakDistance', subject_settings.get('left_legswing_peak_distance_threshold') * sampling_rate_mocap ...
                  );

            % normalize
            [lleg_angle_ap_normalized, lleg_angle_ap_dot_normalized] = normalizePeriodicVariable(left_leg_angle_ap, time_mocap, left_leg_peak_locations);

            % calculate phase
            left_leg_phase = atan2(-lleg_angle_ap_dot_normalized, lleg_angle_ap_normalized);
     
            % set ignore markers for NaNs
            if strcmp(subject_id, '101')
                right_arm_phase = zeros(length(right_arm_angle_ap),1);               
                left_arm_phase = zeros(length(right_arm_angle_ap),1);             
                left_leg_phase = zeros(length(right_arm_angle_ap),1);             
                right_leg_phase = zeros(length(right_arm_angle_ap),1);
            end
            
            % add new variables to be saved
            variables_to_save.left_arm_angle_ap = left_arm_angle_ap;
            variables_to_save.left_arm_phase = left_arm_phase;
            variables_to_save.right_leg_angle_ap = right_leg_angle;
            variables_to_save.right_leg_phase = right_leg_phase;
            addAvailableData_new('left_arm_angle_ap', 'time_mocap', 'sampling_rate_mocap', 'Left Arm Angle (A/P)', {'+', '-'}, save_folder, save_file_name);
            addAvailableData_new('right_leg_angle_ap', 'time_mocap', 'sampling_rate_mocap', 'Right Leg Angle (A/P)', {'+', '-'}, save_folder, save_file_name);
            addAvailableData_new('left_arm_phase', 'time_mocap', 'sampling_rate_mocap', 'Left Arm Phase', {'+', '-'}, save_folder, save_file_name);
            addAvailableData_new('right_leg_phase', 'time_mocap', 'sampling_rate_mocap', 'Right Leg Phase', {'+', '-'}, save_folder, save_file_name);
            variables_to_save.right_arm_angle_ap = right_arm_angle_ap;
            variables_to_save.right_arm_phase = right_arm_phase;
            variables_to_save.left_leg_angle_ap = left_leg_angle_ap;
            variables_to_save.left_leg_phase = left_leg_phase;
            addAvailableData_new('right_arm_angle_ap', 'time_mocap', 'sampling_rate_mocap', 'Right Arm Angle (A/P)', {'+', '-'}, save_folder, save_file_name);
            addAvailableData_new('right_arm_phase', 'time_mocap', 'sampling_rate_mocap', 'Right Arm Phase', {'+', '-'}, save_folder, save_file_name);
            addAvailableData_new('left_leg_angle_ap', 'time_mocap', 'sampling_rate_mocap', 'Left Leg Angle (A/P)', {'+', '-'}, save_folder, save_file_name);
            addAvailableData_new('left_leg_phase', 'time_mocap', 'sampling_rate_mocap', 'Left Leg Phase', {'+', '-'}, save_folder, save_file_name);

            variables_to_save.sampling_rate_mocap = sampling_rate_mocap;
            variables_to_save.time_mocap = time_mocap; 
            saveDataToFile([save_folder filesep save_file_name], variables_to_save);
                    
            if any(isnan(left_leg_phase)) || any(isnan(left_leg_angle_ap)) || any(isnan(left_arm_phase)) || any(isnan(left_arm_angle_ap)) || any(isnan(right_leg_phase)) || any(isnan(right_leg_angle_ap)) || any(isnan(right_arm_phase)) || any(isnan(right_arm_angle_ap))
                left_leg_ignore_indices = find(isnan(left_leg_phase));
                left_arm_ignore_indices = find(isnan(left_arm_phase));
                right_leg_ignore_indices = find(isnan(right_leg_phase));
                right_arm_ignore_indices = find(isnan(right_arm_phase));
                all_ignore_indices = unique([left_leg_ignore_indices; right_leg_ignore_indices; left_arm_ignore_indices; right_arm_ignore_indices]);

                stepEvent_file_name = makeFileName(date, subject_id, condition_list{i_condition}, i_trial, 'stepEvents.mat');
                % just to grab ignore_times.. is there a way to specifiy
                % this variable and no others?
                load(['analysis' filesep stepEvent_file_name]);

                if ~exist('ignore_times', 'var')
                   ignore_times = [];
                end

                stretches_to_ignore = struct;
                this_ignore_times = time_mocap(all_ignore_indices);
                ignore_times = [ignore_times; this_ignore_times];
                stretches_to_ignore.ignore_times = sort(ignore_times);
                saveDataToFile(['analysis' filesep stepEvent_file_name], stretches_to_ignore)
                findRelevantDataStretches('condition',condition_list{i_condition}, 'trial', i_trial)
                clear('ignore_times');
            end
           
           disp(['Finding Phase Trajectories: Trial ', num2str(i_trial), ' completed',  ' relevant stretches, saved as ', save_file_name]);                
        end
    end
end