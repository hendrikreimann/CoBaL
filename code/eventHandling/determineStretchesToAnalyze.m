%     This file is part of the CoBaL code base
%     Copyright (C) 2017 Hendrik Reimann <hendrikreimann@gmail.com>
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% This script finds the stretches of data relevant to the experimental paradigm and determines the condition factors.

% Stretches are further subdivided into bands, e.g. double stance and single stance. Stretches are defined by the times
% where each band starts and ends, i.e. N+1 numbers per stretch, where N is the number of bands per stretch. These are
% stored in stretch_times, an array with one row per stretch and N columns.
% Conditions are saved in a struct conditions_trial, with one field for each factor. Factors are column cell arrays with
% one entry per stretch. Stance foot information is saved in stance_foot_data, an array with one row per stretch and one
% column per band, where each entry specifies which foot is on the ground for that band. 

% pushoff_times is a legacy variable which might disappear soon

% stimulus_condition dictates the appropriate relevant data stretches

% input: 
% - subjectInfo.mat
% - optional: conditions.mat
%
% output:
% file relevantDataStretches.mat, containing
% 
% - stretch_pushoff_times
% - stretch_times
% - stance_foot_data
% - conditions_trial
% - bands_per_stretch


function determineStretchesToAnalyze(varargin)
    % parse arguments
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});

    
    variables_to_save_latency_check = struct;
    latency_value_list = {};
    
    %% prepare
%     load('subjectInfo.mat', 'date', 'subject_id', 'most_affected', 'gender');
    % load settings
    study_settings_file = '';
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    experimental_paradigm = study_settings.get('experimental_paradigm', 1);

    subject_settings = SettingsCustodian('subjectSettings.txt');
    collection_date = num2str(subject_settings.get('collection_date'));
    gender = subject_settings.get('gender');
    subject_id = subject_settings.get('subject_id');
    
    conditions_file_name = [];
    if exist('conditions.csv', 'file')
        conditions_file_name = 'conditions.csv';
    end
    if exist(makeFileName(collection_date, subject_id, 'conditions.csv'), 'file')
        conditions_file_name = makeFileName(collection_date, subject_id, 'conditions.csv');
    end
    
    if strcmp(experimental_paradigm, 'CadenceVision') || strcmp(experimental_paradigm, 'CadenceGVS')
        protocol_data = load('protocolInfo.mat');
    end
    if strcmp(experimental_paradigm, 'FatigueGVS')
        fatigue_trials = subject_settings.get('fatigue_trials');
    end
    if strcmp(experimental_paradigm, 'CognitiveLoadVision') || strcmp(experimental_paradigm, 'CognitiveLoadGvs')
        back_7_trials = subject_settings.get('back_7_trials');
    end
    

    %% process
    time_to_nearest_heelstrike_before_trigger_threshold = 0.10; % a heelstrike should happen less than this long before a trigger
    time_to_nearest_heelstrike_after_trigger_threshold = 0.3; % a heelstrike should happen less than this long after a trigger
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            %% load data
            ignore_times = [];
            load(['analysis' filesep makeFileName(collection_date, subject_id, condition_list{i_condition}, i_trial, 'events')]);
%             load(['analysis' filesep makeFileName(date, subject_id, condition_list{i_condition}, i_trial, 'stepEvents')]);
%             load(['processed' filesep makeFileName(date, subject_id, condition_list{i_condition}, i_trial, 'kinematicTrajectories')]);
            
            % determine experimental condition
            this_trial_type = condition_list{i_condition};
            condition_experimental = study_settings.get('experimental_condition');
            if strcmp(condition_experimental, 'load_from_conditions_file')
                condition_experimental = loadConditionFromFile(conditions_file_name, 'experimental', i_trial);
            end
            if strcmp(condition_experimental, 'determine_from_file_name')
                condition_experimental = condition_list{i_condition};
            end
            if strcmp(condition_experimental, 'determine_from_type_day_combination')
                % do nothing, we'll deal with this depending on experiment type
            end
            if strcmp(experimental_paradigm, 'platformShift')            
                [condition_table, condition_header] = loadConditionTableFromFile(conditions_file_name);
                this_trial_line_indicator = strcmp(num2str(i_trial), condition_table(:, strcmp(condition_header, 'trial')));
                condition_experimental = condition_table{this_trial_line_indicator, strcmp(condition_header, 'condition')};
            end
            
            % determine stimulus type
            condition_stimulus = study_settings.get('stimulus_condition');
            if strcmp(condition_stimulus, 'load_from_conditions_file')
                condition_stimulus = loadConditionFromFile(conditions_file_name, 'stimulus', i_trial);
            end
            if strcmp(condition_stimulus, 'determine_from_file_name')
                condition_stimulus = condition_list{i_condition};
            end
            
            % determine day
            condition_day = study_settings.get('day_condition', 1);
            if strcmp(condition_day, 'load_from_conditions_file')
                condition_day = loadConditionFromFile(conditions_file_name, 'day', i_trial);
            end
            if strcmp(condition_day, 'determine_from_file_name')
                condition_day = condition_list{i_condition};
            end
            if strcmp(condition_day, 'determine_from_subject_settings')
                condition_day = subject_settings.get('session_label');
            end
            
            % marker data
            [marker_trajectories, time_marker, sampling_rate_marker, marker_labels, marker_directions] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'marker_trajectories');
            if study_settings.get('prune_gaps_com', 1)
                [com_trajectories, time_marker, sampling_rate_marker, com_labels, com_directions] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'com_trajectories');
            
%                 [com_trajectories, time_marker, sampling_rate_marker, com_labels, com_directions] = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'com_trajectories_optimized')
            end
            if study_settings.get('prune_gaps_angles', 1)
                [joint_angle_trajectories, time_marker, sampling_rate_marker, joint_labels, joint_directions] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'joint_angle_trajectories');
%                 [joint_angle_trajectories, time_marker, sampling_rate_marker, com_labels, joint_directions] = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'joint_trajectories_optimized')
            
            end
            
            % forceplate data
            [left_forceplate_cop_world_trajectory, time_left_forceplate, ~, ~, ~, left_forceplate_available] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'left_foot_cop_world', 'optional');
            [right_forceplate_cop_world_trajectory, time_right_forceplate, ~, ~, ~, right_forceplate_available] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'right_foot_cop_world', 'optional');
            [cop_world_trajectory, time_forceplate, ~, ~, cop_available] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'total_forceplate_cop_world', 'optional');
            if left_forceplate_available && right_forceplate_available
                left_copx_trajectory = left_forceplate_cop_world_trajectory(:, 1);
                right_copx_trajectory = right_forceplate_cop_world_trajectory(:, 1);
                
                if any(strcmp(condition_list{i_condition}, study_settings.get('trial_types_with_inverted_forceplate_sides', 1)))
                    % this was for the VEPO lab, where the sides are inverted
                    left_copx_trajectory = right_forceplate_cop_world_trajectory(:, 3);
                    right_copx_trajectory = left_forceplate_cop_world_trajectory(:, 3);
                end
            end
            
            % stimulus data
            if strcmp(experimental_paradigm, 'GVS_old')
                load(['processed' filesep makeFileName(collection_date, subject_id, condition_list{i_condition}, i_trial, 'labviewData')]);
%                 GVS_out_trajectory = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'GVS_out_trajectory');
                GVS_stim_trajectory = GVS_out_trajectory + subject_settings.get('gvs_offset');
%                 [stimulus_state_trajectory, time_stimulus] = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'stimulus_state_trajectory');
                time_stimulus = time;
            end
            if strcmp(condition_stimulus, 'VISUAL')
                % this if for TU data
                visual_scene_ml_translation_trajectory = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'visual_scene_ml_translation__trajectory'); %take note of the double "_"
                [stimulus_state_trajectory, time_stimulus] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'stimulus_state_trajectory');
            end
            if strcmp(experimental_paradigm, 'Vision') || strcmp(experimental_paradigm, 'CadenceVision') || strcmp(experimental_paradigm, 'CognitiveLoadVision')
%                 current_rotation_trajectory = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'current_rotation_trajectory');
                current_rotation_trajectory = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'visual_rotation_angle_trajectory');
                current_acceleration_trajectory = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'visual_rotation_acceleration_trajectory');
                [stimulus_state_trajectory, time_stimulus] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'stimulus_state_trajectory');
            end
            if strcmp(experimental_paradigm, 'GVS') || strcmp(experimental_paradigm, 'CadenceGVS') || strcmp(experimental_paradigm, 'FatigueGVS') || strcmp(experimental_paradigm, 'CognitiveLoadGvs')
                gvs_trajectory = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'GVS_current_trajectory');
                [stimulus_state_trajectory, time_stimulus] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'stimulus_state_trajectory');
            end
            if strcmp(experimental_paradigm, 'GvsOverground')
                [analog_trajectories, time_analog, sampling_rate_analog, analog_labels, analog_directions] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'analog_trajectories');
                gvs_trajectory = analog_trajectories(:, strcmp(analog_labels, 'GVS_out'));
            end
            if strcmp(experimental_paradigm, 'OculusLaneRestriction')
                gvs_trajectory = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'GVS_current_trajectory');%'GVS_current_trajectory');
                [stimulus_state_trajectory, time_stimulus] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'stimulus_state_trajectory');
                scene_translation_trajectory = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'SceneTranslation_trajectory');
                load('virtualobjectInfo');
            end


            % determine indices for optional markers
            marker_weight_table = study_settings.get('marker_weights');
            if isempty(marker_weight_table)
                optional_marker_indices = [];
            else
                optional_marker_list = marker_weight_table(:, 1);
                optional_marker_indices = [];
                for i_marker = 1 : length(optional_marker_list)
                    marker_indices = extractMarkerData(marker_trajectories, marker_labels, optional_marker_list{i_marker}, 'indices');
                    optional_marker_indices = [optional_marker_indices marker_indices];
                end
            end
            essential_marker_indicator = ~ismember(1 : size(marker_trajectories, 2), optional_marker_indices);
            
            % determine illusion
            if strcmp(experimental_paradigm, 'GVS_old')
                illusion_trajectory = zeros(size(time_stimulus)); % 1 = RIGHT, -1 = LEFT
                % use GVS_out_trajectory as illusion
                for i_time = 1 : length(time_stimulus)
                    if GVS_stim_trajectory(i_time) > 0
                        % anode is on the right, cathode is on the left, illusory fall is towards the cathode, i.e. LEFT
                        illusion_trajectory(i_time) = -1;
                    end
                    if GVS_stim_trajectory(i_time) < 0
                        % anode is on the left, cathode is on the right, illusory fall is towards the cathode, i.e. RIGHT
                        illusion_trajectory(i_time) = 1;
                    end
                end
            end
            if strcmp(condition_stimulus, 'VISUAL')
                illusion_trajectory = zeros(size(time_stimulus)); % 1 = RIGHT, -1 = LEFT
                for i_time = 1 : length(time_stimulus)
                    if visual_scene_ml_translation_trajectory(i_time) > 0
                        % angle change is positive horizon rotates counter-clockwise, illusion is to the RIGHT
                        illusion_trajectory(i_time) = 1;
                    end
                    if visual_scene_ml_translation_trajectory(i_time) < 0 
                        % angle change is negative, horizon rotates clockwise, illusion is to the LEFT
                        illusion_trajectory(i_time) = -1;
                        
                        
                    end
                end
            end
            if strcmp(experimental_paradigm, 'Vision') || strcmp(experimental_paradigm, 'CadenceVision') || strcmp(experimental_paradigm, 'CognitiveLoadVision')
                illusion_trajectory = zeros(size(time_stimulus)); % -1 = LEFT, 1 = RIGHT
                for i_time = 1 : length(time_stimulus)
                    if stimulus_state_trajectory(i_time) == 3
                        % stimulus is currently active
                        if current_rotation_trajectory(i_time) < 0
                            % angle is negative, horizon rotates clockwise, illusion is to the LEFT
                            illusion_trajectory(i_time) = -1;
                        end
                        if current_rotation_trajectory(i_time) > 0
                            % angle is positive, horizon rotates counter-clockwise, illusion is to the RIGHT
                            illusion_trajectory(i_time) = 1;
                        end
                    end
                end
            end
            if strcmp(experimental_paradigm, 'GVS') || strcmp(experimental_paradigm, 'CadenceGVS') || strcmp(experimental_paradigm, 'FatigueGVS') || strcmp(experimental_paradigm, 'CognitiveLoadGvs')
                illusion_trajectory = zeros(size(time_stimulus)); % -1 = LEFT, 1 = RIGHT
                for i_time = 1 : length(time_stimulus)
                    if stimulus_state_trajectory(i_time) == 3
                        % stimulus is currently active
                        if gvs_trajectory(i_time) < 0
                            % negative current = anode on the left = illusory fall to the right
                            illusion_trajectory(i_time) = 1;
                        end
                        if gvs_trajectory(i_time) > 0
                            % positive current = anode on the right = illusory fall to the left
                            illusion_trajectory(i_time) = -1;
                        end
                    end
                end
            end
            if strcmp(experimental_paradigm, 'GvsOverground')
                gvs_threshold = study_settings.get('gvs_threshold');
                illusion_trajectory = zeros(size(time_analog)); % -1 = LEFT, 1 = RIGHT
                
                % positive current = anode on the right = illusory fall to the left
                illusion_trajectory(gvs_trajectory > gvs_threshold) = -1;
                % negative current = anode on the left = illusory fall to the right
                illusion_trajectory(gvs_trajectory < -gvs_threshold) = 1;
            end
            if strcmp(experimental_paradigm, 'OculusLaneRestriction')
                illusion_trajectory = zeros(size(time_stimulus)); % -1 = LEFT, 1 = RIGHT
                for i_time = 1 : length(time_stimulus)
                    if stimulus_state_trajectory(i_time) == 5
                        % stimulus is currently active
                        if gvs_trajectory(i_time) < 0
                            % negative current = anode on the left = illusory fall to the right
                            illusion_trajectory(i_time) = 1;
                        end
                        if gvs_trajectory(i_time) > 0
                            % positive current = anode on the right = illusory fall to the left
                            illusion_trajectory(i_time) = -1;
                        end
                    end
                end
            end
            
            %% extract events
            if strcmp(condition_stimulus, 'NONE') ...
                    || strcmp(condition_stimulus, 'VISUAL') || strcmp(experimental_paradigm, 'GVS_old') ...
                    || strcmp(experimental_paradigm, 'Vision') || strcmp(experimental_paradigm, 'CadenceVision') || strcmp(experimental_paradigm, 'CognitiveLoadVision') ...
                    || strcmp(experimental_paradigm, 'GVS') || strcmp(experimental_paradigm, 'CadenceGVS') || strcmp(experimental_paradigm, 'FatigueGVS') || strcmp(experimental_paradigm, 'CognitiveLoadGvs')  ...
                    || strcmp(condition_stimulus, 'OBSTACLE') || strcmp(condition_stimulus, 'ARMSENSE') ...
                    || strcmp(experimental_paradigm, 'Vision Stochastic') || strcmp(experimental_paradigm, 'GvsOverground') ...
                    || strcmp(experimental_paradigm, 'OculusLaneRestriction')
                right_pushoff_times = event_data{strcmp(event_labels, 'right_pushoff')};
                right_touchdown_times = event_data{strcmp(event_labels, 'right_touchdown')};
                left_pushoff_times = event_data{strcmp(event_labels, 'left_pushoff')};
                left_touchdown_times = event_data{strcmp(event_labels, 'left_touchdown')};
            end
            if strcmp(experimental_paradigm, 'platformShift')
                perturbation_start_times = event_data{strcmp(event_labels, 'perturbation_start')};
                perturbation_end_times = event_data{strcmp(event_labels, 'perturbation_end')};
                response_end_times = event_data{strcmp(event_labels, 'perturbation_end_plus_two')};
            end
            
            %% find triggers
            %
            % Find the triggering events that indicate a stretch of interest. For perturbation experiments, this is the onset of
            % a perturbation. For unperturbed walking, this is any heelstrike.
            % The result is trigger_indices_labview.
            %
            if strcmp(condition_stimulus, 'NONE')
                % use all touchdown events as triggers
                trigger_times = [left_touchdown_times; right_touchdown_times];
            end
            if strcmp(condition_stimulus, 'VISUAL') || strcmp(experimental_paradigm, 'GVS_old')
                % find the time steps where the stimulus state crosses a threshold
                stimulus_threshold = 0.5;
                trigger_indices_labview = find(diff(sign(stimulus_state_trajectory - stimulus_threshold)) > 0) + 1;

                %
                epsilon = 1e-5;
                % remove weird noise in illusion trajectory (check labview
                % for odd behavior i.e. wait time and illusion_trajectory)
                stim_start_indices_labview = find(diff(sign(abs(illusion_trajectory) - epsilon)) > 0) + 1;
                i_stim = 1;
                while i_stim ~= length(stim_start_indices_labview)
                    if illusion_trajectory(stim_start_indices_labview(i_stim) + 5) == 0
                        stim_start_indices_labview(i_stim) = [];
                        i_stim = i_stim - 1;
                    end
                    i_stim = i_stim + 1;
                end
                trigger_indices_labview = trigger_indices_labview(1 : length(stim_start_indices_labview)); % in case a stim is triggered, but not recorded
                
                trigger_times = time_stimulus(trigger_indices_labview);
            end
            if strcmp(experimental_paradigm, 'Vision') || strcmp(experimental_paradigm, 'CadenceVision') || strcmp(experimental_paradigm, 'CognitiveLoadVision')
                % find the time steps where the stimulus state crosses a threshold
                stimulus_threshold = 1.5;
                trigger_indices_labview = find(diff(sign(stimulus_state_trajectory - stimulus_threshold)) > 0) + 2;
                trigger_times = time_stimulus(trigger_indices_labview);
            end
            if strcmp(experimental_paradigm, 'GVS') || strcmp(experimental_paradigm, 'CadenceGVS') || strcmp(experimental_paradigm, 'FatigueGVS') || strcmp(experimental_paradigm, 'CognitiveLoadGvs')
                % find the time steps where the stimulus state crosses a threshold
                stimulus_threshold = 1.5;
                trigger_indices_labview = find(diff(sign(stimulus_state_trajectory - stimulus_threshold)) > 0) + 2;
                trigger_times = time_stimulus(trigger_indices_labview);
            end
            if strcmp(experimental_paradigm, 'GvsOverground')
                % find the time steps where the first forceplate vertical force crosses a threshold
                stimulus_threshold = 20;
                [frist_forceplate_wrench_trajectory, time_forceplate] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'left_foot_wrench_world');
                vertical_force_trajectory = frist_forceplate_wrench_trajectory(:, 3);
                
                trigger_indices_forceplate = find(diff(sign(-vertical_force_trajectory - stimulus_threshold)) > 0) + 2;
                trigger_times = time_forceplate(trigger_indices_forceplate);
            end
            if strcmp(experimental_paradigm, 'OculusLaneRestriction') 
                % find the time steps where the stimulus state crosses a threshold
                stimulus_threshold = 4.5;
                trigger_indices_labview = find(diff(sign(stimulus_state_trajectory - stimulus_threshold)) > 0) + 2;
                trigger_times = time_stimulus(trigger_indices_labview);
            end
            if strcmp(condition_stimulus, 'OBSTACLE')
                trigger_times = [];
            end
            if strcmp(condition_stimulus, 'ARMSENSE')
                 trigger_times = [];
            end
            if strcmp(experimental_paradigm, 'Vision Stochastic')
                trigger_times = [];
            end
            if strcmp(experimental_paradigm, 'platformShift')
                % use all touchdown events as triggers
                trigger_times = perturbation_start_times;
            end
            
            % calculate indices
            trigger_indices_mocap = zeros(size(trigger_times));
            for i_index = 1 : length(trigger_times)
                [~, index_mocap] = min(abs(time_marker - trigger_times(i_index)));
                trigger_indices_mocap(i_index) = index_mocap;
            end

            % visualize triggers
            if visualize
                LHEE = extractMarkerData(marker_trajectories, marker_labels, 'LHEE');
                RHEE = extractMarkerData(marker_trajectories, marker_labels, 'RHEE');
                [left_forceplate_wrench_world_trajectory, time_left_forceplate] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'left_foot_wrench_world', 'optional');
                [right_forceplate_wrench_world_trajectory, time_right_forceplate] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'right_foot_wrench_world', 'optional');
                
%                 figure; axes; hold on
%                 plot(time_stimulus, stimulus_state_trajectory*0.02);
%                 plot(time_marker, LHEE(:, 3), 'Displayname', 'left heel z')
%                 plot(time_marker, RHEE(:, 3), 'Displayname', 'right heel z')
%                 
%                 plot(time_marker, LHEE(:, 2), 'Displayname', 'left heel y')
%                 plot(time_marker, RHEE(:, 2), 'Displayname', 'right heel y')
                
                
%                 if left_forceplate_available
%                     left_cop_x_trajectory_relevant = left_copx_trajectory; left_cop_x_trajectory_relevant(left_cop_x_trajectory_relevant==0) = NaN;
%                     plot(time_left_forceplate, left_cop_x_trajectory_relevant, 'linewidth', 2, 'Displayname', 'cop left');
%                 end
%                 if right_forceplate_available
%                     right_cop_x_trajectory_relevant = right_copx_trajectory; right_cop_x_trajectory_relevant(right_cop_x_trajectory_relevant==0) = NaN;
%                     plot(time_right_forceplate, right_cop_x_trajectory_relevant, 'linewidth', 2, 'Displayname', 'cop right');
%                 end
%                 plot(time_marker(trigger_indices_mocap), zeros(size(trigger_indices_mocap)), 'x', 'Displayname', 'triggers')
%                 plot(time_stimulus, illusion_trajectory, 'Displayname', 'illusion')
%                 legend('stimulus state', 'left cop', 'right cop', 'left touchdown', 'right touchdown', 'trigger', 'stim start')
%                 legend('toggle')
            end

            %% extract data and determine condition variables

            % For each trigger, determine the conditions and the relevant step events.
            number_of_triggers = length(trigger_times);
            removal_flags = zeros(number_of_triggers, 1);
            event_variables_to_save = struct;
            if strcmp(condition_stimulus, 'NONE')
                % determine start and end
                stance_foot_data = {'STANCE_RIGHT', 'STANCE_BOTH', 'STANCE_LEFT'};
                bands_per_stretch = length(stance_foot_data);
                
                stretch_start_times = zeros(number_of_triggers, 1);
                stretch_end_times = zeros(number_of_triggers, 1);
                stretch_pushoff_times = zeros(number_of_triggers, 1);
                closest_heelstrike_distance_times = zeros(number_of_triggers, 1);
                condition_stance_foot_list = cell(number_of_triggers, 1);
                perturbation_list = cell(number_of_triggers, 1);
                condition_delay_list = cell(number_of_triggers, 1);
                condition_index_list = cell(number_of_triggers, 1);
                condition_experimental_list = cell(number_of_triggers, 1);
                condition_stimulus_list = cell(number_of_triggers, 1);
                condition_day_list = cell(number_of_triggers, 1);
                
                for i_trigger = 1 : number_of_triggers
                    perturbation_list{i_trigger, 1} = 'N/A';
                    condition_delay_list{i_trigger, 1} = 'N/A';
                    condition_experimental_list{i_trigger, 1} = condition_experimental;
                    condition_stimulus_list{i_trigger, 1} = condition_stimulus;
                    condition_day_list{i_trigger, 1} = condition_day;
                    
                
                    % find out which heelstrike triggered
                    % XXX change this to use the interval, same as in the stimulus case
                    [distance_to_trigger_left_time, index_left] = min(abs(left_touchdown_times - trigger_times(i_trigger)));
                    [distance_to_trigger_right_time, index_right] = min(abs(right_touchdown_times - trigger_times(i_trigger)));
                    closest_heelstrike_distance_time = min([distance_to_trigger_left_time distance_to_trigger_right_time]);
                    
                    if distance_to_trigger_left_time < distance_to_trigger_right_time
                        condition_stance_foot_list{i_trigger, 1} = 'STANCE_LEFT';
                        condition_index_list{i_trigger, 1} = 'ONE';
                        stretch_start_times(i_trigger, 1) = trigger_times(i_trigger);
                        stretch_end_time_index = find(right_touchdown_times > trigger_times(i_trigger), 1, 'first');
                        stretch_pushoff_time_index = find(right_pushoff_times > trigger_times(i_trigger), 1, 'first');
                        if isempty(stretch_end_time_index) || isempty(stretch_pushoff_time_index)
                            stretch_end_times(i_trigger, 1) = -1;
                            stretch_pushoff_times(i_trigger, 1) = -1;
                            removal_flags(i_trigger) = 1;
                        else
                            stretch_end_times(i_trigger, 1) = right_touchdown_times(stretch_end_time_index);
                            stretch_pushoff_times(i_trigger, 1) = right_pushoff_times(stretch_pushoff_time_index);
                        end
                    end    
                    if distance_to_trigger_right_time < distance_to_trigger_left_time
                        condition_stance_foot_list{i_trigger, 1} = 'STANCE_RIGHT';
                        condition_index_list{i_trigger, 1} = 'TWO';
                        stretch_start_times(i_trigger, 1) = trigger_times(i_trigger);
                        stretch_end_time_index = find(left_touchdown_times > trigger_times(i_trigger), 1, 'first');
                        stretch_pushoff_time_index = find(left_pushoff_times > trigger_times(i_trigger), 1, 'first');
                        if isempty(stretch_end_time_index) || isempty(stretch_pushoff_time_index)
                            stretch_end_times(i_trigger, 1) = -1;
                            stretch_pushoff_times(i_trigger, 1) = -1;
                            removal_flags(i_trigger) = 1;
                        else
                            stretch_end_times(i_trigger, 1) = left_touchdown_times(stretch_end_time_index);
                            stretch_pushoff_times(i_trigger, 1) = left_pushoff_times(stretch_pushoff_time_index);
                        end
                    end
                   
                end
                
                % remove flagged triggers
                unflagged_indices = ~removal_flags;
                trigger_times = trigger_times(unflagged_indices);
                stretch_start_times = stretch_start_times(unflagged_indices, :);
                stretch_end_times = stretch_end_times(unflagged_indices, :);
                stretch_pushoff_times = stretch_pushoff_times(unflagged_indices, :);
                condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
                perturbation_list = perturbation_list(unflagged_indices, :);
                condition_delay_list = condition_delay_list(unflagged_indices, :);
                condition_index_list = condition_index_list(unflagged_indices, :);
                condition_experimental_list = condition_experimental_list(unflagged_indices, :);
                condition_stimulus_list = condition_stimulus_list(unflagged_indices, :);
                condition_day_list = condition_day_list(unflagged_indices, :);
                closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :); 
                
                if visualize
                    for i_trigger = 1 : length(stretch_start_times)
                        if strcmp(condition_stance_foot_list(i_trigger), 'STANCE_RIGHT')
                            stretch_indicator_height = 0.01;
                        else
                            stretch_indicator_height = -0.01;
                        end
                            
                        plot([stretch_start_times(i_trigger) stretch_end_times(i_trigger)], [1 1]*stretch_indicator_height, 'linewidth', 3);
                    end
                end
                     
                stretch_times = [stretch_start_times stretch_end_times];

                conditions_trial = struct;
                conditions_trial.condition_stance_foot_list = condition_stance_foot_list;
                conditions_trial.condition_perturbation_list = perturbation_list;
                conditions_trial.condition_delay_list = condition_delay_list;
                conditions_trial.condition_index_list = condition_index_list;
                conditions_trial.condition_experimental_list = condition_experimental_list;
                conditions_trial.condition_stimulus_list = condition_stimulus_list;
                conditions_trial.condition_day_list = condition_day_list;
                
                event_variables_to_save.stretch_start_times = stretch_start_times;
                event_variables_to_save.stretch_pushoff_times = stretch_pushoff_times;
                event_variables_to_save.stretch_end_times = stretch_end_times;
                event_variables_to_save.stretch_times = stretch_times;               
            end  
            
            if strcmp(condition_stimulus, 'OBSTACLE')
                % determine start and end
                stance_foot_data = {'STANCE_BOTH', 'STANCE_BOTH', 'STANCE_LEFT'};
                bands_per_stretch = length(stance_foot_data);
                
                init_time = right_pushoff_times(1) - 1; % assume heel-off happened at least one second before toes-off, so start looking at that point
                end_time = right_touchdown_times(1);
                unload_time = right_pushoff_times(1);
                
                % determine unload time as maximal backward-right shift of the CoP, following Halliday et al, Gait and Posture 8 (1998) 8?14
                [~, start_time_index_forceplate] = min(abs(time_forceplate - init_time));
                [~, unload_time_index_forceplate] = min(abs(time_forceplate - unload_time));
                cop_data_relevant = cop_world_trajectory(start_time_index_forceplate : unload_time_index_forceplate, :);
                time_forceplate_relevant = time_forceplate(start_time_index_forceplate : unload_time_index_forceplate);
                [~, release_time_index_forceplate] = max(cop_data_relevant(:, 1));
                release_time = time_forceplate_relevant(release_time_index_forceplate);
                start_time = release_time - 0.5;
                
                
                stretch_start_times = right_pushoff_times(1) - 1; % HR: this is probably not right anymore
                stretch_end_times = right_touchdown_times(1);
                stretch_pushoff_times = 0;
                condition_experimental_list = {condition_experimental};
                stretch_times = [start_time release_time unload_time end_time];
                
                if visualize
                    for i_trigger = 1 : length(stretch_start_times)
                        if strcmp(stance_foot_data(i_trigger), 'STANCE_LEFT')
                            stretch_indicator_height = 0.01;
                        else
                            stretch_indicator_height = -0.01;
                        end
                            
                        plot([stretch_start_times(i_trigger) stretch_end_times(i_trigger)], [1 1]*stretch_indicator_height, 'linewidth', 3);
                    end
                end
                
                % add new variables to be saved
                conditions_trial = struct;
                conditions_trial.condition_experimental_list = condition_experimental_list;
                event_variables_to_save.stretch_pushoff_times = stretch_pushoff_times;
                event_variables_to_save.stretch_times = stretch_times;
                event_variables_to_save.stance_foot_data = stance_foot_data;

            end

            if strcmp(experimental_paradigm, 'GvsOverground')
                % determine start and end
                stance_foot_data = {'STANCE_LEFT', 'STANCE_RIGHT', 'STANCE_LEFT', 'STANCE_RIGHT'};
                trigger_foot_list = cell(number_of_triggers, 1); % triggering foot TRIGGER_LEFT or TRIGGER_RIGHT
                bands_per_stretch = length(stance_foot_data);
                
                % find heelstrike that corresponds to the trigger
                if ~isempty(trigger_times) && length(left_touchdown_times) >= 3 && length(right_touchdown_times) >= 2
                    trigger_time = trigger_times(1);

                    % get closest heelstrike on either side
                    [~, index_left] = min(abs(left_touchdown_times - trigger_time));
                    [~, index_right] = min(abs(right_touchdown_times - trigger_time));

                    % is the closest left heelstrike within the acceptable interval?
                    closest_left_heelstrike = left_touchdown_times(index_left);
                    time_difference_left = closest_left_heelstrike - trigger_time; % where does the closest left heelstrike lie relative to the trigger?
                    if -time_to_nearest_heelstrike_before_trigger_threshold < time_difference_left && time_difference_left < time_to_nearest_heelstrike_after_trigger_threshold
                        % left heelstrike is acceptable
                        left_heelstrike_acceptable = true;
                    else
                        left_heelstrike_acceptable = false;
                    end

                    % is the closest right heelstrike within the acceptable interval?
                    closest_right_heelstrike = right_touchdown_times(index_right);
                    time_difference_right = closest_right_heelstrike - trigger_time; % where does the closest right heelstrike lie relative to the trigger?
                    if -time_to_nearest_heelstrike_before_trigger_threshold < time_difference_right && time_difference_right < time_to_nearest_heelstrike_after_trigger_threshold
                        % right heelstrike is acceptable
                        right_heelstrike_acceptable = true;
                    else
                        right_heelstrike_acceptable = false;
                    end

                    % flag for removal if not triggered by left foot
                    if left_heelstrike_acceptable && ~right_heelstrike_acceptable
                        trigger_foot_list = {'TRIGGER_LEFT'};
                    else
                        % not triggered by left heelstrike
                        removal_flags = 1;
                    end                  

                    % extract needed events
                    left_foot_heelstrike_0  = left_touchdown_times(index_left);
                    left_foot_heelstrike_1  = left_touchdown_times(index_left+1);
                    left_foot_heelstrike_2  = left_touchdown_times(index_left+2);
                    left_foot_pushoff_0     = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_0));
                    left_foot_pushoff_1     = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_1));
%                    left_foot_pushoff_2     = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_2));

                    right_foot_heelstrike_0 = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_0));
                    right_foot_heelstrike_1 = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_1));
%                     right_foot_heelstrike_2 = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_2));
                    right_foot_pushoff_0    = max(right_pushoff_times(right_pushoff_times <= left_foot_pushoff_0));
                    right_foot_pushoff_1    = max(right_pushoff_times(right_pushoff_times <= left_foot_pushoff_1));
%                     right_foot_pushoff_2    = max(right_pushoff_times(right_pushoff_times <= left_foot_pushoff_2));

                    % notify if events are not sorted properly
                    if ~issorted ...
                          ( ...
                            [ ...
                              left_foot_heelstrike_0 right_foot_pushoff_0 right_foot_heelstrike_0 left_foot_pushoff_0 ...
                              left_foot_heelstrike_1 right_foot_pushoff_1 right_foot_heelstrike_1 left_foot_pushoff_1 ...
                              left_foot_heelstrike_2 ...
                            ] ...
                          )
                        disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(trigger_time)]);
                    end
                    % check check
                    if visualize
                        hold on
                        plot([left_foot_heelstrike_0 left_foot_heelstrike_1 left_foot_heelstrike_2], [0 0 0]-0.01, 'v', 'linewidth', 3);
                        plot([ left_foot_pushoff_0 left_foot_pushoff_1 left_foot_pushoff_2], [0 0 0]-0.01, '^', 'linewidth', 3);
                        plot([right_foot_heelstrike_0 right_foot_heelstrike_1 right_foot_heelstrike_2], [0 0 0]+0.01, 'v', 'linewidth', 3);
                        plot([ right_foot_pushoff_0 right_foot_pushoff_1 right_foot_pushoff_2], [0 0 0]+0.01, '^', 'linewidth', 3);
                        % note: this can crash if one of thse events is empty, because we are plotting before we
                        % have checked that
                    end                  

                    stretch_times = [left_foot_heelstrike_0 right_foot_heelstrike_0 left_foot_heelstrike_1 right_foot_heelstrike_1 left_foot_heelstrike_2];


                    % determine stimulus
                    check_time_delay = 0.1;
                    check_time = stretch_times(1) + check_time_delay;
                    [~, check_index] = min(abs(time_analog - check_time));
                    if illusion_trajectory(check_index) > 0
                        stimulus_list = {'STIM_RIGHT'};
                    elseif illusion_trajectory(check_index) < 0
                        stimulus_list = {'STIM_LEFT'};
                    elseif illusion_trajectory(check_index) == 0
                        stimulus_list = {'STIM_NONE'};
                    end

                else
                    stimulus_list = cell(0);
                    trigger_foot_list = cell(0);
                    stretch_times = cell(0);
                    stance_foot_data = cell(0);
                end
                % add new variables to be saved
                conditions_trial = struct;
                conditions_trial.stimulus_list = stimulus_list;
                conditions_trial.trigger_foot_list = trigger_foot_list;
                event_variables_to_save.stretch_times = stretch_times;
                event_variables_to_save.stance_foot_data = stance_foot_data;
            end
            
            if strcmp(condition_stimulus, 'ARMSENSE')
                            
                % sort out type-day-combination
                if strcmp(condition_experimental, 'determine_from_type_day_combination')

                    if strcmp(condition_day, 'day1') && strcmp(this_trial_type, 'preOG')
                        condition_experimental = 'pre';
                    end
                    if strcmp(condition_day, 'day1') && strcmp(this_trial_type, 'postOG')
                        condition_experimental = 'post0';
                    end
                    if strcmp(condition_day, 'day2') && strcmp(this_trial_type, 'preOG')
                        condition_experimental = 'post4';
                    end
                end
                    
                % determine stride identification type
                stride_identification = study_settings.get('stride_identification');
                if strcmp(stride_identification, 'step_4-step_5')
                   if length(left_touchdown_times) < 3 | length(right_touchdown_times) < 3
                       this_stretch_start = 0;
                       this_stretch_end = 0;
                   else                   
                       if left_touchdown_times(1) <= right_touchdown_times(1)
                            this_stretch_start = right_touchdown_times(2);
                            this_stretch_end = right_touchdown_times(3);
                            band_delimiter = min(left_touchdown_times(left_touchdown_times>this_stretch_start));
                            first_stance_foot = 'STANCE_RIGHT';
                            second_stance_foot = 'STANCE_LEFT';
                       else
                            this_stretch_start = left_touchdown_times(2);
                            this_stretch_end = left_touchdown_times(3);
                            band_delimiter = min(right_touchdown_times(right_touchdown_times>this_stretch_start));
                            first_stance_foot = 'STANCE_LEFT';
                            second_stance_foot = 'STANCE_RIGHT';
                       end
                   end
                    % go through events and take stretches
                    stretch_times = [];
                    stance_foot_data = {};
                    condition_experimental_list = {};
                    condition_startfoot_list = {};
                    number_of_stretches = size(stretch_times, 1);
                    
                    % check if we have a valid stretch
                    if ~isempty(band_delimiter) && band_delimiter < this_stretch_end
                        stretch_times = [this_stretch_start, band_delimiter, this_stretch_end];
                        stance_foot_data = {first_stance_foot, second_stance_foot};
                        condition_experimental_list = condition_experimental;
                        condition_startfoot_list = {first_stance_foot; second_stance_foot};
                        bands_per_stretch = 2;
                        
                        if most_affected == 'L'
                            condition_affectedSide_list = 'L';
                            % assign affected_startfoot
                            % change to first_stance_foot..
                            if strcmp(first_stance_foot, 'STANCE_LEFT')
                                condition_affected_stancefoot_list = 'STANCE_AFFECTED';
                            elseif strcmp(first_stance_foot, 'STANCE_RIGHT')
                                condition_affected_stancefoot_list = 'STANCE_UNAFFECTED';
                            end
                        elseif most_affected == 'R'
                            condition_affectedSide_list = 'R';
                            % assign affected_startfoot
                            if strcmp(first_stance_foot, 'STANCE_RIGHT')
                                condition_affected_stancefoot_list = 'STANCE_AFFECTED';
                            elseif strcmp(first_stance_foot, 'STANCE_LEFT')
                                condition_affected_stancefoot_list = 'STANCE_UNAFFECTED';
                            end
                        end
                    end
                    
                    % fill in stuff
                    number_of_stretches = size(stretch_times, 1);
                    stretch_pushoff_times = zeros(size(stretch_times));
                    if isempty(stretch_times)
                        stretch_start_times = [];
                        stretch_end_times = [];
                    else
                        stretch_start_times = stretch_times(:, 1);
                        stretch_end_times = stretch_times(:, end);
                    end
                end
                if strcmp(stride_identification, 'all_but_first') 
                    % determine first stance foot
                    if left_touchdown_times(1) <= right_touchdown_times(1)
                        stretch_starter_events = left_touchdown_times;
                        band_delimiter_events = right_touchdown_times;
                        first_stance_foot = 'STANCE_RIGHT';
                        second_stance_foot = 'STANCE_LEFT';
                    else
                        stretch_starter_events = right_touchdown_times;
                        band_delimiter_events = left_touchdown_times;
                        first_stance_foot = 'STANCE_LEFT';
                        second_stance_foot = 'STANCE_RIGHT';
                    end                    

                    % go through events and take stretches
                    stretch_times = [];
                    stance_foot_data = {};
                    condition_experimental_list = {};
                    condition_startfoot_list = {};
                    condition_affectedSide_list = cell(number_of_stretches, 1);
                    condition_affected_stancefoot_list = cell(number_of_stretches, 1);
                  
                    number_of_stretches = length(stretch_starter_events) - 1;
                    
                    for i_stretch = 1 : number_of_stretches % assign i_stretch to # of stretches
                        this_stretch_start = stretch_starter_events(i_event);
                        this_stretch_end = stretch_starter_events(i_event+1);
                        band_delimiter = min(band_delimiter_events(band_delimiter_events>this_stretch_start));
                        % check if we have a valid stretch
                        if ~isempty(band_delimiter) && band_delimiter < this_stretch_end
                            this_stretch = [this_stretch_start band_delimiter this_stretch_end];
                            stretch_times = [stretch_times; this_stretch];
                            stance_foot_data = [stance_foot_data; {first_stance_foot, second_stance_foot}];
                            
                            if most_affected == 'L'
                                condition_affectedSide_list{i_stretch} = 'L';
                                % assign affected_startfoot
                                % change to first_stance_foot..
                                if strcmp(first_stance_foot, 'STANCE_LEFT')
                                    condition_affected_stancefoot_list{i_stretch} = 'STANCE_AFFECTED';
                                elseif strcmp(first_stance_foot, 'STANCE_RIGHT')
                                    condition_affected_stancefoot_list{i_stretch} = 'STANCE_UNAFFECTED';
                                end
                            elseif most_affected == 'R'
                                condition_affectedSide_list{i_stretch} = 'R';
                                % assign affected_startfoot
                                if strcmp(first_stance_foot, 'STANCE_RIGHT')
                                    condition_affected_stancefoot_list{i_stretch} = 'STANCE_AFFECTED';
                                elseif strcmp(first_stance_foot, 'STANCE_LEFT')
                                    condition_affected_stancefoot_list{i_stretch} = 'STANCE_UNAFFECTED';
                                end
                            end
                            
                            condition_startfoot_list = [condition_startfoot_list; first_stance_foot];
                            condition_experimental_list = [condition_experimental_list; condition_experimental];
                        end

                    end

                    if strcmp(this_trial_type(end-1:end), 'OG')
                        % remove all but first two stretches
                        stretch_times(3:end, :) = [];
                    end

                    % fill in stuff
                    number_of_stretches = size(stretch_times, 1);
                    stretch_pushoff_times = zeros(size(stretch_times));
                    if isempty(stretch_times)
                        stretch_start_times = [];
                        stretch_end_times = [];
                    else
                        stretch_start_times = stretch_times(:, 1);
                        stretch_end_times = stretch_times(:, end);
                    end

                    bands_per_stretch = 2;
                end
                % restructure for saving
                conditions_trial = struct;
                conditions_trial.condition_experimental_list = condition_experimental_list;
                conditions_trial.condition_startfoot_list = condition_startfoot_list;
                conditions_trial.condition_affected_stancefoot_list = condition_affected_stancefoot_list;
                conditions_trial.condition_affectedSide_list = condition_affectedSide_list;
                event_variables_to_save.stretch_pushoff_times = stretch_pushoff_times;
                event_variables_to_save.stretch_times = stretch_times;
                
                event_variables_to_save.stance_foot_data = stance_foot_data;
                       
            end
            
            if strcmp(condition_stimulus, 'GVS')
                bands_per_stretch = 1;
                
                % for each trigger, extract conditions and relevant step events
                number_of_triggers = length(trigger_indices_mocap);
                removal_flags = zeros(number_of_triggers, 1);
                stim_start_indices_labview = [stim_start_indices_labview zeros(number_of_triggers, 5)]; %#ok<AGROW>
                stretch_start_times = zeros(number_of_triggers, 6);
                stretch_end_times = zeros(number_of_triggers, 6);
                stretch_pushoff_times = zeros(number_of_triggers, 6);
                closest_heelstrike_distance_times = zeros(number_of_triggers, 1);
                condition_stance_foot_list = cell(number_of_triggers, 6);
                perturbation_list = cell(number_of_triggers, 6);
                condition_delay_list = cell(number_of_triggers, 6);
                condition_index_list = cell(number_of_triggers, 6);
                condition_experimental_list = cell(number_of_triggers, 6);
                condition_stimulus_list = cell(number_of_triggers, 6);
                condition_day_list = cell(number_of_triggers, 6);
                
                for i_trigger = 1 : number_of_triggers
                   % perturbation condition
                    if illusion_trajectory(stim_start_indices_labview(i_trigger)) > 0
                        perturbation_list{i_trigger, 1} = 'ILLUSION_RIGHT';
                        perturbation_list{i_trigger, 2} = 'ILLUSION_RIGHT';
                        perturbation_list{i_trigger, 3} = 'ILLUSION_RIGHT';
                        perturbation_list{i_trigger, 4} = 'ILLUSION_RIGHT';
                    elseif illusion_trajectory(stim_start_indices_labview(i_trigger)) < 0
                        perturbation_list{i_trigger, 1} = 'ILLUSION_LEFT';
                        perturbation_list{i_trigger, 2} = 'ILLUSION_LEFT';
                        perturbation_list{i_trigger, 3} = 'ILLUSION_LEFT';
                        perturbation_list{i_trigger, 4} = 'ILLUSION_LEFT';
                    else
        %                 disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_stimulus(trigger_indices_labview(i_trigger))) ' - no stim']);
                    end
                    perturbation_list{i_trigger, 5} = 'CONTROL';
                    perturbation_list{i_trigger, 6} = 'CONTROL';
                
                    % delay condition
                    wait_time_stim = time_stimulus(stim_start_indices_labview(i_trigger)) - time_stimulus(trigger_indices_labview(i_trigger));
                    delay_time_labels = study_settings.get('delay_time_labels');
                    [~, wait_condition_index] = min(abs(study_settings.get('delay_times') - wait_time_stim));
                    if iscell(study_settings.get('delay_time_labels'))
                        delay_condition_label = delay_time_labels{wait_condition_index};
                        condition_delay_list{i_trigger, 5} = 'CONTROL';
                        condition_delay_list{i_trigger, 6} = 'CONTROL';
                    else
                        delay_condition_label = study_settings.get('delay_time_labels');
                        condition_delay_list{i_trigger, 5} = delay_condition_label;
                        condition_delay_list{i_trigger, 6} = delay_condition_label;
                    end
                    condition_delay_list{i_trigger, 1} = delay_condition_label;
                    condition_delay_list{i_trigger, 2} = delay_condition_label;
                    condition_delay_list{i_trigger, 3} = delay_condition_label;
                    condition_delay_list{i_trigger, 4} = delay_condition_label;
                    
                    % experimental condition
                    condition_experimental_list{i_trigger, 1} = condition_experimental;
                    condition_experimental_list{i_trigger, 2} = condition_experimental;
                    condition_experimental_list{i_trigger, 3} = condition_experimental;
                    condition_experimental_list{i_trigger, 4} = condition_experimental;
                    condition_experimental_list{i_trigger, 5} = condition_experimental;
                    condition_experimental_list{i_trigger, 6} = condition_experimental;
                    
                    % stimulus condition
                    condition_stimulus_list{i_trigger, 1} = condition_stimulus;
                    condition_stimulus_list{i_trigger, 2} = condition_stimulus;
                    condition_stimulus_list{i_trigger, 3} = condition_stimulus;
                    condition_stimulus_list{i_trigger, 4} = condition_stimulus;
                    condition_stimulus_list{i_trigger, 5} = condition_stimulus;
                    condition_stimulus_list{i_trigger, 6} = condition_stimulus;
                    
                    % day condition
                    condition_day_list{i_trigger, 1} = condition_day;
                    condition_day_list{i_trigger, 2} = condition_day;
                    condition_day_list{i_trigger, 3} = condition_day;
                    condition_day_list{i_trigger, 4} = condition_day;
                    condition_day_list{i_trigger, 5} = condition_day;
                    condition_day_list{i_trigger, 6} = condition_day;
                    
                    % get closest heelstrike
                    [~, index_left] = min(abs(left_touchdown_times - trigger_times(i_trigger)));
                    [~, index_right] = min(abs(right_touchdown_times - trigger_times(i_trigger)));
                    
                    % is the closest left heelstrike within the acceptable interval?
                    closest_left_heelstrike = left_touchdown_times(index_left);
                    time_difference_left = closest_left_heelstrike - trigger_times(i_trigger); % where does the closest left heelstrike lie relative to the trigger?
                    if -time_to_nearest_heelstrike_before_trigger_threshold < time_difference_left && time_difference_left < time_to_nearest_heelstrike_after_trigger_threshold
                    	% left heelstrike is acceptable
                        left_heelstrike_acceptable = true;
                    else
                        left_heelstrike_acceptable = false;
                    end
                    
                    % is the closest right heelstrike within the acceptable interval?
                    closest_right_heelstrike = right_touchdown_times(index_right);
                    time_difference_right = closest_right_heelstrike - trigger_times(i_trigger); % where does the closest right heelstrike lie relative to the trigger?
                    if -time_to_nearest_heelstrike_before_trigger_threshold < time_difference_right && time_difference_right < time_to_nearest_heelstrike_after_trigger_threshold
                    	% right heelstrike is acceptable
                        right_heelstrike_acceptable = true;
                    else
                        right_heelstrike_acceptable = false;
                    end
                    
                    % accept the acceptable one
                    if left_heelstrike_acceptable && ~right_heelstrike_acceptable
                        % triggered by left heelstrike
                        trigger_foot = 'left';
                        closest_heelstrike_distance_times(i_trigger) = time_difference_left;
                    elseif ~left_heelstrike_acceptable && right_heelstrike_acceptable
                        % triggered by right heelstrike
                        trigger_foot = 'right';
                        closest_heelstrike_distance_times(i_trigger) = time_difference_right;
                    elseif left_heelstrike_acceptable && right_heelstrike_acceptable
                        trigger_foot = 'unclear';
                        removal_flags(i_trigger) = 1;
                    elseif ~left_heelstrike_acceptable && ~right_heelstrike_acceptable
                        trigger_foot = 'unclear';
                        removal_flags(i_trigger) = 1;
                    end
                    
%                     closest_heelstrike_distance_time = min([distance_to_trigger_left_time distance_to_trigger_right_time]);

                    % confirm that the trigger happened less than 100ms (or whatever time's set) before the actual heelstrike
%                     if closest_heelstrike_distance_time > time_to_nearest_heelstrike_before_trigger_threshold
%                         removal_flags(i_trigger) = 1;
%         %                 disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_stimulus(trigger_indices_labview(i_trigger))) ' - trigger not in sync with heelstrike']);
%                     end

                    if strcmp(trigger_foot, 'left')
                        % check whether time offset is positive or negative
%                         if left_touchdown_times(index_left) - trigger_times(i_trigger) >= 0
%                             closest_heelstrike_distance_times(i_trigger) = closest_heelstrike_distance_time;
%                         else
%                             closest_heelstrike_distance_times(i_trigger) = -closest_heelstrike_distance_time;
%                         end
                        if index_left == 1 || length(left_touchdown_times) < index_left + 1 || removal_flags(i_trigger) == 1
                            % data doesn't include previous or next step
                            removal_flags(i_trigger) = 1;
                            right_foot_heelstrike_minus_1   = NaN;
                            right_foot_heelstrike_0         = NaN;
                            right_foot_heelstrike_plus_1    = NaN;
                            right_foot_heelstrike_plus_2    = NaN;
                            left_foot_heelstrike_minus_1    = NaN;
                            left_foot_heelstrike_0          = NaN;
                            left_foot_heelstrike_plus_1     = NaN;
                            left_foot_heelstrike_plus_2     = NaN;
                            
                            right_foot_pushoff_minus_1      = NaN;
                            right_foot_pushoff_0            = NaN;
                            right_foot_pushoff_plus_1       = NaN;
                            right_foot_pushoff_plus_2       = NaN;
                            left_foot_pushoff_minus_1       = NaN;
                            left_foot_pushoff_0             = NaN;
                            left_foot_pushoff_plus_1        = NaN;
                            left_foot_pushoff_plus_2        = NaN;
                        else
                            left_foot_heelstrike_minus_1    = left_touchdown_times(index_left-1);
                            left_foot_heelstrike_0          = left_touchdown_times(index_left);
                            left_foot_heelstrike_plus_1     = left_touchdown_times(index_left+1);
                            left_foot_heelstrike_plus_2     = left_touchdown_times(index_left+2);
                            left_foot_pushoff_minus_1       = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_minus_1));
                            left_foot_pushoff_0             = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_0));
                            left_foot_pushoff_plus_1        = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_plus_1));
                            left_foot_pushoff_plus_2        = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_plus_2));

                            right_foot_heelstrike_minus_1   = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_minus_1));
                            right_foot_heelstrike_0         = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_0));
                            right_foot_heelstrike_plus_1    = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_plus_1));
                            right_foot_heelstrike_plus_2    = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_plus_2));
                            right_foot_pushoff_minus_1      = max(right_pushoff_times(right_pushoff_times <= left_foot_pushoff_minus_1));
                            right_foot_pushoff_0            = max(right_pushoff_times(right_pushoff_times <= left_foot_pushoff_0));
                            right_foot_pushoff_plus_1       = max(right_pushoff_times(right_pushoff_times <= left_foot_pushoff_plus_1));
                            right_foot_pushoff_plus_2       = max(right_pushoff_times(right_pushoff_times <= left_foot_pushoff_plus_2));

                            % notify if events are not sorted properly
                            if ~issorted ...
                                  ( ...
                                    [ ...
                                      left_foot_heelstrike_minus_1 ...
                                      right_foot_pushoff_minus_1 ...
                                      right_foot_heelstrike_minus_1 ...
                                      left_foot_pushoff_minus_1 ...
                                      left_foot_heelstrike_0 ...
                                      right_foot_pushoff_0 ...
                                      right_foot_heelstrike_0 ...
                                      left_foot_pushoff_0 ...
                                      left_foot_heelstrike_plus_1 ...
                                      right_foot_pushoff_plus_1 ...
                                      right_foot_heelstrike_plus_1 ...
                                      left_foot_pushoff_plus_1 ...
                                      left_foot_heelstrike_plus_2 ...
                                      right_foot_pushoff_plus_2 ...
                                      right_foot_heelstrike_plus_2 ...
                                      left_foot_pushoff_plus_2 ...
                                    ] ...
                                  );
                                disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(time_stimulus(trigger_indices_labview(i_trigger)))]);
                            end
                            % check check
                            if visualize
                                plot([left_foot_heelstrike_minus_1 left_foot_heelstrike_0 left_foot_heelstrike_plus_1 left_foot_heelstrike_plus_2], [0 0 0 0]-0.01, 'v', 'linewidth', 3);
                                plot([left_foot_pushoff_minus_1  left_foot_pushoff_0 left_foot_pushoff_plus_1 left_foot_pushoff_plus_2], [0 0 0 0]-0.01, '^', 'linewidth', 3);
                                plot([right_foot_heelstrike_minus_1 right_foot_heelstrike_0 right_foot_heelstrike_plus_1 right_foot_heelstrike_plus_2], [0 0 0 0]+0.01, 'v', 'linewidth', 3);
                                plot([right_foot_pushoff_minus_1  right_foot_pushoff_0 right_foot_pushoff_plus_1 right_foot_pushoff_plus_2], [0 0 0 0]+0.01, '^', 'linewidth', 3);
                                % note: this can crash if one of thse events is empty, because we are plotting before we
                                % have checked that
                            end                
                        end
                    elseif strcmp(trigger_foot, 'right')
                        % check whether time offset is positive or negative
%                         if right_touchdown_times(index_right) - trigger_times(i_trigger) > 0
%                             closest_heelstrike_distance_times(i_trigger) = closest_heelstrike_distance_time;
%                         else
%                             closest_heelstrike_distance_times(i_trigger) = -closest_heelstrike_distance_time;
%                         end
                        if index_right == 1 || length(right_touchdown_times) < index_right + 1 || removal_flags(i_trigger) == 1
                            % data doesn't include previous or next step
                            removal_flags(i_trigger) = 1;
                            left_foot_heelstrike_minus_1    = NaN;
                            left_foot_heelstrike_0          = NaN;
                            left_foot_heelstrike_plus_1     = NaN;
                            left_foot_heelstrike_plus_2     = NaN;
                            right_foot_heelstrike_minus_1   = NaN;
                            right_foot_heelstrike_0         = NaN;
                            right_foot_heelstrike_plus_1    = NaN;
                            right_foot_heelstrike_plus_2    = NaN;
                            
                            right_foot_pushoff_minus_1      = NaN;
                            right_foot_pushoff_0            = NaN;
                            right_foot_pushoff_plus_1       = NaN;
                            right_foot_pushoff_plus_2       = NaN;
                            left_foot_pushoff_minus_1       = NaN;
                            left_foot_pushoff_0             = NaN;
                            left_foot_pushoff_plus_1        = NaN;
                            left_foot_pushoff_plus_2        = NaN;
                        else
                            right_foot_heelstrike_minus_1       = right_touchdown_times(index_right-1);
                            right_foot_heelstrike_0             = right_touchdown_times(index_right);
                            right_foot_heelstrike_plus_1        = right_touchdown_times(index_right+1);
                            right_foot_heelstrike_plus_2        = right_touchdown_times(index_right+2);
                            
                            right_foot_pushoff_minus_1          = min(right_pushoff_times(right_pushoff_times >= right_foot_heelstrike_minus_1));
                            right_foot_pushoff_0                = min(right_pushoff_times(right_pushoff_times >= right_foot_heelstrike_0));
                            right_foot_pushoff_plus_1           = min(right_pushoff_times(right_pushoff_times >= right_foot_heelstrike_plus_1));
                            right_foot_pushoff_plus_2           = min(right_pushoff_times(right_pushoff_times >= right_foot_heelstrike_plus_2));

                            left_foot_heelstrike_minus_1        = min(left_touchdown_times(left_touchdown_times >= right_foot_heelstrike_minus_1));
                            left_foot_heelstrike_0              = min(left_touchdown_times(left_touchdown_times >= right_foot_heelstrike_0));
                            left_foot_heelstrike_plus_1         = min(left_touchdown_times(left_touchdown_times >= right_foot_heelstrike_plus_1));
                            left_foot_heelstrike_plus_2         = min(left_touchdown_times(left_touchdown_times >= right_foot_heelstrike_plus_2));
                            left_foot_pushoff_minus_1           = max(left_pushoff_times(left_pushoff_times <= right_foot_pushoff_minus_1));
                            left_foot_pushoff_0                 = max(left_pushoff_times(left_pushoff_times <= right_foot_pushoff_0));
                            left_foot_pushoff_plus_1            = max(left_pushoff_times(left_pushoff_times <= right_foot_pushoff_plus_1));
                            left_foot_pushoff_plus_2            = max(left_pushoff_times(left_pushoff_times <= right_foot_pushoff_plus_2));

                            % notify if events are not sorted properly
                            if ~issorted ...
                                  ( ...
                                    [ ...
                                      right_foot_heelstrike_minus_1 ...
                                      left_foot_pushoff_minus_1 ...
                                      left_foot_heelstrike_minus_1 ...
                                      right_foot_pushoff_minus_1 ...
                                      right_foot_heelstrike_0 ...
                                      left_foot_pushoff_0 ...
                                      left_foot_heelstrike_0 ...
                                      right_foot_pushoff_0 ...
                                      right_foot_heelstrike_plus_1 ...
                                      left_foot_pushoff_plus_1 ...
                                      left_foot_heelstrike_plus_1 ...
                                      right_foot_pushoff_plus_1 ...
                                      right_foot_heelstrike_plus_2 ...
                                      left_foot_pushoff_plus_2 ...
                                      left_foot_heelstrike_plus_2 ...
                                      right_foot_pushoff_plus_2 ...
                                    ] ...
                                  );
                              disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(time_stimulus(trigger_indices_labview(i_trigger)))]);

                            end

                            if visualize
                                plot([left_foot_heelstrike_minus_1 left_foot_heelstrike_0 left_foot_heelstrike_plus_1 left_foot_heelstrike_plus_2], [0 0 0 0]-0.01, 'v', 'linewidth', 2, 'color', 'r');
                                plot([left_foot_pushoff_minus_1  left_foot_pushoff_0 left_foot_pushoff_plus_1 left_foot_pushoff_plus_2], [0 0 0 0]-0.01, '^', 'linewidth', 2, 'color', 'g');
                                plot([right_foot_heelstrike_minus_1 right_foot_heelstrike_0 right_foot_heelstrike_plus_1 right_foot_heelstrike_plus_2], [0 0 0 0]+0.01, 'v', 'linewidth', 2, 'color', 'r');
                                plot([right_foot_pushoff_minus_1  right_foot_pushoff_0 right_foot_pushoff_plus_1 right_foot_pushoff_plus_2], [0 0 0 0]+0.01, '^', 'linewidth', 2, 'color', 'g');
                            end            
                        end            
                    else
                        trigger_foot = 'unclear';
                        disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_stimulus(trigger_indices_labview(i_trigger))) ' - triggering heelstrike unclear']);
                        left_foot_heelstrike_minus_1    = 0;
                        left_foot_heelstrike_0          = 0;
                        left_foot_heelstrike_plus_1     = 0;
                        left_foot_heelstrike_plus_2     = 0;
                        left_foot_pushoff_minus_1       = 0;
                        left_foot_pushoff_0             = 0;
                        left_foot_pushoff_plus_1        = 0;
                        left_foot_pushoff_plus_2        = 0;

                        right_foot_heelstrike_minus_1   = 0;
                        right_foot_heelstrike_0         = 0;
                        right_foot_heelstrike_plus_1    = 0;
                        right_foot_heelstrike_plus_2    = 0;
                        right_foot_pushoff_minus_1      = 0;
                        right_foot_pushoff_0            = 0;
                        right_foot_pushoff_plus_1       = 0;
                        right_foot_pushoff_plus_2       = 0;

                        removal_flags(i_trigger) = 1;
                    end


                    % flag for removal if not all events are present
                    if any( ...
                            [ ...
                              isempty(left_foot_heelstrike_minus_1) ...
                              isempty(left_foot_heelstrike_0) ...
                              isempty(left_foot_heelstrike_plus_1) ...
                              isempty(left_foot_heelstrike_plus_2) ...
                              isempty(right_foot_heelstrike_minus_1) ...
                              isempty(right_foot_heelstrike_0) ...
                              isempty(right_foot_heelstrike_plus_1) ...
                              isempty(right_foot_heelstrike_plus_2) ...
                              isempty(left_foot_pushoff_minus_1) ...
                              isempty(left_foot_pushoff_0) ...
                              isempty(left_foot_pushoff_plus_1) ...
                              isempty(left_foot_pushoff_plus_2) ...
                              isempty(right_foot_pushoff_minus_1) ...
                              isempty(right_foot_pushoff_0) ...
                              isempty(right_foot_pushoff_plus_1) ...
                              isempty(right_foot_pushoff_plus_2) ...
                            ] ...
                          ) || removal_flags(i_trigger) == 1
                        % not all events are present
                        removal_flags(i_trigger) = 1;
                        left_foot_heelstrike_minus_1 = NaN;
                        left_foot_heelstrike_0 = NaN;
                        left_foot_heelstrike_plus_1 = NaN;
                        left_foot_heelstrike_plus_2 = NaN;
                        right_foot_heelstrike_minus_1 = NaN;
                        right_foot_heelstrike_0 = NaN;
                        right_foot_heelstrike_plus_1 = NaN;
                        right_foot_heelstrike_plus_2 = NaN;

                    else % all events are present
                        % mark relevant event delimiters depending on wait time and triggering foot
                        if strcmp(trigger_foot, 'right')
                            % triggered by right foot heelstrike
                            % triggering step
                            condition_index_list{i_trigger, 1} = 'ONE';
                            condition_stance_foot_list{i_trigger, 1} = 'STANCE_RIGHT';
                            % post stimulus steps
                            condition_index_list{i_trigger, 2} = 'TWO';
                            condition_stance_foot_list{i_trigger, 2} = 'STANCE_LEFT';
                            condition_index_list{i_trigger, 3} = 'THREE';
                            condition_stance_foot_list{i_trigger, 3} = 'STANCE_RIGHT';
                            condition_index_list{i_trigger, 4} = 'FOUR';
                            condition_stance_foot_list{i_trigger, 4} = 'STANCE_LEFT';
                            % use two previous steps as control
                            condition_index_list{i_trigger, 5} = 'CONTROL';
                            condition_stance_foot_list{i_trigger, 5} = 'STANCE_RIGHT';
                            condition_index_list{i_trigger, 6} = 'CONTROL';
                            condition_stance_foot_list{i_trigger, 6} = 'STANCE_LEFT';

                            % this step
                            stretch_start_times(i_trigger, 1) = right_foot_heelstrike_0;
                            stretch_pushoff_times(i_trigger, 1) = left_foot_pushoff_0;
                            stretch_end_times(i_trigger, 1) = left_foot_heelstrike_0;
                            % post-stimulus steps
                            stretch_start_times(i_trigger, 2) = left_foot_heelstrike_0;
                            stretch_pushoff_times(i_trigger, 2) = right_foot_pushoff_0;
                            stretch_end_times(i_trigger, 2) = right_foot_heelstrike_plus_1;

                            stretch_start_times(i_trigger, 3) = right_foot_heelstrike_plus_1;
                            stretch_pushoff_times(i_trigger, 3) = left_foot_pushoff_plus_1;
                            stretch_end_times(i_trigger, 3) = left_foot_heelstrike_plus_1;

                            stretch_start_times(i_trigger, 4) = left_foot_heelstrike_plus_1;
                            stretch_pushoff_times(i_trigger, 4) = right_foot_pushoff_plus_1;
                            stretch_end_times(i_trigger, 4) = right_foot_heelstrike_plus_2;
                            % use two previous steps as control
                            stretch_start_times(i_trigger, 5) = right_foot_heelstrike_minus_1;
                            stretch_pushoff_times(i_trigger, 5) = left_foot_pushoff_minus_1;
                            stretch_end_times(i_trigger, 5) = left_foot_heelstrike_minus_1;

                            stretch_start_times(i_trigger, 6) = left_foot_heelstrike_minus_1;
                            stretch_pushoff_times(i_trigger, 6) = right_foot_pushoff_minus_1;
                            stretch_end_times(i_trigger, 6) = right_foot_heelstrike_0;

                            % visualize
                            if visualize
                                plot([stretch_start_times(i_trigger, 1) stretch_end_times(i_trigger, 1)], [-1 1]*-0.01, 'color', 'r', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 2) stretch_end_times(i_trigger, 2)], [-1 1]*+0.01, 'color', 'g', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 3) stretch_end_times(i_trigger, 3)], [-1 1]*-0.01, 'color', 'b', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 4) stretch_end_times(i_trigger, 4)], [-1 1]*+0.01, 'color', 'm', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 5) stretch_end_times(i_trigger, 5)], [-1 1]*-0.01, 'color', 'y', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 6) stretch_end_times(i_trigger, 6)], [-1 1]*+0.01, 'color', 'k', 'linewidth', 3);
                            end

                        elseif strcmp(trigger_foot, 'left')
                            % triggered by left foot heelstrike
                            % this step
                            condition_index_list{i_trigger, 1} = 'ONE';
                            condition_stance_foot_list{i_trigger, 1} = 'STANCE_LEFT';
                            % post stimulus steps
                            condition_index_list{i_trigger, 2} = 'TWO';
                            condition_stance_foot_list{i_trigger, 2} = 'STANCE_RIGHT';
                            condition_index_list{i_trigger, 3} = 'THREE';
                            condition_stance_foot_list{i_trigger, 3} = 'STANCE_LEFT';
                            condition_index_list{i_trigger, 4} = 'FOUR';
                            condition_stance_foot_list{i_trigger, 4} = 'STANCE_RIGHT';
                            % use two previous steps as control
                            condition_index_list{i_trigger, 5} = 'CONTROL';
                            condition_stance_foot_list{i_trigger, 5} = 'STANCE_LEFT';
                            condition_index_list{i_trigger, 6} = 'CONTROL';
                            condition_stance_foot_list{i_trigger, 6} = 'STANCE_RIGHT';

                            % this step
                            stretch_start_times(i_trigger, 1) = left_foot_heelstrike_0;
                            stretch_pushoff_times(i_trigger, 1) = right_foot_pushoff_0;
                            stretch_end_times(i_trigger, 1) = right_foot_heelstrike_0;
                            % post stimulus steps
                            stretch_start_times(i_trigger, 2) = right_foot_heelstrike_0;
                            stretch_pushoff_times(i_trigger, 2) = left_foot_pushoff_0;
                            stretch_end_times(i_trigger, 2) = left_foot_heelstrike_plus_1;

                            stretch_start_times(i_trigger, 3) = left_foot_heelstrike_plus_1;
                            stretch_pushoff_times(i_trigger, 3) = right_foot_pushoff_plus_1;
                            stretch_end_times(i_trigger, 3) = right_foot_heelstrike_plus_1;

                            stretch_start_times(i_trigger, 4) = right_foot_heelstrike_plus_1;
                            stretch_pushoff_times(i_trigger, 4) = left_foot_pushoff_plus_1;
                            stretch_end_times(i_trigger, 4) = left_foot_heelstrike_plus_2;
                            % use two previous steps as control
                            stretch_start_times(i_trigger, 5) = left_foot_heelstrike_minus_1;
                            stretch_pushoff_times(i_trigger, 5) = right_foot_pushoff_minus_1;
                            stretch_end_times(i_trigger, 5) = right_foot_heelstrike_minus_1;

                            stretch_start_times(i_trigger, 6) = right_foot_heelstrike_minus_1;
                            stretch_pushoff_times(i_trigger, 6) = left_foot_pushoff_minus_1;
                            stretch_end_times(i_trigger, 6) = left_foot_heelstrike_0;

                            % visualize
                            if visualize
                                plot([stretch_start_times(i_trigger, 1) stretch_end_times(i_trigger, 1)], [-1 1]*0.01, 'color', 'r', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 2) stretch_end_times(i_trigger, 2)], [-1 1]*-0.01, 'color', 'g', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 3) stretch_end_times(i_trigger, 3)], [-1 1]*0.01, 'color', 'b', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 4) stretch_end_times(i_trigger, 4)], [-1 1]*-0.01, 'color', 'm', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 5) stretch_end_times(i_trigger, 5)], [-1 1]*0.01, 'color', 'y', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 6) stretch_end_times(i_trigger, 6)], [-1 1]*-0.01, 'color', 'k', 'linewidth', 3);
                            end
                        else
                            removal_flags(i_trigger) = 1;
                            condition_stance_foot_list{i_trigger, 1} = 'UNCLEAR';
                            condition_stance_foot_list{i_trigger, 2} = 'UNCLEAR';
                            condition_stance_foot_list{i_trigger, 3} = 'UNCLEAR';
                            condition_stance_foot_list{i_trigger, 4} = 'UNCLEAR';
                            condition_stance_foot_list{i_trigger, 5} = 'UNCLEAR';
                            condition_stance_foot_list{i_trigger, 6} = 'UNCLEAR';
                        end
                    end
                end

                % remove flagged triggers
                % needs to be restructured to only populate what is
                % necessary
                unflagged_indices = ~removal_flags;
                trigger_times = trigger_times(unflagged_indices);
                stim_start_indices_labview = stim_start_indices_labview(unflagged_indices, :);
                stretch_start_times = stretch_start_times(unflagged_indices, :);
                stretch_pushoff_times = stretch_pushoff_times(unflagged_indices, :);
                stretch_end_times = stretch_end_times(unflagged_indices, :);
                condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
                perturbation_list = perturbation_list(unflagged_indices, :);
                condition_delay_list = condition_delay_list(unflagged_indices, :);
                condition_index_list = condition_index_list(unflagged_indices, :);
                condition_experimental_list = condition_experimental_list(unflagged_indices, :); % XXX needs to be updated
%                 condition_startfoot_list = condition_startfoot_list(unflagged_indices, :);
                condition_stimulus_list = condition_stimulus_list(unflagged_indices, :); % XXX needs to be updated
                condition_day_list = condition_day_list(unflagged_indices, :); % XXX needs to be updated
%                 closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :);
                
                % reorder
                stretch_start_times = reshape(stretch_start_times, numel(stretch_start_times), 1);
                stretch_pushoff_times = reshape(stretch_pushoff_times, numel(stretch_pushoff_times), 1);
                stretch_end_times = reshape(stretch_end_times, numel(stretch_end_times), 1);
                
                condition_stance_foot_list = reshape(condition_stance_foot_list, numel(condition_stance_foot_list), 1);
                perturbation_list = reshape(perturbation_list, numel(perturbation_list), 1);
                condition_delay_list = reshape(condition_delay_list, numel(condition_delay_list), 1);
                condition_index_list = reshape(condition_index_list, numel(condition_index_list), 1);
                condition_experimental_list = reshape(condition_experimental_list, numel(condition_experimental_list), 1);
%                 condition_startfoot_list = reshape(condition_startfoot_list, numel(condition_startfoot_list), 1);
                condition_stimulus_list = reshape(condition_stimulus_list, numel(condition_stimulus_list), 1);
                condition_day_list = reshape(condition_day_list, numel(condition_day_list), 1);
                
                
                
                
%                 % we now have a neatly ordered list of stretches which we can prune
%                 % check step times and flag outliers
%                 number_of_stretches = length(stretch_start_times);
%                 stretch_durations = stretch_end_times - stretch_start_times;
%                 stretch_duration_outlier_limits = median(stretch_durations) * [0.8 1.2];
% 
%                 removal_flags = zeros(number_of_stretches, 1);
%                 removal_flags(stretch_durations < stretch_duration_outlier_limits(1)) = 1;
%                 removal_flags(stretch_durations > stretch_duration_outlier_limits(2)) = 1;
% 
%                 % remove flagged triggers
%                 unflagged_indices = ~removal_flags;
%                 stretch_start_times = stretch_start_times(unflagged_indices, :);
%                 stretch_pushoff_times = stretch_pushoff_times(unflagged_indices, :);
%                 stretch_end_times = stretch_end_times(unflagged_indices, :);
%                 condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
%                 condition_perturbation_list = condition_perturbation_list(unflagged_indices, :);
%                 condition_delay_list = condition_delay_list(unflagged_indices, :);
%                 condition_index_list = condition_index_list(unflagged_indices, :);
%                 condition_experimental_list = condition_experimental_list(unflagged_indices, :);
% %                 condition_startfoot_list = condition_startfoot_list(unflagged_indices, :);
%                 condition_stimulus_list = condition_stimulus_list(unflagged_indices, :);
%                 condition_day_list = condition_day_list(unflagged_indices, :);

                % restructure for saving
                stretch_times = [stretch_start_times stretch_end_times];
                conditions_trial = struct;
                conditions_trial.condition_stance_foot_list = condition_stance_foot_list;
                conditions_trial.condition_perturbation_list = perturbation_list;
                conditions_trial.condition_delay_list = condition_delay_list;
                conditions_trial.condition_index_list = condition_index_list;
                conditions_trial.condition_experimental_list = condition_experimental_list;
%                 conditions_trial.condition_startfoot_list = condition_startfoot_list;
                conditions_trial.condition_stimulus_list = condition_stimulus_list;
                conditions_trial.condition_day_list = condition_day_list;
                
%                 event_variables_to_save.stretch_start_times = stretch_start_times;
                event_variables_to_save.stretch_pushoff_times = stretch_pushoff_times;
%                 event_variables_to_save.stretch_end_times = stretch_end_times;
                event_variables_to_save.stretch_times = stretch_times;
                
                event_variables_to_save.stance_foot_data = condition_stance_foot_list; % TODO: haven't tested this yet. Adapted during the stretch rework, which is currently in development for Obstacle data

                
                
                
                % determine trigger foot
                condition_trigger_foot_list = cell(size(condition_stance_foot_list));
                for i_stretch = 1 : length(condition_trigger_foot_list)
                    if strcmp(condition_index_list{i_stretch}, 'ONE') || strcmp(condition_index_list{i_stretch}, 'THREE')
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')
                            condition_trigger_foot_list{i_stretch} = 'TRIGGER_RIGHT';
                        elseif strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')
                            condition_trigger_foot_list{i_stretch} = 'TRIGGER_LEFT';
                        end
                    end
                    if strcmp(condition_index_list{i_stretch}, 'TWO') || strcmp(condition_index_list{i_stretch}, 'FOUR')
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')
                            condition_trigger_foot_list{i_stretch} = 'TRIGGER_RIGHT';
                        elseif strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')
                            condition_trigger_foot_list{i_stretch} = 'TRIGGER_LEFT';
                        end
                    end
                    if strcmp(condition_index_list{i_stretch}, 'CONTROL')
                        condition_trigger_foot_list{i_stretch} = 'CONTROL';
                    end
                end
                conditions_trial.condition_trigger_foot_list = condition_trigger_foot_list;

                
                % determine direction
                condition_direction_list = cell(size(condition_stance_foot_list));
                for i_stretch = 1 : length(condition_direction_list)
                    if strcmp(condition_trigger_foot_list{i_stretch}, 'TRIGGER_RIGHT')
                        if strcmp(perturbation_list{i_stretch}, 'ILLUSION_RIGHT')
                            condition_direction_list{i_stretch} = 'TOWARDS';
                        end
                        if strcmp(perturbation_list{i_stretch}, 'ILLUSION_LEFT')
                            condition_direction_list{i_stretch} = 'AWAY';
                        end
                    end
                    if strcmp(condition_trigger_foot_list{i_stretch}, 'TRIGGER_LEFT')
                        if strcmp(perturbation_list{i_stretch}, 'ILLUSION_RIGHT')
                            condition_direction_list{i_stretch} = 'AWAY';
                        end
                        if strcmp(perturbation_list{i_stretch}, 'ILLUSION_LEFT')
                            condition_direction_list{i_stretch} = 'TOWARDS';
                        end
                    end
                    if strcmp(condition_trigger_foot_list{i_stretch}, 'CONTROL')
                        condition_direction_list{i_stretch} = 'CONTROL';
                    end
                end
                conditions_trial.condition_direction_list = condition_direction_list;
                conditions_trial.condition_stance_foot_list = condition_stance_foot_list;
                
                % put in placeholder for group
                conditions_trial.condition_group_list = condition_direction_list;
                [conditions_trial.condition_group_list{:}] = deal('to be determined');
                
                % make copy of stance foot information
                event_variables_to_save.stance_foot_data = condition_stance_foot_list;
            end
            
            if strcmp(experimental_paradigm, 'Vision') ...
                    || strcmp(experimental_paradigm, 'CadenceVision') ...
                    || strcmp(experimental_paradigm, 'GVS') ...
                    || strcmp(experimental_paradigm, 'CadenceGVS') ...
                    || strcmp(experimental_paradigm, 'FatigueGVS') ...
                    || strcmp(experimental_paradigm, 'OculusLaneRestriction') ...
                    || strcmp(experimental_paradigm, 'CognitiveLoadVision') ...
                    || strcmp(experimental_paradigm, 'CognitiveLoadGvs')
                        
                bands_per_stretch = study_settings.get('number_of_steps_to_analyze');
            
                
                conditions_trial = struct;
                number_of_triggers = length(trigger_indices_mocap);
                removal_flags = zeros(number_of_triggers, 1);
                stretch_times = zeros(number_of_triggers, bands_per_stretch+1);
                stretch_times_minus_1 = zeros(number_of_triggers, bands_per_stretch+2);
                closest_heelstrike_distance_times = zeros(number_of_triggers, 1);
                stance_foot_data = cell(number_of_triggers, bands_per_stretch);
                stance_foot_data_minus_1 = cell(number_of_triggers, bands_per_stretch+1);
                stimulus_list = cell(number_of_triggers, 1); % stimulus STIM_LEFT, STIM_RIGHT or STIM_NONE
                amplitude_list = cell(number_of_triggers, 1); % amplitude 30, 60, 120
                trigger_foot_list = cell(number_of_triggers, 1); % triggering foot TRIGGER_LEFT or TRIGGER_RIGHT
                
                for i_trigger = 1 : number_of_triggers
                    % determine stimulus
                    if trigger_indices_labview(i_trigger) == length(illusion_trajectory)
                        removal_flags(i_trigger) = 1;
                    else
                        if illusion_trajectory(trigger_indices_labview(i_trigger)+1) > 0
                            stimulus_list{i_trigger} = 'STIM_RIGHT';
                        end
                        if illusion_trajectory(trigger_indices_labview(i_trigger)+1) < 0
                            stimulus_list{i_trigger} = 'STIM_LEFT';
                        end
                        if illusion_trajectory(trigger_indices_labview(i_trigger)+1) == 0
                            stimulus_list{i_trigger} = 'STIM_NONE';
                        end
                    end
                    % determine trigger foot
                    
                    % get closest heelstrike on either side
                    [~, index_left] = min(abs(left_touchdown_times - trigger_times(i_trigger)));
                    [~, index_right] = min(abs(right_touchdown_times - trigger_times(i_trigger)));
                    
                    % is the closest left heelstrike within the acceptable interval?
                    closest_left_heelstrike = left_touchdown_times(index_left);
                    time_difference_left = closest_left_heelstrike - trigger_times(i_trigger); % where does the closest left heelstrike lie relative to the trigger?
                    if -time_to_nearest_heelstrike_before_trigger_threshold < time_difference_left && time_difference_left < time_to_nearest_heelstrike_after_trigger_threshold
                    	% left heelstrike is acceptable
                        left_heelstrike_acceptable = true;
                    else
                        left_heelstrike_acceptable = false;
                    end
                    
                    % is the closest right heelstrike within the acceptable interval?
                    closest_right_heelstrike = right_touchdown_times(index_right);
                    time_difference_right = closest_right_heelstrike - trigger_times(i_trigger); % where does the closest right heelstrike lie relative to the trigger?
                    if -time_to_nearest_heelstrike_before_trigger_threshold < time_difference_right && time_difference_right < time_to_nearest_heelstrike_after_trigger_threshold
                    	% right heelstrike is acceptable
                        right_heelstrike_acceptable = true;
                    else
                        right_heelstrike_acceptable = false;
                    end
                    
                    % accept the acceptable one
                    if left_heelstrike_acceptable && ~right_heelstrike_acceptable
                        % triggered by left heelstrike
                        trigger_foot = 'left';
                        closest_heelstrike_distance_times(i_trigger) = time_difference_left;
                    elseif ~left_heelstrike_acceptable && right_heelstrike_acceptable
                        % triggered by right heelstrike
                        trigger_foot = 'right';
                        closest_heelstrike_distance_times(i_trigger) = time_difference_right;
                    elseif left_heelstrike_acceptable && right_heelstrike_acceptable
                        trigger_foot = 'unclear';
                        removal_flags(i_trigger) = 1;
                    elseif ~left_heelstrike_acceptable && ~right_heelstrike_acceptable
                        trigger_foot = 'unclear';
                        removal_flags(i_trigger) = 1;
                    end                    
                    
                    % extract relevant events in order
                    if strcmp(trigger_foot, 'left')
                        if length(left_touchdown_times) < index_left + 2 || removal_flags(i_trigger) == 1
                            % data doesn't include the required number of steps after the trigger
                            removal_flags(i_trigger) = 1;
                            right_foot_heelstrike_0 = NaN;
                            right_foot_heelstrike_plus_1 = NaN;
                            right_foot_heelstrike_plus_2 = NaN;
                            left_foot_heelstrike_0 = NaN;
                            left_foot_heelstrike_plus_1 = NaN;
                            left_foot_heelstrike_plus_2 = NaN;
                            
                        else
                            left_foot_heelstrike_0  = left_touchdown_times(index_left);
                            left_foot_heelstrike_1  = left_touchdown_times(index_left+1);
                            left_foot_heelstrike_2  = left_touchdown_times(index_left+2);
                            left_foot_pushoff_0     = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_0));
                            left_foot_pushoff_1     = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_1));
                            left_foot_pushoff_2     = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_2));

                            right_foot_heelstrike_0 = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_0));
                            right_foot_heelstrike_1 = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_1));
                            right_foot_heelstrike_2 = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_2));
                            right_foot_pushoff_0    = max(right_pushoff_times(right_pushoff_times <= left_foot_pushoff_0));
                            right_foot_pushoff_1    = max(right_pushoff_times(right_pushoff_times <= left_foot_pushoff_1));
                            right_foot_pushoff_2    = max(right_pushoff_times(right_pushoff_times <= left_foot_pushoff_2));

                            % notify if events are not sorted properly
                            event_order = [ ...
                                            left_foot_heelstrike_0 right_foot_pushoff_0 right_foot_heelstrike_0 left_foot_pushoff_0 ...
                                            left_foot_heelstrike_1 right_foot_pushoff_1 right_foot_heelstrike_1 left_foot_pushoff_1 ...
                                            left_foot_heelstrike_2  ...
                                          ];
                            if ~issorted(event_order)
                                disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(time_stimulus(trigger_indices_labview(i_trigger)))]);
                            end
                            % check check
                            if visualize
                                plot([left_foot_heelstrike_minus_1 left_foot_heelstrike_0 left_foot_heelstrike_1 left_foot_heelstrike_2], [0 0 0 0]-0.01, 'v', 'linewidth', 3);
                                plot([left_foot_pushoff_minus_1  left_foot_pushoff_0 left_foot_pushoff_1 left_foot_pushoff_2], [0 0 0 0]-0.01, '^', 'linewidth', 3);
                                plot([right_foot_heelstrike_minus_1 right_foot_heelstrike_0 right_foot_heelstrike_1 right_foot_heelstrike_2], [0 0 0 0]+0.01, 'v', 'linewidth', 3);
                                plot([right_foot_pushoff_minus_1  right_foot_pushoff_0 right_foot_pushoff_1 right_foot_pushoff_2], [0 0 0 0]+0.01, '^', 'linewidth', 3);
                                % note: this can crash if one of thse events is empty, because we are plotting before we
                                % have checked that
                            end                
                        end
                    elseif strcmp(trigger_foot, 'right')
                        if length(right_touchdown_times) < index_right + 2 || removal_flags(i_trigger) == 1
                            % data doesn't include the required number of steps after the trigger
                            removal_flags(i_trigger) = 1;
                            left_foot_heelstrike_0 = NaN;
                            left_foot_heelstrike_1 = NaN;
                            left_foot_heelstrike_2 = NaN;
                            right_foot_heelstrike_0 = NaN;
                            right_foot_heelstrike_1 = NaN;
                            right_foot_heelstrike_2 = NaN;
                            
                        else
                            right_foot_heelstrike_0 = right_touchdown_times(index_right);
                            right_foot_heelstrike_1 = right_touchdown_times(index_right+1);
                            right_foot_heelstrike_2 = right_touchdown_times(index_right+2);
                            
                            right_foot_pushoff_0    = min(right_pushoff_times(right_pushoff_times >= right_foot_heelstrike_0));
                            right_foot_pushoff_1    = min(right_pushoff_times(right_pushoff_times >= right_foot_heelstrike_1));
                            right_foot_pushoff_2    = min(right_pushoff_times(right_pushoff_times >= right_foot_heelstrike_2));

                            left_foot_heelstrike_0  = min(left_touchdown_times(left_touchdown_times >= right_foot_heelstrike_0));
                            left_foot_heelstrike_1  = min(left_touchdown_times(left_touchdown_times >= right_foot_heelstrike_1));
                            left_foot_heelstrike_2  = min(left_touchdown_times(left_touchdown_times >= right_foot_heelstrike_2));
                            left_foot_pushoff_0     = max(left_pushoff_times(left_pushoff_times <= right_foot_pushoff_0));
                            left_foot_pushoff_1     = max(left_pushoff_times(left_pushoff_times <= right_foot_pushoff_1));
                            left_foot_pushoff_2     = max(left_pushoff_times(left_pushoff_times <= right_foot_pushoff_2));

                            % notify if events are not sorted properly
                            event_order = [ ...
                                            right_foot_heelstrike_0 left_foot_pushoff_0 left_foot_heelstrike_0 right_foot_pushoff_0 ...
                                            right_foot_heelstrike_1 left_foot_pushoff_1 left_foot_heelstrike_1 right_foot_pushoff_1 ...
                                            right_foot_heelstrike_2 ...
                                          ];
                            if ~issorted(event_order)
                                disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(time_stimulus(trigger_indices_labview(i_trigger)))]);
                            end

                            if visualize
                                plot([left_foot_heelstrike_0 left_foot_heelstrike_1 left_foot_heelstrike_2], [0 0 0]-0.01, 'v', 'linewidth', 2, 'color', 'r');
                                plot([left_foot_pushoff_0 left_foot_pushoff_1 left_foot_pushoff_2], [0 0 0]-0.01, '^', 'linewidth', 2, 'color', 'g');
                                plot([right_foot_heelstrike_0 right_foot_heelstrike_1 right_foot_heelstrike_2], [0 0 0]+0.01, 'v', 'linewidth', 2, 'color', 'r');
                                plot([right_foot_pushoff_0 right_foot_pushoff_1 right_foot_pushoff_2], [0 0 0]+0.01, '^', 'linewidth', 2, 'color', 'g');
                            end            
                        end            
                    else
                        trigger_foot = 'unclear';
                        disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_stimulus(trigger_indices_labview(i_trigger))) ' - triggering heelstrike unclear']);
                        left_foot_heelstrike_0  = 0;
                        left_foot_heelstrike_1  = 0;
                        left_foot_heelstrike_2  = 0;
                        left_foot_pushoff_0     = 0;
                        left_foot_pushoff_1     = 0;
                        left_foot_pushoff_2     = 0;

                        right_foot_heelstrike_0 = 0;
                        right_foot_heelstrike_1 = 0;
                        right_foot_heelstrike_2 = 0;
                        right_foot_pushoff_0    = 0;
                        right_foot_pushoff_1    = 0;
                        right_foot_pushoff_2    = 0;

                        removal_flags(i_trigger) = 1;
                    end
                    
                    % flag for removal if not all events are present
                    if any ...
                         ( ...
                           [ ...
                             isempty(left_foot_heelstrike_0) isempty(left_foot_heelstrike_1) isempty(left_foot_heelstrike_2) ...
                             isempty(right_foot_heelstrike_0) isempty(right_foot_heelstrike_1) isempty(right_foot_heelstrike_2) ...
                             isempty(left_foot_pushoff_0) isempty(left_foot_pushoff_1) isempty(left_foot_pushoff_2) ...
                             isempty(right_foot_pushoff_0) isempty(right_foot_pushoff_1) isempty(right_foot_pushoff_2) ...
                           ] ...
                         ) ...
                       || removal_flags(i_trigger) == 1
                        % not all events are present, flag for removal
                        removal_flags(i_trigger) = 1;
                        left_foot_heelstrike_0 = NaN;
                        left_foot_heelstrike_1 = NaN;
                        left_foot_heelstrike_2 = NaN;
                        right_foot_heelstrike_0 = NaN;
                        right_foot_heelstrike_1 = NaN;
                        right_foot_heelstrike_2 = NaN;
                    end
                    
                    % collect event times to form stretches
                    if ~removal_flags(i_trigger) == 1
                        if strcmp(trigger_foot, 'right')
                            if bands_per_stretch == 4
                                stretch_times(i_trigger, :) = [right_foot_heelstrike_0 left_foot_heelstrike_0 right_foot_heelstrike_1 left_foot_heelstrike_1 right_foot_heelstrike_2];
                                stance_foot_data(i_trigger, :) = {'STANCE_RIGHT', 'STANCE_LEFT', 'STANCE_RIGHT', 'STANCE_LEFT'};
                            elseif bands_per_stretch == 2
                                stretch_times(i_trigger, :) = [right_foot_heelstrike_0 left_foot_heelstrike_0 right_foot_heelstrike_1];
                                stance_foot_data(i_trigger, :) = {'STANCE_RIGHT', 'STANCE_LEFT'};
                            end
                            
                            trigger_foot_list{i_trigger} = 'TRIGGER_RIGHT';
                        end
                        if strcmp(trigger_foot, 'left')
                            if bands_per_stretch == 4
                                stretch_times(i_trigger, :) = [left_foot_heelstrike_0 right_foot_heelstrike_0  left_foot_heelstrike_1 right_foot_heelstrike_1 left_foot_heelstrike_2];
                                stance_foot_data(i_trigger, :) = {'STANCE_LEFT', 'STANCE_RIGHT', 'STANCE_LEFT', 'STANCE_RIGHT'};
                            elseif bands_per_stretch == 2
                                stretch_times(i_trigger, :) = [left_foot_heelstrike_0 right_foot_heelstrike_0  left_foot_heelstrike_1];
                                stance_foot_data(i_trigger, :) = {'STANCE_LEFT', 'STANCE_RIGHT'};
                            end
                            trigger_foot_list{i_trigger} = 'TRIGGER_LEFT';
                        end
                        % visualize
                        if visualize
                            plot([stretch_times(i_trigger, 1) stretch_times(i_trigger, 2)], [-1 1]*-0.01, 'color', 'r', 'linewidth', 3);
                            plot([stretch_times(i_trigger, 2) stretch_times(i_trigger, 3)], [-1 1]*+0.01, 'color', 'g', 'linewidth', 3);
                            plot([stretch_times(i_trigger, 3) stretch_times(i_trigger, 4)], [-1 1]*-0.01, 'color', 'b', 'linewidth', 3);
                            plot([stretch_times(i_trigger, 4) stretch_times(i_trigger, 5)], [-1 1]*+0.01, 'color', 'm', 'linewidth', 3);
                        end
                    end
                    
                    % determine amplitude of this stimulus
                    if strcmp(experimental_paradigm, 'Vision')
                        amplitude = current_acceleration_trajectory(trigger_indices_labview(i_trigger));
                        amplitude_list{i_trigger} = num2str(amplitude);
                    end
                end
                    
                
                
                % remove flagged triggers
                unflagged_indices = ~removal_flags;
                trigger_times = trigger_times(unflagged_indices);
                trigger_indices_labview = trigger_indices_labview(unflagged_indices, :);
                stretch_times = stretch_times(unflagged_indices, :);
                stance_foot_data = stance_foot_data(unflagged_indices, :);
                stimulus_list = stimulus_list(unflagged_indices, :);
                trigger_foot_list = trigger_foot_list(unflagged_indices, :);
                stretch_times_minus_1 = stretch_times_minus_1(unflagged_indices, :);
                stance_foot_data_minus_1 = stance_foot_data_minus_1(unflagged_indices, :);
                
                % determine direction
                direction_list = cell(size(trigger_foot_list));
                for i_stretch = 1 : length(trigger_foot_list)
                    if strcmp(trigger_foot_list{i_stretch}, 'TRIGGER_RIGHT')
                        if strcmp(stimulus_list{i_stretch}, 'STIM_RIGHT')
                            direction_list{i_stretch} = 'STIM_TOWARDS';
                        end
                        if strcmp(stimulus_list{i_stretch}, 'STIM_LEFT')
                            direction_list{i_stretch} = 'STIM_AWAY';
                        end
                        if strcmp(stimulus_list{i_stretch}, 'STIM_NONE')
                            direction_list{i_stretch} = 'STIM_NONE';
                        end
                    end
                    if strcmp(trigger_foot_list{i_stretch}, 'TRIGGER_LEFT')
                        if strcmp(stimulus_list{i_stretch}, 'STIM_RIGHT')
                            direction_list{i_stretch} = 'STIM_AWAY';
                        end
                        if strcmp(stimulus_list{i_stretch}, 'STIM_LEFT')
                            direction_list{i_stretch} = 'STIM_TOWARDS';
                        end
                        if strcmp(stimulus_list{i_stretch}, 'STIM_NONE')
                            direction_list{i_stretch} = 'STIM_NONE';
                        end
                    end
                end
                
                if strcmp(experimental_paradigm, 'OculusLaneRestriction')
                    % determine where the "no step zone" was at stretch
                    % trigger
                    zone_side_list = cell(size(trigger_foot_list));
                    zone_direction_list = cell(size(trigger_foot_list));
                    scene_translation_mod100 = mod(scene_translation_trajectory + 25, 100); %TO DO the origin of the scene is +25 relative to the end of the virtual objects
                    
                    for i_stretch = 1:length(trigger_indices_labview)
                        VR_trigger_position = scene_translation_mod100(trigger_indices_labview(i_stretch));
                        [~,scene_translation_mod100_index] = min(abs(virtual_object_ap_location - VR_trigger_position));
                       
                        
                        %% CHECK THIS %
                        % 2 = NO STEP ZONE RIGHT
                        % 0 = NO STEP ZONE LEFT
                         if virtual_object_ml_location(scene_translation_mod100_index) == 2
                            zone_side_list{i_stretch} = 'STIM_ZONE_LEFT';
                         elseif virtual_object_ml_location(scene_translation_mod100_index) == 0
                            zone_side_list{i_stretch} = 'STIM_ZONE_RIGHT';
                         end
                        if (virtual_object_ml_location(scene_translation_mod100_index) == 2 && strcmp(trigger_foot_list{i_stretch}, 'TRIGGER_LEFT') && strcmp(direction_list{i_stretch}, 'STIM_TOWARDS')) || ...
                                (virtual_object_ml_location(scene_translation_mod100_index) == 2 && strcmp(trigger_foot_list{i_stretch}, 'TRIGGER_RIGHT') && strcmp(direction_list{i_stretch}, 'STIM_AWAY')) || ...
                                (virtual_object_ml_location(scene_translation_mod100_index) == 0 && strcmp(trigger_foot_list{i_stretch}, 'TRIGGER_RIGHT') && strcmp(direction_list{i_stretch}, 'STIM_TOWARDS')) ||...
                                (virtual_object_ml_location(scene_translation_mod100_index) == 0 && strcmp(trigger_foot_list{i_stretch}, 'TRIGGER_LEFT') && strcmp(direction_list{i_stretch}, 'STIM_AWAY'))
                            zone_direction_list{i_stretch} = 'STIM_ZONE_TOWARDS';
                        elseif (virtual_object_ml_location(scene_translation_mod100_index) == 0 && strcmp(trigger_foot_list{i_stretch}, 'TRIGGER_LEFT') && strcmp(direction_list{i_stretch}, 'STIM_TOWARDS')) || ...
                                (virtual_object_ml_location(scene_translation_mod100_index) == 0 && strcmp(trigger_foot_list{i_stretch}, 'TRIGGER_RIGHT') && strcmp(direction_list{i_stretch}, 'STIM_AWAY')) || ...
                                (virtual_object_ml_location(scene_translation_mod100_index) == 2 && strcmp(trigger_foot_list{i_stretch}, 'TRIGGER_RIGHT') && strcmp(direction_list{i_stretch}, 'STIM_TOWARDS')) || ...
                                (virtual_object_ml_location(scene_translation_mod100_index) == 2 && strcmp(trigger_foot_list{i_stretch}, 'TRIGGER_LEFT') && strcmp(direction_list{i_stretch}, 'STIM_AWAY'))
                            zone_direction_list{i_stretch} = 'STIM_ZONE_AWAY';
                        else
                            zone_direction_list{i_stretch} = 'STIM_NONE';
                        end
                    end
                end
                    
                    
                % put in placeholder for group
                group_list = cell(size(direction_list));
                [group_list{:}] = deal('to be determined');
                
                % add cadence list
                if strcmp(experimental_paradigm, 'CadenceVision') || strcmp(experimental_paradigm, 'CadenceGVS')
                    this_trial_cadence = '~';
                    
                    % determine cadence for this trial
                    this_trial_type = condition_list{i_condition};
                    this_trial_number = i_trial;
                    this_trial_protocol_index = find(strcmp(protocol_data.trial_type, this_trial_type) & protocol_data.trial_number==this_trial_number);
                    if isempty(this_trial_protocol_index)
                        error(['Condition ' this_trial_type ', trial ' num2str(this_trial_number) ' - not found in protocol'])
                    end
                    this_trial_cadence = protocol_data.metronome_cadence(this_trial_protocol_index);
                    
                    cadence_list = cell(size(direction_list));
                    [cadence_list{:}] = deal([num2str(this_trial_cadence) 'BPM']);
                    conditions_trial.cadence_list = cadence_list;
                end
                
                if strcmp(experimental_paradigm, 'FatigueGVS')
                    fatigue_list = cell(size(direction_list));
                    if ismember(i_trial, fatigue_trials)
                        [fatigue_list{:}] = deal('FATIGUED');
                    else
                        [fatigue_list{:}] = deal('UNFATIGUED');
                    end
                    conditions_trial.fatigue_list = fatigue_list;
                end
                if strcmp(experimental_paradigm, 'CognitiveLoadVision') || strcmp(experimental_paradigm, 'CognitiveLoadGvs')
                    cognitive_load_list = cell(size(direction_list));
                    if ismember(i_trial, back_7_trials)
                        [cognitive_load_list{:}] = deal('BACK_7');
                    else
                        [cognitive_load_list{:}] = deal('NO_LOAD');
                    end
                    conditions_trial.cognitive_load_list = cognitive_load_list;
                end
                
                if exist('affected_side')
                    for i_stretch = 1 : length(trigger_foot_list)
                        if strcmp(affected_side, 'Left')
                            if strcmp(trigger_foot_list(i_stretch), 'TRIGGER_LEFT')
                                condition_affected_stancefoot_list{i_stretch} = 'TRIGGER_AFFECTED';
                            elseif strcmp(trigger_foot_list(i_stretch), 'TRIGGER_RIGHT')
                                condition_affected_stancefoot_list{i_stretch} = 'TRIGGER_UNAFFECTED';
                            end
                        elseif strcmp(affected_side, 'Right')
                            if strcmp(trigger_foot_list(i_stretch), 'TRIGGER_RIGHT')
                                condition_affected_stancefoot_list{i_stretch} = 'TRIGGER_AFFECTED';
                            elseif strcmp(trigger_foot_list(i_stretch), 'TRIGGER_LEFT')
                                condition_affected_stancefoot_list{i_stretch} = 'TRIGGER_UNAFFECTED';
                            end
                        end
                    end
                    conditions_trial.affected_stancefoot_list = condition_affected_stancefoot_list';
                end
                
                % restructure for saving
                conditions_trial.stimulus_list = stimulus_list;
                conditions_trial.amplitude_list = amplitude_list;
                conditions_trial.trigger_foot_list = trigger_foot_list;
                conditions_trial.direction_list = direction_list;
                conditions_trial.group_list = group_list;
                if strcmp(experimental_paradigm, 'OculusLaneRestriction')
                    conditions_trial.zone_side_list = zone_side_list;
                    conditions_trial.zone_direction_list = zone_direction_list;
                end
                event_variables_to_save.stretch_times = stretch_times;
                event_variables_to_save.stance_foot_data = stance_foot_data;
                event_variables_to_save.stretch_times_minus_1 = stretch_times_minus_1;
                event_variables_to_save.stance_foot_data_minus_1 = stance_foot_data_minus_1;
            end
            
            if strcmp(condition_stimulus, 'VISUAL')
                bands_per_stretch = 2;
                
                number_of_triggers = length(trigger_indices_mocap);
                closest_heelstrike_distance_times = zeros(number_of_triggers, 1);
                removal_flags = zeros(number_of_triggers, 1);
                
                stretch_times_stim = zeros(number_of_triggers, bands_per_stretch+1);
                stance_foot_data_stim = cell(number_of_triggers, bands_per_stretch);
                stimulus_list_stim = cell(number_of_triggers, 1); % stimulus STIM_LEFT, STIM_RIGHT or STIM_NONE
                trigger_foot_list_stim = cell(number_of_triggers, 1); % triggering foot TRIGGER_LEFT or TRIGGER_RIGHT
%                 delay_list_stim = cell(number_of_triggers, 1);
                
                stretch_times_ctrl = zeros(number_of_triggers, bands_per_stretch+1);
                stance_foot_data_ctrl = cell(number_of_triggers, bands_per_stretch);
                stimulus_list_ctrl = cell(number_of_triggers, 1); % stimulus STIM_LEFT, STIM_RIGHT or STIM_NONE
                trigger_foot_list_ctrl = cell(number_of_triggers, 1); % triggering foot TRIGGER_LEFT or TRIGGER_RIGHT
%                 delay_list_ctrl = cell(number_of_triggers, 1);
                
                for i_trigger = 1 : number_of_triggers
                    % determine stimulus
                    if illusion_trajectory(stim_start_indices_labview(i_trigger)+1) > 0
                        stimulus_list_stim{i_trigger} = 'STIM_RIGHT';
                    end
                    if illusion_trajectory(stim_start_indices_labview(i_trigger)+1) < 0
                        stimulus_list_stim{i_trigger} = 'STIM_LEFT';
                    end
                    if illusion_trajectory(stim_start_indices_labview(i_trigger)+1) == 0
                        stimulus_list_stim{i_trigger} = 'STIM_NONE';
                    end
                    stimulus_list_ctrl{i_trigger} = 'STIM_NONE';
                    
                    % determine delay
%                     wait_time_stim = time_stimulus(stim_start_indices_labview(i_trigger)) - time_stimulus(trigger_indices_labview(i_trigger));
%                     delay_time_labels = study_settings.get('delay_time_labels');
%                     [~, wait_condition_index] = min(abs(study_settings.get('delay_times') - wait_time_stim));
%                     if iscell(study_settings.get('delay_time_labels'))
%                         delay_condition_label = delay_time_labels{wait_condition_index};
%                     else
%                         delay_condition_label = study_settings.get('delay_time_labels');
%                     end
%                     delay_list_stim{i_trigger} = delay_condition_label;
%                     delay_list_ctrl{i_trigger} = 'CONTROL';
                    
                    
                    % get closest heelstrike on either side
                    [~, index_left] = min(abs(left_touchdown_times - trigger_times(i_trigger)));
                    [~, index_right] = min(abs(right_touchdown_times - trigger_times(i_trigger)));
                    
                    % is the closest left heelstrike within the acceptable interval?
                    closest_left_heelstrike = left_touchdown_times(index_left);
                    time_difference_left = closest_left_heelstrike - trigger_times(i_trigger); % where does the closest left heelstrike lie relative to the trigger?
                    if -time_to_nearest_heelstrike_before_trigger_threshold < time_difference_left && time_difference_left < time_to_nearest_heelstrike_after_trigger_threshold
                    	% left heelstrike is acceptable
                        left_heelstrike_acceptable = true;
                    else
                        left_heelstrike_acceptable = false;
                    end
                    
                    % is the closest right heelstrike within the acceptable interval?
                    closest_right_heelstrike = right_touchdown_times(index_right);
                    time_difference_right = closest_right_heelstrike - trigger_times(i_trigger); % where does the closest right heelstrike lie relative to the trigger?
                    if -time_to_nearest_heelstrike_before_trigger_threshold < time_difference_right && time_difference_right < time_to_nearest_heelstrike_after_trigger_threshold
                    	% right heelstrike is acceptable
                        right_heelstrike_acceptable = true;
                    else
                        right_heelstrike_acceptable = false;
                    end
                    
                    % accept the acceptable one
                    if left_heelstrike_acceptable && ~right_heelstrike_acceptable
                        % triggered by left heelstrike
                        trigger_foot = 'left';
                        closest_heelstrike_distance_times(i_trigger) = time_difference_left;
                    elseif ~left_heelstrike_acceptable && right_heelstrike_acceptable
                        % triggered by right heelstrike
                        trigger_foot = 'right';
                        closest_heelstrike_distance_times(i_trigger) = time_difference_right;
                    elseif left_heelstrike_acceptable && right_heelstrike_acceptable
                        trigger_foot = 'unclear';
                        removal_flags(i_trigger) = 1;
                    elseif ~left_heelstrike_acceptable && ~right_heelstrike_acceptable
                        trigger_foot = 'unclear';
                        removal_flags(i_trigger) = 1;
                    end                    
                    
                    % extract relevant events in order
                    if strcmp(trigger_foot, 'left')
                        if length(left_touchdown_times) < index_left + 1 || removal_flags(i_trigger) == 1 || index_left == 1
                            % data doesn't include the required number of steps after the trigger
                            removal_flags(i_trigger) = 1;
                            left_foot_heelstrike_0  = NaN;
                            left_foot_heelstrike_1  = NaN;
                            left_foot_pushoff_0     = NaN;
                            right_foot_heelstrike_0 = NaN;
                            right_foot_pushoff_0    = NaN;
                        else
                            left_foot_heelstrike_pre  = left_touchdown_times(index_left-1);
                            left_foot_heelstrike_0  = left_touchdown_times(index_left);
                            left_foot_heelstrike_1  = left_touchdown_times(index_left+1);
                            
                            left_foot_pushoff_pre     = max(left_pushoff_times(left_pushoff_times < left_foot_heelstrike_0));
                            left_foot_pushoff_0     = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_0));
                            
                            right_foot_heelstrike_pre = max(right_touchdown_times(right_touchdown_times < left_foot_heelstrike_0));
                            right_foot_heelstrike_0 = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_0));
                            right_foot_pushoff_pre    = max(right_pushoff_times(right_pushoff_times <= left_foot_heelstrike_0));
                            right_foot_pushoff_0    = max(right_pushoff_times(right_pushoff_times <= left_foot_pushoff_0));

                            % notify if events are not sorted properly
                            if ~issorted ...
                                  ( ...
                                    [ ...
                                      left_foot_heelstrike_pre right_foot_pushoff_pre right_foot_heelstrike_pre left_foot_pushoff_pre ...
                                      left_foot_heelstrike_0 right_foot_pushoff_0 right_foot_heelstrike_0 left_foot_pushoff_0 ...
                                      left_foot_heelstrike_1 ...
                                    ] ...
                                  )
                                disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(time_stimulus(trigger_indices_labview(i_trigger)))]);
                            end
                        end
                    elseif strcmp(trigger_foot, 'right')
                        if length(right_touchdown_times) < index_right + 1 || removal_flags(i_trigger) == 1 || index_right == 1
                            % data doesn't include the required number of steps after the trigger
                            removal_flags(i_trigger) = 1;
                            right_foot_heelstrike_pre = NaN;
                            right_foot_heelstrike_0 = NaN;
                            right_foot_heelstrike_1 = NaN;
                            
                            right_foot_pushoff_pre  = NaN;
                            right_foot_pushoff_0    = NaN;

                            left_foot_heelstrike_pre  = NaN;
                            left_foot_heelstrike_0  = NaN;
                            left_foot_pushoff_pre     = NaN;
                            left_foot_pushoff_0     = NaN;
                        else
                            right_foot_heelstrike_pre = right_touchdown_times(index_right-1);
                            right_foot_heelstrike_0 = right_touchdown_times(index_right);
                            right_foot_heelstrike_1 = right_touchdown_times(index_right+1);
                            
                            right_foot_pushoff_pre    = max(right_pushoff_times(right_pushoff_times < right_foot_heelstrike_0));
                            right_foot_pushoff_0    = min(right_pushoff_times(right_pushoff_times >= right_foot_heelstrike_0));

                            left_foot_heelstrike_pre  = max(left_touchdown_times(left_touchdown_times < right_foot_heelstrike_0));
                            left_foot_heelstrike_0  = min(left_touchdown_times(left_touchdown_times >= right_foot_heelstrike_0));
                            left_foot_pushoff_pre     = max(left_pushoff_times(left_pushoff_times <= right_foot_heelstrike_0));
                            left_foot_pushoff_0     = max(left_pushoff_times(left_pushoff_times <= right_foot_pushoff_0));

                            % notify if events are not sorted properly
                            if ~issorted ...
                                  ( ...
                                    [ ...
                                      right_foot_heelstrike_pre left_foot_pushoff_pre left_foot_heelstrike_pre right_foot_pushoff_pre ...
                                      right_foot_heelstrike_0 left_foot_pushoff_0 left_foot_heelstrike_0 right_foot_pushoff_0 ...
                                      right_foot_heelstrike_1 ...
                                    ] ...
                                  )
                                disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(time_stimulus(trigger_indices_labview(i_trigger)))]);
                            end

                        end            
                    else
                        trigger_foot = 'unclear';
                        disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_stimulus(trigger_indices_labview(i_trigger))) ' - triggering heelstrike unclear']);
                        left_foot_heelstrike_0  = 0;
                        left_foot_heelstrike_1  = 0;
                        left_foot_pushoff_0     = 0;

                        right_foot_heelstrike_0 = 0;
                        right_foot_heelstrike_1 = 0;
                        right_foot_pushoff_0    = 0;

                        removal_flags(i_trigger) = 1;
                    end
                    
                    % collect event times to form stretches
                    if ~removal_flags(i_trigger) == 1
                        if strcmp(trigger_foot, 'right')
%                             stretch_times_stim(i_trigger, :) = [right_foot_heelstrike_0 left_foot_pushoff_0 left_foot_heelstrike_0 right_foot_pushoff_0 right_foot_heelstrike_1];
%                             stance_foot_data_stim(i_trigger, :) = {'STANCE_BOTH', 'STANCE_RIGHT', 'STANCE_BOTH', 'STANCE_LEFT'};
%                             trigger_foot_list_stim{i_trigger} = 'TRIGGER_RIGHT';
%                             
%                             stretch_times_ctrl(i_trigger, :) = [right_foot_heelstrike_pre left_foot_pushoff_pre left_foot_heelstrike_pre right_foot_pushoff_pre right_foot_heelstrike_0];
%                             stance_foot_data_ctrl(i_trigger, :) = {'STANCE_BOTH', 'STANCE_RIGHT', 'STANCE_BOTH', 'STANCE_LEFT'};
%                             trigger_foot_list_ctrl{i_trigger} = 'TRIGGER_RIGHT';

                            stretch_times_stim(i_trigger, :) = [right_foot_heelstrike_0 left_foot_heelstrike_0 right_foot_heelstrike_1];
                            stance_foot_data_stim(i_trigger, :) = {'STANCE_RIGHT', 'STANCE_LEFT'};
                            trigger_foot_list_stim{i_trigger} = 'TRIGGER_RIGHT';
                            
                            stretch_times_ctrl(i_trigger, :) = [right_foot_heelstrike_pre left_foot_heelstrike_pre right_foot_heelstrike_0];
                            stance_foot_data_ctrl(i_trigger, :) = {'STANCE_RIGHT', 'STANCE_LEFT'};
                            trigger_foot_list_ctrl{i_trigger} = 'TRIGGER_RIGHT';
                            
                        end
                        if strcmp(trigger_foot, 'left')
%                             stretch_times_stim(i_trigger, :) = [left_foot_heelstrike_0 right_foot_pushoff_0 right_foot_heelstrike_0 left_foot_pushoff_0 left_foot_heelstrike_1];
%                             stance_foot_data_stim(i_trigger, :) = {'STANCE_BOTH', 'STANCE_LEFT', 'STANCE_BOTH', 'STANCE_RIGHT'};
%                             trigger_foot_list_stim{i_trigger} = 'TRIGGER_LEFT';
%                             
%                             stretch_times_ctrl(i_trigger, :) = [left_foot_heelstrike_pre right_foot_pushoff_pre right_foot_heelstrike_pre left_foot_pushoff_pre left_foot_heelstrike_0];
%                             stance_foot_data_ctrl(i_trigger, :) = {'STANCE_BOTH', 'STANCE_LEFT', 'STANCE_BOTH', 'STANCE_RIGHT'};
%                             trigger_foot_list_ctrl{i_trigger} = 'TRIGGER_LEFT';
                            
                            stretch_times_stim(i_trigger, :) = [left_foot_heelstrike_0 right_foot_heelstrike_0 left_foot_heelstrike_1];
                            stance_foot_data_stim(i_trigger, :) = {'STANCE_LEFT', 'STANCE_RIGHT'};
                            trigger_foot_list_stim{i_trigger} = 'TRIGGER_LEFT';
                            
                            stretch_times_ctrl(i_trigger, :) = [left_foot_heelstrike_pre right_foot_heelstrike_pre left_foot_heelstrike_0];
                            stance_foot_data_ctrl(i_trigger, :) = {'STANCE_LEFT', 'STANCE_RIGHT'};
                            trigger_foot_list_ctrl{i_trigger} = 'TRIGGER_LEFT';
                            
                        end
                    end
                    
                end
                    
                % remove flagged triggers
                unflagged_indices = ~removal_flags;
                trigger_times = trigger_times(unflagged_indices);
                trigger_indices_labview = trigger_indices_labview(unflagged_indices, :);
                stretch_times_stim = stretch_times_stim(unflagged_indices, :);
                stance_foot_data_stim = stance_foot_data_stim(unflagged_indices, :);
                stimulus_list_stim = stimulus_list_stim(unflagged_indices, :);
%                 delay_list_stim = delay_list_stim(unflagged_indices, :);
                trigger_foot_list_stim = trigger_foot_list_stim(unflagged_indices, :);
                stretch_times_ctrl = stretch_times_ctrl(unflagged_indices, :);
                stance_foot_data_ctrl = stance_foot_data_ctrl(unflagged_indices, :);
                stimulus_list_ctrl = stimulus_list_ctrl(unflagged_indices, :);
%                 delay_list_ctrl = delay_list_ctrl(unflagged_indices, :);
                trigger_foot_list_ctrl = trigger_foot_list_ctrl(unflagged_indices, :);
                
                % merge stim and control
                stretch_times = [stretch_times_stim; stretch_times_ctrl];
                stance_foot_data = [stance_foot_data_stim; stance_foot_data_ctrl];
                stimulus_list = [stimulus_list_stim; stimulus_list_ctrl];
%                 delay_list = [delay_list_stim; delay_list_ctrl];
                trigger_foot_list = [trigger_foot_list_stim; trigger_foot_list_ctrl];
                
                % determine direction
                direction_list = cell(size(trigger_foot_list));
                for i_stretch = 1 : length(trigger_foot_list)
                    if strcmp(trigger_foot_list{i_stretch}, 'TRIGGER_RIGHT')
                        if strcmp(stimulus_list{i_stretch}, 'STIM_RIGHT')
                            direction_list{i_stretch} = 'STIM_TOWARDS';
                        end
                        if strcmp(stimulus_list{i_stretch}, 'STIM_LEFT')
                            direction_list{i_stretch} = 'STIM_AWAY';
                        end
                        if strcmp(stimulus_list{i_stretch}, 'STIM_NONE')
                            direction_list{i_stretch} = 'STIM_NONE';
                        end
                    end
                    if strcmp(trigger_foot_list{i_stretch}, 'TRIGGER_LEFT')
                        if strcmp(stimulus_list{i_stretch}, 'STIM_RIGHT')
                            direction_list{i_stretch} = 'STIM_AWAY';
                        end
                        if strcmp(stimulus_list{i_stretch}, 'STIM_LEFT')
                            direction_list{i_stretch} = 'STIM_TOWARDS';
                        end
                        if strcmp(stimulus_list{i_stretch}, 'STIM_NONE')
                            direction_list{i_stretch} = 'STIM_NONE';
                        end
                    end
                end
                    
                % put in placeholder for group
                group_list = cell(size(direction_list));
                [group_list{:}] = deal('to be determined');
                
                % restructure for saving
                conditions_trial = struct;
                conditions_trial.stimulus_list = stimulus_list;
                conditions_trial.trigger_foot_list = trigger_foot_list;
                conditions_trial.direction_list = direction_list;
                
%                 conditions_trial.group_list = group_list;
%                 conditions_trial.delay_list = delay_list;
                
                event_variables_to_save.stretch_times = stretch_times;
                event_variables_to_save.stance_foot_data = stance_foot_data;  
            end
            
            if strcmp(experimental_paradigm, 'GVS_old')
                bands_per_stretch = study_settings.get('strides_to_process') * 2;
%                 bands_per_stretch = 4;
                
                number_of_triggers = length(trigger_indices_mocap);
                closest_heelstrike_distance_times = zeros(number_of_triggers, 1);
                removal_flags = zeros(number_of_triggers, 1);
                
                stretch_times_stim = zeros(number_of_triggers, bands_per_stretch+1);
                stance_foot_data_stim = cell(number_of_triggers, bands_per_stretch);
                stimulus_list_stim = cell(number_of_triggers, 1); % stimulus STIM_LEFT, STIM_RIGHT or STIM_NONE
                trigger_foot_list_stim = cell(number_of_triggers, 1); % triggering foot TRIGGER_LEFT or TRIGGER_RIGHT
                delay_list_stim = cell(number_of_triggers, 1);
                
                stretch_times_ctrl = zeros(number_of_triggers, bands_per_stretch+1);
                stance_foot_data_ctrl = cell(number_of_triggers, bands_per_stretch);
                stimulus_list_ctrl = cell(number_of_triggers, 1); % stimulus STIM_LEFT, STIM_RIGHT or STIM_NONE
                trigger_foot_list_ctrl = cell(number_of_triggers, 1); % triggering foot TRIGGER_LEFT or TRIGGER_RIGHT
                delay_list_ctrl = cell(number_of_triggers, 1);
                
                for i_trigger = 1 : number_of_triggers
                    % determine stimulus
                    if illusion_trajectory(stim_start_indices_labview(i_trigger)+1) > 0
                        stimulus_list_stim{i_trigger} = 'STIM_RIGHT';
                    end
                    if illusion_trajectory(stim_start_indices_labview(i_trigger)+1) < 0
                        stimulus_list_stim{i_trigger} = 'STIM_LEFT';
                    end
                    if illusion_trajectory(stim_start_indices_labview(i_trigger)+1) == 0
                        stimulus_list_stim{i_trigger} = 'STIM_NONE';
                    end
                    stimulus_list_ctrl{i_trigger} = 'STIM_NONE';
                    
                    % determine delay
                    wait_time_stim = time_stimulus(stim_start_indices_labview(i_trigger)) - time_stimulus(trigger_indices_labview(i_trigger));
                    delay_time_labels = study_settings.get('delay_time_labels');
                    [~, wait_condition_index] = min(abs(study_settings.get('delay_times') - wait_time_stim));
                    if iscell(study_settings.get('delay_time_labels'))
                        delay_condition_label = delay_time_labels{wait_condition_index};
                    else
                        delay_condition_label = study_settings.get('delay_time_labels');
                    end
                    delay_list_stim{i_trigger} = delay_condition_label;
                    delay_list_ctrl{i_trigger} = 'CONTROL';
                    
                    
                    % get closest heelstrike on either side
                    [~, index_left] = min(abs(left_touchdown_times - trigger_times(i_trigger)));
                    [~, index_right] = min(abs(right_touchdown_times - trigger_times(i_trigger)));
                    
                    % is the closest left heelstrike within the acceptable interval?
                    closest_left_heelstrike = left_touchdown_times(index_left);
                    time_difference_left = closest_left_heelstrike - trigger_times(i_trigger); % where does the closest left heelstrike lie relative to the trigger?
                    if -time_to_nearest_heelstrike_before_trigger_threshold < time_difference_left && time_difference_left < time_to_nearest_heelstrike_after_trigger_threshold
                    	% left heelstrike is acceptable
                        left_heelstrike_acceptable = true;
                    else
                        left_heelstrike_acceptable = false;
                    end
                    
                    % is the closest right heelstrike within the acceptable interval?
                    closest_right_heelstrike = right_touchdown_times(index_right);
                    time_difference_right = closest_right_heelstrike - trigger_times(i_trigger); % where does the closest right heelstrike lie relative to the trigger?
                    if -time_to_nearest_heelstrike_before_trigger_threshold < time_difference_right && time_difference_right < time_to_nearest_heelstrike_after_trigger_threshold
                    	% right heelstrike is acceptable
                        right_heelstrike_acceptable = true;
                    else
                        right_heelstrike_acceptable = false;
                    end
                    
                    % accept the acceptable one
                    if left_heelstrike_acceptable && ~right_heelstrike_acceptable
                        % triggered by left heelstrike
                        trigger_foot = 'left';
                        closest_heelstrike_distance_times(i_trigger) = time_difference_left;
                    elseif ~left_heelstrike_acceptable && right_heelstrike_acceptable
                        % triggered by right heelstrike
                        trigger_foot = 'right';
                        closest_heelstrike_distance_times(i_trigger) = time_difference_right;
                    elseif left_heelstrike_acceptable && right_heelstrike_acceptable
                        trigger_foot = 'unclear';
                        removal_flags(i_trigger) = 1;
                    elseif ~left_heelstrike_acceptable && ~right_heelstrike_acceptable
                        trigger_foot = 'unclear';
                        removal_flags(i_trigger) = 1;
                    end                    
                    
                    % extract relevant events in order
                    if strcmp(trigger_foot, 'left')
                        if length(left_touchdown_times) < index_left + 2 || removal_flags(i_trigger) == 1
                            removal_flags(i_trigger) = 1;
                        else
                            left_foot_heelstrike_m2 = left_touchdown_times(index_left-2);
                            left_foot_heelstrike_m1 = left_touchdown_times(index_left-1);
                            left_foot_heelstrike_t0 = left_touchdown_times(index_left);
                            left_foot_heelstrike_p1 = left_touchdown_times(index_left+1);
                            left_foot_heelstrike_p2 = left_touchdown_times(index_left+2);
                            
                            left_foot_pushoff_m2    = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_m2));
                            left_foot_pushoff_m1    = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_m1));
                            left_foot_pushoff_t0    = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_t0));
                            left_foot_pushoff_p1    = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_p1));
                            
                            right_foot_heelstrike_m2  = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_m2));
                            right_foot_heelstrike_m1  = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_m1));
                            right_foot_heelstrike_t0  = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_t0));
                            right_foot_heelstrike_p1  = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_p1));
                            
                            right_foot_pushoff_m2     = min(right_pushoff_times(right_pushoff_times >= left_foot_heelstrike_m2));
                            right_foot_pushoff_m1     = min(right_pushoff_times(right_pushoff_times >= left_foot_heelstrike_m1));
                            right_foot_pushoff_t0     = min(right_pushoff_times(right_pushoff_times >= left_foot_heelstrike_t0));
                            right_foot_pushoff_p1     = min(right_pushoff_times(right_pushoff_times >= left_foot_heelstrike_p1));

                            % notify if events are not sorted properly
                            if ~issorted ...
                                  ( ...
                                    [ ...
                                      left_foot_heelstrike_m2 right_foot_pushoff_m2 right_foot_heelstrike_m2 left_foot_pushoff_m2 ...
                                      left_foot_heelstrike_m1 right_foot_pushoff_m1 right_foot_heelstrike_m1 left_foot_pushoff_m1 ...
                                      left_foot_heelstrike_t0 right_foot_pushoff_t0 right_foot_heelstrike_t0 left_foot_pushoff_t0 ...
                                      left_foot_heelstrike_p1 right_foot_pushoff_p1 right_foot_heelstrike_p1 left_foot_pushoff_p1 ...
                                      left_foot_heelstrike_p2 ...
                                    ] ...
                                  )
                                disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(time_stimulus(trigger_indices_labview(i_trigger)))]);
                            end

                        end
                    elseif strcmp(trigger_foot, 'right')
                        if length(right_touchdown_times) < index_left + 2 || removal_flags(i_trigger) == 1
                            removal_flags(i_trigger) = 1;
                        else
                            right_foot_heelstrike_m2 = right_touchdown_times(index_right-2);
                            right_foot_heelstrike_m1 = right_touchdown_times(index_right-1);
                            right_foot_heelstrike_t0 = right_touchdown_times(index_right);
                            right_foot_heelstrike_p1 = right_touchdown_times(index_right+1);
                            right_foot_heelstrike_p2 = right_touchdown_times(index_right+2);
                            
                            right_foot_pushoff_m2    = min(right_pushoff_times(right_pushoff_times >= right_foot_heelstrike_m2));
                            right_foot_pushoff_m1    = min(right_pushoff_times(right_pushoff_times >= right_foot_heelstrike_m1));
                            right_foot_pushoff_t0    = min(right_pushoff_times(right_pushoff_times >= right_foot_heelstrike_t0));
                            right_foot_pushoff_p1    = min(right_pushoff_times(right_pushoff_times >= right_foot_heelstrike_p1));
                            
                            left_foot_heelstrike_m2  = min(left_touchdown_times(left_touchdown_times >= right_foot_heelstrike_m2));
                            left_foot_heelstrike_m1  = min(left_touchdown_times(left_touchdown_times >= right_foot_heelstrike_m1));
                            left_foot_heelstrike_t0  = min(left_touchdown_times(left_touchdown_times >= right_foot_heelstrike_t0));
                            left_foot_heelstrike_p1  = min(left_touchdown_times(left_touchdown_times >= right_foot_heelstrike_p1));
                            
                            left_foot_pushoff_m2     = min(left_pushoff_times(left_pushoff_times >= right_foot_heelstrike_m2));
                            left_foot_pushoff_m1     = min(left_pushoff_times(left_pushoff_times >= right_foot_heelstrike_m1));
                            left_foot_pushoff_t0     = min(left_pushoff_times(left_pushoff_times >= right_foot_heelstrike_t0));
                            left_foot_pushoff_p1     = min(left_pushoff_times(left_pushoff_times >= right_foot_heelstrike_p1));

                            % notify if events are not sorted properly
                            if ~issorted ...
                                  ( ...
                                    [ ...
                                      right_foot_heelstrike_m2 left_foot_pushoff_m2 left_foot_heelstrike_m2 right_foot_pushoff_m2 ...
                                      right_foot_heelstrike_m1 left_foot_pushoff_m1 left_foot_heelstrike_m1 right_foot_pushoff_m1 ...
                                      right_foot_heelstrike_t0 left_foot_pushoff_t0 left_foot_heelstrike_t0 right_foot_pushoff_t0 ...
                                      right_foot_heelstrike_p1 left_foot_pushoff_p1 left_foot_heelstrike_p1 right_foot_pushoff_p1 ...
                                      right_foot_heelstrike_p2 ...
                                    ] ...
                                  )
                                disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(time_stimulus(trigger_indices_labview(i_trigger)))]);
                            end

                        end            
                    else
                        trigger_foot = 'unclear';
                        disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_stimulus(trigger_indices_labview(i_trigger))) ' - triggering heelstrike unclear']);
                        removal_flags(i_trigger) = 1;
                    end
                    
                    % collect event times to form stretches
                    if ~removal_flags(i_trigger) == 1
                        if strcmp(trigger_foot, 'right') && study_settings.get('strides_to_process') == 1
        
                            stretch_times_stim(i_trigger, :) = [right_foot_heelstrike_t0 left_foot_heelstrike_t0 right_foot_heelstrike_p1];
                            stance_foot_data_stim(i_trigger, :) = {'STANCE_RIGHT', 'STANCE_LEFT'};
                            trigger_foot_list_stim{i_trigger} = 'TRIGGER_RIGHT';
                            
                            stretch_times_ctrl(i_trigger, :) = [right_foot_heelstrike_p0 left_foot_heelstrike_p0 right_foot_heelstrike_t0];
                            stance_foot_data_ctrl(i_trigger, :) = {'STANCE_RIGHT', 'STANCE_LEFT'};
                            trigger_foot_list_ctrl{i_trigger} = 'TRIGGER_RIGHT';
                            
                        end
                        if strcmp(trigger_foot, 'left') && study_settings.get('strides_to_process') == 1
                            stretch_times_stim(i_trigger, :) = [left_foot_heelstrike_t0 right_foot_heelstrike_t0 left_foot_heelstrike_p1];
                            stance_foot_data_stim(i_trigger, :) = {'STANCE_LEFT', 'STANCE_RIGHT'};
                            trigger_foot_list_stim{i_trigger} = 'TRIGGER_LEFT';
                            
                            stretch_times_ctrl(i_trigger, :) = [left_foot_heelstrike_p0 right_foot_heelstrike_p0 left_foot_heelstrike_t0];
                            stance_foot_data_ctrl(i_trigger, :) = {'STANCE_LEFT', 'STANCE_RIGHT'};
                            trigger_foot_list_ctrl{i_trigger} = 'TRIGGER_LEFT';
                            
                        end
                        if strcmp(trigger_foot, 'right') && study_settings.get('strides_to_process') == 2
                            stretch_times_stim(i_trigger, :) = ...
                              [ ...
                                right_foot_heelstrike_t0 left_foot_heelstrike_t0 ...
                                right_foot_heelstrike_p1 left_foot_heelstrike_p1 ...
                                right_foot_heelstrike_p2 ...
                              ];
                            stance_foot_data_stim(i_trigger, :) = {'STANCE_RIGHT', 'STANCE_LEFT', 'STANCE_RIGHT', 'STANCE_LEFT'};
                            trigger_foot_list_stim{i_trigger} = 'TRIGGER_RIGHT';
                            
                            stretch_times_ctrl(i_trigger, :) = ...
                              [ ...
                                right_foot_heelstrike_m2 left_foot_heelstrike_m2 ...
                                right_foot_heelstrike_m1 left_foot_heelstrike_m1 ...
                                right_foot_heelstrike_t0 ...
                              ];
                            stance_foot_data_ctrl(i_trigger, :) = {'STANCE_RIGHT', 'STANCE_LEFT', 'STANCE_RIGHT', 'STANCE_LEFT'};
                            trigger_foot_list_ctrl{i_trigger} = 'TRIGGER_RIGHT';
                            
                        end
                        if strcmp(trigger_foot, 'left') && study_settings.get('strides_to_process') == 2
                            stretch_times_stim(i_trigger, :) = ...
                              [ ...
                                left_foot_heelstrike_t0 right_foot_heelstrike_t0 ...
                                left_foot_heelstrike_p1 right_foot_heelstrike_p1 ...
                                left_foot_heelstrike_p2 ...
                              ];
                            stance_foot_data_stim(i_trigger, :) = {'STANCE_LEFT', 'STANCE_RIGHT', 'STANCE_LEFT', 'STANCE_RIGHT'};
                            trigger_foot_list_stim{i_trigger} = 'TRIGGER_LEFT';
                            
                            stretch_times_ctrl(i_trigger, :) = ...
                              [ ...
                                left_foot_heelstrike_m2 right_foot_heelstrike_m2 ...
                                left_foot_heelstrike_m1 right_foot_heelstrike_m1 ...
                                left_foot_heelstrike_t0 ...
                              ];
                            stance_foot_data_ctrl(i_trigger, :) = {'STANCE_LEFT', 'STANCE_RIGHT', 'STANCE_LEFT', 'STANCE_RIGHT'};
                            trigger_foot_list_ctrl{i_trigger} = 'TRIGGER_LEFT';
                            
                        end
                    end
                    
                end
                    
                % remove flagged triggers
                unflagged_indices = ~removal_flags;
                trigger_times = trigger_times(unflagged_indices);
                trigger_indices_labview = trigger_indices_labview(unflagged_indices, :);
                stretch_times_stim = stretch_times_stim(unflagged_indices, :);
                stance_foot_data_stim = stance_foot_data_stim(unflagged_indices, :);
                stimulus_list_stim = stimulus_list_stim(unflagged_indices, :);
                delay_list_stim = delay_list_stim(unflagged_indices, :);
                trigger_foot_list_stim = trigger_foot_list_stim(unflagged_indices, :);
                stretch_times_ctrl = stretch_times_ctrl(unflagged_indices, :);
                stance_foot_data_ctrl = stance_foot_data_ctrl(unflagged_indices, :);
                stimulus_list_ctrl = stimulus_list_ctrl(unflagged_indices, :);
                delay_list_ctrl = delay_list_ctrl(unflagged_indices, :);
                trigger_foot_list_ctrl = trigger_foot_list_ctrl(unflagged_indices, :);
                
                % merge stim and control
                stretch_times = [stretch_times_stim; stretch_times_ctrl];
                stance_foot_data = [stance_foot_data_stim; stance_foot_data_ctrl];
                stimulus_list = [stimulus_list_stim; stimulus_list_ctrl];
                delay_list = [delay_list_stim; delay_list_ctrl];
                trigger_foot_list = [trigger_foot_list_stim; trigger_foot_list_ctrl];
                
                % determine direction
                direction_list = cell(size(trigger_foot_list));
                for i_stretch = 1 : length(trigger_foot_list)
                    if strcmp(trigger_foot_list{i_stretch}, 'TRIGGER_RIGHT')
                        if strcmp(stimulus_list{i_stretch}, 'STIM_RIGHT')
                            direction_list{i_stretch} = 'STIM_TOWARDS';
                        end
                        if strcmp(stimulus_list{i_stretch}, 'STIM_LEFT')
                            direction_list{i_stretch} = 'STIM_AWAY';
                        end
                        if strcmp(stimulus_list{i_stretch}, 'STIM_NONE')
                            direction_list{i_stretch} = 'STIM_NONE';
                        end
                    end
                    if strcmp(trigger_foot_list{i_stretch}, 'TRIGGER_LEFT')
                        if strcmp(stimulus_list{i_stretch}, 'STIM_RIGHT')
                            direction_list{i_stretch} = 'STIM_AWAY';
                        end
                        if strcmp(stimulus_list{i_stretch}, 'STIM_LEFT')
                            direction_list{i_stretch} = 'STIM_TOWARDS';
                        end
                        if strcmp(stimulus_list{i_stretch}, 'STIM_NONE')
                            direction_list{i_stretch} = 'STIM_NONE';
                        end
                    end
                end
                
                    
                    
                % put in placeholder for group
                group_list = cell(size(direction_list));
                [group_list{:}] = deal('to be determined');
                
                % restructure for saving
                conditions_trial = struct;
                conditions_trial.stimulus_list = stimulus_list;
                conditions_trial.trigger_foot_list = trigger_foot_list;
                conditions_trial.direction_list = direction_list;
                conditions_trial.group_list = group_list;
                conditions_trial.delay_list = delay_list;
                
                event_variables_to_save.stretch_times = stretch_times;
                event_variables_to_save.stance_foot_data = stance_foot_data;                
            end
            
            if strcmp(experimental_paradigm, 'Vision Stochastic')
                stim_frequency = loadConditionFromFile(conditions_file_name, 'frequency', i_trial);
                stim_amplitude = loadConditionFromFile(conditions_file_name, 'SD', i_trial);
                if i_trial < 12
                    block = 'FIRST';
                else
                    block = 'SECOND';
                end
                
                stance_foot_data_stretch = {'STANCE_BOTH', 'STANCE_LEFT', 'STANCE_BOTH', 'STANCE_RIGHT'};
                bands_per_stretch = length(stance_foot_data_stretch);
                
                left_touchdown_times_relevant = ...
                    left_touchdown_times ...
                      ( ...
                        left_touchdown_times > study_settings.get('analysis_start_time') ...
                        & left_touchdown_times < study_settings.get('analysis_end_time') ...
                      );
                stretch_start_times = left_touchdown_times_relevant(1:end-1);
                number_of_stretches = length(stretch_start_times);
                stretch_times = zeros(number_of_stretches, bands_per_stretch+1);
                removal_flags = false(number_of_stretches, 1);
                for i_stretch = 1 : number_of_stretches
                    this_stretch_start = stretch_start_times(i_stretch);
                    this_right_pushoff = min(right_pushoff_times(right_pushoff_times > this_stretch_start));
                    this_right_touchdown = min(right_touchdown_times(right_touchdown_times > this_stretch_start));
                    this_left_pushoff = min(left_pushoff_times(left_pushoff_times > this_stretch_start));
                    this_left_touchdown = min(left_touchdown_times(left_touchdown_times > this_stretch_start));
                    this_stretch_times = [this_stretch_start this_right_pushoff this_right_touchdown this_left_pushoff this_left_touchdown];
                    if ~issorted(this_stretch_times)
                        removal_flags(i_stretch) = 1;
                    end
                    stretch_times(i_stretch, :) = this_stretch_times;
                end
                stretch_times(removal_flags, :) = [];
                
                stance_foot_data = repmat(stance_foot_data_stretch, size(stretch_times, 1), 1);
                event_variables_to_save.stretch_times = stretch_times;
                event_variables_to_save.stance_foot_data = stance_foot_data;

                % conditions
                stim_frequency_list = repmat({['FRQ_' stim_frequency]}, size(stretch_times, 1), 1);
                stim_amplitude_list = repmat({['AMPL_' stim_amplitude]}, size(stretch_times, 1), 1);
                block_list = repmat({block}, size(stretch_times, 1), 1);
                conditions_trial = struct;
                conditions_trial.stim_frequency_list = stim_frequency_list;
                conditions_trial.stim_amplitude_list = stim_amplitude_list;
                conditions_trial.block_list = block_list;
            end
            
            if strcmp(experimental_paradigm, 'platformShift')
                stance_foot_data = {'STANCE_BOTH', 'STANCE_BOTH'};
                bands_per_stretch = 2;
                % determine start and end
                
                perturbation_list = cell(number_of_triggers, 1);
                
                for i_trigger = 1 : number_of_triggers
                    perturbation_list{i_trigger, 1} = condition_experimental;
                    stretch_times(i_trigger, :) = [perturbation_start_times(i_trigger); perturbation_end_times(i_trigger); response_end_times(i_trigger)];
                end
                
                     

                conditions_trial = struct;
                conditions_trial.perturbation_list = perturbation_list;
                
                event_variables_to_save.stretch_start_times = stretch_times(:, 1);
                event_variables_to_save.stretch_end_times = stretch_times(:, end);
                event_variables_to_save.stretch_times = stretch_times;
                event_variables_to_save.stance_foot_data = stance_foot_data;
            end

            % add subject
            condition_subject_list = cell(size(event_variables_to_save.stretch_times, 1), 1);
            for i_stretch = 1 : length(condition_subject_list)
                condition_subject_list{i_stretch} = subject_id;
            end
            conditions_trial.subject_list = condition_subject_list;
            
            % add gender
            condition_gender_list = cell(size(event_variables_to_save.stance_foot_data, 1), 1);
            for i_stretch = 1 : length(condition_gender_list)
                condition_gender_list{i_stretch} = gender;
            end
            conditions_trial.gender_list = condition_gender_list;

            %% remove stretches where important variables are missing

            % calculate variables that depend upon the step events to be identified correctly
            stretch_variables = study_settings.get('stretch_variables');

            % prune
            number_of_stretches = size(stretch_times, 1);
            removal_flags = zeros(number_of_stretches, 1);
            
            if study_settings.get('prune_step_time_outliers')
                    stretch_durations = stretch_times(:,2) - stretch_times(:,1);
                    stretch_duration_outlier_limits = median(stretch_durations) * [.75 1.25];
                    removal_flags(stretch_durations < stretch_duration_outlier_limits(1)) = 1;
                    removal_flags(stretch_durations > stretch_duration_outlier_limits(2)) = 1;
                    if any(removal_flags)
                        disp(['Removing a stretch due to innappropriate step length']);
                    end
            end
            
            if study_settings.get('prune_gaps_com', 1)
                for i_stretch = 1 : number_of_stretches
                    [~, start_index_mocap] = min(abs(time_marker - stretch_times(i_stretch, 1)));
                    [~, end_index_mocap] = min(abs(time_marker - stretch_times(i_stretch, end)));
                    if any(any(isnan(com_trajectories(start_index_mocap : end_index_mocap,:))))
                        removal_flags(i_stretch) = 1;
                        disp('Removing a stretch due to gaps in com trajectories')
                    end
                end
            end
            
            if study_settings.get('prune_gaps_angles', 1)
                for i_stretch = 1 : number_of_stretches
                    [~, start_index_mocap] = min(abs(time_marker - stretch_times(i_stretch, 1)));
                    [~, end_index_mocap] = min(abs(time_marker - stretch_times(i_stretch, end)));
                    if any(any(isnan(joint_angle_trajectories(start_index_mocap : end_index_mocap,:))))
                        removal_flags(i_stretch) = 1;
                        disp('Removing a stretch due to gaps in joint trajectories')
                    end
                end
            end
            
            if strcmp(experimental_paradigm, 'OculusLaneRestriction')

                if study_settings.get('prune_step_placements')
                    step_zone_delinquent_list = cell(size(trigger_foot_list));
                    for i_stretch = 1 : number_of_stretches
                        [~, trigger_start_index_mocap] = min(abs(time_marker - stretch_times(i_stretch, 1)));
                        [~, trigger_end_index_mocap] = min(abs(time_marker - stretch_times(i_stretch, 2)));
                        [~, remainder_start_index_mocap] = min(abs(time_marker - stretch_times(i_stretch, 1)));
                        [~, remainder_end_index_mocap] = min(abs(time_marker - stretch_times(i_stretch, 2)));
                        
                        
                        LHEE_marker_data = extractMarkerData(marker_trajectories, marker_labels, 'LHEE');
                        RHEE_marker_data = extractMarkerData(marker_trajectories, marker_labels, 'RHEE');
                        LTOE_marker_data = extractMarkerData(marker_trajectories, marker_labels, 'LTOE');
                        RTOE_marker_data = extractMarkerData(marker_trajectories, marker_labels, 'RTOE');
                        
                        
                         step_zone_delinquent_list{i_stretch} = 'NONE';
                        if strcmp(zone_side_list{i_stretch}, 'STIM_ZONE_LEFT')
                            threshold = -0.1835; % limit on left belt
                            if any(LHEE_marker_data(trigger_start_index_mocap : trigger_end_index_mocap,1) < threshold) || any(RHEE_marker_data(trigger_start_index_mocap : trigger_end_index_mocap,1) < threshold) ||...
                                    any(RTOE_marker_data(trigger_start_index_mocap : trigger_end_index_mocap,1) < threshold) || any(LTOE_marker_data(trigger_start_index_mocap : trigger_end_index_mocap,1) < threshold)
                                
                                step_zone_delinquent_list{i_stretch} = 'TRIGGER';
                                disp('Stretch flagged due stepping in No Step Zone')
                            end
                            if any(LHEE_marker_data(remainder_start_index_mocap : remainder_end_index_mocap,1) < threshold) || any(RHEE_marker_data(remainder_start_index_mocap : remainder_end_index_mocap,1) < threshold) ||...
                                    any(RTOE_marker_data(remainder_start_index_mocap : remainder_end_index_mocap,1) < threshold) || any(LTOE_marker_data(remainder_start_index_mocap : remainder_end_index_mocap,1) < threshold)
                                
                                if strcmp(step_zone_delinquent_list{i_stretch}, 'TRIGGER')
                                    step_zone_delinquent_list{i_stretch} = 'ALL';
                                else
                                    step_zone_delinquent_list{i_stretch} = 'LATER';
                                end
                                disp('Stretch flagged due stepping in No Step Zone')
                            end
                        end
                        if strcmp(zone_side_list{i_stretch}, 'STIM_ZONE_RIGHT')
                            threshold = 0.1835; % limit on right belt
                            if any(LHEE_marker_data(trigger_start_index_mocap : trigger_end_index_mocap,1) > threshold) || any(RHEE_marker_data(trigger_start_index_mocap : trigger_end_index_mocap,1) > threshold) || ...
                                    any(RTOE_marker_data(trigger_start_index_mocap : trigger_end_index_mocap,1) > threshold) || any(LTOE_marker_data(trigger_start_index_mocap : trigger_end_index_mocap,1) > threshold)
                                
                                step_zone_delinquent_list{i_stretch} = 'TRIGGER';
                                %                                 removal_flags(i_stretch) = 1;
                                disp('Stretch flagged due stepping in No Step Zone')
                            end
                            if any(LHEE_marker_data(remainder_start_index_mocap : remainder_end_index_mocap,1) > threshold) || any(RHEE_marker_data(remainder_start_index_mocap : remainder_end_index_mocap,1) > threshold) ||...
                                    any(RTOE_marker_data(remainder_start_index_mocap : remainder_end_index_mocap,1) > threshold) || any(LTOE_marker_data(remainder_start_index_mocap : remainder_end_index_mocap,1) > threshold)
                                
                                if strcmp(step_zone_delinquent_list(i_stretch), 'TRIGGER')
                                    step_zone_delinquent_list{i_stretch} = 'ALL';
                                else
                                    step_zone_delinquent_list{i_stretch} = 'LATER';
                                end
                                disp('Stretch flagged due stepping in No Step Zone')
                            end
                        end
                            
                        
%                         if strcmp(zone_side_list{i_stretch}, 'STIM_ZONE_RIGHT')
%                             threshold = 0.1835; % limit on right belt
%                             if any(LHEE_marker_data(trigger_start_index_mocap : trigger_end_index_mocap,1) > threshold)
%                                 removal_flags(i_stretch) = 1;
%                                 disp('Stretch flagged due stepping in No Step Zone')
%                             end
%                             if any(RHEE_marker_data(trigger_start_index_mocap : trigger_end_index_mocap,1) > threshold)
%                                 removal_flags(i_stretch) = 1;
%                                 disp('Stretch flagged due stepping in No Step Zone')
%                             end
%                             if any(RTOE_marker_data(trigger_start_index_mocap : trigger_end_index_mocap,1) > threshold)
%                                 removal_flags(i_stretch) = 1;
%                                 disp('Stretch flagged due stepping in No Step Zone')
%                             end
%                             if any(LTOE_marker_data(trigger_start_index_mocap : trigger_end_index_mocap,1) > threshold)
%                                 removal_flags(i_stretch) = 1;
%                                 disp('Stretch flagged due stepping in No Step Zone')
%                             end
%                         end
                        
                    end
                    
                        conditions_trial.step_zone_delinquent_list = step_zone_delinquent_list;
                end
            end
                
            %  check data availability for markers with non-zero weight
            marker_weights = study_settings.get('marker_weights');
            if isempty(marker_weights)
                marker_weights = ones(size(marker_labels));
            end
            for i_marker = 1 : 3 : length(marker_labels)
                this_marker_weight = 1; % 1 is default
                
                this_marker_label = marker_labels{i_marker}(1:end-2);
                if any(strcmp(marker_weights(:, 1), this_marker_label))
                    this_marker_weight = marker_weights{strcmp(marker_weights(:, 1), this_marker_label), 2};
                end
                
                if this_marker_weight == 1
                    this_marker_data = extractMarkerData(marker_trajectories, marker_labels, this_marker_label);
                    for i_stretch = 1 : number_of_stretches
                        [~, start_index] = min(abs(time_marker - stretch_times(i_stretch, 1)));
                        [~, end_index] = min(abs(time_marker - stretch_times(i_stretch, end)));
                        if any(any(isnan(this_marker_data(start_index : end_index, :))))
                            removal_flags(i_stretch) = 1;
                            disp(['Removing a stretch due to gap in marker "', this_marker_label '"']);
                        end 
                    end
                    
                end
                
            end
            
            % check ignore markers
            for i_stretch = 1 : number_of_stretches
                if ~isempty(ignore_times)
                    for i_ignore = 1 : length(ignore_times)
                        if stretch_times(i_stretch, 1) <= ignore_times(i_ignore) && ignore_times(i_ignore) <= stretch_times(i_stretch, end)
                            if ~removal_flags(i_stretch) == 1
                                removal_flags(i_stretch) = 1;
                                disp('Removing a stretch because of manually set ignore marker');
                            end
                        end
                    end
                end         
            end
          
            % remove flagged stretches
            unflagged_indices = ~removal_flags;
            event_variables_to_save_names = fieldnames(event_variables_to_save);
            for i_variable = 1 : length(event_variables_to_save_names)
                this_variable_name = event_variables_to_save_names{i_variable};
                
                evalstring = ['this_variable_data = event_variables_to_save.' this_variable_name ';'];
                eval(evalstring);
                this_variable_data = this_variable_data(unflagged_indices, :);
                
                evalstring = ['event_variables_to_save.' this_variable_name ' = this_variable_data;'];
                eval(evalstring);
            end
            
            conditions_trial_names = fieldnames(conditions_trial);
            for i_label = 1 : length(conditions_trial_names)
                this_condition_name = conditions_trial_names{i_label};
                
                evalstring = ['this_condition_data = conditions_trial.' this_condition_name ';'];
                eval(evalstring);
                this_condition_data = this_condition_data(unflagged_indices, :);
                
                evalstring = ['conditions_trial.' this_condition_name ' = this_condition_data;'];
                eval(evalstring);
            end
            
            %% save
            event_variables_to_save.conditions_trial = conditions_trial;
            event_variables_to_save.bands_per_stretch = bands_per_stretch;
            
            stretches_file_name = ['analysis' filesep makeFileName(collection_date, subject_id, condition_list{i_condition}, i_trial, 'relevantDataStretches')];
            saveDataToFile(stretches_file_name, event_variables_to_save);
            
            
            disp(['Finding Relevant Data Stretches: condition ' condition_list{i_condition} ', Trial ' num2str(i_trial) ' completed, found ' num2str(size(event_variables_to_save.stretch_times, 1)) ' relevant stretches, saved as ' stretches_file_name]);                
            % display average heelstrike distance
            disp(['Finding Relevant Data Stretches: condition ' condition_list{i_condition} ', Trial ' num2str(i_trial) ' completed, average discrepancy between trigger and heelstrike: ' num2str(mean(closest_heelstrike_distance_times))])
            
            % create a storage variable for latency values for each subject
%             trial_type_list = [trial_type_list; condition_list{i_condition}];
%             trial_number_list = {trial_number_list; (i_trial)};
            latency_value_list = [latency_value_list; closest_heelstrike_distance_times'];  
        end
    end
%     variables_to_save_latency_check.condition_list = trial_type_list;
    variables_to_save_latency_check.trial_number_list = trial_number_list;
    variables_to_save_latency_check.latency_value_list = latency_value_list;
    save_file_name = 'latency_check.mat';
    save(save_file_name, '-struct', 'variables_to_save_latency_check');
end




















