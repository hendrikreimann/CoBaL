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

    %% prepare
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
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            % create empty container for relevant data
            trial_data = struct;
            
            %% load data
            ignore_times = [];
            load(['analysis' filesep makeFileName(collection_date, subject_id, condition_list{i_condition}, i_trial, 'events')]);
            
            % determine experimental condition
            this_trial_type = condition_list{i_condition};
            trial_data.trial_type = this_trial_type;
            trial_data.trial_number = i_trial;
            
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
            
            % determine stimulus type
            condition_stimulus = study_settings.get('stimulus_condition', true);
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
            if strcmp(experimental_paradigm, 'GvsOverground')
                [first_forceplate_wrench_trajectory, time_forceplate] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'left_foot_wrench_world');
                trial_data.vertical_force_trajectory = first_forceplate_wrench_trajectory(:, 3);
                trial_data.time_forceplate = time_forceplate;
            end
            
            % stimulus data
            if strcmp(experimental_paradigm, 'GVS_old')
                load(['processed' filesep makeFileName(collection_date, subject_id, condition_list{i_condition}, i_trial, 'labviewData')]);
%                 GVS_out_trajectory = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'GVS_out_trajectory');
                GVS_stim_trajectory = GVS_out_trajectory + subject_settings.get('gvs_offset');
%                 [stimulus_state_trajectory, time_stimulus] = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'stimulus_state_trajectory');
                time_stimulus = time;
                
                trial_data.stimulus_state_trajectory = stimulus_state_trajectory;
                trial_data.time_stimulus = time_stimulus;
            end
            if strcmp(experimental_paradigm, 'Vision_old')
                % this if for TU data
                visual_scene_ml_translation_trajectory = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'visual_scene_ml_translation__trajectory'); %take note of the double "_"
                [stimulus_state_trajectory, time_stimulus] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'stimulus_state_trajectory');
                trial_data.stimulus_state_trajectory = stimulus_state_trajectory;
                trial_data.time_stimulus = time_stimulus;
            end
            if strcmp(experimental_paradigm, 'Vision') || strcmp(experimental_paradigm, 'CadenceVision') || strcmp(experimental_paradigm, 'SR_VisualStim') || strcmp(experimental_paradigm, 'CognitiveLoadVision')
                current_rotation_trajectory = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'visual_rotation_angle_trajectory');
                trial_data.current_acceleration_trajectory = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'visual_rotation_acceleration_trajectory');
                [stimulus_state_trajectory, time_stimulus] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'stimulus_state_trajectory');
                
                trial_data.stimulus_state_trajectory = stimulus_state_trajectory;
                trial_data.time_stimulus = time_stimulus;
            end
            if strcmp(experimental_paradigm, 'GVS') || strcmp(experimental_paradigm, 'CadenceGVS') || strcmp(experimental_paradigm, 'FatigueGVS') || strcmp(experimental_paradigm, 'CognitiveLoadGvs')
                gvs_trajectory = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'GVS_current_trajectory');
                [stimulus_state_trajectory, time_stimulus] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'stimulus_state_trajectory');
                
                trial_data.stimulus_state_trajectory = stimulus_state_trajectory;
                trial_data.time_stimulus = time_stimulus;
            end
            if strcmp(experimental_paradigm, 'GvsOverground')
                [analog_trajectories, time_analog, sampling_rate_analog, analog_labels, analog_directions] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'analog_trajectories');
                gvs_trajectory = analog_trajectories(:, strcmp(analog_labels, 'GVS_out'));
                trial_data.time_analog = time_analog;
            end
            if strcmp(experimental_paradigm, 'OculusLaneRestriction')
                gvs_trajectory = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'GVS_current_trajectory');%'GVS_current_trajectory');
                [stimulus_state_trajectory, time_stimulus] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'stimulus_state_trajectory');
                scene_translation_trajectory = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'SceneTranslation_trajectory');
                load('virtualobjectInfo');
                trial_data.scene_translation_trajectory = scene_translation_trajectory; %%%% added by SD ********
                trial_data.virtual_object_ap_location = virtual_object_ap_location; %%%% added by SD ********
                trial_data.virtual_object_ml_location = virtual_object_ml_location; %%%% added by SD ********
                trial_data.stimulus_state_trajectory = stimulus_state_trajectory;
                trial_data.time_stimulus = time_stimulus;
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
                trial_data.illusion_trajectory = illusion_trajectory;
            end
            if strcmp(experimental_paradigm, 'Vision_old')
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
                trial_data.illusion_trajectory = illusion_trajectory;
            end
            if strcmp(experimental_paradigm, 'Vision') || strcmp(experimental_paradigm, 'CadenceVision') || strcmp(experimental_paradigm, 'SR_VisualStim') || strcmp(experimental_paradigm, 'CognitiveLoadVision')
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
                trial_data.illusion_trajectory = illusion_trajectory;
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
                trial_data.illusion_trajectory = illusion_trajectory;
            end
            if strcmp(experimental_paradigm, 'GvsOverground')
                gvs_threshold = study_settings.get('gvs_threshold');
                trial_data.illusion_trajectory = zeros(size(time_analog)); % -1 = LEFT, 1 = RIGHT
                
                % positive current = anode on the right = illusory fall to the left
                trial_data.illusion_trajectory(gvs_trajectory > gvs_threshold) = -1;
                % negative current = anode on the left = illusory fall to the right
                trial_data.illusion_trajectory(gvs_trajectory < -gvs_threshold) = 1;
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
                trial_data.illusion_trajectory = illusion_trajectory; %% added by SD *********
            end
            
            %% extract events
            if strcmp(condition_stimulus, 'NONE') ...
                    || strcmp(condition_stimulus, 'VISUAL') || strcmp(experimental_paradigm, 'GVS_old') ...
                    || strcmp(experimental_paradigm, 'Vision') || strcmp(experimental_paradigm, 'CadenceVision') || strcmp(experimental_paradigm, 'SR_VisualStim') ||strcmp(experimental_paradigm, 'CognitiveLoadVision') ...
                    || strcmp(experimental_paradigm, 'GVS') || strcmp(experimental_paradigm, 'CadenceGVS') || strcmp(experimental_paradigm, 'FatigueGVS') || strcmp(experimental_paradigm, 'CognitiveLoadGvs')  ...
                    || strcmp(condition_stimulus, 'OBSTACLE') || strcmp(condition_stimulus, 'ARMSENSE') ...
                    || strcmp(experimental_paradigm, 'Stochastic Resonance') || strcmp(experimental_paradigm, 'Vision Stochastic') || strcmp(experimental_paradigm, 'GvsOverground')...
                    || strcmp(experimental_paradigm, 'OculusLaneRestriction') %% added by SD ****** 
                trial_data.right_pushoff_times = event_data{strcmp(event_labels, 'right_pushoff')};
                trial_data.right_touchdown_times = event_data{strcmp(event_labels, 'right_touchdown')};
                trial_data.left_pushoff_times = event_data{strcmp(event_labels, 'left_pushoff')};
                trial_data.left_touchdown_times = event_data{strcmp(event_labels, 'left_touchdown')};
            end
            
            %% find triggers
            %
            % Find the triggering events that indicate a stretch of interest. For perturbation experiments, this is the onset of
            % a perturbation. For unperturbed walking, this is any heelstrike.
            % The result is trigger_indices_labview.
            %
            
            % TODO: move the code below into the new function. Do this bit by bit, when able to test stuff.
            % 
            trial_data = determineTriggerTimes(study_settings, trial_data);
            
            if strcmp(condition_stimulus, 'NONE')
                % use all touchdown events as triggers
                trial_data.trigger_times = [left_touchdown_times; right_touchdown_times];
            end
            if strcmp(experimental_paradigm, 'OculusLaneRestriction') 
                % find the time steps where the stimulus state crosses a threshold
                stimulus_threshold = 4.5;
                trigger_indices_stimulus = find(diff(sign(stimulus_state_trajectory - stimulus_threshold)) > 0) + 2;
                trial_data.trigger_times = time_stimulus(trigger_indices_stimulus);
            end
            if strcmp(condition_stimulus, 'ARMSENSE')
                 trial_data.trigger_times = [];
            end
            if strcmp(experimental_paradigm, 'Vision Stochastic')
                trial_data.trigger_times = [];
            end
            
            % calculate indices
            if exist('time_marker', 'var')
                trigger_indices_mocap = zeros(size(trial_data.trigger_times));
                for i_index = 1 : length(trial_data.trigger_times)
                    [~, index_mocap] = min(abs(time_marker - trial_data.trigger_times(i_index)));
                    trigger_indices_mocap(i_index) = index_mocap;
                end
                trial_data.trigger_indices_mocap = trigger_indices_mocap;
            end
            
            if exist('time_stimulus', 'var')
                trigger_indices_stimulus = zeros(size(trial_data.trigger_times));
                for i_index = 1 : length(trial_data.trigger_times)
                    [~, index_labview] = min(abs(time_stimulus - trial_data.trigger_times(i_index)));
                    trigger_indices_stimulus(i_index) = index_labview;
                end
                trial_data.trigger_indices_stimulus = trigger_indices_stimulus;
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
            number_of_triggers = length(trial_data.trigger_times);
            
            [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels(study_settings, subject_settings, trial_data);
            
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
                removal_flags = 0;
                
                init_time = trial_data.right_pushoff_times(1) - 1; % assume heel-off happened at least one second before toes-off, so start looking at that point
                end_time = trial_data.right_touchdown_times(1);
                unload_time = trial_data.right_pushoff_times(1);
                
                % determine unload time as maximal backward-right shift of the CoP, following Halliday et al, Gait and Posture 8 (1998) 8?14
                [~, start_time_index_forceplate] = min(abs(time_forceplate - init_time));
                [~, unload_time_index_forceplate] = min(abs(time_forceplate - unload_time));
                cop_data_relevant = cop_world_trajectory(start_time_index_forceplate : unload_time_index_forceplate, :);
                time_forceplate_relevant = time_forceplate(start_time_index_forceplate : unload_time_index_forceplate);
                [~, release_time_index_forceplate] = max(cop_data_relevant(:, 1));
                release_time = time_forceplate_relevant(release_time_index_forceplate);
                start_time = release_time - 0.5;
                
                
                stretch_start_times = trial_data.right_pushoff_times(1) - 1; % HR: this is probably not right anymore
                stretch_end_times = trial_data.right_touchdown_times(1);
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
                
                trial_data.left_touchdown_times_relevant = ...
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
                    this_right_pushoff = min(right_pushoff_timesright_pushoff_times(right_pushoff_times > this_stretch_start));
                    this_right_touchdown = min(right_pushoff_timesright_touchdown_times(right_touchdown_times > this_stretch_start));
                    this_left_pushoff = min(right_pushoff_timesleft_pushoff_times(left_pushoff_times > this_stretch_start));
                    this_left_touchdown = min(right_pushoff_timesleft_touchdown_times(left_touchdown_times > this_stretch_start));
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
            number_of_stretches = size(event_variables_to_save.stretch_times, 1);
%             removal_flags = zeros(number_of_stretches, 1);
            
            if study_settings.get('prune_step_time_outliers')
                    stretch_durations = event_variables_to_save.stretch_times(:,2) - event_variables_to_save.stretch_times(:,1);
                    stretch_duration_outlier_limits = median(stretch_durations) * [.75 1.25];
                    removal_flags(stretch_durations < stretch_duration_outlier_limits(1)) = 1;
                    removal_flags(stretch_durations > stretch_duration_outlier_limits(2)) = 1;
                    if any(removal_flags)
                        disp(['Removing a stretch due to innappropriate step length']);
                    end
            end
            
            if study_settings.get('prune_gaps_com', 1)
                for i_stretch = 1 : number_of_stretches
                    [~, start_index_mocap] = min(abs(time_marker - event_variables_to_save.stretch_times(i_stretch, 1)));
                    [~, end_index_mocap] = min(abs(time_marker - event_variables_to_save.stretch_times(i_stretch, end)));
                    if any(any(isnan(com_trajectories(start_index_mocap : end_index_mocap,:))))
                        removal_flags(i_stretch) = 1;
                        disp('Removing a stretch due to gaps in com trajectories')
                    end
                end
            end
            
            if study_settings.get('prune_gaps_angles', 1)
                for i_stretch = 1 : number_of_stretches
                    [~, start_index_mocap] = min(abs(time_marker - event_variables_to_save.stretch_times(i_stretch, 1)));
                    [~, end_index_mocap] = min(abs(time_marker - event_variables_to_save.stretch_times(i_stretch, end)));
                    if any(any(isnan(joint_angle_trajectories(start_index_mocap : end_index_mocap,:))))
                        removal_flags(i_stretch) = 1;
                        disp('Removing a stretch due to gaps in joint trajectories')
                    end
                end
            end
            
            if strcmp(experimental_paradigm, 'OculusLaneRestriction')

                if study_settings.get('prune_step_placements')
                    step_zone_delinquent_list = cell(size(conditions_trial.trigger_foot_list));
                    for i_stretch = 1 : number_of_stretches
                        [~, trigger_start_index_mocap] = min(abs(time_marker - event_variables_to_save.stretch_times(i_stretch, 1)));
                        [~, trigger_end_index_mocap] = min(abs(time_marker - event_variables_to_save.stretch_times(i_stretch, 2)));
                        [~, remainder_start_index_mocap] = min(abs(time_marker - event_variables_to_save.stretch_times(i_stretch, 1)));
                        [~, remainder_end_index_mocap] = min(abs(time_marker - event_variables_to_save.stretch_times(i_stretch, 2)));
                        
                        
                        LHEE_marker_data = extractMarkerData(marker_trajectories, marker_labels, 'LHEE');
                        RHEE_marker_data = extractMarkerData(marker_trajectories, marker_labels, 'RHEE');
                        LTOE_marker_data = extractMarkerData(marker_trajectories, marker_labels, 'LTOE');
                        RTOE_marker_data = extractMarkerData(marker_trajectories, marker_labels, 'RTOE');
                        
                        
                        step_zone_delinquent_list{i_stretch} = 'NONE';
                        if strcmp(conditions_trial.zone_side_list{i_stretch}, 'STIM_ZONE_LEFT')
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
                        if strcmp(conditions_trial.zone_side_list{i_stretch}, 'STIM_ZONE_RIGHT')
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
                        [~, start_index] = min(abs(time_marker - event_variables_to_save.stretch_times(i_stretch, 1)));
                        [~, end_index] = min(abs(time_marker - event_variables_to_save.stretch_times(i_stretch, end)));
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
                        if event_variables_to_save.stretch_times(i_stretch, 1) <= ignore_times(i_ignore) && ignore_times(i_ignore) <= event_variables_to_save.stretch_times(i_stretch, end)
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
            event_variables_to_save.bands_per_stretch = size(event_variables_to_save.stretch_times, 2) - 1;
            
            stretches_file_name = ['analysis' filesep makeFileName(collection_date, subject_id, condition_list{i_condition}, i_trial, 'relevantDataStretches')];
            saveDataToFile(stretches_file_name, event_variables_to_save);
            
            
            disp(['Finding Relevant Data Stretches: condition ' condition_list{i_condition} ', Trial ' num2str(i_trial) ' completed, found ' num2str(size(event_variables_to_save.stretch_times, 1)) ' relevant stretches, saved as ' stretches_file_name]);                
        end
    end
end




















