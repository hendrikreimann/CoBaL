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
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    trials_to_process_table = makeTrialsToProcessTable(condition_list, trial_number_list);

    %% prepare
    study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');

    %% process
    
    for i_trial = 1 : height(trials_to_process_table)
        % create container for data
        trial_data = struct;
        
        % fill in basics
        trial_data.trial_type = trials_to_process_table{i_trial, 'trial type'};
        trial_data.trial_number = trials_to_process_table{i_trial, 'trial number'};

        % load data
        trial_data.loaded_events_data = load(['analysis' filesep makeFileName(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'events')]);
        trial_data = loadKinematicData(study_settings, subject_settings, trial_data);
        trial_data = loadForceplateData(study_settings, subject_settings, trial_data);
        trial_data = loadStimulusData(study_settings, subject_settings, trial_data);

        %  determine auxiliary information
        trial_data.condition_experimental = getExperimentalCondition(study_settings, subject_settings, trial_data); % the "experimental condition" is a specially highlighted condition that will be used to make some decisions lateron

        % determine illusion
        trial_data = determineIllusion(study_settings, trial_data);

        % extract events
        trial_data = extractEvents(study_settings, trial_data);

        % find triggers that indicate that a stretch of interest starts
        trial_data = determineTriggerTimes(study_settings, trial_data);

        % for each trigger, determine the condition levels
        [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels(study_settings, subject_settings, trial_data);

        % remove stretches where important variables are missing
        [conditions_trial, event_variables_to_save] = removeProblematicData(study_settings, trial_data, conditions_trial, event_variables_to_save, removal_flags);

        % remove stretches with levels listed in the settings to be ignored
        [conditions_trial, event_variables_to_save] = removeLevelsSpecifiedInSettings(study_settings, conditions_trial, event_variables_to_save);

        % save
        saveDeterminedStretchData(conditions_trial, event_variables_to_save, subject_settings, trial_data)
    end
end

function conditions_file_name = getConditionsFileName(subject_settings)
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');

    conditions_file_name = [];
    if exist('conditions.csv', 'file')
        conditions_file_name = 'conditions.csv';
    end
    if exist(makeFileName(collection_date, subject_id, 'conditions.csv'), 'file')
        conditions_file_name = makeFileName(collection_date, subject_id, 'conditions.csv');
    end    

end

function condition_experimental = getExperimentalCondition(study_settings, subject_settings, trial_data)
    condition_experimental = study_settings.get('experimental_condition');
    if strcmp(condition_experimental, 'load_from_conditions_file')
        conditions_file_name = getConditionsFileName(subject_settings);
        condition_experimental = loadConditionFromFile(conditions_file_name, 'experimental', trial_data.trial_number);
    end
    if strcmp(condition_experimental, 'determine_from_file_name')
        condition_experimental = trial_data.trial_type;
    end
    if strcmp(condition_experimental, 'determine_from_type_day_combination')
        % do nothing, we'll deal with this depending on experiment type
    end
end

function trial_data = loadKinematicData(study_settings, subject_settings, trial_data)
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
                
    [trial_data.marker_trajectories, trial_data.time_marker, ~, trial_data.marker_labels] = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'marker_trajectories');
    if study_settings.get('prune_gaps_com', 1)
        trial_data.com_trajectories = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'com_trajectories');
    end
    if study_settings.get('prune_gaps_angles', 1)
        trial_data.joint_angle_trajectories = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'joint_angle_trajectories');
    end

end

function trial_data = loadForceplateData(study_settings, subject_settings, trial_data)
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
    experimental_paradigm = study_settings.get('experimental_paradigm', 1);   

    % forceplate data
    trial_data.left_forceplate_cop_world_trajectory = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'left_foot_cop_world', 'optional');
    trial_data.right_forceplate_cop_world_trajectory = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'right_foot_cop_world', 'optional');
    [trial_data.cop_world_trajectory, trial_data.time_forceplate] = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'total_forceplate_cop_world', 'optional');
    if strcmp(experimental_paradigm, 'GvsOverground')
        first_forceplate_wrench_trajectory = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'left_foot_wrench_world');
        trial_data.vertical_force_trajectory = first_forceplate_wrench_trajectory(:, 3);
    end
end

function trial_data = loadStimulusData(study_settings, subject_settings, trial_data)
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
    experimental_paradigm = study_settings.get('experimental_paradigm', 1);   

    % stimulus data
    if strcmp(experimental_paradigm, 'GVS_old')
        [trial_data.stimulus_state_trajectory, trial_data.time_stimulus] = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'stimulus_state_trajectory');
        loaded_labview_data = load(['processed' filesep makeFileName(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'labviewData')]);
        trial_data.GVS_stim_trajectory = loaded_labview_data.GVS_out_trajectory + subject_settings.get('gvs_offset');
    end
    if strcmp(experimental_paradigm, 'Vision_old')
        % this if for TU data
        trial_data.visual_scene_ml_translation_trajectory = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'visual_scene_ml_translation__trajectory'); %take note of the double "_"
        [stimulus_state_trajectory, time_stimulus] = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'stimulus_state_trajectory');
        trial_data.stimulus_state_trajectory = stimulus_state_trajectory;
        trial_data.time_stimulus = time_stimulus;
    end
    if strcmp(experimental_paradigm, 'Vision') || strcmp(experimental_paradigm, 'CadenceVision') || strcmp(experimental_paradigm, 'SR_VisualStim') || strcmp(experimental_paradigm, 'CognitiveLoadVision')
%         current_rotation_trajectory = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'visual_rotation_angle_trajectory');
        trial_data.current_acceleration_trajectory = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'visual_rotation_acceleration_trajectory');
        [stimulus_state_trajectory, time_stimulus] = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'stimulus_state_trajectory');

        trial_data.stimulus_state_trajectory = stimulus_state_trajectory;
        trial_data.time_stimulus = time_stimulus;
    end
    if strcmp(experimental_paradigm, 'GVS') || strcmp(experimental_paradigm, 'CadenceGVS') || strcmp(experimental_paradigm, 'FatigueGVS') || strcmp(experimental_paradigm, 'CognitiveLoadGvs')
        trial_data.gvs_trajectory = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'GVS_current_trajectory');
        [stimulus_state_trajectory, time_stimulus] = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'stimulus_state_trajectory');

        trial_data.stimulus_state_trajectory = stimulus_state_trajectory;
        trial_data.time_stimulus = time_stimulus;
    end
    if strcmp(experimental_paradigm, 'GvsOverground')
        [analog_trajectories, time_analog, ~, analog_labels, ~] = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'analog_trajectories');
        trial_data.gvs_trajectory = analog_trajectories(:, strcmp(analog_labels, 'GVS_out'));
        trial_data.time_analog = time_analog;
    end
    if strcmp(experimental_paradigm, 'OculusLaneRestriction')
        trial_data.gvs_trajectory = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'GVS_current_trajectory');%'GVS_current_trajectory');
        [stimulus_state_trajectory, time_stimulus] = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'stimulus_state_trajectory');
        scene_translation_trajectory = loadData(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'SceneTranslation_trajectory');
        virtual_object_info = load('virtualobjectInfo');
        trial_data.scene_translation_trajectory = scene_translation_trajectory; %%%% added by SD ********
        trial_data.virtual_object_ap_location = virtual_object_info.virtual_object_ap_location; %%%% added by SD ********
        trial_data.virtual_object_ml_location = virtual_object_info.virtual_object_ml_location; %%%% added by SD ********
        trial_data.stimulus_state_trajectory = stimulus_state_trajectory;
        trial_data.time_stimulus = time_stimulus;
    end

end

function trial_data = determineIllusion(study_settings, trial_data)
    experimental_paradigm = study_settings.get('experimental_paradigm', 1);   
    if strcmp(experimental_paradigm, 'GVS_old')
        illusion_trajectory = zeros(size(trial_data.time_stimulus)); % 1 = RIGHT, -1 = LEFT
        % use GVS_out_trajectory as illusion
        for i_time = 1 : length(trial_data.time_stimulus)
            if trial_data.GVS_stim_trajectory(i_time) > 0
                % anode is on the right, cathode is on the left, illusory fall is towards the cathode, i.e. LEFT
                illusion_trajectory(i_time) = -1;
            end
            if trial_data.GVS_stim_trajectory(i_time) < 0
                % anode is on the left, cathode is on the right, illusory fall is towards the cathode, i.e. RIGHT
                illusion_trajectory(i_time) = 1;
            end
        end
        trial_data.illusion_trajectory = illusion_trajectory;
    end
    if strcmp(experimental_paradigm, 'Vision_old')
        illusion_trajectory = zeros(size(trial_data.time_stimulus)); % 1 = RIGHT, -1 = LEFT
        for i_time = 1 : length(trial_data.time_stimulus)
            if trial_data.visual_scene_ml_translation_trajectory(i_time) > 0
                % angle change is positive horizon rotates counter-clockwise, illusion is to the RIGHT
                illusion_trajectory(i_time) = 1;
            end
            if trial_data.visual_scene_ml_translation_trajectory(i_time) < 0 
                % angle change is negative, horizon rotates clockwise, illusion is to the LEFT
                illusion_trajectory(i_time) = -1;
            end
        end
        trial_data.illusion_trajectory = illusion_trajectory;
    end
    if strcmp(experimental_paradigm, 'Vision') || strcmp(experimental_paradigm, 'CadenceVision') || strcmp(experimental_paradigm, 'SR_VisualStim') || strcmp(experimental_paradigm, 'CognitiveLoadVision')
        illusion_trajectory = zeros(size(trial_data.time_stimulus)); % -1 = LEFT, 1 = RIGHT
        for i_time = 1 : length(trial_data.time_stimulus)
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
        illusion_trajectory = zeros(size(trial_data.time_stimulus)); % -1 = LEFT, 1 = RIGHT
        for i_time = 1 : length(trial_data.time_stimulus)
            if stimulus_state_trajectory(i_time) == 3
                % stimulus is currently active
                if trial_data.gvs_trajectory(i_time) < 0
                    % negative current = anode on the left = illusory fall to the right
                    illusion_trajectory(i_time) = 1;
                end
                if trial_data.gvs_trajectory(i_time) > 0
                    % positive current = anode on the right = illusory fall to the left
                    illusion_trajectory(i_time) = -1;
                end
            end
        end
        trial_data.illusion_trajectory = illusion_trajectory;
    end
    if strcmp(experimental_paradigm, 'GvsOverground')
        gvs_threshold = study_settings.get('gvs_threshold');
        trial_data.illusion_trajectory = zeros(size(trial_data.time_analog)); % -1 = LEFT, 1 = RIGHT

        % positive current = anode on the right = illusory fall to the left
        trial_data.illusion_trajectory(trial_data.gvs_trajectory > gvs_threshold) = -1;
        % negative current = anode on the left = illusory fall to the right
        trial_data.illusion_trajectory(trial_data.gvs_trajectory < -gvs_threshold) = 1;
    end
    if strcmp(experimental_paradigm, 'OculusLaneRestriction')
        illusion_trajectory = zeros(size(trial_data.time_stimulus)); % -1 = LEFT, 1 = RIGHT
        for i_time = 1 : length(trial_data.time_stimulus)
            if trial_data.stimulus_state_trajectory(i_time) == 5
                % stimulus is currently active
                if trial_data.gvs_trajectory(i_time) < 0
                    % negative current = anode on the left = illusory fall to the right
                    illusion_trajectory(i_time) = 1;
                end
                if trial_data.gvs_trajectory(i_time) > 0
                    % positive current = anode on the right = illusory fall to the left
                    illusion_trajectory(i_time) = -1;
                end
            end
        end
        trial_data.illusion_trajectory = illusion_trajectory; %% added by SD *********
    end
end

function trial_data = extractEvents(study_settings, trial_data)
    experimental_paradigm = study_settings.get('experimental_paradigm', 1);   
    condition_stimulus = study_settings.get('condition_stimulus', 1);   
    if strcmp(condition_stimulus, 'NONE') ...
            || strcmp(condition_stimulus, 'VISUAL') || strcmp(experimental_paradigm, 'GVS_old') || strcmp(experimental_paradigm, 'Vision_old') ...
            || strcmp(experimental_paradigm, 'Vision') || strcmp(experimental_paradigm, 'CadenceVision') || strcmp(experimental_paradigm, 'SR_VisualStim') ||strcmp(experimental_paradigm, 'CognitiveLoadVision') ...
            || strcmp(experimental_paradigm, 'GVS') || strcmp(experimental_paradigm, 'CadenceGVS') || strcmp(experimental_paradigm, 'FatigueGVS') || strcmp(experimental_paradigm, 'CognitiveLoadGvs')  ...
            || strcmp(experimental_paradigm, 'OBSTACLE') || strcmp(condition_stimulus, 'ARMSENSE') ...
            || strcmp(experimental_paradigm, 'GaitInitiationObstacle') || strcmp(experimental_paradigm, 'Vision Stochastic') || strcmp(experimental_paradigm, 'GvsOverground')...
            || strcmp(experimental_paradigm, 'OculusLaneRestriction') %% added by SD ****** 
        trial_data.right_pushoff_times = trial_data.loaded_events_data.event_data{strcmp(trial_data.loaded_events_data.event_labels, 'right_pushoff')};
        trial_data.right_touchdown_times = trial_data.loaded_events_data.event_data{strcmp(trial_data.loaded_events_data.event_labels, 'right_touchdown')};
        trial_data.left_pushoff_times = trial_data.loaded_events_data.event_data{strcmp(trial_data.loaded_events_data.event_labels, 'left_pushoff')};
        trial_data.left_touchdown_times = trial_data.loaded_events_data.event_data{strcmp(trial_data.loaded_events_data.event_labels, 'left_touchdown')};
    end
end

function [conditions_trial, event_variables_to_save] ...
    = removeProblematicData ...
      ( ...
        study_settings, ...
        trial_data, ...
        conditions_trial, ...
        event_variables_to_save, ...
        removal_flags ...
      )

    experimental_paradigm = study_settings.get('experimental_paradigm', 1);   
    number_of_stretches = size(event_variables_to_save.stretch_times, 1);

    if study_settings.get('prune_step_time_outliers')
        stretch_durations = event_variables_to_save.stretch_times(:,2) - event_variables_to_save.stretch_times(:,1);
        stretch_duration_outlier_limits = median(stretch_durations) * [.75 1.25];
        removal_flags(stretch_durations < stretch_duration_outlier_limits(1)) = 1;
        removal_flags(stretch_durations > stretch_duration_outlier_limits(2)) = 1;
        if any(removal_flags)
            disp('Removing a stretch due to innappropriate step length');
        end
    end

    if study_settings.get('prune_gaps_com', 1)
        for i_stretch = 1 : number_of_stretches
            [~, start_index_mocap] = min(abs(trial_data.time_marker - event_variables_to_save.stretch_times(i_stretch, 1)));
            [~, end_index_mocap] = min(abs(trial_data.time_marker - event_variables_to_save.stretch_times(i_stretch, end)));
            if any(any(isnan(trial_data.com_trajectories(start_index_mocap : end_index_mocap,:))))
                removal_flags(i_stretch) = 1;
                disp('Removing a stretch due to gaps in com trajectories')
            end
        end
    end

    if study_settings.get('prune_gaps_angles', 1)
        for i_stretch = 1 : number_of_stretches
            [~, start_index_mocap] = min(abs(trial_data.time_marker - event_variables_to_save.stretch_times(i_stretch, 1)));
            [~, end_index_mocap] = min(abs(trial_data.time_marker - event_variables_to_save.stretch_times(i_stretch, end)));
            if any(any(isnan(trial_data.joint_angle_trajectories(start_index_mocap : end_index_mocap,:))))
                removal_flags(i_stretch) = 1;
                disp('Removing a stretch due to gaps in joint trajectories')
            end
        end
    end

    if strcmp(experimental_paradigm, 'OculusLaneRestriction')

        if study_settings.get('prune_step_placements')
            step_zone_delinquent_list = cell(size(conditions_trial.trigger_foot_list));
            for i_stretch = 1 : number_of_stretches
                [~, trigger_start_index_mocap] = min(abs(trial_data.time_marker - event_variables_to_save.stretch_times(i_stretch, 1)));
                [~, trigger_end_index_mocap] = min(abs(trial_data.time_marker - event_variables_to_save.stretch_times(i_stretch, 2)));
                [~, remainder_start_index_mocap] = min(abs(trial_data.time_marker - event_variables_to_save.stretch_times(i_stretch, 1)));
                [~, remainder_end_index_mocap] = min(abs(trial_data.time_marker - event_variables_to_save.stretch_times(i_stretch, 2)));


                LHEE_marker_data = extractMarkerData(trial_data.marker_trajectories, trial_data.marker_labels, 'LHEE');
                RHEE_marker_data = extractMarkerData(trial_data.marker_trajectories, trial_data.marker_labels, 'RHEE');
                LTOE_marker_data = extractMarkerData(trial_data.marker_trajectories, trial_data.marker_labels, 'LTOE');
                RTOE_marker_data = extractMarkerData(trial_data.marker_trajectories, trial_data.marker_labels, 'RTOE');


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
        marker_weights = ones(size(trial_data.marker_labels));
    end
    for i_marker = 1 : 3 : length(trial_data.marker_labels)
        this_marker_weight = 1; % 1 is default

        this_marker_label = trial_data.marker_labels{i_marker}(1:end-2);
        if any(strcmp(marker_weights(:, 1), this_marker_label))
            this_marker_weight = marker_weights{strcmp(marker_weights(:, 1), this_marker_label), 2};
        end

        if this_marker_weight == 1
            this_marker_data = extractMarkerData(trial_data.marker_trajectories, trial_data.marker_labels, this_marker_label);
            for i_stretch = 1 : number_of_stretches
                [~, start_index] = min(abs(trial_data.time_marker - event_variables_to_save.stretch_times(i_stretch, 1)));
                [~, end_index] = min(abs(trial_data.time_marker - event_variables_to_save.stretch_times(i_stretch, end)));
                if any(any(isnan(this_marker_data(start_index : end_index, :))))
                    removal_flags(i_stretch) = 1;
                    disp(['Removing a stretch due to gap in marker "', this_marker_label '"']);
                end 
            end

        end

    end

    % check ignore markers
    for i_stretch = 1 : number_of_stretches
        if isfield(trial_data.loaded_events_data, 'ignore_times')
            for i_ignore = 1 : length(trial_data.loaded_events_data.ignore_times)
                if event_variables_to_save.stretch_times(i_stretch, 1) <= trial_data.loaded_events_data.ignore_times(i_ignore) && trial_data.loaded_events_data.ignore_times(i_ignore) <= event_variables_to_save.stretch_times(i_stretch, end)
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

end

function [conditions_trial, event_variables_to_save] ...
    = removeLevelsSpecifiedInSettings ...
      ( ...
        study_settings, ...
        conditions_trial, ...
        event_variables_to_save ...
      )
    % get info
    levels_to_remove = study_settings.get('levels_to_remove', 1);

    % go through list of levels to remove
    for i_entry = 1 : size(levels_to_remove, 1)
        % get info
        this_factor = levels_to_remove(i_entry, 1);
        this_level = levels_to_remove(i_entry, 2);
        conditions_table = study_settings.get('conditions');
        
        % find and flag conditions in the list matching this
        this_factor_list_label = conditions_table{strcmp(conditions_table(:, 1), this_factor), 2};
        this_factor_list = conditions_trial.(this_factor_list_label);
        this_factor_match = strcmp(this_factor_list, this_level);
        
        % remove flagged stretches from data
        unflagged_indices = ~this_factor_match;
        event_variables_to_save_names = fieldnames(event_variables_to_save);
        for i_variable = 1 : length(event_variables_to_save_names)
            this_variable_name = event_variables_to_save_names{i_variable};

            evalstring = ['this_variable_data = event_variables_to_save.' this_variable_name ';'];
            eval(evalstring);
            this_variable_data = this_variable_data(unflagged_indices, :);

            evalstring = ['event_variables_to_save.' this_variable_name ' = this_variable_data;'];
            eval(evalstring);
        end
        
        % remove flagged stretches from conditions
        conditions_trial_names = fieldnames(conditions_trial);
        for i_label = 1 : length(conditions_trial_names)
            this_condition_name = conditions_trial_names{i_label};

            evalstring = ['this_condition_data = conditions_trial.' this_condition_name ';'];
            eval(evalstring);
            this_condition_data = this_condition_data(unflagged_indices, :);

            evalstring = ['conditions_trial.' this_condition_name ' = this_condition_data;'];
            eval(evalstring);
        end
    end





end

function saveDeterminedStretchData(conditions_trial, event_variables_to_save, subject_settings, trial_data)
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');

    event_variables_to_save.conditions_trial = conditions_trial;
    event_variables_to_save.bands_per_stretch = size(event_variables_to_save.stretch_times, 2) - 1;

    stretches_file_name = ['analysis' filesep makeFileName(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'relevantDataStretches')];
    saveDataToFile(stretches_file_name, event_variables_to_save);


    disp(['Finding Relevant Data Stretches: condition ' char(trial_data.trial_type) ', Trial ' num2str(trial_data.trial_number) ' completed, found ' num2str(size(event_variables_to_save.stretch_times, 1)) ' relevant stretches, saved as ' stretches_file_name]);                
end

