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

% This script finds the stretches of data relevant to the experimental paradigm.

% We determine the start and end times of the relevant stretches and store them as a list, accompanied by a list of
% condition indicators. There are eight condition types. Each condition type has a list variable that is of the 
% same length as the list of stretch start/end times and gives the condition label for this stretch.

% The condition types are:
% stance_foot:      this is the foot that is in stance throughout the whole stretch, will be STANCE_LEFT, STANCE_RIGHT or STANCE_BOTH
% perturbation:     direction of the perturbation in stimulus paradigms, will be ILLUSION_LEFT, ILLUSION_RIGHT or CONTROL (should be ILLUSION_NONE instead of CONTROL)
% delay:            delay between trigger and stimulus in stimulus paradigms, options are specified in studySettings.txt
% index:            number of steps after the stimulus in stimulus paradigms, will be STEP_ONE, STEP_TWO, STEP_THREE, STEP_FOUR or CONTROL
% experimental:     condition affecting a whole trial instead of stretches of time within a trial. This will be loaded from
%                   conditions.mat if available, otherwise it will be set to the trial_type as determined from the file name.
% stimulus:         the type of stimulus used to generate the perturbation, e.g. STIM_GVS, STIM_VISUAL, STIM_NONE
% day:              the day of the data collection, e.g. DAY_PRE, DAY_POST_0, DAY_POST_24H, DAY_POST_1W

% input: 
% - subjectInfo.mat
% - optional: conditions.mat
%
% output:
% file relevantDataStretches.mat, containing
% - condition_stance_foot_list
% - condition_perturbation_list
% - condition_delay_list
% - condition_index_list
% - condition_experimental_list
% - condition_stimulus_list
% - condition_day_list
% - stretch_pushoff_times
% - stretch_start_times
% - stretch_end_times
                    
function findRelevantDataStretches(varargin)
    % parse arguments
    parser = inputParser;
    parser.KeepUnmatched = true;
%     default_stimulus_type = 'none';
%     valid_stimulus_types = {'none', 'GVS', 'VISUAL', 'obstacle'};
%     check_stimulus_type = @(x) any(validatestring(x,valid_stimulus_types));
%     addParameter(parser, 'stimulus', default_stimulus_type, check_stimulus_type)
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
%     stimulus_type = parser.Results.stimulus;
    visualize = parser.Results.visualize;
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});

    
    % stretches_to_analyze = 'single stance';
    stretches_to_analyze = 'full stride';

    %% prepare
    load('subjectInfo.mat', 'date', 'subject_id');
    % load settings
    study_settings_file = '';
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = loadSettingsFile(study_settings_file);
    
    conditions_file_name = [];
    if exist('conditions.csv', 'file')
        conditions_file_name = 'conditions.csv';
    end
    if exist(makeFileName(date, subject_id, 'conditions.csv'), 'file')
        conditions_file_name = makeFileName(date, subject_id, 'conditions.csv');
    end

    subject_settings = loadSettingsFile('subjectSettings.txt');

    time_to_nearest_heelstrike_before_trigger_threshold = 0.10; % a heelstrike should happen less than this long after a trigger
    time_to_nearest_heelstrike_after_trigger_threshold = 0.3; % a heelstrike should happen less than this long after a trigger
    
    
    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            %% load data
            load(['analysis' filesep makeFileName(date, subject_id, condition_list{i_condition}, i_trial, 'stepEvents')]);
            
            % determine experimental condition
            condition_experimental = study_settings.experimental_condition;
            if strcmp(condition_experimental, 'load_from_conditions_file')
                condition_experimental = loadConditionFromFile(conditions_file_name, 'experimental', i_trial);
            end
            if strcmp(condition_experimental, 'determine_from_file_name')
                condition_experimental = condition_list{i_condition};
            end
            
            % determine stimulus type
            condition_stimulus = study_settings.stimulus_condition;
            if strcmp(condition_stimulus, 'load_from_conditions_file')
                condition_stimulus = loadConditionFromFile(conditions_file_name, 'stimulus', i_trial);
            end
            if strcmp(condition_stimulus, 'determine_from_file_name')
                condition_stimulus = condition_list{i_condition};
            end
            
            % determine day
            condition_day = study_settings.day_condition;
            if strcmp(condition_day, 'load_from_conditions_file')
                condition_day = loadConditionFromFile(conditions_file_name, 'day', i_trial);
            end
            if strcmp(condition_day, 'determine_from_file_name')
                condition_day = condition_list{i_condition};
            end
            
            
            
            % marker data
            [marker_trajectories, time_marker, sampling_rate_marker, marker_labels] = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'marker_trajectories');
            
            % forceplate data
            [left_forceplate_cop_world_trajectory, time_left_forceplate, ~, ~, left_forceplate_available] = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'left_forceplate_cop_world', 'optional');
            [right_forceplate_cop_world_trajectory, time_right_forceplate, ~, ~, right_forceplate_available] = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'right_forceplate_cop_world', 'optional');
            if left_forceplate_available && right_forceplate_available
                left_copx_trajectory = left_forceplate_cop_world_trajectory(:, 1);
                right_copx_trajectory = right_forceplate_cop_world_trajectory(:, 1);
                
                if any(strcmp(condition_list{i_condition}, study_settings.trial_types_with_inverted_forceplate_sides))
                    % this was for the VEPO lab, where the sides are inverted
                    left_copx_trajectory = right_forceplate_cop_world_trajectory(:, 3);
                    right_copx_trajectory = left_forceplate_cop_world_trajectory(:, 3);
                end
            end
            
            % stimulus data
            if strcmp(condition_stimulus, 'GVS')
                GVS_out_trajectory = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'GVS_out_trajectory');
                [stimulus_state_trajectory, time_stimulus] = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'stimulus_state_trajectory');
            end
            if strcmp(condition_stimulus, 'VISUAL')
                visual_scene_ml_translation__trajectory = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'visual_scene_ml_translation__trajectory');
                [stimulus_state_trajectory, time_stimulus] = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'stimulus_state_trajectory');
            end
            
            % determine indices for optional markers
            optional_marker_indices = [];
            for i_marker = 1 : length(study_settings.optional_markers)
                marker = find(strcmp(marker_labels, study_settings.optional_markers{i_marker}));
                marker_indices = reshape([(marker - 1) * 3 + 1; (marker - 1) * 3 + 2; (marker - 1) * 3 + 3], 1, length(marker)*3);
                optional_marker_indices = [optional_marker_indices marker_indices];
            end
            essential_marker_indicator = ~ismember(1 : size(marker_trajectories, 2), optional_marker_indices);
            
            % determine illusion
            if strcmp(condition_stimulus, 'GVS')
                illusion_trajectory = zeros(size(time_stimulus)); % 1 = RIGHT, -1 = LEFT
                % use GVS_out_trajectory as illusion
                for i_time = 1 : length(time_stimulus)
                    if GVS_out_trajectory(i_time) > 0
                        % anode is on the right, cathode is on the left, illusory fall is towards the cathode, i.e. LEFT
                        illusion_trajectory(i_time) = -1;
                    end
                    if GVS_out_trajectory(i_time) < 0
                        % anode is on the left, cathode is on the right, illusory fall is towards the cathode, i.e. RIGHT
                        illusion_trajectory(i_time) = 1;
                    end
                end
            end
            if strcmp(condition_stimulus, 'VISUAL')
                illusion_trajectory = zeros(size(time_stimulus)); % 1 = RIGHT, -1 = LEFT
                for i_time = 1 : length(time_stimulus)
                    if visual_scene_ml_translation__trajectory(i_time) > 0
                        % angle change is positive, horizon rotates counter-clockwise, illusion is to the RIGHT
                        illusion_trajectory(i_time) = 1;
                    end
                    if visual_scene_ml_translation__trajectory(i_time) < 0
                        % angle change is negative, horizon rotates clockwise, illusion is to the LEFT
                        illusion_trajectory(i_time) = -1;
                    end
                end
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
            if strcmp(condition_stimulus, 'VISUAL') || strcmp(condition_stimulus, 'GVS')
                % find the time steps where the stimulus state crosses a threshold
                stimulus_threshold = 0.5;
                trigger_indices_labview = find(diff(sign(stimulus_state_trajectory - stimulus_threshold)) > 0) + 1;

                %
                epsilon = 1e-5;
                stim_start_indices_labview = find(diff(sign(abs(illusion_trajectory) - epsilon)) > 0) + 1;
                trigger_indices_labview = trigger_indices_labview(1 : length(stim_start_indices_labview)); % in case a stim is triggered, but not recorded
                
                trigger_times = time_stimulus(trigger_indices_labview);
            end
            if strcmp(condition_stimulus, 'OBSTACLE')
                trigger_times = [];
            end
            
            % calculate indices
            trigger_indices_mocap = zeros(size(trigger_times));
            for i_index = 1 : length(trigger_times)
                [~, index_mocap] = min(abs(time_marker - trigger_times(i_index)));
                trigger_indices_mocap(i_index) = index_mocap;
            end

            % visualize triggers
            if visualize

                figure; axes; hold on
%                 plot(time_stimulus, stimulus_state_trajectory*0.02);
                if left_forceplate_available
                    left_cop_x_trajectory_relevant = left_copx_trajectory; left_cop_x_trajectory_relevant(left_cop_x_trajectory_relevant==0) = NaN;
                    plot(time_left_forceplate, left_cop_x_trajectory_relevant, 'linewidth', 2, 'Displayname', 'cop left');
                end
                
                if right_forceplate_available
                    right_cop_x_trajectory_relevant = right_copx_trajectory; right_cop_x_trajectory_relevant(right_cop_x_trajectory_relevant==0) = NaN;
                    plot(time_right_forceplate, right_cop_x_trajectory_relevant, 'linewidth', 2, 'Displayname', 'cop right');
                end
                plot(time_marker(trigger_indices_mocap), zeros(size(trigger_indices_mocap)), 'x', 'Displayname', 'triggers')
%                 plot(time_stimulus, illusion_trajectory, 'Displayname', 'illusion')
%                 legend('stimulus state', 'left cop', 'right cop', 'left touchdown', 'right touchdown', 'trigger', 'stim start')
                legend('toggle')
            end

            %% extract event data
            %
            % For each trigger, determine the conditions and the relevant step events.
            % The result is 
            % stretch_start_indices_forceplate, stretch_end_indices_forceplate
            % stretch_start_times, stretch_end_times
            % condition_stance_foot_list, condition_perturbation_list, condition_delay_list, condition_index_list

            number_of_triggers = length(trigger_times);
            removal_flags = zeros(number_of_triggers, 1);
            event_variables_to_save = struct;
            if strcmp(condition_stimulus, 'NONE')
                stretch_start_times = zeros(number_of_triggers, 1);
                stretch_end_times = zeros(number_of_triggers, 1);
                stretch_pushoff_times = zeros(number_of_triggers, 1);
                closest_heelstrike_distance_times = zeros(number_of_triggers, 1);
                condition_stance_foot_list_trial = cell(number_of_triggers, 1);
                condition_perturbation_list_trial = cell(number_of_triggers, 1);
                condition_delay_list_trial = cell(number_of_triggers, 1);
                condition_index_list_trial = cell(number_of_triggers, 1);
                condition_experimental_list_trial = cell(number_of_triggers, 1);
                condition_stimulus_list_trial = cell(number_of_triggers, 1);
                condition_day_list_trial = cell(number_of_triggers, 1);
                
                for i_trigger = 1 : number_of_triggers
                    condition_perturbation_list_trial{i_trigger, 1} = 'N/A';
                    condition_delay_list_trial{i_trigger, 1} = 'N/A';
                    condition_index_list_trial{i_trigger, 1} = 'N/A';
                    condition_experimental_list_trial{i_trigger, 1} = condition_experimental;
                    condition_stimulus_list_trial{i_trigger, 1} = condition_stimulus;
                    condition_day_list_trial{i_trigger, 1} = condition_day;
                    
                
                    % find out which heelstrike triggered
                    % XXX change this to use the interval, same as in the stimulus case
                    [distance_to_trigger_left_time, index_left] = min(abs(left_touchdown_times - trigger_times(i_trigger)));
                    [distance_to_trigger_right_time, index_right] = min(abs(right_touchdown_times - trigger_times(i_trigger)));
                    closest_heelstrike_distance_time = min([distance_to_trigger_left_time distance_to_trigger_right_time]);
                    
                    if distance_to_trigger_left_time < distance_to_trigger_right_time
                        condition_stance_foot_list_trial{i_trigger, 1} = 'STANCE_LEFT';
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
                        condition_stance_foot_list_trial{i_trigger, 1} = 'STANCE_RIGHT';
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
                condition_stance_foot_list_trial = condition_stance_foot_list_trial(unflagged_indices, :);
                condition_perturbation_list_trial = condition_perturbation_list_trial(unflagged_indices, :);
                condition_delay_list_trial = condition_delay_list_trial(unflagged_indices, :);
                condition_index_list_trial = condition_index_list_trial(unflagged_indices, :);
                condition_experimental_list_trial = condition_experimental_list_trial(unflagged_indices, :);
                condition_stimulus_list_trial = condition_stimulus_list_trial(unflagged_indices, :);
                condition_day_list_trial = condition_day_list_trial(unflagged_indices, :);
                closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :); 
                
                if visualize
                    for i_trigger = 1 : length(stretch_start_times)
                        if strcmp(condition_stance_foot_list_trial(i_trigger), 'STANCE_RIGHT')
                            stretch_indicator_height = 0.01;
                        else
                            stretch_indicator_height = -0.01;
                        end
                            
                        plot([stretch_start_times(i_trigger) stretch_end_times(i_trigger)], [1 1]*stretch_indicator_height, 'linewidth', 3);
                    end
                end
                
%                 % check step times and flag outliers
%                 number_of_stretches = length(stretch_start_times);
%                 stretch_times = stretch_end_times - stretch_start_times;
%                 stretch_time_outlier_limits = median(stretch_times) * [0.8 1.2];
%                 
%                 removal_flags = zeros(number_of_stretches, 1);
%                 removal_flags(stretch_times < stretch_time_outlier_limits(1)) = 1;
%                 removal_flags(stretch_times > stretch_time_outlier_limits(2)) = 1;
%                 
%                 % check data availability and flag stretches with gaps
%                 for i_stretch = 1 : number_of_stretches
%                     [~, start_index_mocap] = min(abs(time_marker - stretch_start_times(i_stretch)));
%                     [~, end_index_mocap] = min(abs(time_marker - stretch_end_times(i_stretch)));
%                     if any(any(isnan(marker_trajectories(start_index_mocap : end_index_mocap, essential_marker_indicator))))
%                         removal_flags(i_stretch) = 1;
%                     end
%                 end
%                 
%                 
%                 % remove flagged triggers
%                 unflagged_indices = ~removal_flags;
%                 trigger_times = trigger_times(unflagged_indices);
%                 stretch_start_times = stretch_start_times(unflagged_indices, :);
%                 stretch_end_times = stretch_end_times(unflagged_indices, :);
%                 stretch_pushoff_times = stretch_pushoff_times(unflagged_indices, :);
%                 condition_stance_foot_list_trial = condition_stance_foot_list_trial(unflagged_indices, :);
%                 condition_perturbation_list_trial = condition_perturbation_list_trial(unflagged_indices, :);
%                 condition_delay_list_trial = condition_delay_list_trial(unflagged_indices, :);
%                 condition_index_list_trial = condition_index_list_trial(unflagged_indices, :);
%                 condition_experimental_list_trial = condition_experimental_list_trial(unflagged_indices, :);
%                 condition_stimulus_list_trial = condition_stimulus_list_trial(unflagged_indices, :);
%                 condition_day_list_trial = condition_day_list_trial(unflagged_indices, :);
%                 closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :);
                
                % add new variables to be saved
                event_variables_to_save.condition_stance_foot_list_trial = condition_stance_foot_list_trial;
                event_variables_to_save.condition_perturbation_list_trial = condition_perturbation_list_trial;
                event_variables_to_save.condition_delay_list_trial = condition_delay_list_trial;
                event_variables_to_save.condition_index_list_trial = condition_index_list_trial;
                event_variables_to_save.condition_experimental_list_trial = condition_experimental_list_trial;
                event_variables_to_save.condition_stimulus_list_trial = condition_stimulus_list_trial;
                event_variables_to_save.condition_day_list_trial = condition_day_list_trial;
                event_variables_to_save.stretch_start_times = stretch_start_times;
                event_variables_to_save.stretch_pushoff_times = stretch_pushoff_times;
                event_variables_to_save.stretch_end_times = stretch_end_times;
                
                
%                 % save data
%                 data_stretches_file_name = ['analysis' filesep makeFileName(date, subject_id, condition_list{i_condition}, i_trial, 'relevantDataStretches')];
%                 save ...
%                   ( ...
%                     data_stretches_file_name, ...
%                     'condition_stance_foot_list_trial', ...
%                     'condition_perturbation_list_trial', ...
%                     'condition_delay_list_trial', ...
%                     'condition_index_list_trial', ...
%                     'condition_experimental_list_trial', ...
%                     'condition_stimulus_list_trial', ...
%                     'condition_day_list_trial', ...
%                     'stretch_start_times', ...
%                     'stretch_pushoff_times', ...
%                     'stretch_end_times' ...
%                   )
% 
%                 disp(['Finding Relevant Data Stretches: condition ' condition_list{i_condition} ', Trial ' num2str(i_trial) ' completed, found ' num2str(length(stretch_start_times)) ' relevant stretches, saved as ' data_stretches_file_name]);                
               
                
                
                
            end  
            if strcmp(condition_stimulus, 'OBSTACLE')
                stretch_start_times = zeros(2, 1);
                stretch_end_times = zeros(2, 1);
                stretch_pushoff_times = zeros(2, 1);
                
                condition_stance_foot_list_trial = cell(2, 1);
                condition_perturbation_list_trial = cell(2, 1);
                condition_delay_list_trial = cell(2, 1);
                condition_index_list_trial = cell(2, 1);
                condition_experimental_list_trial = cell(2, 1);
                condition_stimulus_list_trial = cell(2, 1);
                condition_day_list_trial = cell(2, 1);
                
                % second stretch goes from first right pushoff to first right heelstrike
                stretch_start_times(2, 1) = right_pushoff_times(1);
                stretch_end_times(2, 1) = right_touchdown_times(1);
                condition_perturbation_list_trial{2, 1} = 'N/A';
                condition_delay_list_trial{2, 1} = 'N/A';
                condition_index_list_trial{2, 1} = 'TWO';
                condition_experimental_list_trial{2, 1} = condition_experimental;
                condition_stance_foot_list_trial{2, 1} = 'STANCE_LEFT';
                condition_stimulus_list_trial{2, 1} = condition_stimulus;
                condition_day_list_trial{2, 1} = condition_day;
                
%                 % third stretch goes from first right heelstrike to first left heelstrike
%                 stretch_start_times(3, 1) = right_touchdown_times(1);
%                 stretch_end_times(3, 1) = left_touchdown_times(1);
%                 condition_perturbation_list_trial{3, 1} = 'N/A';
%                 condition_delay_list_trial{3, 1} = 'N/A';
%                 condition_index_list_trial{3, 1} = 'THREE';
%                 condition_experimental_list_trial{3, 1} = condition_experimental;
%                 condition_stance_foot_list_trial{3, 1} = 'STANCE_RIGHT';
%                 condition_stimulus_list_trial{3, 1} = condition_stimulus;
%                 condition_day_list_trial{3, 1} = condition_day;
                
                % first stretch goes from first right pushoff minus 1 second to first right pushoff
                stretch_start_times(1, 1) = right_pushoff_times(1) - 1;
                stretch_end_times(1, 1) = right_pushoff_times(1);
                condition_perturbation_list_trial{1, 1} = 'N/A';
                condition_delay_list_trial{1, 1} = 'N/A';
                condition_index_list_trial{1, 1} = 'ONE';
                condition_experimental_list_trial{1, 1} = condition_experimental;
                condition_stance_foot_list_trial{1, 1} = 'STANCE_BOTH';
                condition_stimulus_list_trial{1, 1} = condition_stimulus;
                condition_day_list_trial{1, 1} = condition_day;
                
%                 % remove flagged triggers
%                 unflagged_indices = ~removal_flags;
%                 trigger_times = trigger_times(unflagged_indices);
%                 stretch_start_times = stretch_start_times(unflagged_indices, :);
%                 stretch_end_times = stretch_end_times(unflagged_indices, :);
%                 stretch_pushoff_times = stretch_pushoff_times(unflagged_indices, :);
%                 condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
%                 condition_perturbation_list = condition_perturbation_list(unflagged_indices, :);
%                 condition_delay_list = condition_delay_list(unflagged_indices, :);
%                 condition_index_list = condition_index_list(unflagged_indices, :);
%                 condition_experimental_list = condition_experimental_list(unflagged_indices, :);
%                 closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :); 
                
                if visualize
                    for i_trigger = 1 : length(stretch_start_times)
                        if strcmp(condition_stance_foot_list_trial(i_trigger), 'STANCE_LEFT')
                            stretch_indicator_height = 0.01;
                        else
                            stretch_indicator_height = -0.01;
                        end
                            
                        plot([stretch_start_times(i_trigger) stretch_end_times(i_trigger)], [1 1]*stretch_indicator_height, 'linewidth', 3);
                    end
                end
                
%                 % check step times and flag outliers
%                 number_of_stretches = length(stretch_start_times);
%                 stretch_times = stretch_end_times - stretch_start_times;
%                 stretch_time_outlier_limits = median(stretch_times) * [0.8 1.2];
%                 
%                 removal_flags = zeros(number_of_stretches, 1);
%                 removal_flags(stretch_times < stretch_time_outlier_limits(1)) = 1;
%                 removal_flags(stretch_times > stretch_time_outlier_limits(2)) = 1;
                
%                 % check data availability and flag stretches with gaps
%                 for i_stretch = 1 : number_of_stretches
%                     [~, start_index_mocap] = min(abs(time_marker - stretch_start_times(i_stretch)));
%                     [~, end_index_mocap] = min(abs(time_marker - stretch_end_times(i_stretch)));
%                     if any(any(isnan(marker_trajectories(start_index_mocap : end_index_mocap, essential_marker_indicator))))
%                         removal_flags(i_stretch) = 1;
%                     end
%                 end
                
                
%                 % remove flagged triggers
%                 unflagged_indices = ~removal_flags;
%                 trigger_times = trigger_times(unflagged_indices);
%                 stretch_start_times = stretch_start_times(unflagged_indices, :);
%                 stretch_end_times = stretch_end_times(unflagged_indices, :);
%                 stretch_pushoff_times = stretch_pushoff_times(unflagged_indices, :);
%                 condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
%                 condition_perturbation_list = condition_perturbation_list(unflagged_indices, :);
%                 condition_delay_list = condition_delay_list(unflagged_indices, :);
%                 condition_index_list = condition_index_list(unflagged_indices, :);
%                 condition_experimental_list = condition_experimental_list(unflagged_indices, :);
%                 closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :);
                

                % add new variables to be saved
                event_variables_to_save.condition_stance_foot_list_trial = condition_stance_foot_list_trial;
                event_variables_to_save.condition_perturbation_list_trial = condition_perturbation_list_trial;
                event_variables_to_save.condition_delay_list_trial = condition_delay_list_trial;
                event_variables_to_save.condition_index_list_trial = condition_index_list_trial;
                event_variables_to_save.condition_experimental_list_trial = condition_experimental_list_trial;
                event_variables_to_save.condition_stimulus_list_trial = condition_stimulus_list_trial;
                event_variables_to_save.condition_day_list_trial = condition_day_list_trial;
                event_variables_to_save.stretch_start_times = stretch_start_times;
                event_variables_to_save.stretch_pushoff_times = stretch_pushoff_times;
                event_variables_to_save.stretch_end_times = stretch_end_times;

%                 % save data
%                 data_stretches_file_name = ['analysis' filesep makeFileName(date, subject_id, condition_list{i_condition}, i_trial, 'relevantDataStretches')];
%                 save ...
%                   ( ...
%                     data_stretches_file_name, ...
%                     'condition_stance_foot_list_trial', ...
%                     'condition_perturbation_list_trial', ...
%                     'condition_delay_list_trial', ...
%                     'condition_index_list_trial', ...
%                     'condition_experimental_list_trial', ...
%                     'condition_stimulus_list_trial', ...
%                     'condition_day_list_trial', ...
%                     'stretch_pushoff_times', ...
%                     'stretch_start_times', ...
%                     'stretch_end_times' ...
%                   )
% 
%                 disp(['Finding Relevant Data Stretches: condition ' condition_list{i_condition} ', Trial ' num2str(i_trial) ' completed, found ' num2str(length(stretch_start_times)) ' relevant stretches, saved as ' data_stretches_file_name]);                

            end
            if strcmp(condition_stimulus, 'VISUAL') || strcmp(condition_stimulus, 'GVS')
                % for each trigger, extract conditions and relevant step events
                number_of_triggers = length(trigger_indices_mocap);
                removal_flags = zeros(number_of_triggers, 1);
                stim_start_indices_labview = [stim_start_indices_labview zeros(number_of_triggers, 5)]; %#ok<AGROW>
                stretch_start_times = zeros(number_of_triggers, 6);
                stretch_end_times = zeros(number_of_triggers, 6);
                stretch_pushoff_times = zeros(number_of_triggers, 6);
                closest_heelstrike_distance_times = zeros(number_of_triggers, 1);
                condition_stance_foot_list_trial = cell(number_of_triggers, 6);
                condition_perturbation_list_trial = cell(number_of_triggers, 6);
                condition_delay_list_trial = cell(number_of_triggers, 6);
                condition_index_list_trial = cell(number_of_triggers, 6);
                condition_experimental_list_trial = cell(number_of_triggers, 6);
                condition_stimulus_list_trial = cell(number_of_triggers, 6);
                condition_day_list_trial = cell(number_of_triggers, 6);
                
                for i_trigger = 1 : number_of_triggers
                   % perturbation condition
                    if illusion_trajectory(stim_start_indices_labview(i_trigger)) > 0
                        condition_perturbation_list_trial{i_trigger, 1} = 'ILLUSION_RIGHT';
                        condition_perturbation_list_trial{i_trigger, 2} = 'ILLUSION_RIGHT';
                        condition_perturbation_list_trial{i_trigger, 3} = 'ILLUSION_RIGHT';
                        condition_perturbation_list_trial{i_trigger, 4} = 'ILLUSION_RIGHT';
                    elseif illusion_trajectory(stim_start_indices_labview(i_trigger)) < 0
                        condition_perturbation_list_trial{i_trigger, 1} = 'ILLUSION_LEFT';
                        condition_perturbation_list_trial{i_trigger, 2} = 'ILLUSION_LEFT';
                        condition_perturbation_list_trial{i_trigger, 3} = 'ILLUSION_LEFT';
                        condition_perturbation_list_trial{i_trigger, 4} = 'ILLUSION_LEFT';
                    else
        %                 disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_stimulus(trigger_indices_labview(i_trigger))) ' - no stim']);
                    end
                    condition_perturbation_list_trial{i_trigger, 5} = 'CONTROL';
                    condition_perturbation_list_trial{i_trigger, 6} = 'CONTROL';
                
                    % delay condition
                    wait_time_stim = time_stimulus(stim_start_indices_labview(i_trigger)) - time_stimulus(trigger_indices_labview(i_trigger));
                    [~, wait_condition_index] = min(abs(study_settings.delay_times - wait_time_stim));
                    if iscell(study_settings.delay_time_labels)
                        delay_condition_label = study_settings.delay_time_labels{wait_condition_index};
                        condition_delay_list_trial{i_trigger, 5} = 'CONTROL';
                        condition_delay_list_trial{i_trigger, 6} = 'CONTROL';
                    else
                        delay_condition_label = study_settings.delay_time_labels;
                        condition_delay_list_trial{i_trigger, 5} = delay_condition_label;
                        condition_delay_list_trial{i_trigger, 6} = delay_condition_label;
                    end
                    condition_delay_list_trial{i_trigger, 1} = delay_condition_label;
                    condition_delay_list_trial{i_trigger, 2} = delay_condition_label;
                    condition_delay_list_trial{i_trigger, 3} = delay_condition_label;
                    condition_delay_list_trial{i_trigger, 4} = delay_condition_label;
                    
                    % experimental condition
                    condition_experimental_list_trial{i_trigger, 1} = condition_experimental;
                    condition_experimental_list_trial{i_trigger, 2} = condition_experimental;
                    condition_experimental_list_trial{i_trigger, 3} = condition_experimental;
                    condition_experimental_list_trial{i_trigger, 4} = condition_experimental;
                    condition_experimental_list_trial{i_trigger, 5} = condition_experimental;
                    condition_experimental_list_trial{i_trigger, 6} = condition_experimental;
                    
                    % stimulus condition
                    condition_stimulus_list_trial{i_trigger, 1} = condition_stimulus;
                    condition_stimulus_list_trial{i_trigger, 2} = condition_stimulus;
                    condition_stimulus_list_trial{i_trigger, 3} = condition_stimulus;
                    condition_stimulus_list_trial{i_trigger, 4} = condition_stimulus;
                    condition_stimulus_list_trial{i_trigger, 5} = condition_stimulus;
                    condition_stimulus_list_trial{i_trigger, 6} = condition_stimulus;
                    
                    % day condition
                    condition_day_list_trial{i_trigger, 1} = condition_day;
                    condition_day_list_trial{i_trigger, 2} = condition_day;
                    condition_day_list_trial{i_trigger, 3} = condition_day;
                    condition_day_list_trial{i_trigger, 4} = condition_day;
                    condition_day_list_trial{i_trigger, 5} = condition_day;
                    condition_day_list_trial{i_trigger, 6} = condition_day;
                    
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
                            right_foot_heelstrike_minus_1 = NaN;
                            right_foot_heelstrike_0 = NaN;
                            right_foot_heelstrike_plus_1 = NaN;
                            right_foot_heelstrike_plus_2 = NaN;
                            left_foot_heelstrike_minus_1 = NaN;
                            left_foot_heelstrike_0 = NaN;
                            left_foot_heelstrike_plus_1 = NaN;
                            left_foot_heelstrike_plus_2 = NaN;
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
                            left_foot_heelstrike_minus_1 = NaN;
                            left_foot_heelstrike_0 = NaN;
                            left_foot_heelstrike_plus_1 = NaN;
                            left_foot_heelstrike_plus_2 = NaN;
                            right_foot_heelstrike_minus_1 = NaN;
                            right_foot_heelstrike_0 = NaN;
                            right_foot_heelstrike_plus_1 = NaN;
                            right_foot_heelstrike_plus_2 = NaN;
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
                                plot([left_foot_heelstrike_minus_1 left_foot_heelstrike_0 left_foot_heelstrike_plus_1 left_foot_heelstrike_plus_2], [0 0 0 0]-0.01, 'v', 'linewidth', 3);
                                plot([left_foot_pushoff_minus_1  left_foot_pushoff_0 left_foot_pushoff_plus_1 left_foot_pushoff_plus_2], [0 0 0 0]-0.01, '^', 'linewidth', 3);
                                plot([right_foot_heelstrike_minus_1 right_foot_heelstrike_0 right_foot_heelstrike_plus_1 right_foot_heelstrike_plus_2], [0 0 0 0]+0.01, 'v', 'linewidth', 3);
                                plot([right_foot_pushoff_minus_1  right_foot_pushoff_0 right_foot_pushoff_plus_1 right_foot_pushoff_plus_2], [0 0 0 0]+0.01, '^', 'linewidth', 3);
                            end            
                        end            
                    else
                        trigger_foot = 'unclear';
                        disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_stimulus(trigger_indices_labview(i_trigger))) ' - trigger exactly between two heelstrikes']);
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
                            condition_index_list_trial{i_trigger, 1} = 'ONE';
                            condition_stance_foot_list_trial{i_trigger, 1} = 'STANCE_RIGHT';
                            % post stimulus steps
                            condition_index_list_trial{i_trigger, 2} = 'TWO';
                            condition_stance_foot_list_trial{i_trigger, 2} = 'STANCE_LEFT';
                            condition_index_list_trial{i_trigger, 3} = 'THREE';
                            condition_stance_foot_list_trial{i_trigger, 3} = 'STANCE_RIGHT';
                            condition_index_list_trial{i_trigger, 4} = 'FOUR';
                            condition_stance_foot_list_trial{i_trigger, 4} = 'STANCE_LEFT';
                            % use two previous steps as control
                            condition_index_list_trial{i_trigger, 5} = 'CONTROL';
                            condition_stance_foot_list_trial{i_trigger, 5} = 'STANCE_RIGHT';
                            condition_index_list_trial{i_trigger, 6} = 'CONTROL';
                            condition_stance_foot_list_trial{i_trigger, 6} = 'STANCE_LEFT';

                            if strcmp(stretches_to_analyze, 'single stance')
                                % ACHTUNG: this is not up to date
                                
                                % this step
                                stretch_start_times(i_trigger, 1) = left_foot_pushoff_0;
                                stretch_end_times(i_trigger, 1) = left_foot_heelstrike_0;

                                % post-stimulus steps
                                stretch_start_times(i_trigger, 2) = right_foot_pushoff_0;
                                stretch_end_times(i_trigger, 2) = right_foot_heelstrike_plus_1;

                                stretch_start_times(i_trigger, 3) = left_foot_pushoff_plus_1;
                                stretch_end_times(i_trigger, 3) = left_foot_heelstrike_plus_1;
                                stretch_start_times(i_trigger, 4) = right_foot_pushoff_plus_1;
                                stretch_end_times(i_trigger, 4) = right_foot_heelstrike_plus_2;

                                % use two previous steps as control
                                stretch_start_times(i_trigger, 5) = left_foot_pushoff_minus_1;
                                stretch_end_times(i_trigger, 5) = left_foot_heelstrike_minus_1;
                                stretch_start_times(i_trigger, 6) = right_foot_pushoff_minus_1;
                                stretch_end_times(i_trigger, 6) = right_foot_heelstrike_0;

                            elseif strcmp(stretches_to_analyze, 'full stride')
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
                            end

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
                            condition_index_list_trial{i_trigger, 1} = 'ONE';
                            condition_stance_foot_list_trial{i_trigger, 1} = 'STANCE_LEFT';
                            % post stimulus steps
                            condition_index_list_trial{i_trigger, 2} = 'TWO';
                            condition_stance_foot_list_trial{i_trigger, 2} = 'STANCE_RIGHT';
                            condition_index_list_trial{i_trigger, 3} = 'THREE';
                            condition_stance_foot_list_trial{i_trigger, 3} = 'STANCE_LEFT';
                            condition_index_list_trial{i_trigger, 4} = 'FOUR';
                            condition_stance_foot_list_trial{i_trigger, 4} = 'STANCE_RIGHT';
                            % use two previous steps as control
                            condition_index_list_trial{i_trigger, 5} = 'CONTROL';
                            condition_stance_foot_list_trial{i_trigger, 5} = 'STANCE_LEFT';
                            condition_index_list_trial{i_trigger, 6} = 'CONTROL';
                            condition_stance_foot_list_trial{i_trigger, 6} = 'STANCE_RIGHT';
                            if strcmp(stretches_to_analyze, 'single stance')
                                % this step
                                stretch_start_times(i_trigger, 1) = right_foot_pushoff_0;
                                stretch_end_times(i_trigger, 1) = right_foot_heelstrike_0;
                                % post stimulus steps
                                stretch_start_times(i_trigger, 2) = left_foot_pushoff_0;
                                stretch_end_times(i_trigger, 2) = left_foot_heelstrike_plus_1;
                                stretch_start_times(i_trigger, 3) = right_foot_pushoff_plus_1;
                                stretch_end_times(i_trigger, 3) = right_foot_heelstrike_plus_1;
                                stretch_start_times(i_trigger, 4) = left_foot_pushoff_plus_1;
                                stretch_end_times(i_trigger, 4) = left_foot_heelstrike_plus_2;
                                % use two previous steps as control
                                stretch_start_times(i_trigger, 3) = right_foot_pushoff_minus_1;
                                stretch_end_times(i_trigger, 3) = right_foot_heelstrike_minus_1;
                                stretch_start_times(i_trigger, 4) = left_foot_pushoff_minus_1;
                                stretch_end_times(i_trigger, 4) = left_foot_heelstrike_0;

                            elseif strcmp(stretches_to_analyze, 'full stride')
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
                            end


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
                            condition_stance_foot_list_trial{i_trigger, 1} = 'UNCLEAR';
                            condition_stance_foot_list_trial{i_trigger, 2} = 'UNCLEAR';
                            condition_stance_foot_list_trial{i_trigger, 3} = 'UNCLEAR';
                            condition_stance_foot_list_trial{i_trigger, 4} = 'UNCLEAR';
                            condition_stance_foot_list_trial{i_trigger, 5} = 'UNCLEAR';
                            condition_stance_foot_list_trial{i_trigger, 6} = 'UNCLEAR';
                        end
                    end
                end

                % remove flagged triggers
                unflagged_indices = ~removal_flags;
                trigger_times = trigger_times(unflagged_indices);
                stim_start_indices_labview = stim_start_indices_labview(unflagged_indices, :);
                stretch_start_times = stretch_start_times(unflagged_indices, :);
                stretch_pushoff_times = stretch_pushoff_times(unflagged_indices, :);
                stretch_end_times = stretch_end_times(unflagged_indices, :);
                condition_stance_foot_list_trial = condition_stance_foot_list_trial(unflagged_indices, :);
                condition_perturbation_list_trial = condition_perturbation_list_trial(unflagged_indices, :);
                condition_delay_list_trial = condition_delay_list_trial(unflagged_indices, :);
                condition_index_list_trial = condition_index_list_trial(unflagged_indices, :);
                condition_experimental_list_trial = condition_experimental_list_trial(unflagged_indices, :); % XXX needs to be updated
                condition_stimulus_list_trial = condition_stimulus_list_trial(unflagged_indices, :); % XXX needs to be updated
                condition_day_list_trial = condition_day_list_trial(unflagged_indices, :); % XXX needs to be updated
                closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :);
                
                % reorder
                % XXX this should be replaced with a reshape
                stretch_start_times = reshape(stretch_start_times, numel(stretch_start_times), 1);
                stretch_pushoff_times = reshape(stretch_pushoff_times, numel(stretch_pushoff_times), 1);
                stretch_end_times = reshape(stretch_end_times, numel(stretch_end_times), 1);
                
                condition_stance_foot_list_trial = reshape(condition_stance_foot_list_trial, numel(condition_stance_foot_list_trial), 1);
                condition_perturbation_list_trial = reshape(condition_perturbation_list_trial, numel(condition_perturbation_list_trial), 1);
                condition_delay_list_trial = reshape(condition_delay_list_trial, numel(condition_delay_list_trial), 1);
                condition_index_list_trial = reshape(condition_index_list_trial, numel(condition_index_list_trial), 1);
                condition_experimental_list_trial = reshape(condition_experimental_list_trial, numel(condition_experimental_list_trial), 1);
                condition_stimulus_list_trial = reshape(condition_stimulus_list_trial, numel(condition_stimulus_list_trial), 1);
                condition_day_list_trial = reshape(condition_day_list_trial, numel(condition_day_list_trial), 1);
                
                % we now have a neatly ordered list of stretches which we can prune
                % check step times and flag outliers
                number_of_stretches = length(stretch_start_times);
                stretch_times = stretch_end_times - stretch_start_times;
                stretch_time_outlier_limits = median(stretch_times) * [0.8 1.2];

                removal_flags = zeros(number_of_stretches, 1);
                removal_flags(stretch_times < stretch_time_outlier_limits(1)) = 1;
                removal_flags(stretch_times > stretch_time_outlier_limits(2)) = 1;

                % check data availability and flag stretches with gaps
%                 for i_stretch = 1 : number_of_stretches
%                     [~, start_index_mocap] = min(abs(time_marker - stretch_start_times(i_stretch)));
%                     [~, end_index_mocap] = min(abs(time_marker - stretch_end_times(i_stretch)));
%                     if any(any(isnan(marker_trajectories(start_index_mocap : end_index_mocap, :))))
%                         removal_flags(i_stretch) = 1;
%                     end
%                 end
% XXX removed for now. Only those stretches with gaps in markers that actually interest me should be flagged.
% information about markers of interest should go to the experiment file

                % remove flagged triggers
                unflagged_indices = ~removal_flags;
                stretch_start_times = stretch_start_times(unflagged_indices, :);
                stretch_pushoff_times = stretch_pushoff_times(unflagged_indices, :);
                stretch_end_times = stretch_end_times(unflagged_indices, :);
                condition_stance_foot_list_trial = condition_stance_foot_list_trial(unflagged_indices, :);
                condition_perturbation_list_trial = condition_perturbation_list_trial(unflagged_indices, :);
                condition_delay_list_trial = condition_delay_list_trial(unflagged_indices, :);
                condition_index_list_trial = condition_index_list_trial(unflagged_indices, :);
                condition_experimental_list_trial = condition_experimental_list_trial(unflagged_indices, :);
                condition_stimulus_list_trial = condition_stimulus_list_trial(unflagged_indices, :);
                condition_day_list_trial = condition_day_list_trial(unflagged_indices, :);

                
                % add new variables to be saved
                event_variables_to_save.condition_stance_foot_list_trial = condition_stance_foot_list_trial;
                event_variables_to_save.condition_perturbation_list_trial = condition_perturbation_list_trial;
                event_variables_to_save.condition_delay_list_trial = condition_delay_list_trial;
                event_variables_to_save.condition_index_list_trial = condition_index_list_trial;
                event_variables_to_save.condition_experimental_list_trial = condition_experimental_list_trial;
                event_variables_to_save.condition_stimulus_list_trial = condition_stimulus_list_trial;
                event_variables_to_save.condition_day_list_trial = condition_day_list_trial;
                event_variables_to_save.stretch_start_times = stretch_start_times;
                event_variables_to_save.stretch_pushoff_times = stretch_pushoff_times;
                event_variables_to_save.stretch_end_times = stretch_end_times;
%                 event_variables_to_save.trigger_times = trigger_times;
%                 event_variables_to_save.closest_heelstrike_distance_times = closest_heelstrike_distance_times;
                
%                 % save data
%                 data_stretches_file_name = ['analysis' filesep makeFileName(date, subject_id, condition_list{i_condition}, i_trial, 'relevantDataStretches')];
%                 save ...
%                   ( ...
%                     data_stretches_file_name, ...
%                     'condition_stance_foot_list_trial', ...
%                     'condition_perturbation_list_trial', ...
%                     'condition_delay_list_trial', ...
%                     'condition_index_list_trial', ...
%                     'condition_experimental_list_trial', ...
%                     'condition_stimulus_list_trial', ...
%                     'condition_day_list_trial', ...
%                     'trigger_times', ...
%                     'stretch_start_times', ...
%                     'stretch_pushoff_times', ...
%                     'stretch_end_times', ...
%                     'closest_heelstrike_distance_times' ...
%                   )
% 
%                 disp(['Finding Relevant Data Stretches: condition ' condition_list{i_condition} ', Trial ' num2str(i_trial) ' completed, found ' num2str(length(stretch_start_times)) ' relevant stretches, saved as ' data_stretches_file_name]);                
                
                
            end            

            %% remove stretches where important variables are missing
            %
            % calculate variables that depend upon the step events to be identified correctly
            variables_to_save = struct;
            variables_to_prune_for = {};
            save_folder = 'processed';
            save_file_name = makeFileName(date, subject_id, condition_list{i_condition}, i_trial, 'kinematicTrajectories.mat');
            if any(strcmp(study_settings.variables_to_analyze(:, 1), 'left_arm_phase'))
                LELB_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'LELB');
                LWRA_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'LWRA');
                LWRB_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'LWRB');
                
                % calculate vectors
                left_wrist_center_trajectory = (LWRA_trajectory + LWRB_trajectory) * 0.5;
                left_arm_vector_trajectory = LELB_trajectory - left_wrist_center_trajectory;

                % calculate angles
                left_arm_angle = rad2deg(atan2(-left_arm_vector_trajectory(:, 2), left_arm_vector_trajectory(:, 3)));

                % find negative peaks
                [~, left_arm_peak_locations] = findpeaks(-left_arm_angle, 'MinPeakProminence', subject_settings.left_armswing_peak_prominence_threshold, 'MinPeakDistance', subject_settings.left_armswing_peak_distance_threshold * sampling_rate_marker);
                
                % normalize
                [larm_angle_normalized, larm_angle_dot_normalized] = normalizePeriodicVariable(left_arm_angle, time_marker, left_arm_peak_locations);

                % calculate phase
                left_arm_phase = atan2(larm_angle_dot_normalized, larm_angle_normalized);
                
%                 % XXX plot some stuff to check
%                 figure; hold on
%                 plot(left_arm_phase)
%                 plot(left_arm_phase_atan2)
%                 
                
                % add new variables to be saved
                variables_to_save.left_arm_angle = left_arm_angle;
                variables_to_save.left_arm_phase = left_arm_phase;
                variables_to_save.sampling_rate_marker = sampling_rate_marker;
                variables_to_save.time_marker = time_marker;
                saveDataToFile([save_folder filesep save_file_name], variables_to_save);
                addAvailableData('left_arm_angle', 'time_marker', 'sampling_rate_marker', '', save_folder, save_file_name);
                addAvailableData('left_arm_phase', 'time_marker', 'sampling_rate_marker', '', save_folder, save_file_name);
                variables_to_prune_for = [variables_to_prune_for; 'left_arm_angle']; %#ok<AGROW>
                variables_to_prune_for = [variables_to_prune_for; 'left_arm_phase']; %#ok<AGROW>
            end
            if any(strcmp(study_settings.variables_to_analyze(:, 1), 'right_arm_phase'))
                RELB_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'RELB');
                RWRA_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'RWRA');
                RWRB_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'RWRB');
                
                % calculate vectors
                right_wrist_center_trajectory = (RWRA_trajectory + RWRB_trajectory) * 0.5;
                right_arm_vector_trajectory = RELB_trajectory - right_wrist_center_trajectory;
                
                % calculate angles
                right_arm_angle = rad2deg(atan2(-right_arm_vector_trajectory(:, 2), right_arm_vector_trajectory(:, 3)));

                % find negative peaks
                [~, right_arm_peak_locations] = findpeaks(-right_arm_angle, 'MinPeakProminence', subject_settings.right_armswing_peak_prominence_threshold, 'MinPeakDistance', subject_settings.right_armswing_peak_distance_threshold * sampling_rate_marker);
                
                % normalize
                [larm_angle_normalized, larm_angle_dot_normalized] = normalizePeriodicVariable(right_arm_angle, time_marker, right_arm_peak_locations);

                % calculate phase
                right_arm_phase = atan2(larm_angle_dot_normalized, larm_angle_normalized);
                
                % add new variables to be saved
                variables_to_save.right_arm_angle = right_arm_angle;
                variables_to_save.right_arm_phase = right_arm_phase;
                saveDataToFile([save_folder filesep save_file_name], variables_to_save);
                addAvailableData('right_arm_angle', 'time_marker', 'sampling_rate_marker', '', save_folder, save_file_name);
                addAvailableData('right_arm_phase', 'time_marker', 'sampling_rate_marker', '', save_folder, save_file_name);
                variables_to_prune_for = [variables_to_prune_for; 'right_arm_angle']; %#ok<AGROW>
                variables_to_prune_for = [variables_to_prune_for; 'right_arm_phase']; %#ok<AGROW>
            end
            if any(strcmp(study_settings.variables_to_analyze(:, 1), 'left_leg_phase'))
                LANK_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'LANK');
                LPSI_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'LPSI');
                LASI_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'LASI');
                
                % calculate vectors
                left_pelvis_center_trajectory = (LPSI_trajectory + LASI_trajectory) * 0.5;
                left_leg_vector_trajectory = left_pelvis_center_trajectory - LANK_trajectory;
                
                % calculate angles
                left_leg_angle = rad2deg(atan2(-left_leg_vector_trajectory(:, 2), left_leg_vector_trajectory(:, 3)));

                % find negative peaks
                [~, left_leg_peak_locations] = findpeaks(-left_leg_angle, 'MinPeakProminence', subject_settings.left_legswing_peak_prominence_threshold, 'MinPeakDistance', subject_settings.left_legswing_peak_distance_threshold * sampling_rate_marker);
                
                % normalize
                [lleg_angle_normalized, lleg_angle_dot_normalized] = normalizePeriodicVariable(left_leg_angle, time_marker, left_leg_peak_locations);

                % calculate phase
                left_leg_phase = atan2(lleg_angle_dot_normalized, lleg_angle_normalized);
                
                % add new variables to be saved
                variables_to_save.left_leg_angle = left_leg_angle;
                variables_to_save.left_leg_phase = left_leg_phase;
                saveDataToFile([save_folder filesep save_file_name], variables_to_save);
                addAvailableData('left_leg_angle', 'time_marker', 'sampling_rate_marker', '', save_folder, save_file_name);
                addAvailableData('left_leg_phase', 'time_marker', 'sampling_rate_marker', '', save_folder, save_file_name);
                variables_to_prune_for = [variables_to_prune_for; 'left_leg_angle']; %#ok<AGROW>
                variables_to_prune_for = [variables_to_prune_for; 'left_leg_phase']; %#ok<AGROW>
            end
            if any(strcmp(study_settings.variables_to_analyze(:, 1), 'right_leg_phase'))
                RANK_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'RANK');
                RPSI_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'RPSI');
                RASI_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'RASI');
                
                % calculate vectors
                right_pelvis_center_trajectory = (RPSI_trajectory + RASI_trajectory) * 0.5;
                right_leg_vector_trajectory = right_pelvis_center_trajectory - RANK_trajectory;
                
                % calculate angles
                right_leg_angle = rad2deg(atan2(-right_leg_vector_trajectory(:, 2), right_leg_vector_trajectory(:, 3)));

                % find negative peaks
                [~, right_leg_peak_locations] = findpeaks(-right_leg_angle, 'MinPeakProminence', subject_settings.right_legswing_peak_prominence_threshold, 'MinPeakDistance', subject_settings.right_legswing_peak_distance_threshold * sampling_rate_marker);
                
                % normalize
                [lleg_angle_normalized, lleg_angle_dot_normalized] = normalizePeriodicVariable(right_leg_angle, time_marker, right_leg_peak_locations);

                % calculate phase
                right_leg_phase = atan2(lleg_angle_dot_normalized, lleg_angle_normalized);
                
                % add new variables to be saved
                variables_to_save.right_leg_angle = right_leg_angle;
                variables_to_save.right_leg_phase = right_leg_phase;
                saveDataToFile([save_folder filesep save_file_name], variables_to_save);
                addAvailableData('right_leg_angle', 'time_marker', 'sampling_rate_marker', '', save_folder, save_file_name);
                addAvailableData('right_leg_phase', 'time_marker', 'sampling_rate_marker', '', save_folder, save_file_name);
                variables_to_prune_for = [variables_to_prune_for; 'right_leg_angle']; %#ok<AGROW>
                variables_to_prune_for = [variables_to_prune_for; 'right_leg_phase']; %#ok<AGROW>
            end
            
            
            
            
            
            
            
            
            % prune
            number_of_stretches = length(stretch_start_times);
            removal_flags = zeros(number_of_stretches, 1);
            
            % take care of steps with very large or small step time
            stretch_times = stretch_end_times - stretch_start_times;
            stretch_time_outlier_limits = median(stretch_times) * [0.8 1.2];
            removal_flags(stretch_times < stretch_time_outlier_limits(1)) = 1;
            removal_flags(stretch_times > stretch_time_outlier_limits(2)) = 1;

            % check data availability for markers and flag stretches with gaps
            for i_stretch = 1 : number_of_stretches
                [~, start_index_mocap] = min(abs(time_marker - stretch_start_times(i_stretch)));
                [~, end_index_mocap] = min(abs(time_marker - stretch_end_times(i_stretch)));
                if any(any(isnan(marker_trajectories(start_index_mocap : end_index_mocap, essential_marker_indicator))))
                    removal_flags(i_stretch) = 1;
                end
            end
            
            % check data availability for variables just calculated here
            for i_variable = 1 : length(variables_to_prune_for)
                [data, time] = loadData(date, subject_id, condition_list{i_condition}, i_trial, variables_to_prune_for{i_variable});
                
                for i_stretch = 1 : number_of_stretches
                    [~, start_index] = min(abs(time - stretch_start_times(i_stretch)));
                    [~, end_index] = min(abs(time - stretch_end_times(i_stretch)));
                    if any(isnan(data(start_index : end_index)))
                        removal_flags(i_stretch) = 1;
                    end
                end
                
                
            end


            % remove flagged triggers
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
            
            
            
            %% save
            %
            stretches_file_name = ['analysis' filesep makeFileName(date, subject_id, condition_list{i_condition}, i_trial, 'relevantDataStretches')];
            save(stretches_file_name, '-struct', 'event_variables_to_save');
            
            
            disp(['Finding Relevant Data Stretches: condition ' condition_list{i_condition} ', Trial ' num2str(i_trial) ' completed, found ' num2str(length(stretch_start_times)) ' relevant stretches, saved as ' stretches_file_name]);                
        end

    end
end




















