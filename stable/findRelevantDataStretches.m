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
% condition indicators. There are five condition types: stance_foot, perturbation, delay, index and experimental. Each
% condition type has an indicator variable that is of the same length as the list of stretch start/end times. For each
% pair of start/end-times, the condition indicator will be 1 if the data stretch between these two points is of this
% condition, and 0 otherwise.

% The condition types are:
% stance_foot: this is the foot that is in stance throughout the whole stretch, will be STANCE_LEFT, STANCE_RIGHT or STANCE_BOTH
% perturbation: direction of the perturbation in stimulus paradigms, will be ILLUSION_LEFT, ILLUSION_RIGHT or CONTROL (should be ILLUSION_NONE instead of CONTROL)
% delay: delay between trigger and stimulus in stimulus paradigms, options are specified in studySettings.txt
% index: number of steps after the stimulus in stimulus paradigms, will be STEP_ONE, STEP_TWO, STEP_THREE, STEP_FOUR or CONTROL
% experimental: condition affecting a whole trial instead of stretches of time within a trial. This will be loaded from
%               conditions.mat if available, otherwise it will be set to the trial_type as determined from the file name.

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
% - stretch_pushoff_times
% - stretch_start_times
% - stretch_end_times
                    
function findRelevantDataStretches(varargin)
    % parse arguments
    parser = inputParser;
    parser.KeepUnmatched = true;
    default_stimulus_type = 'none';
    valid_stimulus_types = {'none', 'gvs', 'visual', 'obstacle'};
    check_stimulus_type = @(x) any(validatestring(x,valid_stimulus_types));
    addParameter(parser, 'stimulus', default_stimulus_type, check_stimulus_type)
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    stimulus_type = parser.Results.stimulus;
    visualize = parser.Results.visualize;
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});

    % stretches_to_analyze = 'single stance';
    stretches_to_analyze = 'full stride';

    %% prepare
    load('subjectInfo.mat', 'date', 'subject_id');
    study_settings = loadSettingsFile(['..' filesep 'studySettings.txt']);

    time_to_nearest_heelstrike_before_trigger_threshold = 0.10; % a heelstrike should happen less than this long after a trigger
    time_to_nearest_heelstrike_after_trigger_threshold = 0.3; % a heelstrike should happen less than this long after a trigger


    %% extract data
    if strcmp(study_settings.experimental_condition_source, 'conditions_file')
        load('conditions.mat');
    end    
    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            %% load data
            condition = condition_list{i_condition};
            load(['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'stepEvents')]);
            
            % marker data
            [marker_trajectories, time_marker, ~, marker_labels] = loadData(date, subject_id, condition, i_trial, 'marker_trajectories');
            
            % forceplate data
            [left_forceplate_cop_world_trajectory, time_left_forceplate, ~, ~, left_forceplate_available] = loadData(date, subject_id, condition, i_trial, 'left_forceplate_cop_world', 'optional');
            [right_forceplate_cop_world_trajectory, time_right_forceplate, ~, ~, right_forceplate_available] = loadData(date, subject_id, condition, i_trial, 'right_forceplate_cop_world', 'optional');
            if left_forceplate_available && right_forceplate_available
                left_copx_trajectory = left_forceplate_cop_world_trajectory(:, 1);
                right_copx_trajectory = right_forceplate_cop_world_trajectory(:, 1);
                
                if any(strcmp(condition, study_settings.trial_types_with_inverted_forceplate_sides))
                    % this was for the VEPO lab, where the sides are inverted
                    left_copx_trajectory = right_forceplate_cop_world_trajectory(:, 3);
                    right_copx_trajectory = left_forceplate_cop_world_trajectory(:, 3);
                end
            end
            
            % stimulus data
            if strcmp(stimulus_type, 'gvs')
                GVS_out_trajectory = loadData(date, subject_id, condition, i_trial, 'GVS_out_trajectory');
                [stimulus_state_trajectory, time_stimulus] = loadData(date, subject_id, condition, i_trial, 'stimulus_state_trajectory');
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
            if strcmp(stimulus_type, 'gvs')
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
            if strcmp(stimulus_type, 'visual')
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
            if strcmp(stimulus_type, 'none')
                % use all touchdown events as triggers
                trigger_times = [left_touchdown_times; right_touchdown_times];
            elseif strcmp(stimulus_type, 'visual') || strcmp(stimulus_type, 'gvs')
                % find the time steps where the stimulus state crosses a threshold
                stimulus_threshold = 0.5;
                trigger_indices_labview = find(diff(sign(stimulus_state_trajectory - stimulus_threshold)) > 0) + 1;

                %
                epsilon = 1e-5;
                stim_start_indices_labview = find(diff(sign(abs(illusion_trajectory) - epsilon)) > 0) + 1;
                trigger_indices_labview = trigger_indices_labview(1 : length(stim_start_indices_labview)); % in case a stim is triggered, but not recorded
                
                trigger_times = time_stimulus(trigger_indices_labview);
            elseif strcmp(stimulus_type, 'obstacle')
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
            if strcmp(stimulus_type, 'none')
                stretch_start_times = zeros(number_of_triggers, 1);
                stretch_end_times = zeros(number_of_triggers, 1);
                stretch_pushoff_times = zeros(number_of_triggers, 1);
                closest_heelstrike_distance_times = zeros(number_of_triggers, 1);
                condition_stance_foot_list = cell(number_of_triggers, 1);
                condition_perturbation_list = cell(number_of_triggers, 1);
                condition_delay_list = cell(number_of_triggers, 1);
                condition_index_list = cell(number_of_triggers, 1);
                condition_experimental_list = cell(number_of_triggers, 1);
                
                for i_trigger = 1 : number_of_triggers
                    condition_perturbation_list{i_trigger, 1} = 'CONTROL';
                    condition_delay_list{i_trigger, 1} = 'CONTROL';
                    condition_index_list{i_trigger, 1} = 'CONTROL';
                    condition_experimental_list{i_trigger, 1} = condition;
                
                    % find out which heelstrike triggered
                    % XXX change this to use the interval, same as in the stimulus case
                    [distance_to_trigger_left_time, index_left] = min(abs(left_touchdown_times - trigger_times(i_trigger)));
                    [distance_to_trigger_right_time, index_right] = min(abs(right_touchdown_times - trigger_times(i_trigger)));
                    closest_heelstrike_distance_time = min([distance_to_trigger_left_time distance_to_trigger_right_time]);
                    
                    if distance_to_trigger_left_time < distance_to_trigger_right_time
                        condition_stance_foot_list{i_trigger, 1} = 'STANCE_LEFT';
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
                    condition_experimental_list{i_trigger, 1} = condition;
                   
                end
                
                % remove flagged triggers
                unflagged_indices = ~removal_flags;
                trigger_times = trigger_times(unflagged_indices);
                stretch_start_times = stretch_start_times(unflagged_indices, :);
                stretch_end_times = stretch_end_times(unflagged_indices, :);
                stretch_pushoff_times = stretch_pushoff_times(unflagged_indices, :);
                condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
                condition_perturbation_list = condition_perturbation_list(unflagged_indices, :);
                condition_delay_list = condition_delay_list(unflagged_indices, :);
                condition_index_list = condition_index_list(unflagged_indices, :);
                condition_experimental_list = condition_experimental_list(unflagged_indices, :);
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
                
                % check step times and flag outliers
                number_of_stretches = length(stretch_start_times);
                stretch_times = stretch_end_times - stretch_start_times;
                stretch_time_outlier_limits = median(stretch_times) * [0.8 1.2];
                
                removal_flags = zeros(number_of_stretches, 1);
                removal_flags(stretch_times < stretch_time_outlier_limits(1)) = 1;
                removal_flags(stretch_times > stretch_time_outlier_limits(2)) = 1;
                
                % check data availability and flag stretches with gaps
                for i_stretch = 1 : number_of_stretches
                    [~, start_index_mocap] = min(abs(time_marker - stretch_start_times(i_stretch)));
                    [~, end_index_mocap] = min(abs(time_marker - stretch_end_times(i_stretch)));
                    if any(any(isnan(marker_trajectories(start_index_mocap : end_index_mocap, essential_marker_indicator))))
                        removal_flags(i_stretch) = 1;
                    end
                end
                
                
                % remove flagged triggers
                unflagged_indices = ~removal_flags;
                trigger_times = trigger_times(unflagged_indices);
                stretch_start_times = stretch_start_times(unflagged_indices, :);
                stretch_end_times = stretch_end_times(unflagged_indices, :);
                stretch_pushoff_times = stretch_pushoff_times(unflagged_indices, :);
                condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
                condition_perturbation_list = condition_perturbation_list(unflagged_indices, :);
                condition_delay_list = condition_delay_list(unflagged_indices, :);
                condition_index_list = condition_index_list(unflagged_indices, :);
                condition_experimental_list = condition_experimental_list(unflagged_indices, :);
                closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :);
                
                % save data
                data_stretches_file_name = ['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'relevantDataStretches')];
                save ...
                  ( ...
                    data_stretches_file_name, ...
                    'condition_stance_foot_list', ...
                    'condition_perturbation_list', ...
                    'condition_delay_list', ...
                    'condition_index_list', ...
                    'condition_experimental_list', ...
                    'stretch_start_times', ...
                    'stretch_pushoff_times', ...
                    'stretch_end_times' ...
                  )

                disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, found ' num2str(length(stretch_start_times)) ' relevant stretches, saved as ' data_stretches_file_name]);                
                
                
                
                
            end  
            if strcmp(stimulus_type, 'obstacle')
                stretch_start_times = zeros(3, 1);
                stretch_end_times = zeros(3, 1);
                stretch_pushoff_times = zeros(3, 1);
                
                condition_stance_foot_list = cell(3, 1);
                condition_perturbation_list = cell(3, 1);
                condition_delay_list = cell(3, 1);
                condition_index_list = cell(3, 1);
                condition_experimental_list = cell(3, 1);
                
                if ismember(i_trial, no_obstacle)
                    condition_experimental = 'NO';
                elseif ismember(i_trial, near_obstacle)
                    condition_experimental = 'NEAR';
                elseif ismember(i_trial, far_obstacle)
                    condition_experimental = 'FAR';
                else
                    error(['No condition specified for trial no. ' num2str(i_trial)]);
                end
                
                % second stretch goes from first right pushoff to first right heelstrike
                stretch_start_times(2, 1) = right_pushoff_times(1);
                stretch_end_times(2, 1) = right_touchdown_times(1);
                condition_perturbation_list{2, 1} = 'CONTROL';
                condition_delay_list{2, 1} = 'CONTROL';
                condition_index_list{2, 1} = 'TWO';
                condition_experimental_list{2, 1} = condition_experimental;
                condition_stance_foot_list{2, 1} = 'STANCE_LEFT';
                
                % third stretch goes from first right heelstrike to first left heelstrike
                stretch_start_times(3, 1) = right_touchdown_times(1);
                stretch_end_times(3, 1) = left_touchdown_times(1);
                condition_perturbation_list{3, 1} = 'CONTROL';
                condition_delay_list{3, 1} = 'CONTROL';
                condition_index_list{3, 1} = 'THREE';
                condition_experimental_list{3, 1} = condition_experimental;
                condition_stance_foot_list{3, 1} = 'STANCE_RIGHT';
                
                % first stretch goes from first right pushoff minus 1 second to first right pushoff
                stretch_start_times(1, 1) = right_pushoff_times(1) - 1;
                stretch_end_times(1, 1) = right_pushoff_times(1);
                condition_perturbation_list{1, 1} = 'CONTROL';
                condition_delay_list{1, 1} = 'CONTROL';
                condition_index_list{1, 1} = 'ONE';
                condition_experimental_list{1, 1} = condition_experimental;
                condition_stance_foot_list{1, 1} = 'STANCE_BOTH';
                
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
                        if strcmp(condition_stance_foot_list(i_trigger), 'STANCE_LEFT')
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
                
                % save data
                data_stretches_file_name = ['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'relevantDataStretches')];
                save ...
                  ( ...
                    data_stretches_file_name, ...
                    'condition_stance_foot_list', ...
                    'condition_perturbation_list', ...
                    'condition_delay_list', ...
                    'condition_index_list', ...
                    'condition_experimental_list', ...
                    'stretch_pushoff_times', ...
                    'stretch_start_times', ...
                    'stretch_end_times' ...
                  )

                disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, found ' num2str(length(stretch_start_times)) ' relevant stretches, saved as ' data_stretches_file_name]);                

            end
            if strcmp(stimulus_type, 'visual') || strcmp(stimulus_type, 'gvs')
                % for each trigger, extract conditions and relevant step events
                number_of_triggers = length(trigger_indices_mocap);
                removal_flags = zeros(number_of_triggers, 1);
                stim_start_indices_labview = [stim_start_indices_labview zeros(number_of_triggers, 5)]; %#ok<AGROW>
                stretch_start_times = zeros(number_of_triggers, 6);
                stretch_end_times = zeros(number_of_triggers, 6);
                stretch_pushoff_times = zeros(number_of_triggers, 6);
                closest_heelstrike_distance_times = zeros(number_of_triggers, 1);
                condition_stance_foot_list = cell(number_of_triggers, 6);
                condition_perturbation_list = cell(number_of_triggers, 6);
                condition_delay_list = cell(number_of_triggers, 6);
                condition_index_list = cell(number_of_triggers, 6);
                condition_experimental_list = cell(number_of_triggers, 6);
                
                for i_trigger = 1 : number_of_triggers
                   % perturbation condition
                    if illusion_trajectory(stim_start_indices_labview(i_trigger)) > 0
                        condition_perturbation_list{i_trigger, 1} = 'ILLUSION_RIGHT';
                        condition_perturbation_list{i_trigger, 2} = 'ILLUSION_RIGHT';
                        condition_perturbation_list{i_trigger, 3} = 'ILLUSION_RIGHT';
                        condition_perturbation_list{i_trigger, 4} = 'ILLUSION_RIGHT';
                    elseif illusion_trajectory(stim_start_indices_labview(i_trigger)) < 0
                        condition_perturbation_list{i_trigger, 1} = 'ILLUSION_LEFT';
                        condition_perturbation_list{i_trigger, 2} = 'ILLUSION_LEFT';
                        condition_perturbation_list{i_trigger, 3} = 'ILLUSION_LEFT';
                        condition_perturbation_list{i_trigger, 4} = 'ILLUSION_LEFT';
                    else
        %                 disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_stimulus(trigger_indices_labview(i_trigger))) ' - no stim']);
                    end
                    condition_perturbation_list{i_trigger, 5} = 'CONTROL';
                    condition_perturbation_list{i_trigger, 6} = 'CONTROL';
                
                    % delay condition
                    wait_time_stim = time_stimulus(stim_start_indices_labview(i_trigger)) - time_stimulus(trigger_indices_labview(i_trigger));
                    [~, wait_condition_index] = min(abs(study_settings.delay_times - wait_time_stim));
                    condition_delay_list{i_trigger, 1} = num2str(study_settings.delay_times(wait_condition_index));
                    condition_delay_list{i_trigger, 2} = num2str(study_settings.delay_times(wait_condition_index));
                    condition_delay_list{i_trigger, 3} = num2str(study_settings.delay_times(wait_condition_index));
                    condition_delay_list{i_trigger, 4} = num2str(study_settings.delay_times(wait_condition_index));
                    condition_delay_list{i_trigger, 5} = 'CONTROL';
                    condition_delay_list{i_trigger, 6} = 'CONTROL';
                    
                    % experimental condition
                    condition_experimental_list{i_trigger, 1} = condition;
                    condition_experimental_list{i_trigger, 2} = condition;
                    condition_experimental_list{i_trigger, 3} = condition;
                    condition_experimental_list{i_trigger, 4} = condition;
                    condition_experimental_list{i_trigger, 5} = condition;
                    condition_experimental_list{i_trigger, 6} = condition;
                
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
                unflagged_indices = ~removal_flags;
                trigger_times = trigger_times(unflagged_indices);
                stim_start_indices_labview = stim_start_indices_labview(unflagged_indices, :);
                stretch_start_times = stretch_start_times(unflagged_indices, :);
                stretch_pushoff_times = stretch_pushoff_times(unflagged_indices, :);
                stretch_end_times = stretch_end_times(unflagged_indices, :);
                condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
                condition_perturbation_list = condition_perturbation_list(unflagged_indices, :);
                condition_delay_list = condition_delay_list(unflagged_indices, :);
                condition_index_list = condition_index_list(unflagged_indices, :);
                condition_experimental_list = condition_experimental_list(unflagged_indices, :); % XXX needs to be updated
                closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :);
                
                % reorder
                % XXX this should be replaced with a reshape
                stretch_start_times = reshape(stretch_start_times, numel(stretch_start_times), 1);
                stretch_pushoff_times = reshape(stretch_pushoff_times, numel(stretch_pushoff_times), 1);
                stretch_end_times = reshape(stretch_end_times, numel(stretch_end_times), 1);
                
                condition_stance_foot_list = reshape(condition_stance_foot_list, numel(condition_stance_foot_list), 1);
                condition_perturbation_list = reshape(condition_perturbation_list, numel(condition_perturbation_list), 1);
                condition_delay_list = reshape(condition_delay_list, numel(condition_delay_list), 1);
                condition_index_list = reshape(condition_index_list, numel(condition_index_list), 1);
                condition_experimental_list = reshape(condition_experimental_list, numel(condition_experimental_list), 1);
                
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
                condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
                condition_perturbation_list = condition_perturbation_list(unflagged_indices, :);
                condition_delay_list = condition_delay_list(unflagged_indices, :);
                condition_index_list = condition_index_list(unflagged_indices, :);
                condition_experimental_list = condition_experimental_list(unflagged_indices, :);
                                
                % save data
                data_stretches_file_name = ['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'relevantDataStretches')];
                save ...
                  ( ...
                    data_stretches_file_name, ...
                    'condition_stance_foot_list', ...
                    'condition_perturbation_list', ...
                    'condition_delay_list', ...
                    'condition_index_list', ...
                    'condition_experimental_list', ...
                    'trigger_times', ...
                    'stretch_start_times', ...
                    'stretch_pushoff_times', ...
                    'stretch_end_times', ...
                    'closest_heelstrike_distance_times' ...
                  )

                disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, found ' num2str(length(stretch_start_times)) ' relevant stretches, saved as ' data_stretches_file_name]);                
                
                
            end            



        end

    end
end




















