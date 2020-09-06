%     This file is part of the CoBaL code base
%     Copyright (C) 2019 Hendrik Reimann <hendrikreimann@gmail.com>
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

function [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels(study_settings, subject_settings, trial_data)
    experimental_paradigm = study_settings.get('experimental_paradigm');
    subject_id = subject_settings.get('subject_id');
    gender = subject_settings.get('gender');
    
    % allocate
    conditions_trial = struct;
    event_variables_to_save = struct;
    removal_flags = zeros(size(trial_data.trigger_times));
    
    % determine levels particular to the experimental paradigm
    paradigms_with_intermittent_perturbation = {'Vision', 'CadenceVision', 'GVS', 'CadenceGVS', 'FatigueGVS', 'OculusLaneRestriction', 'CognitiveLoadVision', 'CognitiveLoadGvs', 'SR_VisualStim'};
    if any(strcmp(experimental_paradigm, paradigms_with_intermittent_perturbation))
        [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_intermittentPerturbations(study_settings, subject_settings, trial_data);
    end
    if strcmp(experimental_paradigm, 'Normal Walking')
        [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_normalWalking(subject_settings, trial_data);
    end    
    if strcmp(experimental_paradigm, 'Stochastic Resonance')
        [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_stochasticResonance(subject_settings, trial_data);
    end
    if strcmp(experimental_paradigm, 'Vision_old')
        [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_visionOld(study_settings, subject_settings, trial_data);
    end
    if strcmp(experimental_paradigm, 'GVS_old')
        [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_gvsOld(study_settings, subject_settings, trial_data);
    end
    if strcmp(experimental_paradigm, 'GvsOverground')
        [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_gvsOverground(study_settings, subject_settings, trial_data);
    end
    if strcmp(experimental_paradigm, 'GaitInitiationObstacle')
        [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_gaitInitiationObstacle(study_settings, subject_settings, trial_data);
    end
    if strcmp(experimental_paradigm, 'Vision Stochastic')
        [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_visionStochastic(study_settings, subject_settings, trial_data);
    end
    
    % add levels that are the same for all experimental paradigms
    condition_subject_list = cell(size(event_variables_to_save.stretch_times, 1), 1);
    for i_stretch = 1 : length(condition_subject_list)
        condition_subject_list{i_stretch} = subject_id;
    end
    conditions_trial.subject_list = condition_subject_list;
    
    condition_gender_list = cell(size(event_variables_to_save.stance_foot_data, 1), 1);
    for i_stretch = 1 : length(condition_gender_list)
        condition_gender_list{i_stretch} = gender;
    end
    conditions_trial.gender_list = condition_gender_list;
    
    % 2020-APR-02 HR: cleaning up, but these two are so old that there's no
    % data to test the code on. Moving code from determineStretchesToAnalyze
    % into sub-functions, I'll leave this here for historic reasons for now
%     if strcmp(condition_stimulus, 'NONE')
%         % determine start and end
%         stance_foot_data = {'STANCE_RIGHT', 'STANCE_BOTH', 'STANCE_LEFT'};
%         bands_per_stretch = length(stance_foot_data);
% 
%         stretch_start_times = zeros(number_of_triggers, 1);
%         stretch_end_times = zeros(number_of_triggers, 1);
%         stretch_pushoff_times = zeros(number_of_triggers, 1);
%         closest_heelstrike_distance_times = zeros(number_of_triggers, 1);
%         condition_stance_foot_list = cell(number_of_triggers, 1);
%         perturbation_list = cell(number_of_triggers, 1);
%         condition_delay_list = cell(number_of_triggers, 1);
%         condition_index_list = cell(number_of_triggers, 1);
%         condition_experimental_list = cell(number_of_triggers, 1);
%         condition_stimulus_list = cell(number_of_triggers, 1);
%         condition_day_list = cell(number_of_triggers, 1);
% 
%         for i_trigger = 1 : number_of_triggers
%             perturbation_list{i_trigger, 1} = 'N/A';
%             condition_delay_list{i_trigger, 1} = 'N/A';
%             condition_experimental_list{i_trigger, 1} = trial_data.condition_experimental;
%             condition_stimulus_list{i_trigger, 1} = condition_stimulus;
%             condition_day_list{i_trigger, 1} = condition_day;
% 
% 
%             % find out which heelstrike triggered
%             % XXX change this to use the interval, same as in the stimulus case
%             [distance_to_trigger_left_time, index_left] = min(abs(left_touchdown_times - trigger_times(i_trigger)));
%             [distance_to_trigger_right_time, index_right] = min(abs(right_touchdown_times - trigger_times(i_trigger)));
%             closest_heelstrike_distance_time = min([distance_to_trigger_left_time distance_to_trigger_right_time]);
% 
%             if distance_to_trigger_left_time < distance_to_trigger_right_time
%                 condition_stance_foot_list{i_trigger, 1} = 'STANCE_LEFT';
%                 condition_index_list{i_trigger, 1} = 'ONE';
%                 stretch_start_times(i_trigger, 1) = trigger_times(i_trigger);
%                 stretch_end_time_index = find(right_touchdown_times > trigger_times(i_trigger), 1, 'first');
%                 stretch_pushoff_time_index = find(right_pushoff_times > trigger_times(i_trigger), 1, 'first');
%                 if isempty(stretch_end_time_index) || isempty(stretch_pushoff_time_index)
%                     stretch_end_times(i_trigger, 1) = -1;
%                     stretch_pushoff_times(i_trigger, 1) = -1;
%                     removal_flags(i_trigger) = 1;
%                 else
%                     stretch_end_times(i_trigger, 1) = right_touchdown_times(stretch_end_time_index);
%                     stretch_pushoff_times(i_trigger, 1) = right_pushoff_times(stretch_pushoff_time_index);
%                 end
%             end    
%             if distance_to_trigger_right_time < distance_to_trigger_left_time
%                 condition_stance_foot_list{i_trigger, 1} = 'STANCE_RIGHT';
%                 condition_index_list{i_trigger, 1} = 'TWO';
%                 stretch_start_times(i_trigger, 1) = trigger_times(i_trigger);
%                 stretch_end_time_index = find(left_touchdown_times > trigger_times(i_trigger), 1, 'first');
%                 stretch_pushoff_time_index = find(left_pushoff_times > trigger_times(i_trigger), 1, 'first');
%                 if isempty(stretch_end_time_index) || isempty(stretch_pushoff_time_index)
%                     stretch_end_times(i_trigger, 1) = -1;
%                     stretch_pushoff_times(i_trigger, 1) = -1;
%                     removal_flags(i_trigger) = 1;
%                 else
%                     stretch_end_times(i_trigger, 1) = left_touchdown_times(stretch_end_time_index);
%                     stretch_pushoff_times(i_trigger, 1) = left_pushoff_times(stretch_pushoff_time_index);
%                 end
%             end
% 
%         end
% 
%         % remove flagged triggers
%         unflagged_indices = ~removal_flags;
%         trigger_times = trigger_times(unflagged_indices);
%         stretch_start_times = stretch_start_times(unflagged_indices, :);
%         stretch_end_times = stretch_end_times(unflagged_indices, :);
%         stretch_pushoff_times = stretch_pushoff_times(unflagged_indices, :);
%         condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
%         perturbation_list = perturbation_list(unflagged_indices, :);
%         condition_delay_list = condition_delay_list(unflagged_indices, :);
%         condition_index_list = condition_index_list(unflagged_indices, :);
%         condition_experimental_list = condition_experimental_list(unflagged_indices, :);
%         condition_stimulus_list = condition_stimulus_list(unflagged_indices, :);
%         condition_day_list = condition_day_list(unflagged_indices, :);
%         closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :); 
% 
%         if visualize
%             for i_trigger = 1 : length(stretch_start_times)
%                 if strcmp(condition_stance_foot_list(i_trigger), 'STANCE_RIGHT')
%                     stretch_indicator_height = 0.01;
%                 else
%                     stretch_indicator_height = -0.01;
%                 end
% 
%                 plot([stretch_start_times(i_trigger) stretch_end_times(i_trigger)], [1 1]*stretch_indicator_height, 'linewidth', 3);
%             end
%         end
% 
%         stretch_times = [stretch_start_times stretch_end_times];
% 
%         conditions_trial = struct;
%         conditions_trial.condition_stance_foot_list = condition_stance_foot_list;
%         conditions_trial.condition_perturbation_list = perturbation_list;
%         conditions_trial.condition_delay_list = condition_delay_list;
%         conditions_trial.condition_index_list = condition_index_list;
%         conditions_trial.condition_experimental_list = condition_experimental_list;
%         conditions_trial.condition_stimulus_list = condition_stimulus_list;
%         conditions_trial.condition_day_list = condition_day_list;
% 
%         event_variables_to_save.stretch_start_times = stretch_start_times;
%         event_variables_to_save.stretch_pushoff_times = stretch_pushoff_times;
%         event_variables_to_save.stretch_end_times = stretch_end_times;
%         event_variables_to_save.stretch_times = stretch_times;               
%     end  
% 
%     if strcmp(condition_stimulus, 'ARMSENSE')
% 
%         % sort out type-day-combination
%         if strcmp(trial_data.condition_experimental, 'determine_from_type_day_combination')
% 
%             if strcmp(condition_day, 'day1') && strcmp(this_trial_type, 'preOG')
%                 trial_data.condition_experimental = 'pre';
%             end
%             if strcmp(condition_day, 'day1') && strcmp(this_trial_type, 'postOG')
%                 trial_data.condition_experimental = 'post0';
%             end
%             if strcmp(condition_day, 'day2') && strcmp(this_trial_type, 'preOG')
%                 trial_data.condition_experimental = 'post4';
%             end
%         end
% 
%         % determine stride identification type
%         stride_identification = study_settings.get('stride_identification');
%         if strcmp(stride_identification, 'step_4-step_5')
%            if length(left_touchdown_times) < 3 | length(right_touchdown_times) < 3
%                this_stretch_start = 0;
%                this_stretch_end = 0;
%            else                   
%                if left_touchdown_times(1) <= right_touchdown_times(1)
%                     this_stretch_start = right_touchdown_times(2);
%                     this_stretch_end = right_touchdown_times(3);
%                     band_delimiter = min(left_touchdown_times(left_touchdown_times>this_stretch_start));
%                     first_stance_foot = 'STANCE_RIGHT';
%                     second_stance_foot = 'STANCE_LEFT';
%                else
%                     this_stretch_start = left_touchdown_times(2);
%                     this_stretch_end = left_touchdown_times(3);
%                     band_delimiter = min(right_touchdown_times(right_touchdown_times>this_stretch_start));
%                     first_stance_foot = 'STANCE_LEFT';
%                     second_stance_foot = 'STANCE_RIGHT';
%                end
%            end
%             % go through events and take stretches
%             stretch_times = [];
%             stance_foot_data = {};
%             condition_experimental_list = {};
%             condition_startfoot_list = {};
%             number_of_stretches = size(stretch_times, 1);
% 
%             % check if we have a valid stretch
%             if ~isempty(band_delimiter) && band_delimiter < this_stretch_end
%                 stretch_times = [this_stretch_start, band_delimiter, this_stretch_end];
%                 stance_foot_data = {first_stance_foot, second_stance_foot};
%                 condition_experimental_list = trial_data.condition_experimental;
%                 condition_startfoot_list = {first_stance_foot; second_stance_foot};
%                 bands_per_stretch = 2;
% 
%                 if most_affected == 'L'
%                     condition_affectedSide_list = 'L';
%                     % assign affected_startfoot
%                     % change to first_stance_foot..
%                     if strcmp(first_stance_foot, 'STANCE_LEFT')
%                         condition_affected_stancefoot_list = 'STANCE_AFFECTED';
%                     elseif strcmp(first_stance_foot, 'STANCE_RIGHT')
%                         condition_affected_stancefoot_list = 'STANCE_UNAFFECTED';
%                     end
%                 elseif most_affected == 'R'
%                     condition_affectedSide_list = 'R';
%                     % assign affected_startfoot
%                     if strcmp(first_stance_foot, 'STANCE_RIGHT')
%                         condition_affected_stancefoot_list = 'STANCE_AFFECTED';
%                     elseif strcmp(first_stance_foot, 'STANCE_LEFT')
%                         condition_affected_stancefoot_list = 'STANCE_UNAFFECTED';
%                     end
%                 end
%             end
% 
%             % fill in stuff
%             number_of_stretches = size(stretch_times, 1);
%             stretch_pushoff_times = zeros(size(stretch_times));
%             if isempty(stretch_times)
%                 stretch_start_times = [];
%                 stretch_end_times = [];
%             else
%                 stretch_start_times = stretch_times(:, 1);
%                 stretch_end_times = stretch_times(:, end);
%             end
%         end
%         if strcmp(stride_identification, 'all_but_first') 
%             % determine first stance foot
%             if left_touchdown_times(1) <= right_touchdown_times(1)
%                 stretch_starter_events = left_touchdown_times;
%                 band_delimiter_events = right_touchdown_times;
%                 first_stance_foot = 'STANCE_RIGHT';
%                 second_stance_foot = 'STANCE_LEFT';
%             else
%                 stretch_starter_events = right_touchdown_times;
%                 band_delimiter_events = left_touchdown_times;
%                 first_stance_foot = 'STANCE_LEFT';
%                 second_stance_foot = 'STANCE_RIGHT';
%             end                    
% 
%             % go through events and take stretches
%             stretch_times = [];
%             stance_foot_data = {};
%             condition_experimental_list = {};
%             condition_startfoot_list = {};
%             condition_affectedSide_list = cell(number_of_stretches, 1);
%             condition_affected_stancefoot_list = cell(number_of_stretches, 1);
% 
%             number_of_stretches = length(stretch_starter_events) - 1;
% 
%             for i_stretch = 1 : number_of_stretches % assign i_stretch to # of stretches
%                 this_stretch_start = stretch_starter_events(i_event);
%                 this_stretch_end = stretch_starter_events(i_event+1);
%                 band_delimiter = min(band_delimiter_events(band_delimiter_events>this_stretch_start));
%                 % check if we have a valid stretch
%                 if ~isempty(band_delimiter) && band_delimiter < this_stretch_end
%                     this_stretch = [this_stretch_start band_delimiter this_stretch_end];
%                     stretch_times = [stretch_times; this_stretch];
%                     stance_foot_data = [stance_foot_data; {first_stance_foot, second_stance_foot}];
% 
%                     if most_affected == 'L'
%                         condition_affectedSide_list{i_stretch} = 'L';
%                         % assign affected_startfoot
%                         % change to first_stance_foot..
%                         if strcmp(first_stance_foot, 'STANCE_LEFT')
%                             condition_affected_stancefoot_list{i_stretch} = 'STANCE_AFFECTED';
%                         elseif strcmp(first_stance_foot, 'STANCE_RIGHT')
%                             condition_affected_stancefoot_list{i_stretch} = 'STANCE_UNAFFECTED';
%                         end
%                     elseif most_affected == 'R'
%                         condition_affectedSide_list{i_stretch} = 'R';
%                         % assign affected_startfoot
%                         if strcmp(first_stance_foot, 'STANCE_RIGHT')
%                             condition_affected_stancefoot_list{i_stretch} = 'STANCE_AFFECTED';
%                         elseif strcmp(first_stance_foot, 'STANCE_LEFT')
%                             condition_affected_stancefoot_list{i_stretch} = 'STANCE_UNAFFECTED';
%                         end
%                     end
% 
%                     condition_startfoot_list = [condition_startfoot_list; first_stance_foot];
%                     condition_experimental_list = [condition_experimental_list; trial_data.condition_experimental];
%                 end
% 
%             end
% 
%             if strcmp(this_trial_type(end-1:end), 'OG')
%                 % remove all but first two stretches
%                 stretch_times(3:end, :) = [];
%             end
% 
%             % fill in stuff
%             number_of_stretches = size(stretch_times, 1);
%             stretch_pushoff_times = zeros(size(stretch_times));
%             if isempty(stretch_times)
%                 stretch_start_times = [];
%                 stretch_end_times = [];
%             else
%                 stretch_start_times = stretch_times(:, 1);
%                 stretch_end_times = stretch_times(:, end);
%             end
% 
%             bands_per_stretch = 2;
%         end
%         % restructure for saving
%         conditions_trial = struct;
%         conditions_trial.condition_experimental_list = condition_experimental_list;
%         conditions_trial.condition_startfoot_list = condition_startfoot_list;
%         conditions_trial.condition_affected_stancefoot_list = condition_affected_stancefoot_list;
%         conditions_trial.condition_affectedSide_list = condition_affectedSide_list;
%         event_variables_to_save.stretch_pushoff_times = stretch_pushoff_times;
%         event_variables_to_save.stretch_times = stretch_times;
% 
%         event_variables_to_save.stance_foot_data = stance_foot_data;
% 
%     end


end