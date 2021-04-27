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

function [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_gvsOverground(study_settings, subject_settings, trial_data)

    % get parameters from settings
    experimental_paradigm = study_settings.get('experimental_paradigm');
    time_to_nearest_heelstrike_before_trigger_threshold = study_settings.get('time_to_nearest_heelstrike_before_trigger_threshold');
    time_to_nearest_heelstrike_after_trigger_threshold = study_settings.get('time_to_nearest_heelstrike_after_trigger_threshold');
    
    % allocate output variables
    number_of_triggers = length(trial_data.trigger_indices_mocap);
    conditions_trial = struct;
    event_variables_to_save = struct;
    removal_flags = false(number_of_triggers, 1);
    
    % determine start and end
    stance_foot_data = {'STANCE_LEFT', 'STANCE_RIGHT', 'STANCE_LEFT', 'STANCE_RIGHT'};
    trigger_foot_list = cell(number_of_triggers, 1); % triggering foot TRIGGER_LEFT or TRIGGER_RIGHT
    bands_per_stretch = length(stance_foot_data);

    % find heelstrike that corresponds to the trigger
    if ~isempty(trial_data.trigger_times) && length(trial_data.left_touchdown_times) >= 3 && length(trial_data.right_touchdown_times) >= 2
        trigger_time = trial_data.trigger_times(1);

        % get closest heelstrike on either side
        [~, index_left] = min(abs(trial_data.left_touchdown_times - trigger_time));
        [~, index_right] = min(abs(trial_data.right_touchdown_times - trigger_time));

        % is the closest left heelstrike within the acceptable interval?
        closest_left_heelstrike = trial_data.left_touchdown_times(index_left);
        time_difference_left = closest_left_heelstrike - trigger_time; % where does the closest left heelstrike lie relative to the trigger?
        if -time_to_nearest_heelstrike_before_trigger_threshold < time_difference_left && time_difference_left < time_to_nearest_heelstrike_after_trigger_threshold
            % left heelstrike is acceptable
            left_heelstrike_acceptable = true;
        else
            left_heelstrike_acceptable = false;
        end

        % is the closest right heelstrike within the acceptable interval?
        closest_right_heelstrike = trial_data.right_touchdown_times(index_right);
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
        left_foot_heelstrike_0  = trial_data.left_touchdown_times(index_left);
        left_foot_heelstrike_1  = trial_data.left_touchdown_times(index_left+1);
        left_foot_heelstrike_2  = trial_data.left_touchdown_times(index_left+2);
        left_foot_pushoff_0     = min(trial_data.left_pushoff_times(trial_data.left_pushoff_times >= left_foot_heelstrike_0));
        left_foot_pushoff_1     = min(trial_data.left_pushoff_times(trial_data.left_pushoff_times >= left_foot_heelstrike_1));
%                    left_foot_pushoff_2     = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_2));

        right_foot_heelstrike_0 = min(trial_data.right_touchdown_times(trial_data.right_touchdown_times >= left_foot_heelstrike_0));
        right_foot_heelstrike_1 = min(trial_data.right_touchdown_times(trial_data.right_touchdown_times >= left_foot_heelstrike_1));
%                     right_foot_heelstrike_2 = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_2));
        right_foot_pushoff_0    = max(trial_data.right_pushoff_times(trial_data.right_pushoff_times <= left_foot_pushoff_0));
        right_foot_pushoff_1    = max(trial_data.right_pushoff_times(trial_data.right_pushoff_times <= left_foot_pushoff_1));
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
%         if visualize
%             hold on
%             plot([left_foot_heelstrike_0 left_foot_heelstrike_1 left_foot_heelstrike_2], [0 0 0]-0.01, 'v', 'linewidth', 3);
%             plot([ left_foot_pushoff_0 left_foot_pushoff_1 left_foot_pushoff_2], [0 0 0]-0.01, '^', 'linewidth', 3);
%             plot([right_foot_heelstrike_0 right_foot_heelstrike_1 right_foot_heelstrike_2], [0 0 0]+0.01, 'v', 'linewidth', 3);
%             plot([ right_foot_pushoff_0 right_foot_pushoff_1 right_foot_pushoff_2], [0 0 0]+0.01, '^', 'linewidth', 3);
%             % note: this can crash if one of thse events is empty, because we are plotting before we
%             % have checked that
%         end                  

        stretch_times = [left_foot_heelstrike_0 right_foot_heelstrike_0 left_foot_heelstrike_1 right_foot_heelstrike_1 left_foot_heelstrike_2];


        % determine stimulus
        check_time_delay = 0.1;
        check_time = stretch_times(1) + check_time_delay;
        [~, check_index] = min(abs(trial_data.time_analog - check_time));
        if trial_data.illusion_trajectory(check_index) > 0
            stimulus_list = {'STIM_RIGHT'};
        elseif trial_data.illusion_trajectory(check_index) < 0
            stimulus_list = {'STIM_LEFT'};
        elseif trial_data.illusion_trajectory(check_index) == 0
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

