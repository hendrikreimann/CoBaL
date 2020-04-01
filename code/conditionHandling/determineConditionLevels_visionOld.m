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

function [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_visionOld(study_settings, subject_settings, trial_data)

    % get parameters from settings
    experimental_paradigm = study_settings.get('experimental_paradigm');
    time_to_nearest_heelstrike_before_trigger_threshold = study_settings.get('time_to_nearest_heelstrike_before_trigger_threshold');
    time_to_nearest_heelstrike_after_trigger_threshold = study_settings.get('time_to_nearest_heelstrike_after_trigger_threshold');
    
    % allocate output variables
    number_of_triggers = length(trial_data.trigger_indices_mocap);
    conditions_trial = struct;
    event_variables_to_save = struct;
    removal_flags = false(number_of_triggers, 1);
    

    
    bands_per_stretch = 2;

    number_of_triggers = length(trial_data.trigger_indices_mocap);
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
        if trial_data.illusion_trajectory(trial_data.stim_start_indices_stimulus(i_trigger)+1) > 0
            stimulus_list_stim{i_trigger} = 'STIM_RIGHT';
        end
        if trial_data.illusion_trajectory(trial_data.stim_start_indices_stimulus(i_trigger)+1) < 0
            stimulus_list_stim{i_trigger} = 'STIM_LEFT';
        end
        if trial_data.illusion_trajectory(trial_data.stim_start_indices_stimulus(i_trigger)+1) == 0
            stimulus_list_stim{i_trigger} = 'STIM_NONE';
        end
        stimulus_list_ctrl{i_trigger} = 'STIM_NONE';

        % get closest heelstrike on either side
        [~, index_left] = min(abs(trial_data.left_touchdown_times - trial_data.trigger_times(i_trigger)));
        [~, index_right] = min(abs(trial_data.right_touchdown_times - trial_data.trigger_times(i_trigger)));

        % is the closest left heelstrike within the acceptable interval?
        closest_left_heelstrike = trial_data.left_touchdown_times(index_left);
        time_difference_left = closest_left_heelstrike - trial_data.trigger_times(i_trigger); % where does the closest left heelstrike lie relative to the trigger?
        if -time_to_nearest_heelstrike_before_trigger_threshold < time_difference_left && time_difference_left < time_to_nearest_heelstrike_after_trigger_threshold
            % left heelstrike is acceptable
            left_heelstrike_acceptable = true;
        else
            left_heelstrike_acceptable = false;
        end

        % is the closest right heelstrike within the acceptable interval?
        closest_right_heelstrike = trial_data.right_touchdown_times(index_right);
        time_difference_right = closest_right_heelstrike - trial_data.trigger_times(i_trigger); % where does the closest right heelstrike lie relative to the trigger?
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
            if length(trial_data.left_touchdown_times) < index_left + 1 || removal_flags(i_trigger) == 1 || index_left == 1
                % data doesn't include the required number of steps after the trigger
                removal_flags(i_trigger) = 1;
                left_foot_heelstrike_0  = NaN;
                left_foot_heelstrike_1  = NaN;
                left_foot_pushoff_0     = NaN;
                right_foot_heelstrike_0 = NaN;
                right_foot_pushoff_0    = NaN;
            else
                left_foot_heelstrike_pre  = trial_data.left_touchdown_times(index_left-1);
                left_foot_heelstrike_0  = trial_data.left_touchdown_times(index_left);
                left_foot_heelstrike_1  = trial_data.left_touchdown_times(index_left+1);

                left_foot_pushoff_pre     = max(trial_data.left_pushoff_times(trial_data.left_pushoff_times < left_foot_heelstrike_0));
                left_foot_pushoff_0     = min(trial_data.left_pushoff_times(trial_data.left_pushoff_times >= left_foot_heelstrike_0));

                right_foot_heelstrike_pre = max(trial_data.right_touchdown_times(trial_data.right_touchdown_times < left_foot_heelstrike_0));
                right_foot_heelstrike_0 = min(trial_data.right_touchdown_times(trial_data.right_touchdown_times >= left_foot_heelstrike_0));
                right_foot_pushoff_pre    = max(trial_data.right_pushoff_times(trial_data.right_pushoff_times <= left_foot_heelstrike_0));
                right_foot_pushoff_0    = max(trial_data.right_pushoff_times(trial_data.right_pushoff_times <= left_foot_pushoff_0));

                % notify if events are not sorted properly
                if ~issorted ...
                      ( ...
                        [ ...
                          left_foot_heelstrike_pre right_foot_pushoff_pre right_foot_heelstrike_pre left_foot_pushoff_pre ...
                          left_foot_heelstrike_0 right_foot_pushoff_0 right_foot_heelstrike_0 left_foot_pushoff_0 ...
                          left_foot_heelstrike_1 ...
                        ] ...
                      )
                    disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(time_stimulus(trigger_indices_stimulus(i_trigger)))]);
                end
            end
        elseif strcmp(trigger_foot, 'right')
            if length(trial_data.right_touchdown_times) < index_right + 1 || removal_flags(i_trigger) == 1 || index_right == 1
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
                right_foot_heelstrike_pre = trial_data.right_touchdown_times(index_right-1);
                right_foot_heelstrike_0 = trial_data.right_touchdown_times(index_right);
                right_foot_heelstrike_1 = trial_data.right_touchdown_times(index_right+1);

                right_foot_pushoff_pre    = max(trial_data.right_pushoff_times(trial_data.right_pushoff_times < right_foot_heelstrike_0));
                right_foot_pushoff_0    = min(trial_data.right_pushoff_times(trial_data.right_pushoff_times >= right_foot_heelstrike_0));

                left_foot_heelstrike_pre  = max(trial_data.left_touchdown_times(trial_data.left_touchdown_times < right_foot_heelstrike_0));
                left_foot_heelstrike_0  = min(trial_data.left_touchdown_times(trial_data.left_touchdown_times >= right_foot_heelstrike_0));
                left_foot_pushoff_pre     = max(trial_data.left_pushoff_times(trial_data.left_pushoff_times <= right_foot_heelstrike_0));
                left_foot_pushoff_0     = max(trial_data.left_pushoff_times(trial_data.left_pushoff_times <= right_foot_pushoff_0));

                % notify if events are not sorted properly
                if ~issorted ...
                      ( ...
                        [ ...
                          right_foot_heelstrike_pre left_foot_pushoff_pre left_foot_heelstrike_pre right_foot_pushoff_pre ...
                          right_foot_heelstrike_0 left_foot_pushoff_0 left_foot_heelstrike_0 right_foot_pushoff_0 ...
                          right_foot_heelstrike_1 ...
                        ] ...
                      )
                    disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(time_stimulus(trigger_indices_stimulus(i_trigger)))]);
                end

            end            
        else
            trigger_foot = 'unclear';
            disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_stimulus(trigger_indices_stimulus(i_trigger))) ' - triggering heelstrike unclear']);
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

%     % remove flagged triggers
%     unflagged_indices = ~removal_flags;
%     trial_data.trigger_times = trial_data.trigger_times(unflagged_indices);
%     trigger_indices_stimulus = trigger_indices_stimulus(unflagged_indices, :);
%     stretch_times_stim = stretch_times_stim(unflagged_indices, :);
%     stance_foot_data_stim = stance_foot_data_stim(unflagged_indices, :);
%     stimulus_list_stim = stimulus_list_stim(unflagged_indices, :);
% %                 delay_list_stim = delay_list_stim(unflagged_indices, :);
%     trigger_foot_list_stim = trigger_foot_list_stim(unflagged_indices, :);
%     stretch_times_ctrl = stretch_times_ctrl(unflagged_indices, :);
%     stance_foot_data_ctrl = stance_foot_data_ctrl(unflagged_indices, :);
%     stimulus_list_ctrl = stimulus_list_ctrl(unflagged_indices, :);
% %                 delay_list_ctrl = delay_list_ctrl(unflagged_indices, :);
%     trigger_foot_list_ctrl = trigger_foot_list_ctrl(unflagged_indices, :);

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

