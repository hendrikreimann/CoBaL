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

function [conditions_trial, event_variables_to_save] = determineConditionLevels_intermittentPerturbations(study_settings, trial_data)
    experimental_paradigm = study_settings.get('experimental_paradigm');

    bands_per_stretch = study_settings.get('number_of_steps_to_analyze');
    time_to_nearest_heelstrike_before_trigger_threshold = 0.10; % a heelstrike should happen less than this long before a trigger
    time_to_nearest_heelstrike_after_trigger_threshold = 0.3; % a heelstrike should happen less than this long after a trigger

    conditions_trial = struct;
    number_of_triggers = length(trial_data.trigger_indices_mocap);
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
        if trial_data.trigger_indices_labview(i_trigger) == length(trial_data.illusion_trajectory)
            removal_flags(i_trigger) = 1;
        else
            if trial_data.illusion_trajectory(trial_data.trigger_indices_labview(i_trigger)+1) > 0
                stimulus_list{i_trigger} = 'STIM_RIGHT';
            end
            if trial_data.illusion_trajectory(trial_data.trigger_indices_labview(i_trigger)+1) < 0
                stimulus_list{i_trigger} = 'STIM_LEFT';
            end
            if trial_data.illusion_trajectory(trial_data.trigger_indices_labview(i_trigger)+1) == 0
                stimulus_list{i_trigger} = 'STIM_NONE';
            end
        end
        % determine trigger foot

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
            if length(trial_data.left_touchdown_times) < index_left + 2 || removal_flags(i_trigger) == 1
                % data doesn't include the required number of steps after the trigger
                removal_flags(i_trigger) = 1;
                right_foot_heelstrike_0 = NaN;
                right_foot_heelstrike_plus_1 = NaN;
                right_foot_heelstrike_plus_2 = NaN;
                left_foot_heelstrike_0 = NaN;
                left_foot_heelstrike_plus_1 = NaN;
                left_foot_heelstrike_plus_2 = NaN;

            else
                left_foot_heelstrike_0  = trial_data.left_touchdown_times(index_left);
                left_foot_heelstrike_1  = trial_data.left_touchdown_times(index_left+1);
                left_foot_heelstrike_2  = trial_data.left_touchdown_times(index_left+2);
                left_foot_pushoff_0     = min(trial_data.left_pushoff_times(trial_data.left_pushoff_times >= left_foot_heelstrike_0));
                left_foot_pushoff_1     = min(trial_data.left_pushoff_times(trial_data.left_pushoff_times >= left_foot_heelstrike_1));
                left_foot_pushoff_2     = min(trial_data.left_pushoff_times(trial_data.left_pushoff_times >= left_foot_heelstrike_2));

                right_foot_heelstrike_0 = min(trial_data.right_touchdown_times(trial_data.right_touchdown_times >= left_foot_heelstrike_0));
                right_foot_heelstrike_1 = min(trial_data.right_touchdown_times(trial_data.right_touchdown_times >= left_foot_heelstrike_1));
                right_foot_heelstrike_2 = min(trial_data.right_touchdown_times(trial_data.right_touchdown_times >= left_foot_heelstrike_2));
                right_foot_pushoff_0    = max(trial_data.right_pushoff_times(trial_data.right_pushoff_times <= left_foot_pushoff_0));
                right_foot_pushoff_1    = max(trial_data.right_pushoff_times(trial_data.right_pushoff_times <= left_foot_pushoff_1));
                right_foot_pushoff_2    = max(trial_data.right_pushoff_times(trial_data.right_pushoff_times <= left_foot_pushoff_2));

                % notify if events are not sorted properly
                event_order = [ ...
                                left_foot_heelstrike_0 right_foot_pushoff_0 right_foot_heelstrike_0 left_foot_pushoff_0 ...
                                left_foot_heelstrike_1 right_foot_pushoff_1 right_foot_heelstrike_1 left_foot_pushoff_1 ...
                                left_foot_heelstrike_2  ...
                              ];
                if ~issorted(event_order)
                    disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(time_stimulus(trial_data.trigger_indices_labview(i_trigger)))]);
                end
                % check check
%                     if visualize
%                         plot([left_foot_heelstrike_minus_1 left_foot_heelstrike_0 left_foot_heelstrike_1 left_foot_heelstrike_2], [0 0 0 0]-0.01, 'v', 'linewidth', 3);
%                         plot([left_foot_pushoff_minus_1  left_foot_pushoff_0 left_foot_pushoff_1 left_foot_pushoff_2], [0 0 0 0]-0.01, '^', 'linewidth', 3);
%                         plot([right_foot_heelstrike_minus_1 right_foot_heelstrike_0 right_foot_heelstrike_1 right_foot_heelstrike_2], [0 0 0 0]+0.01, 'v', 'linewidth', 3);
%                         plot([right_foot_pushoff_minus_1  right_foot_pushoff_0 right_foot_pushoff_1 right_foot_pushoff_2], [0 0 0 0]+0.01, '^', 'linewidth', 3);
%                         % note: this can crash if one of thse events is empty, because we are plotting before we
%                         % have checked that
%                     end                
            end
        elseif strcmp(trigger_foot, 'right')
            if length(trial_data.right_touchdown_times) < index_right + 2 || removal_flags(i_trigger) == 1
                % data doesn't include the required number of steps after the trigger
                removal_flags(i_trigger) = 1;
                left_foot_heelstrike_0 = NaN;
                left_foot_heelstrike_1 = NaN;
                left_foot_heelstrike_2 = NaN;
                right_foot_heelstrike_0 = NaN;
                right_foot_heelstrike_1 = NaN;
                right_foot_heelstrike_2 = NaN;

            else
                right_foot_heelstrike_0 = trial_data.right_touchdown_times(index_right);
                right_foot_heelstrike_1 = trial_data.right_touchdown_times(index_right+1);
                right_foot_heelstrike_2 = trial_data.right_touchdown_times(index_right+2);

                right_foot_pushoff_0    = min(trial_data.right_pushoff_times(trial_data.right_pushoff_times >= right_foot_heelstrike_0));
                right_foot_pushoff_1    = min(trial_data.right_pushoff_times(trial_data.right_pushoff_times >= right_foot_heelstrike_1));
                right_foot_pushoff_2    = min(trial_data.right_pushoff_times(trial_data.right_pushoff_times >= right_foot_heelstrike_2));

                left_foot_heelstrike_0  = min(trial_data.left_touchdown_times(trial_data.left_touchdown_times >= right_foot_heelstrike_0));
                left_foot_heelstrike_1  = min(trial_data.left_touchdown_times(trial_data.left_touchdown_times >= right_foot_heelstrike_1));
                left_foot_heelstrike_2  = min(trial_data.left_touchdown_times(trial_data.left_touchdown_times >= right_foot_heelstrike_2));
                left_foot_pushoff_0     = max(trial_data.left_pushoff_times(trial_data.left_pushoff_times <= right_foot_pushoff_0));
                left_foot_pushoff_1     = max(trial_data.left_pushoff_times(trial_data.left_pushoff_times <= right_foot_pushoff_1));
                left_foot_pushoff_2     = max(trial_data.left_pushoff_times(trial_data.left_pushoff_times <= right_foot_pushoff_2));

                % notify if events are not sorted properly
                event_order = [ ...
                                right_foot_heelstrike_0 left_foot_pushoff_0 left_foot_heelstrike_0 right_foot_pushoff_0 ...
                                right_foot_heelstrike_1 left_foot_pushoff_1 left_foot_heelstrike_1 right_foot_pushoff_1 ...
                                right_foot_heelstrike_2 ...
                              ];
                if ~issorted(event_order)
                    disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(time_stimulus(trial_data.trigger_indices_labview(i_trigger)))]);
                end

%                     if visualize
%                         plot([left_foot_heelstrike_0 left_foot_heelstrike_1 left_foot_heelstrike_2], [0 0 0]-0.01, 'v', 'linewidth', 2, 'color', 'r');
%                         plot([left_foot_pushoff_0 left_foot_pushoff_1 left_foot_pushoff_2], [0 0 0]-0.01, '^', 'linewidth', 2, 'color', 'g');
%                         plot([right_foot_heelstrike_0 right_foot_heelstrike_1 right_foot_heelstrike_2], [0 0 0]+0.01, 'v', 'linewidth', 2, 'color', 'r');
%                         plot([right_foot_pushoff_0 right_foot_pushoff_1 right_foot_pushoff_2], [0 0 0]+0.01, '^', 'linewidth', 2, 'color', 'g');
%                     end            
            end            
        else
            trigger_foot = 'unclear';
            disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_stimulus(trial_data.trigger_indices_labview(i_trigger))) ' - triggering heelstrike unclear']);
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
%                 if visualize
%                     plot([stretch_times(i_trigger, 1) stretch_times(i_trigger, 2)], [-1 1]*-0.01, 'color', 'r', 'linewidth', 3);
%                     plot([stretch_times(i_trigger, 2) stretch_times(i_trigger, 3)], [-1 1]*+0.01, 'color', 'g', 'linewidth', 3);
%                     plot([stretch_times(i_trigger, 3) stretch_times(i_trigger, 4)], [-1 1]*-0.01, 'color', 'b', 'linewidth', 3);
%                     plot([stretch_times(i_trigger, 4) stretch_times(i_trigger, 5)], [-1 1]*+0.01, 'color', 'm', 'linewidth', 3);
%                 end
        end

        % determine amplitude of this stimulus
        if strcmp(experimental_paradigm, 'Vision')
            amplitude = trial_data.current_acceleration_trajectory(trial_data.trigger_indices_labview(i_trigger));
            amplitude_list{i_trigger} = num2str(amplitude);
        end
    end

    % remove flagged triggers
    unflagged_indices = ~removal_flags;
    trial_data.trigger_times = trial_data.trigger_times(unflagged_indices);
    trial_data.trigger_indices_labview = trial_data.trigger_indices_labview(unflagged_indices, :);
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

        for i_stretch = 1:length(trial_data.trigger_indices_labview)
            VR_trigger_position = scene_translation_mod100(trial_data.trigger_indices_labview(i_stretch));
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