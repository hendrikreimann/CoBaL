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

function [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_intermittentPerturbations_flex(study_settings, trial_data)

    % get parameters from settings
    experimental_paradigm = study_settings.get('experimental_paradigm');
    bands_per_stretch = study_settings.get('number_of_steps_to_analyze');
    
    time_to_nearest_heelstrike_before_trigger_threshold = study_settings.get('time_to_nearest_heelstrike_before_trigger_threshold', 1);
    time_to_nearest_heelstrike_after_trigger_threshold = study_settings.get('time_to_nearest_heelstrike_after_trigger_threshold', 1);
    
%     time_to_nearest_heelstrike_before_trigger_threshold = 0.10; % a heelstrike should happen less than this long before a trigger
%     time_to_nearest_heelstrike_after_trigger_threshold = 0.3; % a heelstrike should happen less than this long after a trigger

    % allocate output variables
    number_of_triggers = length(trial_data.trigger_indices_mocap);
    conditions_trial = struct;
    event_variables_to_save = struct;
    removal_flags = zeros(number_of_triggers, 1);
    

    stretch_times = zeros(number_of_triggers, bands_per_stretch+1);
    closest_heelstrike_distance_times = zeros(number_of_triggers, 1);
    stance_foot_data = cell(number_of_triggers, bands_per_stretch);
    stimulus_list = cell(number_of_triggers, 1); % stimulus STIM_LEFT, STIM_RIGHT or STIM_NONE
    amplitude_list = cell(number_of_triggers, 1); % amplitude of the visual stim, e.g. 30, 60, 120 deg/sec^2
    trigger_foot_list = cell(number_of_triggers, 1); % triggering foot TRIGGER_LEFT or TRIGGER_RIGHT

    for i_trigger = 1 : number_of_triggers
        % determine stimulus
%         if trial_data.trigger_indices_stimulus(i_trigger) == length(trial_data.illusion_trajectory)
%             removal_flags(i_trigger) = 1;
%         else
%             if trial_data.illusion_trajectory(trial_data.trigger_indices_stimulus(i_trigger)+1) > 0
%                 stimulus_list{i_trigger} = 'STIM_RIGHT';
%             end
%             if trial_data.illusion_trajectory(trial_data.trigger_indices_stimulus(i_trigger)+1) < 0
%                 stimulus_list{i_trigger} = 'STIM_LEFT';
%             end
%             if trial_data.illusion_trajectory(trial_data.trigger_indices_stimulus(i_trigger)+1) == 0
%                 stimulus_list{i_trigger} = 'STIM_NONE';
%             end
%         end
        
        stimulus_list{i_trigger} = determineStimulus(i_trigger);
        
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
        
        % I have successfully determined the trigger foot, now let's extract trigger foot and contra foot events
        if strcmp(trigger_foot, 'left')
            trigger_foot_touchdown_times = trial_data.left_touchdown_times;
            contra_foot_touchdown_times = trial_data.right_touchdown_times;
            stretch_times(i_trigger, 1) = trigger_foot_touchdown_times(index_left);
            trigger_foot_stance_label = 'STANCE_LEFT';
            contra_foot_stance_label = 'STANCE_RIGHT';
            trigger_foot_list{i_trigger} = 'TRIGGER_LEFT';
        end
        if strcmp(trigger_foot, 'right')
            trigger_foot_touchdown_times = trial_data.right_touchdown_times;
            contra_foot_touchdown_times = trial_data.left_touchdown_times;
            stretch_times(i_trigger, 1) = trigger_foot_touchdown_times(index_right);
            trigger_foot_stance_label = 'STANCE_RIGHT';
            contra_foot_stance_label = 'STANCE_LEFT';
            trigger_foot_list{i_trigger} = 'TRIGGER_RIGHT';
        end
        
        % now fill the other times
        for i_band = 1 : bands_per_stretch
            band_start_time = stretch_times(i_trigger, i_band);
            % find push-off time within the band
            if mod(i_band, 2) == 0
                swing_foot_touchdown_times = trigger_foot_touchdown_times;
                stance_foot_data{i_trigger, i_band} = contra_foot_stance_label;
            end
            if mod(i_band, 2) == 1
                swing_foot_touchdown_times = contra_foot_touchdown_times;
                stance_foot_data{i_trigger, i_band} = trigger_foot_stance_label;
            end
            this_band_end_time = min(swing_foot_touchdown_times(swing_foot_touchdown_times > band_start_time));
            stretch_times(i_trigger, i_band+1) = this_band_end_time;
        end

        % determine amplitude of this stimulus
        if strcmp(experimental_paradigm, 'Vision')
            amplitude = trial_data.current_acceleration_trajectory(trial_data.trigger_indices_stimulus(i_trigger));
            amplitude_list{i_trigger} = num2str(amplitude);
        end
    end

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

        for i_stretch = 1:length(trial_data.trigger_indices_stimulus)
            VR_trigger_position = scene_translation_mod100(trial_data.trigger_indices_stimulus(i_stretch));
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

%     if exist('affected_side')
%         for i_stretch = 1 : length(trigger_foot_list)
%             if strcmp(affected_side, 'Left')
%                 if strcmp(trigger_foot_list(i_stretch), 'TRIGGER_LEFT')
%                     condition_affected_stancefoot_list{i_stretch} = 'TRIGGER_AFFECTED';
%                 elseif strcmp(trigger_foot_list(i_stretch), 'TRIGGER_RIGHT')
%                     condition_affected_stancefoot_list{i_stretch} = 'TRIGGER_UNAFFECTED';
%                 end
%             elseif strcmp(affected_side, 'Right')
%                 if strcmp(trigger_foot_list(i_stretch), 'TRIGGER_RIGHT')
%                     condition_affected_stancefoot_list{i_stretch} = 'TRIGGER_AFFECTED';
%                 elseif strcmp(trigger_foot_list(i_stretch), 'TRIGGER_LEFT')
%                     condition_affected_stancefoot_list{i_stretch} = 'TRIGGER_UNAFFECTED';
%                 end
%             end
%         end
%         conditions_trial.affected_stancefoot_list = condition_affected_stancefoot_list';
%     end

    % restructure for saving
    conditions_trial.stimulus_list = stimulus_list;
    conditions_trial.amplitude_list = amplitude_list;
    conditions_trial.trigger_foot_list = trigger_foot_list;
    conditions_trial.direction_list = direction_list;
    if strcmp(experimental_paradigm, 'OculusLaneRestriction')
        conditions_trial.zone_side_list = zone_side_list;
        conditions_trial.zone_direction_list = zone_direction_list;
    end
    event_variables_to_save.stretch_times = stretch_times;
    event_variables_to_save.stance_foot_data = stance_foot_data;
    
function stimulus_label = determineStimulus(trigger_index)
    % we know that the stimulus started on the time step given by trigger_indices_stimulus(trigger_index)
    % so the next time step has the relevant information
    if trial_data.trigger_indices_stimulus(trigger_index) == length(trial_data.illusion_trajectory)
        % for some reason the trigger is the last time step - that can't be right, so flag this trigger for removal
        removal_flags(trigger_index) = 1;
        stimulus_label = 'N/A';
    else
        first_stimulus_index = trial_data.trigger_indices_stimulus(trigger_index) + 1;
        if trial_data.illusion_trajectory(first_stimulus_index) > 0
            stimulus_label = 'STIM_RIGHT';
        end
        if trial_data.illusion_trajectory(first_stimulus_index) < 0
            stimulus_label = 'STIM_LEFT';
        end
        if trial_data.illusion_trajectory(first_stimulus_index) == 0
            stimulus_label = 'STIM_NONE';
        end
    end
    
    
%     if trial_data.trigger_indices_stimulus(trigger_index) == length(trial_data.illusion_trajectory)
%         removal_flags(trigger_index) = 1;
%     else
%         if trial_data.illusion_trajectory(trial_data.trigger_indices_stimulus(trigger_index)+1) > 0
%             stimulus_label = 'STIM_RIGHT';
%         end
%         if trial_data.illusion_trajectory(trial_data.trigger_indices_stimulus(trigger_index)+1) < 0
%             stimulus_label = 'STIM_LEFT';
%         end
%         if trial_data.illusion_trajectory(trial_data.trigger_indices_stimulus(trigger_index)+1) == 0
%             stimulus_label = 'STIM_NONE';
%         end
%     end

end    
    
end

