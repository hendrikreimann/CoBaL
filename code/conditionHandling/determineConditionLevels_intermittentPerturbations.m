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

function [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_intermittentPerturbations(study_settings, subject_settings, trial_data)

    % get parameters from settings
    experimental_paradigm = study_settings.get('experimental_paradigm');
    bands_per_stretch = study_settings.get('number_of_steps_to_analyze');
    
    % allocate output variables
    number_of_triggers = length(trial_data.trigger_indices_mocap);
    conditions_trial = struct;
    event_variables_to_save = struct;
    removal_flags = false(number_of_triggers, 1);
    

    stretch_times = zeros(number_of_triggers, bands_per_stretch+1);
    stance_foot_data = cell(number_of_triggers, bands_per_stretch);
    stimulus_list = cell(number_of_triggers, 1); % stimulus STIM_LEFT, STIM_RIGHT or STIM_NONE
    amplitude_list = cell(number_of_triggers, 1); % amplitude of the visual stim, e.g. 30, 60, 120 deg/sec^2
    trigger_foot_list = cell(number_of_triggers, 1); % triggering foot TRIGGER_LEFT or TRIGGER_RIGHT
    direction_list = cell(number_of_triggers, 1); % direction relative to body, STIM_TOWARDS or STIM_AWAY

    for i_trigger = 1 : number_of_triggers
        % determine stimulus
        stimulus_list{i_trigger} = determineStimulus(i_trigger);
        
        % determine trigger foot
        [trigger_foot, trigger_time] = determineTriggerFoot(i_trigger);
        if strcmp(trigger_foot, 'left')
            trigger_foot_list{i_trigger} = 'TRIGGER_LEFT';
        end        
        if strcmp(trigger_foot, 'right')
            trigger_foot_list{i_trigger} = 'TRIGGER_RIGHT';
        end
        
        % determine stretch times
        [stretch_times_this_trigger, stance_foot_data_this_trigger, removal_flag_this_trigger] = determineStretchTimes(trigger_time, trigger_foot);
        stretch_times(i_trigger, :) = stretch_times_this_trigger;
        stance_foot_data(i_trigger, :) = stance_foot_data_this_trigger;
        removal_flags(i_trigger) = removal_flag_this_trigger;

        % determine stimulus direction
        direction_list{i_trigger} = determineStimulusDirection(trigger_foot_list{i_trigger}, stimulus_list{i_trigger});
        
        % determine stimulus amplitude
        if strcmp(experimental_paradigm, 'Vision')
            amplitude = trial_data.current_acceleration_trajectory(trial_data.trigger_indices_stimulus(i_trigger));
            amplitude_list{i_trigger} = num2str(amplitude);
        end
    end

    % add location of the no-step zone for lane restriction paradigm
    if strcmp(experimental_paradigm, 'OculusLaneRestriction')
        % TODO: this part probably broke when moving the code into a function, fix later
        
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
        conditions_trial.zone_side_list = zone_side_list;
        conditions_trial.zone_direction_list = zone_direction_list;
    end

    % add cadence information for cadence paradigm
    if strcmp(experimental_paradigm, 'CadenceVision') || strcmp(experimental_paradigm, 'CadenceGVS')
        % TODO: this part probably broke when moving the code into a function, fix when necessary
        
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
    
    % add fatigue information for fatigue paradigm
    if strcmp(experimental_paradigm, 'FatigueGVS')
        % TODO: this part probably broke when moving the code into a function, fix when necessary
        fatigue_list = cell(size(direction_list));
        if ismember(i_trial, fatigue_trials)
            [fatigue_list{:}] = deal('FATIGUED');
        else
            [fatigue_list{:}] = deal('UNFATIGUED');
        end
        conditions_trial.fatigue_list = fatigue_list;
    end
    
    % add cognitive load information for cognitive load paradigm
    if strcmp(experimental_paradigm, 'CognitiveLoadVision') || strcmp(experimental_paradigm, 'CognitiveLoadGvs')
        % TODO: this part probably broke when moving the code into a function, fix when necessary
        cognitive_load_list = cell(size(direction_list));
        if ismember(i_trial, back_7_trials)
            [cognitive_load_list{:}] = deal('BACK_7');
        else
            [cognitive_load_list{:}] = deal('NO_LOAD');
        end
        conditions_trial.cognitive_load_list = cognitive_load_list;
    end

    if strcmp(experimental_paradigm, 'SR_VisualStim')
        % get stimulus strength
        stimulus_strength = determineStochasticResonanceStrength(subject_settings, trial_data.trial_type, trial_data.trial_number);
        stim_amplitude_list = repmat({stimulus_strength}, number_of_triggers, 1);
        conditions_trial.stim_amplitude_list = stim_amplitude_list;
    end
    
    % check if 'group' is listed as a condition and extract it from the subject settings if needed
    conditions_table = study_settings.get('conditions');
    if any(strcmp(conditions_table(:, 1), 'group'))
        group = subject_settings.get('group');
        condition_group_list = cell(number_of_triggers, 1);
        for i_stretch = 1 : number_of_triggers
            condition_group_list{i_stretch} = group;
        end
        conditions_trial.group_list = condition_group_list;
    end
    
%     %Ash %add affected side
    conditions_table = study_settings.get('conditions');
    if any(strcmp(conditions_table(:, 1), 'affected_side'))
        affected_side = subject_settings.get('affected_side');
        condition_affected_side_list = cell(number_of_triggers, 1);
        for i_stretch = 1 : length(trigger_foot_list)
            if strcmp(affected_side, 'Left')
                if strcmp(trigger_foot_list(i_stretch), 'TRIGGER_LEFT')
                    condition_affected_side_list{i_stretch} = 'TRIGGER_AFFECTED';
                elseif strcmp(trigger_foot_list(i_stretch), 'TRIGGER_RIGHT')
                    condition_affected_side_list{i_stretch} = 'TRIGGER_UNAFFECTED';
                end
            elseif strcmp(affected_side, 'Right')
                if strcmp(trigger_foot_list(i_stretch), 'TRIGGER_RIGHT')
                    condition_affected_side_list{i_stretch} = 'TRIGGER_AFFECTED';
                elseif strcmp(trigger_foot_list(i_stretch), 'TRIGGER_LEFT')
                    condition_affected_side_list{i_stretch} = 'TRIGGER_UNAFFECTED';
                end
            end
        end
        conditions_trial.affected_side_list = condition_affected_side_list;
    end
    
%Ash %           conditions_table = study_settings.get('conditions');
%     if any(strcmp(conditions_table(:, 1), 'affected_side'))
%         affected_side = subject_settings.get('affected_side');
%         condition_affected_side_list = cell(number_of_triggers, 1);
%         for i_stretch = 1 : number_of_triggers
%             condition_affected_side_list{i_stretch} = affected_side;
%         end
%         conditions_trial.affected_side_list = condition_affected_side_list;
%     end
%     
    
    % add information about trigger relative to more affected side 
    % HR: I commented this out since it is legacy code and in a messy format 
    %     if it is still needed, e.g. for the Parkinson's or CP study, I'll look at getting it back in
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
end    
function [trigger_foot, trigger_time] = determineTriggerFoot(trigger_index)
        % get thresholds for what time difference is acceptable
        time_to_nearest_heelstrike_before_trigger_threshold = study_settings.get('time_to_nearest_heelstrike_before_trigger_threshold', 1);
        time_to_nearest_heelstrike_after_trigger_threshold = study_settings.get('time_to_nearest_heelstrike_after_trigger_threshold', 1);
    
        % get closest heelstrike on either side
        [~, index_candidate_left] = min(abs(trial_data.left_touchdown_times - trial_data.trigger_times(trigger_index)));
        [~, index_candidate_right] = min(abs(trial_data.right_touchdown_times - trial_data.trigger_times(trigger_index)));

        % is the closest left heelstrike within the acceptable interval?
        time_candidate_left = trial_data.left_touchdown_times(index_candidate_left);
        time_difference_left = time_candidate_left - trial_data.trigger_times(trigger_index); % where does the closest left heelstrike lie relative to the trigger?
        if -time_to_nearest_heelstrike_before_trigger_threshold < time_difference_left && time_difference_left < time_to_nearest_heelstrike_after_trigger_threshold
            % left heelstrike is acceptable
            candidate_left_acceptable = true;
        else
            candidate_left_acceptable = false;
        end

        % is the closest right heelstrike within the acceptable interval?
        time_candidate_right = trial_data.right_touchdown_times(index_candidate_right);
        time_difference_right = time_candidate_right - trial_data.trigger_times(trigger_index); % where does the closest right heelstrike lie relative to the trigger?
        if -time_to_nearest_heelstrike_before_trigger_threshold < time_difference_right && time_difference_right < time_to_nearest_heelstrike_after_trigger_threshold
            % right heelstrike is acceptable
            candidate_right_acceptable = true;
        else
            candidate_right_acceptable = false;
        end

        % accept the acceptable one
        if candidate_left_acceptable && ~candidate_right_acceptable
            % triggered by left heelstrike
            trigger_foot = 'left';
            trigger_time = time_candidate_left;
        elseif ~candidate_left_acceptable && candidate_right_acceptable
            % triggered by right heelstrike
            trigger_foot = 'right';
            trigger_time = time_candidate_right;
        elseif candidate_left_acceptable && candidate_right_acceptable
            trigger_foot = 'unclear';
            trigger_time = 0;
            removal_flags(trigger_index) = 1;
        elseif ~candidate_left_acceptable && ~candidate_right_acceptable
            trigger_foot = 'unclear';
            trigger_time = 0;
            removal_flags(trigger_index) = 1;
        end

end
function [stretch_times_this_trigger, stance_foot_data_this_trigger, removal_flag_this_trigger] = determineStretchTimes(trigger_time, trigger_foot)
    stretch_times_this_trigger = zeros(1, bands_per_stretch+1);
    stance_foot_data_this_trigger = cell(1, bands_per_stretch);
    removal_flag_this_trigger = false;

    % re-label time arrays as trigger and contralateral
    if strcmp(trigger_foot, 'left')
        trigger_foot_touchdown_times = trial_data.left_touchdown_times;
        contra_foot_touchdown_times = trial_data.right_touchdown_times;
        trigger_foot_stance_label = 'STANCE_LEFT';
        contra_foot_stance_label = 'STANCE_RIGHT';
        stretch_times_this_trigger(1) = trigger_time;
    end
    if strcmp(trigger_foot, 'right')
        trigger_foot_touchdown_times = trial_data.right_touchdown_times;
        contra_foot_touchdown_times = trial_data.left_touchdown_times;
        trigger_foot_stance_label = 'STANCE_RIGHT';
        contra_foot_stance_label = 'STANCE_LEFT';
        stretch_times_this_trigger(1) = trigger_time;
    end

    % now fill the other times
    for i_band = 1 : bands_per_stretch
        band_start_time = stretch_times_this_trigger(i_band);
        % find push-off time within the band
        if mod(i_band, 2) == 0
            swing_foot_touchdown_times = trigger_foot_touchdown_times;
            stance_foot_data_this_trigger{i_band} = contra_foot_stance_label;
        end
        if mod(i_band, 2) == 1
            swing_foot_touchdown_times = contra_foot_touchdown_times;
            stance_foot_data_this_trigger{i_band} = trigger_foot_stance_label;
        end
        this_band_end_time = min(swing_foot_touchdown_times(swing_foot_touchdown_times > band_start_time));
        if isempty(this_band_end_time)
            removal_flag_this_trigger = true;
        else
            stretch_times_this_trigger(i_band+1) = this_band_end_time;
        end
    end    
end
function stimulus_direction = determineStimulusDirection(trigger_foot, stimulus)
    if strcmp(trigger_foot, 'TRIGGER_RIGHT')
        if strcmp(stimulus, 'STIM_RIGHT')
            stimulus_direction = 'STIM_TOWARDS';
        end
        if strcmp(stimulus, 'STIM_LEFT')
            stimulus_direction = 'STIM_AWAY';
        end
        if strcmp(stimulus, 'STIM_NONE')
            stimulus_direction = 'STIM_NONE';
        end
    end
    if strcmp(trigger_foot, 'TRIGGER_LEFT')
        if strcmp(stimulus, 'STIM_RIGHT')
            stimulus_direction = 'STIM_AWAY';
        end
        if strcmp(stimulus, 'STIM_LEFT')
            stimulus_direction = 'STIM_TOWARDS';
        end
        if strcmp(stimulus, 'STIM_NONE')
            stimulus_direction = 'STIM_NONE';
        end
    end
end
end

