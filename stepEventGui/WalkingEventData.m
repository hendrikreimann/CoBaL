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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.% compare the kinematic tree against the kinematic chain

classdef WalkingEventData < handle
    properties
        trial_data;
        
        left_pushoff;
        left_touchdown;
        left_fullstance_times;
        right_pushoff;
        right_touchdown;
        right_fullstance_times
        
        left_arm_swing_onset_times;
        right_arm_swing_onset_times;
        left_leg_swing_onset_times;
        right_leg_swing_onset_times;
        
        ignore_times;
        
        selected_event_label;
        selected_event_time;
        
        stretch_start_times;
        stretch_end_times;
        
    end
    methods
        function this = WalkingEventData(trialData)
            this.trial_data = trialData;
            this.loadEvents();
            this.loadStretches();
        end
        function loadEvents(this)
            step_events_file_name = [this.trial_data.data_directory filesep 'analysis' filesep makeFileName(this.trial_data.date, this.trial_data.subject_id, this.trial_data.condition, this.trial_data.trial_number, 'stepEvents.mat')];
            if exist(step_events_file_name, 'file')
                loaded_event_data = load(step_events_file_name);
                this.left_pushoff = loaded_event_data.left_pushoff_times;
                this.left_touchdown = loaded_event_data.left_touchdown_times;
                this.right_pushoff = loaded_event_data.right_pushoff_times;
                this.right_touchdown = loaded_event_data.right_touchdown_times;
            else
                loaded_event_data = struct;
            end
            if isfield(loaded_event_data, 'left_fullstance_times')
                this.left_fullstance_times = loaded_event_data.left_fullstance_times;
            end
            if isfield(loaded_event_data, 'right_fullstance_times')
                this.right_fullstance_times = loaded_event_data.right_fullstance_times;
            end
            if isfield(loaded_event_data, 'left_arm_swing_onset_times')
                this.left_arm_swing_onset_times = loaded_event_data.left_arm_swing_onset_times;
            end
            if isfield(loaded_event_data, 'right_arm_swing_onset_times')
                this.right_arm_swing_onset_times = loaded_event_data.right_arm_swing_onset_times;
            end
            if isfield(loaded_event_data, 'left_leg_swing_onset_times')
                this.left_leg_swing_onset_times = loaded_event_data.left_leg_swing_onset_times;
            end
            if isfield(loaded_event_data, 'right_leg_swing_onset_times')
                this.right_leg_swing_onset_times = loaded_event_data.right_leg_swing_onset_times;
            end
            if isfield(loaded_event_data, 'ignore_times')
                this.ignore_times = loaded_event_data.ignore_times;
            else
                this.ignore_times = [];
            end
            
            this.removeDuplicates();
        end
        function loadStretches(this)
            relevant_stretches_file_name = [this.trial_data.data_directory filesep 'analysis' filesep makeFileName(this.trial_data.date, this.trial_data.subject_id, this.trial_data.condition, this.trial_data.trial_number, 'relevantDataStretches.mat')];
            if exist(relevant_stretches_file_name, 'file')
                loaded_stretch_data = load(relevant_stretches_file_name);
                this.stretch_start_times = loaded_stretch_data.stretch_start_times;
                this.stretch_end_times = loaded_stretch_data.stretch_end_times;
            end
        end
        function saveEvents(this, sender, eventdata) %#ok<INUSD>
            % prepare data
            this.removeDuplicates();
            variables_to_save = struct;
            
            variables_to_save.left_pushoff_times = this.left_pushoff;
            variables_to_save.left_touchdown_times = this.left_touchdown;
            variables_to_save.right_pushoff_times = this.right_pushoff;
            variables_to_save.right_touchdown_times = this.right_touchdown;
            
            if ~isempty(this.left_fullstance_times)
                variables_to_save.left_fullstance_times = this.left_fullstance_times;
            end
            if ~isempty(this.right_fullstance_times)
                variables_to_save.right_fullstance_times = this.right_fullstance_times;
            end
            if ~isempty(this.left_arm_swing_onset_times)
                variables_to_save.left_arm_swing_onset_times = this.left_arm_swing_onset_times;
            end
            if ~isempty(this.right_arm_swing_onset_times)
                variables_to_save.right_arm_swing_onset_times = this.right_arm_swing_onset_times;
            end
            if ~isempty(this.left_leg_swing_onset_times)
                variables_to_save.left_leg_swing_onset_times = this.left_leg_swing_onset_times;
            end
            if ~isempty(this.right_leg_swing_onset_times)
                variables_to_save.right_leg_swing_onset_times = this.right_leg_swing_onset_times;
            end
            if ~isempty(this.ignore_times)
                variables_to_save.ignore_times = this.ignore_times;
            end
            
            step_events_file_name = [this.trial_data.data_directory filesep 'analysis' filesep makeFileName(this.trial_data.date, this.trial_data.subject_id, this.trial_data.condition, this.trial_data.trial_number, 'stepEvents')];
            saveDataToFile(step_events_file_name, variables_to_save);


            disp(['Step events saved as "' step_events_file_name '"']);
        end
        
        function setEventTimes(this, event_times, event_label) %#ok<INUSL>
            eval(['this.' event_label ' = event_times;']);
        end
        function addEventTime(this, event_time, event_label)
            event_times = this.getEventTimes(event_label);
            event_times = [event_times; event_time];
            event_times = sort(event_times);
            this.setEventTimes(sort(event_times), event_label);
            this.selected_event_time = event_time;
            this.selected_event_label = event_label;
        end
        function event_times = getEventTimes(this, event_label) %#ok<INUSL,STOUT>
            eval(['event_times = this.' event_label ';']);
        end
        function event_index = getEventIndex(this, event_label, event_time)
            event_times = this.getEventTimes(event_label);
            [~, event_index] = min(abs(event_times - event_time));
        end
        function new_event_time = updateEventTime(this, event_label, current_event_time, new_event_time)
            event_index = this.getEventIndex(event_label, current_event_time);
            event_times = this.getEventTimes(event_label);
            if isnan(new_event_time)
                % delete event
                event_times(event_index) = [];
                this.setEventTimes(event_times, event_label);
                this.selectNextEvent();
                
            else
                % clamp new time to limits
                new_event_time = max([new_event_time, 0]);
                new_event_time = min([new_event_time, this.trial_data.recording_time]);
                
                % update
                event_times(event_index) = new_event_time;
                
                % sort and store
                this.setEventTimes(sort(event_times), event_label);
            end
            
            
            
            
        end
        function selectNextEvent(this)
            % find index of currently selected event
            currently_selected_event_time = this.selected_event_time;
            event_data_of_current_type = this.getEventTimes(this.selected_event_label);
            event_data_of_current_type_after_currently_selected = event_data_of_current_type(event_data_of_current_type > currently_selected_event_time);
            
            if isempty(event_data_of_current_type_after_currently_selected)
                % the selected event was the last of this type, so go to next type
                this.selected_event_label = this.getNextEventTypeLabel(this.selected_event_label);
                event_data_of_current_type = this.getEventTimes(this.selected_event_label);
                this.selected_event_time = event_data_of_current_type(1);
            else
                % select next one
                this.selected_event_time = event_data_of_current_type_after_currently_selected(1);
            end
            this.trial_data.selected_time = this.selected_event_time;
        end
        function selectPreviousEvent(this)
            % find index of currently selected event
            currently_selected_event_time = this.selected_event_time;
            event_data_of_current_type = this.getEventTimes(this.selected_event_label);
            event_data_of_current_type_before_currently_selected = event_data_of_current_type(event_data_of_current_type < currently_selected_event_time);
            
            if isempty(event_data_of_current_type_before_currently_selected)
                % the selected event was the last of this type, so go to next type
                this.selected_event_label = this.getPreviousEventTypeLabel(this.selected_event_label);
                event_data_of_current_type = this.getEventTimes(this.selected_event_label);
                this.selected_event_time = event_data_of_current_type(end);
            else
                % select next one
                this.selected_event_time = event_data_of_current_type_before_currently_selected(end);
            end
            this.trial_data.selected_time = this.selected_event_time;
        end
        function next_event_label = getNextEventTypeLabel(this, current_event_type_label) %#ok<INUSL>
            next_event_label = 'left_pushoff';
            if strcmp(current_event_type_label, 'left_touchdown')
                next_event_label = 'left_pushoff';
            end
            if strcmp(current_event_type_label, 'left_pushoff')
                next_event_label = 'right_touchdown';
            end
            if strcmp(current_event_type_label, 'right_touchdown')
                next_event_label = 'right_pushoff';
            end
                

        end
        function previous_event_label = getPreviousEventTypeLabel(this, current_event_type_label) %#ok<INUSL>
            previous_event_label = 'left_pushoff';
            if strcmp(current_event_type_label, 'right_pushoff')
                previous_event_label = 'right_touchdown';
            end
            if strcmp(current_event_type_label, 'left_pushoff')
                previous_event_label = 'left_touchdown';
            end
            if strcmp(current_event_type_label, 'left_touchdown')
                previous_event_label = 'right_pushoff';
            end
        end
        function removeDuplicates(this)
            this.left_pushoff = unique(this.left_pushoff);
            this.left_touchdown = unique(this.left_touchdown);
            this.right_pushoff = unique(this.right_pushoff);
            this.right_touchdown = unique(this.right_touchdown);
        end
    end
end
