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

classdef eventData < handle
    properties
        data_custodian;
        
        event_data;
        event_labels;
        
        ignore_times;
        
        selected_event_label;
        selected_event_time;
        selected_time;
        
        stretch_start_times;
        stretch_end_times;
        
    end
    methods
        function this = eventData(data_custodian)
            this.data_custodian = data_custodian;
            this.loadEvents();
            this.loadStretches();
        end
        function loadEvents(this)
            step_events_file_name = [this.data_custodian.data_directory filesep 'analysis' filesep makeFileName(this.data_custodian.date, this.data_custodian.subject_id, this.data_custodian.trial_type, this.data_custodian.trial_number, 'events.mat')];
            if exist(step_events_file_name, 'file')
                loaded_data = load(step_events_file_name);
                this.event_data = loaded_data.event_data;
                this.event_labels = loaded_data.event_labels;
                
            else
                loaded_data = struct;
                this.event_data = {};
                this.event_labels = {};
            end
            
            if isfield(loaded_data, 'ignore_times')
                this.ignore_times = loaded_data.ignore_times;
            else
                this.ignore_times = [];
            end
            
            this.removeDuplicates();
        end
        function loadStretches(this)
            relevant_stretches_file_name = [this.data_custodian.data_directory filesep 'analysis' filesep makeFileName(this.data_custodian.date, this.data_custodian.subject_id, this.data_custodian.trial_type, this.data_custodian.trial_number, 'relevantDataStretches.mat')];
            if exist(relevant_stretches_file_name, 'file')
                loaded_stretch_data = load(relevant_stretches_file_name);
                if isempty(loaded_stretch_data.stretch_times)
                    this.stretch_start_times = [];
                    this.stretch_end_times = [];
                else
                    this.stretch_start_times = loaded_stretch_data.stretch_times(:, 1);
                    this.stretch_end_times = loaded_stretch_data.stretch_times(:, end);
                end
            end
        end
        function saveEvents(this, sender, eventdata) %#ok<INUSD>
            % prepare data
            this.removeDuplicates();
            variables_to_save = struct;
            
%             variables_to_save.left_pushoff_times = this.left_pushoff;
%             variables_to_save.left_touchdown_times = this.left_touchdown;
%             variables_to_save.right_pushoff_times = this.right_pushoff;
%             variables_to_save.right_touchdown_times = this.right_touchdown;
%             
%             if ~isempty(this.left_fullstance_times)
%                 variables_to_save.left_fullstance_times = this.left_fullstance_times;
%             end
%             if ~isempty(this.right_fullstance_times)
%                 variables_to_save.right_fullstance_times = this.right_fullstance_times;
%             end
%             if ~isempty(this.left_arm_swing_onset_times)
%                 variables_to_save.left_arm_swing_onset_times = this.left_arm_swing_onset_times;
%             end
%             if ~isempty(this.right_arm_swing_onset_times)
%                 variables_to_save.right_arm_swing_onset_times = this.right_arm_swing_onset_times;
%             end
%             if ~isempty(this.left_leg_swing_onset_times)
%                 variables_to_save.left_leg_swing_onset_times = this.left_leg_swing_onset_times;
%             end
%             if ~isempty(this.right_leg_swing_onset_times)
%                 variables_to_save.right_leg_swing_onset_times = this.right_leg_swing_onset_times;
%             end
            
% check this out later, get things to work first
            variables_to_save.event_data = this.event_data;
            variables_to_save.event_labels = this.event_labels;
            variables_to_save.ignore_times = this.ignore_times;
            
            step_events_file_name = [this.data_custodian.data_directory filesep 'analysis' filesep makeFileName(this.data_custodian.date, this.data_custodian.subject_id, this.data_custodian.trial_type, this.data_custodian.trial_number, 'events.mat')];
            saveDataToFile(step_events_file_name, variables_to_save);
            
            disp(['Step events saved as "' step_events_file_name '"']);
        end
        
        function setEventTimes(this, event_times, event_label)
            if strcmp(event_label, 'ignore_times')
                this.ignore_times = event_times;
            else
                event_index = strcmp(this.event_labels, event_label);
                this.event_data{event_index} = event_times;
            end
        end
        function addEventTime(this, event_time, event_label)
            if strcmp(event_label, 'ignore_times')
                ignore_times_current = this.getEventTimes(event_label);
                ignore_times_new = [ignore_times_current; event_time];
                this.ignore_times = sort(ignore_times_new);
            else
                event_times = this.getEventTimes(event_label);
                event_times = [event_times; event_time];
                event_times = sort(event_times);
                this.setEventTimes(sort(event_times), event_label);
            end
            this.selected_event_time = event_time;
            this.selected_event_label = event_label;
        end
        function event_times = getEventTimes(this, event_label)
            if strcmp(event_label, 'ignore_times')
                event_times = this.ignore_times;
            else
                event_index = strcmp(this.event_labels, event_label);
                if any(event_index)
                    event_times = this.event_data{event_index};
                else
                    event_times = [];
                end
            end
            
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
                this.selectClosestEvent();
                
            else
                % clamp new time to limits
                new_event_time = max([new_event_time, this.data_custodian.getRecordingTimeStart]);
                new_event_time = min([new_event_time, this.data_custodian.getRecordingTimeEnd]);
                
                % update
                event_times(event_index) = new_event_time;
                
                % sort and store
                this.setEventTimes(sort(event_times), event_label);
            end
        end
        function selectNextEvent(this)
            % find index of currently selected event
            currently_selected_event_time = this.selected_event_time;
            if ~isempty(currently_selected_event_time)
                event_data_of_current_type = this.getEventTimes(this.selected_event_label);
                if isempty(event_data_of_current_type)
                    this.selected_event_time = [];
                end
                if ~isempty(event_data_of_current_type)
                    event_data_of_current_type_after_currently_selected = event_data_of_current_type(event_data_of_current_type > currently_selected_event_time);
                    if ~isempty(event_data_of_current_type_after_currently_selected)
                        this.selected_event_time = event_data_of_current_type_after_currently_selected(1);
                    end
                end
                this.selected_time = this.selected_event_time;
                
            end
        end
        function selectPreviousEvent(this)
            % find index of currently selected event
            currently_selected_event_time = this.selected_event_time;
            if ~isempty(currently_selected_event_time)
                event_data_of_current_type = this.getEventTimes(this.selected_event_label);
                if isempty(event_data_of_current_type)
                    this.selected_event_time = [];
                end
                if ~isempty(event_data_of_current_type)
                    event_data_of_current_type_before_currently_selected = event_data_of_current_type(event_data_of_current_type < currently_selected_event_time);
                    if ~isempty(event_data_of_current_type_before_currently_selected)
                        this.selected_event_time = event_data_of_current_type_before_currently_selected(end);
                    end
                end
                this.selected_time = this.selected_event_time;
            end
        end
        function setSelectedTime(this, new_time)
            this.selected_time = new_time;
        end
        function selectClosestEvent(this)
            % find index of currently selected event
            currently_selected_event_time = this.selected_event_time;
            if ~isempty(currently_selected_event_time)
                event_data_of_current_type = this.getEventTimes(this.selected_event_label);
                if isempty(event_data_of_current_type)
                    this.selected_event_time = [];
                end
                if ~isempty(event_data_of_current_type)
                    if ~isempty(event_data_of_current_type)
                        [~, candidate_index] = min(abs(currently_selected_event_time - event_data_of_current_type));
                        
                        this.selected_event_time = event_data_of_current_type(candidate_index);
                    end
                end
                this.selected_time = this.selected_event_time;
            end
            
        end
        
        function stepSelectedTime(this, direction, stepsize)
            warning('Stepping through time currently not implemented, wait until scene figures are back!')
%             if nargin < 3
%                 stepsize = 1;
%             end
%             
%             % get current time step
%             [~, time_index_mocap] = min(abs(this.time_marker - this.selected_time));
%             
%             if strcmp(direction, 'back')
%                 new_time_index_mocap = time_index_mocap - stepsize;
%             elseif strcmp(direction, 'forward')
%                 new_time_index_mocap = time_index_mocap + stepsize;
%             else
%                 error('Direction must be either "back" or "forward"');
%             end
%             
%             % enforce limits
%             if new_time_index_mocap < 1
%                 new_time_index_mocap = 1;
%             elseif new_time_index_mocap > length(this.time_marker)
%                 new_time_index_mocap = length(this.time_marker);
%             end
%             
%             % set result
%             this.selected_time = this.time_marker(new_time_index_mocap);
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
            for i_type = 1 : length(this.event_data)
                this.event_data{i_type} = unique(this.event_data{i_type});
            end
            
        end
    end
end
