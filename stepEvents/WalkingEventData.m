classdef WalkingEventData < handle
    properties
        trial_data;
        
        left_pushoff;
        left_touchdown;
        right_pushoff;
        right_touchdown;
    end
    methods
        function this = WalkingEventData(trialData)
            this.trial_data = trialData;
            this.loadEvents();
        end
        function loadEvents(this)
            loaded_event_data = load([this.trial_data.data_directory filesep makeFileName(this.trial_data.date, this.trial_data.subject_id, 'walking', this.trial_data.trial_number, 'stepEvents')]);
            this.left_pushoff = loaded_event_data.left_pushoff_times;
            this.left_touchdown = loaded_event_data.left_touchdown_times;
            this.right_pushoff = loaded_event_data.right_pushoff_times;
            this.right_touchdown = loaded_event_data.right_touchdown_times;
        end
        function saveEvents(this, sender, eventdata)
            % prepare data
            left_pushoff_times = this.left_pushoff;
            left_touchdown_times = this.left_touchdown;
            right_pushoff_times = this.right_pushoff;
            right_touchdown_times = this.right_touchdown;
            
            % find corresponding indices - mocap
            left_pushoff_indices_mocap = zeros(size(this.left_pushoff));
            for i_event = 1 : length(this.left_pushoff)
                [~, left_pushoff_indices_mocap(i_event)] = min(abs(this.trial_data.time_mocap - this.left_pushoff(i_event)));
            end
            left_touchdown_indices_mocap = zeros(size(this.left_touchdown));
            for i_event = 1 : length(this.left_touchdown)
                [~, left_touchdown_indices_mocap(i_event)] = min(abs(this.trial_data.time_mocap - this.left_touchdown(i_event)));
            end
            right_pushoff_indices_mocap = zeros(size(this.right_pushoff));
            for i_event = 1 : length(this.right_pushoff)
                [~, right_pushoff_indices_mocap(i_event)] = min(abs(this.trial_data.time_mocap - this.right_pushoff(i_event)));
            end
            right_touchdown_indices_mocap = zeros(size(this.right_touchdown));
            for i_event = 1 : length(this.right_touchdown)
                [~, right_touchdown_indices_mocap(i_event)] = min(abs(this.trial_data.time_mocap - this.right_touchdown(i_event)));
            end

            % find corresponding indices - forceplate
            left_pushoff_indices_forceplate = zeros(size(this.left_pushoff));
            for i_event = 1 : length(this.left_pushoff)
                [~, left_pushoff_indices_forceplate(i_event)] = min(abs(this.trial_data.time_forceplate - this.left_pushoff(i_event)));
            end
            left_touchdown_indices_forceplate = zeros(size(this.left_touchdown));
            for i_event = 1 : length(this.left_touchdown)
                [~, left_touchdown_indices_forceplate(i_event)] = min(abs(this.trial_data.time_forceplate - this.left_touchdown(i_event)));
            end
            right_pushoff_indices_forceplate = zeros(size(this.right_pushoff));
            for i_event = 1 : length(this.right_pushoff)
                [~, right_pushoff_indices_forceplate(i_event)] = min(abs(this.trial_data.time_forceplate - this.right_pushoff(i_event)));
            end
            right_touchdown_indices_forceplate = zeros(size(this.right_touchdown));
            for i_event = 1 : length(this.right_touchdown)
                [~, right_touchdown_indices_forceplate(i_event)] = min(abs(this.trial_data.time_forceplate - this.right_touchdown(i_event)));
            end
            
            % form contact indicators
            left_contact_indicators_mocap = formContactIndicatorTrajectory(left_pushoff_indices_mocap, left_touchdown_indices_mocap, length(this.trial_data.time_mocap));
            right_contact_indicators_mocap = formContactIndicatorTrajectory(right_pushoff_indices_mocap, right_touchdown_indices_mocap, length(this.trial_data.time_mocap));
            left_contact_indicators_forceplate = formContactIndicatorTrajectory(left_pushoff_indices_forceplate, left_touchdown_indices_forceplate, length(this.trial_data.time_forceplate));
            right_contact_indicators_forceplate = formContactIndicatorTrajectory(right_pushoff_indices_forceplate, right_touchdown_indices_forceplate, length(this.trial_data.time_forceplate));

            step_events_file_name = makeFileName(this.trial_data.date, this.trial_data.subject_id, 'walking', this.trial_data.trial_number, 'stepEvents');
            save ...
              ( ...
                [this.trial_data.data_directory filesep step_events_file_name], ...
                'left_pushoff_indices_mocap', ...
                'left_touchdown_indices_mocap', ...
                'right_pushoff_indices_mocap', ...
                'right_touchdown_indices_mocap', ...
                'left_contact_indicators_mocap', ...
                'right_contact_indicators_mocap', ...
                'left_pushoff_indices_forceplate', ...
                'left_touchdown_indices_forceplate', ...
                'right_pushoff_indices_forceplate', ...
                'right_touchdown_indices_forceplate', ...
                'left_contact_indicators_forceplate', ...
                'right_contact_indicators_forceplate', ...
                'left_pushoff_times', ...
                'left_touchdown_times', ...
                'right_pushoff_times', ...
                'right_touchdown_times' ...
              );
        end
        
        function setEventTimes(this, event_times, event_label)
            eval(['this.' event_label ' = event_times;']);
        end
        function event_times = getEventTimes(this, event_label)
            eval(['event_times = this.' event_label ';']);
        end
        function event_index = getEventIndex(this, event_label, event_time)
            event_times = this.getEventTimes(event_label);
            [~, event_index] = min(abs(event_times - event_time));
        end
        function new_event_time = updateEventTime(this, event_label, current_event_time, new_event_time)
            % clamp to limits
            new_event_time = max([new_event_time, 0]);
            new_event_time = min([new_event_time, this.trial_data.recording_time]);
            
            % update
            event_index = this.getEventIndex(event_label, current_event_time);
            event_times = this.getEventTimes(event_label);
            event_times(event_index) = new_event_time;
            
            % sort and store
            this.setEventTimes(sort(event_times), event_label);
            
        end
        function next_event_label = getNextEventTypeLabel(this, current_event_type_label)
            if strcmp(current_event_type_label, 'left_touchdown')
                next_event_label = 'left_pushoff';
            elseif strcmp(current_event_type_label, 'left_pushoff')
                next_event_label = 'right_touchdown';
            elseif strcmp(current_event_type_label, 'right_touchdown')
                next_event_label = 'right_pushoff';
            elseif strcmp(current_event_type_label, 'right_pushoff')
                next_event_label = 'left_touchdown';
            end
                

        end
        function previous_event_label = getPreviousEventTypeLabel(this, current_event_type_label)
            if strcmp(current_event_type_label, 'right_pushoff')
                previous_event_label = 'right_touchdown';
            elseif strcmp(current_event_type_label, 'right_touchdown')
                previous_event_label = 'left_pushoff';
            elseif strcmp(current_event_type_label, 'left_pushoff')
                previous_event_label = 'left_touchdown';
            elseif strcmp(current_event_type_label, 'left_touchdown')
                previous_event_label = 'right_pushoff';
            end
                

        end
    end
end
