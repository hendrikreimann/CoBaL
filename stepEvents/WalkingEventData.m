classdef WalkingEventData < handle
    properties
        trial_data;
        
        left_touchdown;
        left_pushoff;
        right_touchdown;
        right_pushoff;
    end
    methods
        function object = WalkingEventData(trialData)
            object.trial_data = trialData;
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
