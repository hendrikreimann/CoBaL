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
        function setLeftTouchdownTimes(this, left_touchdown)
            this.left_touchdown = left_touchdown;
        end
        function setLeftPushoffTimes(this, left_pushoff)
            this.left_pushoff = left_pushoff;
        end
        function setRightTouchdownTimes(this, right_touchdown)
            this.right_touchdown = right_touchdown;
        end
        function setRightPushoffTimes(this, right_pushoff)
            this.right_pushoff = right_pushoff;
        end
        function event_times = getEventTimes(this, event_label)
            eval(['event_times = this.' event_label ';']);
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
        function next_event_label = getPreviousEventTypeLabel(this, current_event_type_label)
            if strcmp(current_event_type_label, 'right_pushoff')
                next_event_label = 'right_touchdown';
            elseif strcmp(current_event_type_label, 'right_touchdown')
                next_event_label = 'left_pushoff';
            elseif strcmp(current_event_type_label, 'left_pushoff')
                next_event_label = 'left_touchdown';
            elseif strcmp(current_event_type_label, 'left_touchdown')
                next_event_label = 'right_pushoff';
            end
                

        end
    end
end
