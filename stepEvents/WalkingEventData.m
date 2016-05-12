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
        function data = getData(this, data_label)
            eval(['data = this.' data_label ';']);
        end
    end
end
