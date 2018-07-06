function sawtooth = generateSawtooth(time, period, offset_time, offset_value, amplitude)
    
    slope = amplitude / (period/2);
    sawtooth = zeros(size(time));
    for i_time = 1 : length(time)
        this_time_original = time(i_time);
        
        % remove offset
        this_time_offset_removed = this_time_original - offset_time;
        
        % remove period
        this_time_current = this_time_offset_removed;
        while this_time_current > period
            this_time_current = this_time_current - period;
        end
        while this_time_current < 0
            this_time_current = this_time_current + period;
        end
        this_time_period_removed = this_time_current;
        
        % apply slope
        if this_time_period_removed < period/2
            this_time_slope_applied = this_time_period_removed * slope;
        end
        if this_time_period_removed >= period/2
            this_time_slope_applied = - this_time_period_removed * slope + period * slope;
        end
        
        % apply value offset
        this_time_value_offset_applied = this_time_slope_applied + offset_value;
        
        sawtooth(i_time) = this_time_value_offset_applied;
    end
end