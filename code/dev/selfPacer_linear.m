function [pos_paced, pace] = selfPacer_linear(pos, time)
    p_fit = polyfit(time, pos, 1);
    belt_position_timeseries = p_fit(1) * time + p_fit(2);
    
    pos_paced = pos - belt_position_timeseries;
    pace = p_fit(1) * ones(size(time));
    
end