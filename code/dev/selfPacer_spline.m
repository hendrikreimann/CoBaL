function [pos_paced, pace] = selfPacer_spline(pos, time, heelstrike_indices)
    pelvis_pos_belt_splined = spline(heelstrike_indices, pos(heelstrike_indices), (1 : length(time))');

    belt_position_timeseries = pelvis_pos_belt_splined;
    
    pos_paced = pos - belt_position_timeseries;
    pace = deriveByTime(belt_position_timeseries, time);
end