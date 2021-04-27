function [pos_paced, pace] = selfPacer_pid(pos, time, gain_p, gain_d, v_belt_to_world_init)
    number_of_time_steps = length(time);
    pelvis_velocity_belt_timeseries = deriveByTime(pos, time);
   
    p_treadmill = 0;
    v_treadmill = v_belt_to_world_init;
    
    belt_position_timeseries = zeros(number_of_time_steps, 1);
    belt_position_timeseries(1) = 0;
    belt_velocity_timeseries = zeros(number_of_time_steps, 1);
    belt_velocity_timeseries(1) = v_belt_to_world_init;
    for i_time = 2 : number_of_time_steps
        p_pelvis_belt = pos(i_time);
        v_pelvis_belt = pelvis_velocity_belt_timeseries(i_time);
        
        % calculate state
        pelvis_pos_world = p_pelvis_belt - p_treadmill;
        pelvis_vel_world = v_pelvis_belt - v_treadmill;
        
        % calculate controller signal
        a_treadmill = gain_p * pelvis_pos_world + gain_d * pelvis_vel_world;
        
        % time step
        dt = time(i_time) - time(i_time-1);
        v_treadmill = v_treadmill + dt * a_treadmill;
        p_treadmill = p_treadmill + dt * v_treadmill + dt^2 * a_treadmill;
        
        % store
        belt_position_timeseries(i_time) = p_treadmill;
        belt_velocity_timeseries(i_time) = v_treadmill;
    end
    
    pos_paced = pos - belt_position_timeseries;
    pace = belt_velocity_timeseries;
    
    
    
    
    
    
    
end