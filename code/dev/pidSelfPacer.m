function belt_velocity_timeseries = pidSelfPacer(pelvis_position_timeseries, time, gain_p, gain_i, gain_d, v_treadmill_init)
    number_of_time_steps = length(time);
    belt_velocity_timeseries = zeros(1, number_of_time_steps);
    belt_velocity_timeseries(1) = v_treadmill_init;
    pelvis_velocity_belt_timeseries = deriveByTime(pelvis_position_timeseries, time);
   
    p_treadmill = 0;
    v_treadmill = v_treadmill_init;
    
    for i_time = 2 : number_of_time_steps
        p_pelvis_belt = pelvis_position_timeseries(i_time);
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
        belt_velocity_timeseries(i_time) = v_treadmill;
    end
    
end