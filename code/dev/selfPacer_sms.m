function [pos_paced, pace] = selfPacer_sms(pos, time, v_belt_to_world_init, g_p, g_v, heelstrike_indices)
    number_of_time_steps = length(time);
    pelvis_velocity_belt_timeseries = deriveByTime(pos, time);
   
    p_treadmill = 0;
    v_treadmill = v_belt_to_world_init;
    v_treadmill_target = v_treadmill;
    
    belt_position_timeseries = zeros(number_of_time_steps, 1);
    belt_position_timeseries(1) = 0;
    belt_velocity_timeseries = zeros(number_of_time_steps, 1);
    belt_velocity_timeseries(1) = v_belt_to_world_init;
    
    pelvis_pos_world_timeseries = zeros(number_of_time_steps, 1);
    pelvis_pos_world_timeseries(1) = pos(1) - p_treadmill;
    pelvis_vel_world_timeseries = zeros(number_of_time_steps, 1);
    pelvis_vel_world_timeseries(1) = pelvis_velocity_belt_timeseries(1) - v_belt_to_world_init;
    belt_velocity_avg = v_belt_to_world_init;
    
    for i_time = 2 : number_of_time_steps
        p_pelvis_belt = pos(i_time);
        v_pelvis_belt = pelvis_velocity_belt_timeseries(i_time);
        
        % calculate state
        pelvis_pos_world = p_pelvis_belt - p_treadmill;
        pelvis_vel_world = v_pelvis_belt - v_treadmill;
        pelvis_pos_world_timeseries(i_time) = pelvis_pos_world;
        pelvis_vel_world_timeseries(i_time) = pelvis_vel_world;
        
        % set new target velocity on heel-strikes
        if ismember(i_time, heelstrike_indices) && i_time > heelstrike_indices(1)
            index_in_heelstrikes = find(heelstrike_indices==i_time);
            indices_per_step = heelstrike_indices(index_in_heelstrikes) - heelstrike_indices(index_in_heelstrikes-1);
            
            pelvis_pos_world_avg = mean(pelvis_pos_world_timeseries(i_time - indices_per_step : i_time-1));
            pelvis_vel_world_avg = mean(pelvis_vel_world_timeseries(i_time - indices_per_step : i_time-1));
            belt_velocity_avg = mean(belt_velocity_timeseries(i_time - indices_per_step : i_time-1));
            
            v_treadmill_target = belt_velocity_avg + g_v * pelvis_vel_world_avg + g_p * pelvis_pos_world_avg;
            
        end
        
        delta_treadmill_target = 0.5;
%         a_treadmill = - (v_treadmill - v_treadmill_target) / delta_treadmill_target;
        a_treadmill = - (belt_velocity_avg - v_treadmill_target) / delta_treadmill_target;
        
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