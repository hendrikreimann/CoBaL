function [pos_paced, pace] = selfPacer_pidWithAveraging(pos, time, gain_p, gain_d, v_belt_to_world_init, heelstrike_indices)
    number_of_time_steps = length(time);
    pelvis_velocity_belt_timeseries = deriveByTime(pos, time);
   
    p_treadmill = 0;
    v_treadmill = v_belt_to_world_init;
    
    indices_per_step = heelstrike_indices(2) - heelstrike_indices(1);
    
    belt_position_timeseries = zeros(number_of_time_steps, 1);
    belt_position_timeseries(1) = 0;
    belt_velocity_timeseries = zeros(number_of_time_steps, 1);
    belt_velocity_timeseries(1) = v_belt_to_world_init;
    
    pelvis_pos_world_timeseries = zeros(number_of_time_steps, 1);
    pelvis_pos_world_timeseries(1) = pos(1) - p_treadmill;
    pelvis_vel_world_timeseries = zeros(number_of_time_steps, 1);
    pelvis_vel_world_timeseries(1) = pelvis_velocity_belt_timeseries(1) - v_belt_to_world_init;
    
    pelvis_pos_world_avg_timeseries = zeros(number_of_time_steps, 1);
    pelvis_pos_world_avg_timeseries(1) = pos(1) - p_treadmill;
    pelvis_vel_world_avg_timeseries = zeros(number_of_time_steps, 1);
    pelvis_vel_world_avg_timeseries(1) = pelvis_velocity_belt_timeseries(1) - v_belt_to_world_init;
    
    step_memory_length = 3;
    step_memory_weights = (1 : step_memory_length) * 1 / mean(1 : step_memory_length);
    step_duration_update_mode = 'single';
    step_duration_update_mode = 'weighted_average';
    
    for i_time = 2 : number_of_time_steps
        p_pelvis_belt = pos(i_time);
        v_pelvis_belt = pelvis_velocity_belt_timeseries(i_time);
        
        % calculate state
        pelvis_pos_world = p_pelvis_belt - p_treadmill;
        pelvis_vel_world = v_pelvis_belt - v_treadmill;
        pelvis_pos_world_timeseries(i_time) = pelvis_pos_world;
        pelvis_vel_world_timeseries(i_time) = pelvis_vel_world;
        
        % average state across previous step
        if i_time >= indices_per_step
            pelvis_pos_world_avg = mean(pelvis_pos_world_timeseries(i_time - indices_per_step + 1 : i_time));
            pelvis_vel_world_avg = mean(pelvis_vel_world_timeseries(i_time - indices_per_step + 1 : i_time));
        else
            pelvis_pos_world_avg = pelvis_pos_world;
            pelvis_vel_world_avg = pelvis_vel_world;
        end
        pelvis_pos_world_avg_timeseries(i_time) = pelvis_pos_world_avg;
        pelvis_vel_world_avg_timeseries(i_time) = pelvis_vel_world_avg;
        
        % calculate controller signal
%         a_treadmill = gain_p * pelvis_pos_world + gain_d * pelvis_vel_world;
        a_treadmill = gain_p * pelvis_pos_world_avg + gain_d * pelvis_vel_world_avg;
        
        % time step
        dt = time(i_time) - time(i_time-1);
        v_treadmill = v_treadmill + dt * a_treadmill;
        p_treadmill = p_treadmill + dt * v_treadmill + dt^2 * a_treadmill;
        
        % store
        belt_position_timeseries(i_time) = p_treadmill;
        belt_velocity_timeseries(i_time) = v_treadmill;
        
        % update step duration - weighted average
        if strcmp(step_duration_update_mode, 'single')
            if ismember(i_time, heelstrike_indices) && i_time > heelstrike_indices(1)
                index_in_heelstrikes = find(heelstrike_indices==i_time);
                indices_per_step = heelstrike_indices(index_in_heelstrikes) - heelstrike_indices(index_in_heelstrikes-1);
            end
        end
        if strcmp(step_duration_update_mode, 'weighted average')
            if ismember(i_time, heelstrike_indices) && i_time > heelstrike_indices(step_memory_length)
                index_in_heelstrikes = find(heelstrike_indices==i_time);

                step_memory_indices_in_heelstrikes = index_in_heelstrikes - step_memory_length : index_in_heelstrikes;
                step_memory_indices_in_time = heelstrike_indices(step_memory_indices_in_heelstrikes);
                step_memory_differences = diff(step_memory_indices_in_time);
                step_memory_differences_weighted = step_memory_weights .* step_memory_differences';
                indices_per_step = round(mean(step_memory_differences_weighted));
            end
        end
    end
    
    pos_paced = pos - belt_position_timeseries;
    pace = belt_velocity_timeseries;
end