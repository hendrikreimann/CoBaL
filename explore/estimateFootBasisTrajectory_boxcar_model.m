function foot_basis_trajectory = estimateFootBasisTrajectory_boxcar_model ...
  ( ...
    time, ...
    contact_trajectory, ...
    left_step_times, ...
    right_step_times, ...
    sampling_rate, ...
    metronome_frequency, ...
    visualize ...
  )

    if nargin < 7
        visualize = 0;
    end

    left_foot_x = contact_trajectory;
    right_foot_x = contact_trajectory;

    % remove swing time
    left_heel_strikes = left_step_times;
    right_heel_strikes = right_step_times;
    left_pushoffs = right_step_times;
    right_pushoffs = left_step_times;
    left_heel_strike_indices = findClosestIndex(left_heel_strikes, time);
    right_heel_strike_indices = findClosestIndex(right_heel_strikes, time);
    left_pushoff_indices = findClosestIndex(left_pushoffs, time);
    right_pushoff_indices = findClosestIndex(right_pushoffs, time);

    left_foot_x_stance = left_foot_x;
    for i_heelstrike = 1 : length(right_heel_strike_indices)
        this_heelstrike_time_index = right_heel_strike_indices(i_heelstrike);
        next_pushoff_time_index = min(right_pushoff_indices(right_pushoff_indices>this_heelstrike_time_index));
        left_foot_x_stance(this_heelstrike_time_index : next_pushoff_time_index) = NaN;
    end
    right_foot_x_stance = right_foot_x;
    for i_heelstrike = 1 : length(left_heel_strike_indices)
        this_heelstrike_time_index = left_heel_strike_indices(i_heelstrike);
        next_pushoff_time_index = min(left_pushoff_indices(left_pushoff_indices>this_heelstrike_time_index));
        right_foot_x_stance(this_heelstrike_time_index : next_pushoff_time_index) = NaN;
    end

    % boxcar
    boxcar_width_time = 1 / metronome_frequency;
    boxcar_width = ceil(boxcar_width_time * sampling_rate);
    left_foot_x_boxcar = movmean(left_foot_x_stance, boxcar_width, 'omitnan');
    right_foot_x_boxcar = movmean(right_foot_x_stance, boxcar_width, 'omitnan');
    
    % butterworth
    [b_filter, a_filter] = butter(2, 0.5/(sampling_rate/2));
    left_foot_x_butter = nanfiltfilt(b_filter, a_filter, left_foot_x_boxcar);
    right_foot_x_butter = nanfiltfilt(b_filter, a_filter, right_foot_x_boxcar);
    
    % average
    foot_basis_trajectory = (left_foot_x_butter + right_foot_x_butter) * 0.5;
    
    if visualize
        colors = lines(4);
        figure; hold on;
        plot(time, left_foot_x, ':', 'color', colors(1, :), 'DisplayName', 'left foot position');
        plot(time, left_foot_x_stance, 'linewidth', 4, 'color', colors(1, :), 'DisplayName', 'left foot single stance');
        plot(time, left_foot_x_boxcar, 'color', colors(1, :), 'linewidth', 1, 'DisplayName', 'left foot boxcar');
        plot(time, left_foot_x_butter, 'color', colors(3, :), 'linewidth', 1, 'DisplayName', 'left foot filtered');
        plot(time, right_foot_x, ':', 'color', colors(2, :), 'DisplayName', 'right foot position');
        plot(time, right_foot_x_stance, 'linewidth', 4, 'color', colors(2, :), 'DisplayName', 'right foot single stance');
        plot(time, right_foot_x_boxcar, 'color', colors(2, :), 'linewidth', 1, 'DisplayName', 'right foot boxcar');
        plot(time, right_foot_x_butter, 'color', colors(3, :), 'linewidth', 1, 'DisplayName', 'right foot filtered');

        plot(time, foot_basis_trajectory, 'linewidth', 4, 'color', colors(4, :), 'DisplayName', 'foot basis');
        legend('Location', 'best')
    end
end
