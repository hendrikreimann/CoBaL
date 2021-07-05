function foot_basis_trajectory = estimateFootBasisTrajectory_boxcar ...
  ( ...
    time, ...
    marker_trajectories, ...
    marker_labels, ...
    events_data, ...
    sampling_rate, ...
    metronome_frequency ...
  )

    % extract foot markers for each foot
    LHEE = extractMarkerData(marker_trajectories, marker_labels, 'LHEE');
    LANK = extractMarkerData(marker_trajectories, marker_labels, 'LANK');
    LTOE = extractMarkerData(marker_trajectories, marker_labels, 'LTOE');
    LTOEL = extractMarkerData(marker_trajectories, marker_labels, 'LTOEL');

    RHEE = extractMarkerData(marker_trajectories, marker_labels, 'RHEE');
    RANK = extractMarkerData(marker_trajectories, marker_labels, 'RANK');
    RTOE = extractMarkerData(marker_trajectories, marker_labels, 'RTOE');
    RTOEL = extractMarkerData(marker_trajectories, marker_labels, 'RTOEL');

    left_foot = (LHEE + LANK + LTOE + LTOEL) * 0.25;
%     left_foot = (LHEE + LANK) * 0.25;
%     left_foot = (LTOE + LTOEL) * 0.5;
    right_foot = (RHEE + RANK + RTOE + RTOEL) * 0.25;
    left_foot_x = left_foot(:, 1);
    right_foot_x = right_foot(:, 1);

    % remove swing time
    left_heel_strikes = events_data.event_data{strcmp(events_data.event_labels, 'left_touchdown')};
    right_heel_strikes = events_data.event_data{strcmp(events_data.event_labels, 'right_touchdown')};
    left_pushoffs = events_data.event_data{strcmp(events_data.event_labels, 'left_pushoff')};
    right_pushoffs = events_data.event_data{strcmp(events_data.event_labels, 'right_pushoff')};
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
    
    colors = lines(4);
    figure; hold on;
    plot(time, left_foot_x, ':', 'color', colors(1, :));
    plot(time, left_foot_x_stance, 'linewidth', 4, 'color', colors(1, :));
    plot(time, left_foot_x_boxcar, 'color', colors(1, :), 'linewidth', 1);
    plot(time, left_foot_x_butter, 'color', colors(3, :), 'linewidth', 1);
    plot(time, right_foot_x, ':', 'color', colors(2, :));
    plot(time, right_foot_x_stance, 'linewidth', 4, 'color', colors(2, :));
    plot(time, right_foot_x_boxcar, 'color', colors(2, :), 'linewidth', 1);
    plot(time, right_foot_x_butter, 'color', colors(3, :), 'linewidth', 1);

    plot(time, foot_basis_trajectory, 'linewidth', 4, 'color', colors(4, :));
end
