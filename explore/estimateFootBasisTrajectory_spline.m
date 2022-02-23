function foot_basis_trajectory = estimateFootBasisTrajectory_spline ...
  ( ...
    time, ...
    marker_trajectories, ...
    marker_labels, ...
    events_data, ...
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
    left_foot = (LHEE + LANK) * 0.25;
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
    for i_pushoff = 1 : length(left_pushoff_indices)
        this_pushoff_time_index = left_pushoff_indices(i_pushoff);
        next_heelstrike_time_index = min(left_heel_strike_indices(left_heel_strike_indices>this_pushoff_time_index));
        left_foot_x_stance(this_pushoff_time_index : next_heelstrike_time_index) = NaN;
    end


    figure; hold on;
    plot(time, left_foot_x);
    plot(time, left_foot_x_stance, 'linewidth', 3);

end
