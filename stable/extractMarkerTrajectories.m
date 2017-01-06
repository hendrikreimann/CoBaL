function extracted_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, specified_marker)
    marker_number = find(strcmp(marker_labels, specified_marker));
    markers_indices = reshape([(marker_number - 1) * 3 + 1; (marker_number - 1) * 3 + 2; (marker_number - 1) * 3 + 3], 1, length(marker_number)*3);
    extracted_trajectory = marker_trajectories(:, markers_indices);
end