
transform_force_plate_data_to_Vcw = 1;
plot_force_vs_position  = 1;
plot_cop_vs_position    = 0;




load subjectInfo.mat;

load(makeFileName(date, subject_id, 'walking', 4, 'markerTrajectories'));
load(makeFileName(date, subject_id, 'walking', 4, 'forcePlateData'));

left_heel_marker = 34;
left_toes_marker = 35;
right_heel_marker = 42;
right_toes_marker = 43;
left_heel_marker_indices = reshape([(left_heel_marker - 1) * 3 + 1; (left_heel_marker - 1) * 3 + 2; (left_heel_marker - 1) * 3 + 3], 1, length(left_heel_marker)*3);
left_toes_marker_indices = reshape([(left_toes_marker - 1) * 3 + 1; (left_toes_marker - 1) * 3 + 2; (left_toes_marker - 1) * 3 + 3], 1, length(left_toes_marker)*3);
right_heel_marker_indices = reshape([(right_heel_marker - 1) * 3 + 1; (right_heel_marker - 1) * 3 + 2; (right_heel_marker - 1) * 3 + 3], 1, length(right_heel_marker)*3);
right_toes_marker_indices = reshape([(right_toes_marker - 1) * 3 + 1; (right_toes_marker - 1) * 3 + 2; (right_toes_marker - 1) * 3 + 3], 1, length(right_toes_marker)*3);
left_heel_marker_trajectory = marker_trajectories(:, left_heel_marker_indices);
left_toes_marker_trajectory = marker_trajectories(:, left_toes_marker_indices);
right_heel_marker_trajectory = marker_trajectories(:, right_heel_marker_indices);
right_toes_marker_trajectory = marker_trajectories(:, right_toes_marker_indices);

force_scaler = 3e-4;
copxl_trajectory_relevant = copxl_trajectory; copxl_trajectory_relevant(copxl_trajectory_relevant==0) = NaN;
copyl_trajectory_relevant = copyl_trajectory; copyl_trajectory_relevant(copyl_trajectory_relevant==0) = NaN;
copxr_trajectory_relevant = copxr_trajectory; copxr_trajectory_relevant(copxr_trajectory_relevant==0) = NaN;
copyr_trajectory_relevant = copyr_trajectory; copyr_trajectory_relevant(copyr_trajectory_relevant==0) = NaN;


if transform_force_plate_data_to_Vcw
    % extract and process data
    left_force_plate_wrench = [fxl_trajectory fyl_trajectory fzl_trajectory mxl_trajectory myl_trajectory mzl_trajectory];
    left_force_plate_cop = [copxl_trajectory_relevant copyl_trajectory_relevant zeros(size(copxl_trajectory_relevant))];
    right_force_plate_wrench = [fxr_trajectory fyr_trajectory fzr_trajectory mxr_trajectory myr_trajectory mzr_trajectory];
    right_force_plate_cop = [copxr_trajectory_relevant copyr_trajectory_relevant zeros(size(copxr_trajectory_relevant))];
    
    % apply forceplate rotation
    right_forceplate_to_world_rotation = [-1 0 0; 0 1 0; 0 0 -1];
    right_force_plate_to_world_translation = [0.5588; 0; 0];
    right_force_plate_world_trafo = [right_forceplate_to_world_rotation right_force_plate_to_world_translation; 0 0 0 1];
    left_forceplate_to_world_rotation = [-1 0 0; 0 1 0; 0 0 -1];
    left_force_plate_to_world_translation = [-0.5588; 0; 0];
    left_force_plate_world_trafo = [left_forceplate_to_world_rotation left_force_plate_to_world_translation; 0 0 0 1];

    % transform
    right_force_plate_to_origin_adjoint = rigidToAdjointTransformation(right_force_plate_world_trafo);
    right_force_plate_wrench_origin = (right_force_plate_to_origin_adjoint' * right_force_plate_wrench')';
    right_force_plate_cop_origin = (eye(3, 4) * right_force_plate_world_trafo * [right_force_plate_cop ones(size(right_force_plate_cop, 1), 1)]')';
    left_force_plate_to_origin_adjoint = rigidToAdjointTransformation(left_force_plate_world_trafo);
    left_force_plate_wrench_origin = (left_force_plate_to_origin_adjoint' * left_force_plate_wrench')';
    left_force_plate_cop_origin = (eye(3, 4) * left_force_plate_world_trafo * [left_force_plate_cop ones(size(left_force_plate_cop, 1), 1)]')';
end


if plot_force_vs_position
    % left foot force vs. heel/toes
    figure; axes; hold on
    plot(time_force_plate, fzl_trajectory*force_scaler)
    plot(time_mocap, left_heel_marker_trajectory(:, 3));
    plot(time_mocap, left_toes_marker_trajectory(:, 3));

    % right foot force vs. heel
    figure; axes; hold on
    plot(time_force_plate, fzr_trajectory*force_scaler)
    plot(time_mocap, right_heel_marker_trajectory(:, 3));
    plot(time_mocap, right_toes_marker_trajectory(:, 3));
end

if plot_cop_vs_position
    % left foot cop vs. heel/toes
    figure; axes; hold on
    plot(time_mocap, left_heel_marker_trajectory(:, 1));
    plot(time_force_plate, copxl_trajectory_relevant)
    plot(time_force_plate, left_force_plate_cop_origin(:, 1))
    
    figure; axes; hold on
    plot(time_mocap, left_heel_marker_trajectory(:, 2));
    plot(time_force_plate, copyl_trajectory_relevant)
    plot(time_force_plate, left_force_plate_cop_origin(:, 2))
    
    % left foot cop vs. heel/toes
    figure; axes; hold on
    plot(time_mocap, right_heel_marker_trajectory(:, 1));
    plot(time_force_plate, copxr_trajectory_relevant)
    plot(time_force_plate, right_force_plate_cop_origin(:, 1))
    
    figure; axes; hold on
    plot(time_mocap, right_heel_marker_trajectory(:, 2));
    plot(time_force_plate, copyr_trajectory_relevant)
    plot(time_force_plate, right_force_plate_cop_origin(:, 2))
end



















