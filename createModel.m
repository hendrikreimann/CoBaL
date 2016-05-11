
estimate_joint_cors_with_score          = 1;
estimate_hip_centers_from_landmarks     = 1;
create_plant                            = 1;
show_visualization                      = 1;

static_calibration_file_index   = 1;
left_hip_calibration_file_index = 2;
right_hip_calibration_file_index = 3;
left_knee_calibration_file_index = 4;
right_knee_calibration_file_index = 5;


load subjectInfo.mat;

%% create static reference

% load static reference file
static_reference_file_name = makeFileName(date, subject_id, 'calibration', static_calibration_file_index, 'markerTrajectories');
load(static_reference_file_name);

i_time = 1;
while any(isnan(marker_trajectories(i_time, :)))
    i_time = i_time + 1;
end
marker_reference = marker_trajectories(i_time, :);

head_markers = [1 2 3 4];
trunk_markers = [5 6 7 8 9];
pelvis_markers = [24 25 26 27];
left_thigh_markers = [28 29 30];
left_shank_markers = [31 32 33];
left_foot_markers = [34 35];
right_thigh_markers = [36 37 38];
right_shank_markers = [39 40 41];
right_foot_markers = [42 43];
LASIS_marker = 24;
RASIS_marker = 25;
LPSIS_marker = 26;
RPSIS_marker = 27;

% marker indices
head_markers_indices = reshape([(head_markers - 1) * 3 + 1; (head_markers - 1) * 3 + 2; (head_markers - 1) * 3 + 3], 1, length(head_markers)*3);
trunk_markers_indices = reshape([(trunk_markers - 1) * 3 + 1; (trunk_markers - 1) * 3 + 2; (trunk_markers - 1) * 3 + 3], 1, length(trunk_markers)*3);
pelvis_markers_indices = reshape([(pelvis_markers - 1) * 3 + 1; (pelvis_markers - 1) * 3 + 2; (pelvis_markers - 1) * 3 + 3], 1, length(pelvis_markers)*3);
left_thigh_markers_indices = reshape([(left_thigh_markers - 1) * 3 + 1; (left_thigh_markers - 1) * 3 + 2; (left_thigh_markers - 1) * 3 + 3], 1, length(left_thigh_markers)*3);
left_shank_markers_indices = reshape([(left_shank_markers - 1) * 3 + 1; (left_shank_markers - 1) * 3 + 2; (left_shank_markers - 1) * 3 + 3], 1, length(left_shank_markers)*3);
left_foot_markers_indices = reshape([(left_foot_markers - 1) * 3 + 1; (left_foot_markers - 1) * 3 + 2; (left_foot_markers - 1) * 3 + 3], 1, length(left_foot_markers)*3);
right_thigh_markers_indices = reshape([(right_thigh_markers - 1) * 3 + 1; (right_thigh_markers - 1) * 3 + 2; (right_thigh_markers - 1) * 3 + 3], 1, length(right_thigh_markers)*3);
right_shank_markers_indices = reshape([(right_shank_markers - 1) * 3 + 1; (right_shank_markers - 1) * 3 + 2; (right_shank_markers - 1) * 3 + 3], 1, length(right_shank_markers)*3);
right_foot_markers_indices = reshape([(right_foot_markers - 1) * 3 + 1; (right_foot_markers - 1) * 3 + 2; (right_foot_markers - 1) * 3 + 3], 1, length(right_foot_markers)*3);
LASIS_markers_indices = reshape([(LASIS_marker - 1) * 3 + 1; (LASIS_marker - 1) * 3 + 2; (LASIS_marker - 1) * 3 + 3], 1, length(LASIS_marker)*3);
RASIS_markers_indices = reshape([(RASIS_marker - 1) * 3 + 1; (RASIS_marker - 1) * 3 + 2; (RASIS_marker - 1) * 3 + 3], 1, length(RASIS_marker)*3);
LPSIS_markers_indices = reshape([(LPSIS_marker - 1) * 3 + 1; (LPSIS_marker - 1) * 3 + 2; (LPSIS_marker - 1) * 3 + 3], 1, length(LPSIS_marker)*3);
RPSIS_markers_indices = reshape([(RPSIS_marker - 1) * 3 + 1; (RPSIS_marker - 1) * 3 + 2; (RPSIS_marker - 1) * 3 + 3], 1, length(RPSIS_marker)*3);

% marker references
head_markers_reference = marker_reference(head_markers_indices);
trunk_markers_reference = marker_reference(trunk_markers_indices);
pelvis_markers_reference = marker_reference(pelvis_markers_indices);
left_thigh_markers_reference = marker_reference(left_thigh_markers_indices);
left_shank_markers_reference = marker_reference(left_shank_markers_indices);
left_foot_markers_reference = marker_reference(left_foot_markers_indices);
right_thigh_markers_reference = marker_reference(right_thigh_markers_indices);
right_shank_markers_reference = marker_reference(right_shank_markers_indices);
right_foot_markers_reference = marker_reference(right_foot_markers_indices);
LASIS_marker_reference = marker_reference(LASIS_markers_indices);
RASIS_marker_reference = marker_reference(RASIS_markers_indices);
LPSIS_marker_reference = marker_reference(LPSIS_markers_indices);
RPSIS_marker_reference = marker_reference(RPSIS_markers_indices);

%% estimate_joint_cors
if estimate_joint_cors_with_score
    % find CoRs
    pelvis_center_reference = mean(reshape(pelvis_markers_reference, 3, size(pelvis_markers_reference, 2)/3), 2);

%     % find left hip CoR
%     left_hip_reference_file_name = makeFileName(date, subject_id, 'calibration', left_hip_calibration_file_index, 'markerTrajectories');
%     disp(['Left hip reference file name: ' left_hip_reference_file_name]);
%     load(left_hip_reference_file_name);
%     hip_reference = marker_trajectories;
% 
%     pelvis_markers_trajectory = hip_reference(:, pelvis_markers_indices);
%     left_thigh_markers_trajectory = hip_reference(:, left_thigh_markers_indices);
%     [left_hip_cor, left_hip_cor_error] = estimateJointKinematics ...
%       ( ...
%         pelvis_markers_reference, ...
%         left_thigh_markers_reference, ...
%         pelvis_markers_trajectory, ...
%         left_thigh_markers_trajectory ...
%       );
% 
%     % find right hip CoR
%     right_hip_reference_file_name = makeFileName(date, subject_id, 'calibration', right_hip_calibration_file_index, 'markerTrajectories');
%     disp(['Right hip reference file name: ' right_hip_reference_file_name]);
%     load(right_hip_reference_file_name);
%     hip_reference = marker_trajectories;
% 
%     pelvis_markers_trajectory = hip_reference(:, pelvis_markers_indices);
%     right_thigh_markers_trajectory = hip_reference(:, right_thigh_markers_indices);
%     [right_hip_cor, right_hip_cor_error] = estimateJointKinematics ...
%       ( ...
%         pelvis_markers_reference, ...
%         right_thigh_markers_reference, ...
%         pelvis_markers_trajectory, ...
%         right_thigh_markers_trajectory ...
%       );

    % find left knee CoR
    left_knee_reference_file_name = makeFileName(date, subject_id, 'calibration', left_knee_calibration_file_index, 'markerTrajectories');
    disp(['Left knee reference file name: ' left_knee_reference_file_name]);
    load(left_knee_reference_file_name);
    knee_reference = marker_trajectories;

    left_thigh_markers_trajectory = knee_reference(:, left_thigh_markers_indices);
    left_shank_markers_trajectory = knee_reference(:, left_shank_markers_indices);
    [left_knee_point, left_knee_cor_error, left_knee_aor] = estimateJointKinematics ...
      ( ...
        left_thigh_markers_reference, ...
        left_shank_markers_reference, ...
        left_thigh_markers_trajectory, ...
        left_shank_markers_trajectory, ...
        1 ...
      );

    % find right knee CoR
    right_knee_reference_file_name = makeFileName(date, subject_id, 'calibration', right_knee_calibration_file_index, 'markerTrajectories');
    disp(['Left knee reference file name: ' right_knee_reference_file_name]);
    load(right_knee_reference_file_name);
    knee_reference = marker_trajectories;

    right_thigh_markers_trajectory = knee_reference(:, right_thigh_markers_indices);
    right_shank_markers_trajectory = knee_reference(:, right_shank_markers_indices);
    [right_knee_point, right_knee_cor_error, right_knee_aor] = estimateJointKinematics ...
      ( ...
        right_thigh_markers_reference, ...
        right_shank_markers_reference, ...
        right_thigh_markers_trajectory, ...
        right_shank_markers_trajectory, ...
        1 ...
      );

    left_ankle_cor = marker_reference(left_shank_markers_indices(7:9))';
    right_ankle_cor = marker_reference(right_shank_markers_indices(7:9))';
    
    left_knee_cor = left_knee_point;
    right_knee_cor = right_knee_point;
end

if estimate_hip_centers_from_landmarks
    
end

%% estimate_hip_centers_from_landmarks
if estimate_hip_centers_from_landmarks
    anterior = [0; 1; 0];
    posterior = - anterior;
    leftward = [-1; 0; 0];
    rightward = - leftward;
    proximal = [0; 0; 1];
    distal = - proximal;
    centroid_to_skin_correction = 0.0152; % in meters
    skin_to_bone_correction = 0.01; % in meters
    centroid_to_bone_correction = centroid_to_skin_correction + skin_to_bone_correction;
    LASIS_position_bone = LASIS_marker_reference' + skin_to_bone_correction * posterior;
    RASIS_position_bone = RASIS_marker_reference' + skin_to_bone_correction * posterior;
    inter_ASIS_distance = norm(LASIS_position_bone - RASIS_position_bone);
    left_hip_cor = LASIS_position_bone ...
                    + 0.11 * inter_ASIS_distance * rightward ...
                    + 0.12 * inter_ASIS_distance * distal ...
                    + 0.21 * inter_ASIS_distance * posterior;
    right_hip_cor = RASIS_position_bone ...
                    + 0.11 * inter_ASIS_distance * leftward ...
                    + 0.12 * inter_ASIS_distance * distal ...
                    + 0.21 * inter_ASIS_distance * posterior;
                
end

%% create geometric model
if create_plant
    plant = walkingModel ...
      ( ...
        marker_reference, ...
        {left_knee_aor, right_knee_aor}, ...
        weight, ...
        gender ...
      );
    plant.updateInternals();
    
    model_file_name = makeFileName(date, subject_id, 'model');
    save( ...
         model_file_name, ...
         'plant', ...
         'marker_reference', ...
         'left_hip_cor', ...
         'left_knee_cor', ...
         'left_ankle_cor', ...
         'right_hip_cor', ...
         'right_knee_cor', ...
         'right_ankle_cor', ...
         'left_knee_aor', ...
         'right_knee_aor' ...
        );    
end

%% show visualization
if show_visualization

    % show stick figure of geometric model
    hip_center = (right_hip_cor + left_hip_cor) * 0.5;
    scene_bound = repmat(hip_center, 1, 2) + 2*[-0.5 0.5; -0.5 0.5; -1 1];
    stick_figure = KinematicTreeController(plant, scene_bound, 'ellipsoid');
end







return


red = [1 0 0];
green = [0 1 0];
blue = [0 0 1];
figure; axes; hold on; axis equal
xlabel('x'); ylabel('y'); zlabel('z');
plot3(marker_reference(01), marker_reference(02), marker_reference(03), 'o', 'color', blue)
plot3(marker_reference(04), marker_reference(05), marker_reference(06), 'o', 'color', blue)
plot3(marker_reference(07), marker_reference(08), marker_reference(09), 'o', 'color', blue)
plot3(marker_reference(10), marker_reference(11), marker_reference(12), 'o', 'color', blue)
plot3(marker_reference(13), marker_reference(14), marker_reference(15), 'o', 'color', blue)
plot3(marker_reference(16), marker_reference(17), marker_reference(18), 'o', 'color', blue)
plot3(marker_reference(19), marker_reference(20), marker_reference(21), 'o', 'color', blue)
plot3(marker_reference(22), marker_reference(23), marker_reference(24), 'o', 'color', blue)
plot3(marker_reference(25), marker_reference(26), marker_reference(27), 'o', 'color', blue)

plot3(marker_reference(28), marker_reference(29), marker_reference(30), 'o', 'color', red)
plot3(marker_reference(31), marker_reference(32), marker_reference(33), 'o', 'color', red)
plot3(marker_reference(34), marker_reference(35), marker_reference(36), 'o', 'color', red)
plot3(marker_reference(37), marker_reference(38), marker_reference(39), 'o', 'color', red)
plot3(marker_reference(40), marker_reference(41), marker_reference(42), 'o', 'color', red)
plot3(marker_reference(43), marker_reference(44), marker_reference(45), 'o', 'color', red)
plot3(marker_reference(46), marker_reference(47), marker_reference(48), 'o', 'color', red)

plot3(marker_reference(49), marker_reference(50), marker_reference(51), 'o', 'color', green)
plot3(marker_reference(52), marker_reference(53), marker_reference(54), 'o', 'color', green)
plot3(marker_reference(55), marker_reference(56), marker_reference(57), 'o', 'color', green)
plot3(marker_reference(58), marker_reference(59), marker_reference(60), 'o', 'color', green)
plot3(marker_reference(61), marker_reference(62), marker_reference(63), 'o', 'color', green)
plot3(marker_reference(64), marker_reference(65), marker_reference(66), 'o', 'color', green)
plot3(marker_reference(67), marker_reference(68), marker_reference(69), 'o', 'color', green)

plot3(marker_reference(70), marker_reference(71), marker_reference(72), 'o', 'color', red)
plot3(marker_reference(73), marker_reference(74), marker_reference(75), 'o', 'color', green)
plot3(marker_reference(76), marker_reference(77), marker_reference(78), 'o', 'color', red)
plot3(marker_reference(79), marker_reference(80), marker_reference(81), 'o', 'color', green)

plot3(marker_reference(82), marker_reference(83), marker_reference(84), 'o', 'color', red)
plot3(marker_reference(85), marker_reference(86), marker_reference(87), 'o', 'color', red)
plot3(marker_reference(88), marker_reference(89), marker_reference(90), 'o', 'color', red)
plot3(marker_reference(91), marker_reference(92), marker_reference(93), 'o', 'color', red)
plot3(marker_reference(94), marker_reference(95), marker_reference(96), 'o', 'color', red)
plot3(marker_reference(97), marker_reference(98), marker_reference(99), 'o', 'color', red)
plot3(marker_reference(100), marker_reference(101), marker_reference(102), 'o', 'color', red)
plot3(marker_reference(103), marker_reference(104), marker_reference(105), 'o', 'color', red)

plot3(marker_reference(106), marker_reference(107), marker_reference(108), 'o', 'color', green)
plot3(marker_reference(109), marker_reference(110), marker_reference(111), 'o', 'color', green)
plot3(marker_reference(112), marker_reference(113), marker_reference(114), 'o', 'color', green)
plot3(marker_reference(115), marker_reference(116), marker_reference(117), 'o', 'color', green)
plot3(marker_reference(118), marker_reference(119), marker_reference(120), 'o', 'color', green)
plot3(marker_reference(121), marker_reference(122), marker_reference(123), 'o', 'color', green)
plot3(marker_reference(124), marker_reference(125), marker_reference(126), 'o', 'color', green)
plot3(marker_reference(127), marker_reference(128), marker_reference(129), 'o', 'color', green)

knee_axis_length = [-0.1 0.1];
plot3(left_hip_cor(1), left_hip_cor(2), left_hip_cor(3), 'x', 'color', 'k')
plot3(right_hip_cor(1), right_hip_cor(2), right_hip_cor(3), 'x', 'color', 'k')
plot3(left_knee_point(1), left_knee_point(2), left_knee_point(3), 'x', 'color', 'k')
plot3(right_knee_point(1), right_knee_point(2), right_knee_point(3), 'x', 'color', 'k')
plot3(left_knee_point(1) + left_knee_aor(1)*knee_axis_length, left_knee_point(2) + left_knee_aor(2)*knee_axis_length, left_knee_point(3) + left_knee_aor(3)*knee_axis_length, '-', 'color', 'k')
plot3(right_knee_point(1) + right_knee_aor(1)*knee_axis_length, right_knee_point(2) + right_knee_aor(2)*knee_axis_length, right_knee_point(3) + right_knee_aor(3)*knee_axis_length, '-', 'color', 'k')





















