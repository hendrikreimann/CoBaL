%     This file is part of the CoBaL code base
%     Copyright (C) 2017 Hendrik Reimann <hendrikreimann@gmail.com>
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.% compare the kinematic tree against the kinematic chain



function joint_angles = ...
    markerToAngles ...
      ( ...
        marker_reference, ...
        marker_current, ...
        marker_headers, ...
        joint_center_reference, ...
        joint_center_current, ...
        joint_center_headers, ...
        body_direction_matrix, ...
        left_hip_direction_matrix, ...
        right_hip_direction_matrix, ...
        left_knee_direction_matrix, ...
        right_knee_direction_matrix ...
      );
    
    joint_angles = zeros(40, 1);
    
    body_direction_matrix_inverse = body_direction_matrix^(-1);
    left_hip_direction_matrix_inverse = left_hip_direction_matrix^(-1);
    right_hip_direction_matrix_inverse = right_hip_direction_matrix^(-1);
    left_knee_direction_matrix_inverse = left_knee_direction_matrix^(-1);
    right_knee_direction_matrix_inverse = right_knee_direction_matrix^(-1);
    
    % extract positions
    LASI_current = extractMarkerTrajectories(marker_current, marker_headers, 'LASI')';
    RASI_current = extractMarkerTrajectories(marker_current, marker_headers, 'RASI')';
    LPSI_current = extractMarkerTrajectories(marker_current, marker_headers, 'LPSI')';
    RPSI_current = extractMarkerTrajectories(marker_current, marker_headers, 'RPSI')';
    LKNE_current = extractMarkerTrajectories(marker_current, marker_headers, 'LKNE')';
    RKNE_current = extractMarkerTrajectories(marker_current, marker_headers, 'RKNE')';
    left_hip_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LHIPCOR')';
    left_knee_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LKNEECOR')';
    left_ankle_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LANKLECOR')';
    right_hip_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'RHIPCOR')';
    right_knee_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'RKNEECOR')';
    right_ankle_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'RANKLECOR')';
    lumbar_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LUMBARCOR')';
    cervix_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'CERVIXCOR')';

    LASI_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LASI')';
    RASI_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RASI')';
    LPSI_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LPSI')';
    RPSI_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RPSI')';
    LKNE_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LKNE')';
    RKNE_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RKNE')';
    left_hip_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LHIPCOR')';
    left_knee_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LKNEECOR')';
    left_ankle_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LANKLECOR')';
    right_hip_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RHIPCOR')';
    right_knee_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RKNEECOR')';
    right_ankle_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RANKLECOR')';
    lumbar_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LUMBARCOR')';
    cervix_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'CERVIXCOR')';
    
    
    % calculate pelvis translation and rotation
    pelvis_transformation_reference_cell = calculateMcsToWcsTransformations_detailed(marker_reference, marker_headers, {'PELVIS'});
    pelvis_transformation_reference = pelvis_transformation_reference_cell{1};
    pelvis_rotation_reference = pelvis_transformation_reference(1:3, 1:3);
    pelvis_translation_reference = pelvis_transformation_reference(1:3, 4);
    
    pelvis_transformation_current_cell = calculateMcsToWcsTransformations_detailed(marker_current, marker_headers, {'PELVIS'});    
    pelvis_transformation_current = pelvis_transformation_current_cell{1};
    pelvis_rotation_current = pelvis_transformation_current(1:3, 1:3);
    pelvis_translation_current = pelvis_transformation_current(1:3, 4);
    
    pelvis_transformation_reference_to_current = pelvis_transformation_reference^(-1) * pelvis_transformation_current;
    pelvis_rotation_reference_to_current = pelvis_transformation_reference_to_current(1:3, 1:3);
    pelvis_translation_reference_to_current = pelvis_translation_current - pelvis_translation_reference;
    
    joint_angles(1:3) = pelvis_translation_reference_to_current;
    joint_angles(4:6) = eulerAnglesFromRotationMatrixZXY(pelvis_rotation_reference_to_current);
    
    % calculate left hip flexion-extension
    left_hip_cor_to_left_knee_cor_current_world = left_knee_cor_current - left_hip_cor_current;
    left_hip_cor_to_left_knee_cor_reference_world = left_knee_cor_reference - left_hip_cor_reference;
    left_hip_cor_to_left_knee_cor_current_hip = left_hip_direction_matrix_inverse * pelvis_rotation_current^(-1) * left_hip_cor_to_left_knee_cor_current_world;
    left_hip_cor_to_left_knee_cor_reference_hip = left_hip_direction_matrix_inverse * pelvis_rotation_reference^(-1) * left_hip_cor_to_left_knee_cor_reference_world;

    left_hip_rotation_phi_current = atan2(left_hip_cor_to_left_knee_cor_current_hip(2), -left_hip_cor_to_left_knee_cor_current_hip(3));
    left_hip_rotation_phi_reference = atan2(left_hip_cor_to_left_knee_cor_reference_hip(2), -left_hip_cor_to_left_knee_cor_reference_hip(3));
    left_hip_rotation_phi_reference_to_current = left_hip_rotation_phi_current - left_hip_rotation_phi_reference;
    joint_angles(7) = left_hip_rotation_phi_reference_to_current;
    
    % left hip abduction-adduction
    R_hip_phi_current = expAxis([1; 0; 0], -left_hip_rotation_phi_current);
    R_hip_phi_reference = expAxis([1; 0; 0], -left_hip_rotation_phi_reference);
    left_hip_cor_to_left_knee_cor_current_phi = R_hip_phi_current * left_hip_cor_to_left_knee_cor_current_hip;
    left_hip_cor_to_left_knee_cor_reference_phi = R_hip_phi_reference * left_hip_cor_to_left_knee_cor_reference_hip;
    left_hip_rotation_psi_current = atan2(-left_hip_cor_to_left_knee_cor_current_phi(1), -left_hip_cor_to_left_knee_cor_current_phi(3));
    left_hip_rotation_psi_reference = atan2(-left_hip_cor_to_left_knee_cor_reference_phi(1), -left_hip_cor_to_left_knee_cor_reference_phi(3));
    left_hip_rotation_psi_reference_to_current = left_hip_rotation_psi_current - left_hip_rotation_psi_reference;
    joint_angles(8) = left_hip_rotation_psi_reference_to_current;
    
    % left hip internal rotation
    R_psi_current = expAxis([0; 1; 0], -left_hip_rotation_psi_current);
    R_psi_reference = expAxis([0; 1; 0], -left_hip_rotation_psi_reference);
    left_knee_cor_to_LKNE_current_world = LKNE_current - left_knee_cor_current;
    left_knee_cor_to_LKNE_reference_world = LKNE_reference - left_knee_cor_reference;
    left_knee_cor_to_LKNE_current_phipsi = R_psi_current * R_hip_phi_current * left_hip_direction_matrix_inverse * pelvis_rotation_current^(-1) * left_knee_cor_to_LKNE_current_world;
    left_knee_cor_to_LKNE_reference_phipsi = R_psi_reference * R_hip_phi_reference * left_hip_direction_matrix_inverse * pelvis_rotation_reference^(-1) * left_knee_cor_to_LKNE_reference_world;
    left_hip_rotation_gamma_current = atan2(-left_knee_cor_to_LKNE_current_phipsi(2), -left_knee_cor_to_LKNE_current_phipsi(1));
    left_hip_rotation_gamma_reference = atan2(-left_knee_cor_to_LKNE_reference_phipsi(2), -left_knee_cor_to_LKNE_reference_phipsi(1));
    left_hip_rotation_gamma_reference_to_current_1 = left_hip_rotation_gamma_current - left_hip_rotation_gamma_reference;
    joint_angles(9) = -left_hip_rotation_gamma_reference_to_current_1;
    
    % left knee flexion-extension
    R_gamma_current = expAxis([0; 0; 1], -left_hip_rotation_gamma_current);
    R_gamma_reference = expAxis([0; 0; 1], -left_hip_rotation_gamma_reference);
    R_hip_current = R_gamma_current * R_psi_current * R_hip_phi_current * left_hip_direction_matrix_inverse;
    R_hip_reference = R_gamma_reference * R_psi_reference * R_hip_phi_reference * left_hip_direction_matrix_inverse;
    
    left_knee_cor_to_left_ankle_cor_current_world = left_ankle_cor_current - left_knee_cor_current;
    left_knee_cor_to_left_ankle_cor_reference_world = left_ankle_cor_reference - left_knee_cor_reference;
    left_knee_cor_to_left_ankle_cor_current_knee = left_knee_direction_matrix_inverse * R_hip_current * pelvis_rotation_current^(-1) * left_knee_cor_to_left_ankle_cor_current_world;
    left_knee_cor_to_left_ankle_cor_reference_knee = left_knee_direction_matrix_inverse * R_hip_reference * pelvis_rotation_reference^(-1) * left_knee_cor_to_left_ankle_cor_reference_world;
    left_knee_rotation_current = atan2(-left_knee_cor_to_left_ankle_cor_current_knee(2), -left_knee_cor_to_left_ankle_cor_current_knee(3));
    left_knee_rotation_reference = atan2(-left_knee_cor_to_left_ankle_cor_reference_knee(2), -left_knee_cor_to_left_ankle_cor_reference_knee(3));
    left_knee_rotation_reference_to_current = left_knee_rotation_current - left_knee_rotation_reference;
    joint_angles(10) = left_knee_rotation_reference_to_current;
    
    % experimental: correct internal rotation angle to bring the ankle of the kinematic chain CoR closer to the previously estimated one
    R_knee_phi_current = expAxis([-1; 0; 0], -left_knee_rotation_current);
    left_knee_cor_to_left_hip_cor_current_knee_phi = - R_knee_phi_current * left_knee_direction_matrix_inverse * R_hip_current * pelvis_rotation_current^(-1) * left_hip_cor_to_left_knee_cor_current_world;
    % this is the direction from the knee to the hip, in current knee coordinates
    
    left_knee_cor_to_left_ankle_cor_current_knee_phi = R_knee_phi_current * left_knee_direction_matrix_inverse * R_hip_current * pelvis_rotation_current^(-1) * left_knee_cor_to_left_ankle_cor_current_world;
    % this is how the ankle looks like from the knee
    
    left_knee_cor_to_left_ankle_cor_reference_knee_phi = left_knee_direction_matrix_inverse * left_knee_cor_to_left_ankle_cor_reference_world;
    % this is how the ankle is supposed to look like from the knee
    
    left_knee_cor_to_left_hip_cor_reference_knee_phi = - left_knee_direction_matrix_inverse * left_hip_cor_to_left_knee_cor_current_world;
    % this is how the hip looks like from the knee
    
    
    
    
    
    
    
    left_knee_cor_to_left_hip_cor_current_knee = - R_knee_phi_current * left_knee_direction_matrix_inverse * R_hip_current * pelvis_rotation_current^(-1) * left_hip_cor_to_left_knee_cor_current_world;
    % this is the direction from the knee to the hip, in current flexed knee coordinates
    
    left_knee_cor_to_left_ankle_cor_current_knee_phi = R_knee_phi_current * left_knee_direction_matrix_inverse * R_hip_current * pelvis_rotation_current^(-1) * left_knee_cor_to_left_ankle_cor_current_world;
    % this is the direction from the knee to the ankle, in current flexed knee coordinates
    
    left_knee_cor_to_left_hip_cor_reference_knee_phi = - left_knee_direction_matrix_inverse * left_hip_cor_to_left_knee_cor_reference_world; % not needed
    left_knee_cor_to_left_ankle_cor_reference_knee_phi = left_knee_direction_matrix_inverse * left_knee_cor_to_left_ankle_cor_reference_world;
    % this is how the ankle is supposed to look like from the knee
    
    
    
    
    hip_knee_ankle_current_normal = normVector(cross(left_knee_cor_to_left_hip_cor_current_knee, left_knee_cor_to_left_ankle_cor_current_knee_phi));
    hip_knee_ankle_reference_normal = normVector(cross(left_knee_cor_to_left_hip_cor_current_knee, left_knee_cor_to_left_ankle_cor_reference_knee_phi));
    
    
    
    
    
    
    
    % new approach, try thinking from the hip
    left_knee_cor_to_left_hip_cor_current_hip = - R_hip_current * pelvis_rotation_current^(-1) * left_hip_cor_to_left_knee_cor_current_world;
    left_knee_cor_to_left_ankle_cor_current_hip = - R_hip_current * pelvis_rotation_current^(-1) * left_knee_cor_to_left_ankle_cor_current_world;
    
    % the above two lines tell me where the estimated ankle CoR is, seen from the hip. Now I need to find out where the
    % one in the kinematic chain is
    load subjectModel.mat;
    kinematic_tree.jointAngles = joint_angles;
    kinematic_tree.updateConfiguration;

    % try to calcualate left_knee_cor_to_left_ankle_cor_current_world_from_tree as just a rotation of a vector
    left_knee_cor_to_left_ankle_cor_reference_world_from_tree = kinematic_tree.referenceJointTransformations{13}(:, 4) - kinematic_tree.referenceJointTransformations{11}(:, 4);
    T_knee_from_tree = kinematic_tree.productsOfExponentials{11};
    left_knee_cor_to_left_ankle_cor_current_world_from_tree = T_knee_from_tree * left_knee_cor_to_left_ankle_cor_reference_world_from_tree;
    
%     left_knee_cor_from_tree_current = kinematic_tree.jointTransformations{10}(1:3, 4);
%     left_ankle_cor_from_tree_current = kinematic_tree.jointTransformations{12}(1:3, 4);
%     left_knee_cor_to_left_ankle_cor_current_world_from_tree_check = left_ankle_cor_from_tree_current - left_knee_cor_from_tree_current;
    
    left_knee_cor_to_left_ankle_cor_reference_world = left_ankle_cor_reference - left_knee_cor_reference;
    R_knee_from_tree = kinematic_tree.productsOfExponentials{11}(1:3, 1:3);
    left_knee_cor_to_left_ankle_cor_current_world_from_tree_2 = R_knee_from_tree * left_knee_cor_to_left_ankle_cor_reference_world;
    
    % how do we get R_knee without using the tree?
    pelvis_rotation_from_tree = kinematic_tree.productsOfExponentials{6};
    hip_rotation_from_tree = kinematic_tree.twistExponentials{7} * kinematic_tree.twistExponentials{8} * kinematic_tree.twistExponentials{9};
    knee_rotation_from_tree = kinematic_tree.twistExponentials{10} * kinematic_tree.twistExponentials{11};
    R_knee_from_tree_2 = pelvis_rotation_from_tree * hip_rotation_from_tree * knee_rotation_from_tree;
    
    pelvis_rotation = pelvis_transformation_reference_to_current;
    
    kinematic_tree.twistExponentials{7}
%     kinematic_tree.twistExponentials{8}
%     kinematic_tree.twistExponentials{9}
    left_hip_direction_matrix * R_hip_phi_current * left_hip_direction_matrix_inverse
    (R_hip_phi_current * left_hip_direction_matrix_inverse)^(-1)
    
%     hip_rotation_from_tree
%     R_hip_current
%     R_hip_reference
    
    
    
    % take the knee-to-ankle in world frame and transform it into hip frame
    left_knee_cor_to_left_ankle_cor_current_hip_from_tree = - R_hip_current * pelvis_rotation_current^(-1) * left_knee_cor_to_left_ankle_cor_current_world_from_tree;
    
    % now, by how much do I have to rotate around the z-axis of the hip frame, i.e. internal rotation, to get
    % left_knee_cor_to_left_ankle_cor_current_hip_from_tree onto left_knee_cor_to_left_ankle_cor_current_hip
    hip_internal_rotation_correction = subspace([[0; 0; 1] left_knee_cor_to_left_ankle_cor_current_hip], [[0; 0; 1] left_knee_cor_to_left_ankle_cor_current_hip_from_tree]);
    
%     hip_internal_rotation_correction = atan2(left_knee_cor_to_left_ankle_cor_current_knee_phi(1), -left_knee_cor_to_left_ankle_cor_current_knee_phi(3));
    joint_angles(9) = joint_angles(9) + hip_internal_rotation_correction;
    
    
    % this works! now, figure out how to calculate the above without having to update the kinematic tree object
    
    return
    
    
    
    
    
    
    
%     R_z = expAxis([0; 0; 1], theta_456(1));
%     R_x = expAxis([1; 0; 0], theta_456(2));
%     R_y = expAxis([0; 1; 0], theta_456(3));
%     
%     R_zxy = R_y * R_x * R_z;
%     
%     
%     load('subjectModel.mat')
%     pelvis_reference_trafo = kinematic_tree.referenceJointTransformations{6};
%     kinematic_tree.jointAngles(4:6) = theta_456;
%     kinematic_tree.updateConfiguration;
%     
%     
%     T_456 = kinematic_tree.twistExponentials{4} * kinematic_tree.twistExponentials{5} * kinematic_tree.twistExponentials{6};
%     T_123 = pelvis_transformation_reference_to_current * (T_456 * pelvis_reference_trafo)^(-1);
%     theta_123 = T_123(1:3, 4);
    
    % Probe
%     T_123 * T_456 * pelvis_reference_trafo
    
%     kinematic_tree.jointAngles(1:6) = [theta_123; theta_456];
%     kinematic_tree.updateConfiguration;
%     kinematic_tree.jointTransformations{6}
%     pelvis_transformation_reference_to_current
    
%     joint_angles(1:6) = [theta_123; theta_456];
%     joint_angles(1:6) = [pelvis_translation_reference_to_current; theta_456];
    
    
    
    
    
    
    

    % calculate pelvis translation
    pelvis_base_point_current = mean([LASI_current, RASI_current, LPSI_current, RPSI_current], 2);
    pelvis_base_point_reference = mean([LASI_reference, RASI_reference, LPSI_reference, RPSI_reference], 2);
    pelvis_base_point_translation = pelvis_base_point_current - pelvis_base_point_reference;
    joint_angles(1:3) = pelvis_base_point_translation;
    
    % pelvis rotation in horizontal plane
    MASIS_current = (LASI_current + RASI_current) / 2;
    MASIS_reference = (LASI_reference + RASI_reference) / 2;
    MPSIS_current = (LPSI_current + RPSI_current) / 2;
    MPSIS_reference = (LPSI_reference + RPSI_reference) / 2;
    MPSIS_to_MASIS_current = MASIS_current - MPSIS_current;
    MPSIS_to_MPSIS_reference = MASIS_reference - MPSIS_reference;
    pelvis_rotation_phi_current = -atan2(MPSIS_to_MASIS_current(1), MPSIS_to_MASIS_current(2));
    pelvis_rotation_phi_reference = -atan2(MPSIS_to_MPSIS_reference(1), MPSIS_to_MPSIS_reference(2));
    pelvis_rotation_phi_reference_to_current = pelvis_rotation_phi_current - pelvis_rotation_phi_reference;
    
    % pelvis rotation around anterior-posterior axis
    R_phi_current = expAxis([0; 0; 1], -pelvis_rotation_phi_current);
    R_phi_reference = expAxis([0; 0; 1], -pelvis_rotation_phi_reference);
    MPSIS_to_MASIS_phi_current = R_phi_current * MPSIS_to_MASIS_current;
    MPSIS_to_MASIS_phi_reference = R_phi_reference * MPSIS_to_MPSIS_reference;
    pelvis_rotation_psi_current = atan2(MPSIS_to_MASIS_phi_current(3), MPSIS_to_MASIS_phi_current(2));
    pelvis_rotation_psi_reference = atan2(MPSIS_to_MASIS_phi_reference(3), MPSIS_to_MASIS_phi_reference(2));
    pelvis_rotation_psi_reference_to_current = pelvis_rotation_psi_current - pelvis_rotation_psi_reference;
    
    % pelvis rotation around medial-lateral axis
    pelvis_left_current = (LASI_current + LPSI_current) / 2;
    pelvis_right_current = (RASI_current + RPSI_current) / 2;
    pelvis_left_to_right_current = pelvis_right_current - pelvis_left_current;
    pelvis_left_reference = (LASI_reference + LPSI_reference) / 2;
    pelvis_right_reference = (RASI_reference + RPSI_reference) / 2;
    pelvis_left_to_right_reference = pelvis_right_reference - pelvis_left_reference;
    
    R_psi_current = expAxis([1; 0; 0], -pelvis_rotation_psi_current);
    R_psi_reference = expAxis([1; 0; 0], -pelvis_rotation_psi_reference);
    pelvis_left_to_right_phipsi_current = R_psi_current * R_phi_current * pelvis_left_to_right_current;
    pelvis_left_to_right_phipsi_reference = R_psi_reference * R_phi_reference * pelvis_left_to_right_reference;
    pelvis_rotation_theta_current = atan2(pelvis_left_to_right_phipsi_current(3), pelvis_left_to_right_phipsi_current(1));
    pelvis_rotation_theta_reference = atan2(pelvis_left_to_right_phipsi_reference(3), pelvis_left_to_right_phipsi_reference(1));
    pelvis_rotation_theta_reference_to_current = pelvis_rotation_theta_current - pelvis_rotation_theta_reference;    
    
%     % store calculated angles
%     joint_angles(4) = pelvis_rotation_phi_reference_to_current;
%     joint_angles(5) = pelvis_rotation_psi_reference_to_current;
%     joint_angles(6) = pelvis_rotation_theta_reference_to_current;

    % calculate pelvis rotation matrix
    R_z = expAxis([0; 0; 1], theta_456(1));
    R_x = expAxis([1; 0; 0], theta_456(2));
    R_y = expAxis([0; 1; 0], theta_456(3));
    R_zxy = R_y * R_x * R_z;

%     R_theta_current = expAxis([1; 0; 0], -pelvis_rotation_theta_current);
%     R_theta_reference = expAxis([1; 0; 0], -pelvis_rotation_theta_reference);
%     R_pelvis_current = R_theta_current * R_psi_current * R_phi_current;
%     R_pelvis_reference = R_theta_reference * R_psi_reference * R_phi_reference;


%     % calculate lumbar flexion
%     lumbar_cor_to_cervix_cor_current_world = cervix_cor_current - lumbar_cor_current;
%     lumbar_cor_to_cervix_cor_reference_world = cervix_cor_reference - lumbar_cor_reference;
%     lumbar_cor_to_cervix_cor_current_pelvis = pelvis_rotation_current * lumbar_cor_to_cervix_cor_current_world;
%     lumbar_cor_to_cervix_cor_reference_pelvis = pelvis_rotation_reference * lumbar_cor_to_cervix_cor_reference_world;
% 
%     lumbar_rotation_phi_current = -atan2(lumbar_cor_to_cervix_cor_current_pelvis(3), lumbar_cor_to_cervix_cor_current_pelvis(2));
%     lumbar_rotation_phi_reference = -atan2(lumbar_cor_to_cervix_cor_reference_pelvis(3), lumbar_cor_to_cervix_cor_reference_pelvis(2));
%     lumbar_rotation_phi_reference_to_current = lumbar_rotation_phi_current - lumbar_rotation_phi_reference;
%     joint_angles(19) = lumbar_rotation_phi_reference_to_current;


    






    
    
    
    
    
    
    return
    
 
    
    
    
    
    
    
    
    
    
    
    
    %% below is old code that used center between the hip CoRs as reference
    if false
        % extract positions
        left_hip_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LHIPCOR')';
        left_knee_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LKNEECOR')';
        right_hip_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'RHIPCOR')';
        right_knee_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'RKNEECOR')';
        lumbar_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LUMBARCOR')';
        left_hip_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LHIPCOR')';
        left_knee_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LKNEECOR')';
        right_hip_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RHIPCOR')';
        right_knee_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RKNEECOR')';
        lumbar_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LUMBARCOR')';

        % calculate pelvis translation
        pelvis_base_point_current = (left_hip_cor_current + right_hip_cor_current) / 2;
        pelvis_base_point_reference = (left_hip_cor_reference + right_hip_cor_reference) / 2;
        pelvis_base_point_translation = pelvis_base_point_current - pelvis_base_point_reference;
        joint_angles(1:3) = pelvis_base_point_translation;

        % pelvis rotation in horizontal plane
        pelvis_base_point_to_right_hip_cor_current = right_hip_cor_current - pelvis_base_point_current;
        pelvis_base_point_to_right_hip_cor_reference = right_hip_cor_reference - pelvis_base_point_reference;
        pelvis_rotation_phi_current = atan2(pelvis_base_point_to_right_hip_cor_current(2), pelvis_base_point_to_right_hip_cor_current(1));
        pelvis_rotation_phi_reference = atan2(pelvis_base_point_to_right_hip_cor_reference(2), pelvis_base_point_to_right_hip_cor_reference(1));
        pelvis_rotation_phi_reference_to_current = pelvis_rotation_phi_current - pelvis_rotation_phi_reference;

        % pelvis rotation around anterior-posterior axis
        R_phi_current = expAxis([0; 0; 1], -pelvis_rotation_phi_current);
        R_phi_reference = expAxis([0; 0; 1], -pelvis_rotation_phi_reference);
        pelvis_base_point_to_right_hip_cor_phi_current = R_phi_current * pelvis_base_point_to_right_hip_cor_current;
        pelvis_base_point_to_right_hip_cor_phi_reference = R_phi_reference * pelvis_base_point_to_right_hip_cor_reference;
        pelvis_rotation_psi_current = -atan2(pelvis_base_point_to_right_hip_cor_phi_current(3), pelvis_base_point_to_right_hip_cor_phi_current(1));
        pelvis_rotation_psi_reference = -atan2(pelvis_base_point_to_right_hip_cor_phi_reference(3), pelvis_base_point_to_right_hip_cor_phi_reference(1));
        pelvis_rotation_psi_reference_to_current = pelvis_rotation_psi_current - pelvis_rotation_psi_reference;

        % pelvis rotation around medial-lateral axis
        pelvis_base_point_to_lumbar_cor_current = lumbar_cor_current - pelvis_base_point_current;
        pelvis_base_point_to_lumbar_cor_reference = lumbar_cor_reference - pelvis_base_point_reference;
        R_psi_current = expAxis([0; 1; 0], -pelvis_rotation_psi_current);
        R_psi_reference = expAxis([0; 1; 0], -pelvis_rotation_psi_reference);
        pelvis_base_point_to_lumbar_cor_phipsi_current = R_psi_current * R_phi_current * pelvis_base_point_to_lumbar_cor_current;
        pelvis_base_point_to_lumbar_cor_phipsi_reference = R_psi_reference * R_phi_reference * pelvis_base_point_to_lumbar_cor_reference;
        pelvis_rotation_theta_current = atan2(pelvis_base_point_to_lumbar_cor_phipsi_current(3), pelvis_base_point_to_lumbar_cor_phipsi_current(2));
        pelvis_rotation_theta_reference = atan2(pelvis_base_point_to_lumbar_cor_phipsi_reference(3), pelvis_base_point_to_lumbar_cor_phipsi_reference(2));
        pelvis_rotation_theta_reference_to_current = pelvis_rotation_theta_current - pelvis_rotation_theta_reference;    

        joint_angles(4) = pelvis_rotation_phi_reference_to_current;
        joint_angles(5) = pelvis_rotation_psi_reference_to_current;
        joint_angles(6) = pelvis_rotation_theta_reference_to_current;

        % calculate pelvis rotation matrix
        R_theta_current = expAxis([1; 0; 0], -pelvis_rotation_theta_current);
        R_theta_reference = expAxis([1; 0; 0], -pelvis_rotation_theta_reference);
        R_pelvis_current = R_theta_current * R_psi_current * R_phi_current;
        R_pelvis_reference = R_theta_reference * R_psi_reference * R_phi_reference;

        % calculate left hip flexion
        left_hip_cor_to_left_knee_cor_current_world = left_knee_cor_current - left_hip_cor_current;
        left_hip_cor_to_left_knee_cor_reference_world = left_knee_cor_reference - left_hip_cor_reference;
        left_hip_cor_to_left_knee_cor_current_pelvis = R_pelvis_current * left_hip_cor_to_left_knee_cor_current_world;
        left_hip_cor_to_left_knee_cor_reference_pelvis = R_pelvis_reference * left_hip_cor_to_left_knee_cor_reference_world;

        left_hip_rotation_phi_current = atan2(left_hip_cor_to_left_knee_cor_current_pelvis(3), left_hip_cor_to_left_knee_cor_current_pelvis(2));
        left_hip_rotation_phi_reference = atan2(left_hip_cor_to_left_knee_cor_reference_pelvis(3), left_hip_cor_to_left_knee_cor_reference_pelvis(2));
        left_hip_rotation_phi_reference_to_current = left_hip_rotation_phi_current - left_hip_rotation_phi_reference;
        joint_angles(7) = left_hip_rotation_phi_reference_to_current;
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    %% below is old code that used the MPSIS as a reference
    if false
        % extract needed marker positions
        LASI_current = extractMarkerTrajectories(marker_current, marker_headers, 'LASI')';
        RASI_current = extractMarkerTrajectories(marker_current, marker_headers, 'RASI')';
        LPSI_current = extractMarkerTrajectories(marker_current, marker_headers, 'LPSI')';
        RPSI_current = extractMarkerTrajectories(marker_current, marker_headers, 'RPSI')';

        LASI_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LASI')';
        RASI_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RASI')';
        LPSI_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LPSI')';
        RPSI_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RPSI')';

        % pelvis translation
        MPSIS_current = (LPSI_current + RPSI_current) / 2;
        MPSIS_reference = (LPSI_reference + RPSI_reference) / 2;
        MPSIS_translation = MPSIS_current - MPSIS_reference;
        joint_angles(1:3) = MPSIS_translation;

        % pelvis rotation
        MPSIS_to_RPSIS_current = RPSI_current - MPSIS_current;
        MPSIS_to_RPSIS_reference = RPSI_reference - MPSIS_reference;
        pelvis_rotation_phi_current = atan2(MPSIS_to_RPSIS_current(2), MPSIS_to_RPSIS_current(1));
        pelvis_rotation_phi_reference = atan2(MPSIS_to_RPSIS_reference(2), MPSIS_to_RPSIS_reference(1));
        pelvis_rotation_phi_reference_to_current = pelvis_rotation_phi_current - pelvis_rotation_phi_reference;

        R_phi_current = expAxis([0; 0; 1], -pelvis_rotation_phi_current);
        R_phi_reference = expAxis([0; 0; 1], -pelvis_rotation_phi_reference);
        MPSIS_to_RPSIS_phi_current = R_phi_current * MPSIS_to_RPSIS_current;
        MPSIS_to_RPSIS_phi_reference = R_phi_reference * MPSIS_to_RPSIS_reference;
        pelvis_rotation_psi_current = -atan2(MPSIS_to_RPSIS_phi_current(3), MPSIS_to_RPSIS_phi_current(1));
        pelvis_rotation_psi_reference = -atan2(MPSIS_to_RPSIS_phi_reference(3), MPSIS_to_RPSIS_phi_reference(1));
        pelvis_rotation_psi_reference_to_current = pelvis_rotation_psi_current - pelvis_rotation_psi_reference;

        MASIS_current = (LASI_current + RASI_current) / 2;
        MASIS_reference = (LASI_reference + RASI_reference) / 2;
        MPSIS_to_MASIS_current = MASIS_current - MPSIS_current;
        MPSIS_to_MASIS_reference = MASIS_reference - MPSIS_reference;
        R_psi_current = expAxis([0; 1; 0], -pelvis_rotation_psi_current);
        R_psi_reference = expAxis([0; 1; 0], -pelvis_rotation_psi_reference);
        MPSIS_to_MASIS_phipsi_current = R_psi_current * R_phi_current * MPSIS_to_MASIS_current;
        MPSIS_to_MASIS_phipsi_reference = R_psi_reference * R_phi_reference * MPSIS_to_MASIS_reference;

        pelvis_rotation_theta_current = atan2(MPSIS_to_MASIS_phipsi_current(3), MPSIS_to_MASIS_phipsi_current(2));
        pelvis_rotation_theta_reference = atan2(MPSIS_to_MASIS_phipsi_reference(3), MPSIS_to_MASIS_phipsi_reference(2));
        pelvis_rotation_theta_reference_to_current = pelvis_rotation_theta_current - pelvis_rotation_theta_reference;

        joint_angles(4) = pelvis_rotation_phi_reference_to_current;
        joint_angles(5) = pelvis_rotation_psi_reference_to_current;
        joint_angles(6) = pelvis_rotation_theta_reference_to_current;
    end
    

end

