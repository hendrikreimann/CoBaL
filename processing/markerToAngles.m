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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.



function joint_angles = ...
    markerToAngles ...
      ( ...
        marker_reference, ...
        marker_current, ...
        marker_labels, ...
        joint_center_reference, ...
        joint_center_current, ...
        joint_center_headers, ...
        direction_matrices, ...
        direction_matrix_labels ...
      )
    
    % get matrices for comparison
%     load('/Users/reimajbi/Box Sync/inverseKinematics/BRC/processed/00000000_XXX_simulation_001_markerTrajectories.mat')
%     load('subjectModel.mat')
%     knee_cor_reference = kinematic_tree.jointPositions{10};
%     kinematic_tree.jointAngles = joint_angle_trajectories_simulated(1, :)';
%     kinematic_tree.updateConfiguration;
    
%     left_knee_cor_ground_truth = kinematic_tree.jointPositions{10}
%     left_ankle_cor_ground_truth = kinematic_tree.jointPositions{12}
%     left_toes_eef_current_ground_truth = kinematic_tree.endEffectorPositions{2}
%     lumbar_cor_ground_truth = kinematic_tree.jointPositions{23};
%     cervix_cor_ground_truth = kinematic_tree.jointPositions{26};
% 
%     pelvis_transformation_from_simulation = kinematic_tree.productsOfExponentials{6};
%     lthigh_transformation_from_simulation = kinematic_tree.productsOfExponentials{9};
%     exp1 = kinematic_tree.twistExponentials{1};
%     exp2 = kinematic_tree.twistExponentials{2};
%     exp3 = kinematic_tree.twistExponentials{3};
%     exp4 = kinematic_tree.twistExponentials{4};
%     exp5 = kinematic_tree.twistExponentials{5};
%     exp6 = kinematic_tree.twistExponentials{6};
%     exp7 = kinematic_tree.twistExponentials{7};
%     exp8 = kinematic_tree.twistExponentials{8};
%     exp9 = kinematic_tree.twistExponentials{9};
%     exp10 = kinematic_tree.twistExponentials{10};
%     exp11 = kinematic_tree.twistExponentials{11};
%     exp12 = kinematic_tree.twistExponentials{12};
%     exp13 = kinematic_tree.twistExponentials{13};
%     exp14 = kinematic_tree.twistExponentials{14};
%     exp15 = kinematic_tree.twistExponentials{15};
%     exp16 = kinematic_tree.twistExponentials{16};
%     exp17 = kinematic_tree.twistExponentials{17};
%     exp18 = kinematic_tree.twistExponentials{18};
%     exp19 = kinematic_tree.twistExponentials{19};
%     exp20 = kinematic_tree.twistExponentials{20};
%     exp21 = kinematic_tree.twistExponentials{21};
%     exp22 = kinematic_tree.twistExponentials{22};
%     exp23 = kinematic_tree.twistExponentials{23};
%     exp24 = kinematic_tree.twistExponentials{24};
%     exp25 = kinematic_tree.twistExponentials{25};
%     exp26 = kinematic_tree.twistExponentials{26};
%     exp27 = kinematic_tree.twistExponentials{27};
%     exp28 = kinematic_tree.twistExponentials{28};
%     exp29 = kinematic_tree.twistExponentials{29};
%     exp30 = kinematic_tree.twistExponentials{30};
%     exp31 = kinematic_tree.twistExponentials{31};
%     exp32 = kinematic_tree.twistExponentials{32};
%     exp33 = kinematic_tree.twistExponentials{33};
%     exp34 = kinematic_tree.twistExponentials{34};
%     exp35 = kinematic_tree.twistExponentials{35};
%     exp36 = kinematic_tree.twistExponentials{36};
%     exp37 = kinematic_tree.twistExponentials{37};
%     exp38 = kinematic_tree.twistExponentials{38};
  
  
    %% preparations
    
    body_direction_matrix = direction_matrices{strcmp(direction_matrix_labels, 'body')};
    left_hip_direction_matrix = direction_matrices{strcmp(direction_matrix_labels, 'left_hip')};
    right_hip_direction_matrix = direction_matrices{strcmp(direction_matrix_labels, 'right_hip')};
    left_knee_direction_matrix = direction_matrices{strcmp(direction_matrix_labels, 'left_knee')};
    right_knee_direction_matrix = direction_matrices{strcmp(direction_matrix_labels, 'right_knee')};
    left_ankle_direction_matrix = direction_matrices{strcmp(direction_matrix_labels, 'left_ankle')};
    right_ankle_direction_matrix = direction_matrices{strcmp(direction_matrix_labels, 'right_ankle')};
    lumbar_direction_matrix = direction_matrices{strcmp(direction_matrix_labels, 'lumbar')};
    cervix_direction_matrix = direction_matrices{strcmp(direction_matrix_labels, 'cervix')};
    left_shoulder_direction_matrix = direction_matrices{strcmp(direction_matrix_labels, 'left_shoulder')};
    right_shoulder_direction_matrix = direction_matrices{strcmp(direction_matrix_labels, 'right_shoulder')};
    left_elbow_direction_matrix = direction_matrices{strcmp(direction_matrix_labels, 'left_elbow')};
    right_elbow_direction_matrix = direction_matrices{strcmp(direction_matrix_labels, 'right_elbow')};
    left_wrist_direction_matrix = direction_matrices{strcmp(direction_matrix_labels, 'left_wrist')};
    right_wrist_direction_matrix = direction_matrices{strcmp(direction_matrix_labels, 'right_wrist')};
    
    joint_angles = zeros(38, 1);
    
    body_direction_matrix_inverse = body_direction_matrix^(-1);
    left_hip_direction_matrix_inverse = left_hip_direction_matrix^(-1);
    left_knee_direction_matrix_inverse = left_knee_direction_matrix^(-1);
    left_ankle_direction_matrix_inverse = left_ankle_direction_matrix^(-1);
    right_hip_direction_matrix_inverse = right_hip_direction_matrix^(-1);
    right_knee_direction_matrix_inverse = right_knee_direction_matrix^(-1);
    right_ankle_direction_matrix_inverse = right_ankle_direction_matrix^(-1);
    lumbar_direction_matrix_inverse = lumbar_direction_matrix^(-1);
    cervix_direction_matrix_inverse = cervix_direction_matrix^(-1);
    left_shoulder_direction_matrix_inverse = left_shoulder_direction_matrix^(-1);
    left_elbow_direction_matrix_inverse = left_elbow_direction_matrix^(-1);
    left_wrist_direction_matrix_inverse = left_wrist_direction_matrix^(-1);
    right_shoulder_direction_matrix_inverse = right_shoulder_direction_matrix^(-1);
    right_elbow_direction_matrix_inverse = right_elbow_direction_matrix^(-1);
    right_wrist_direction_matrix_inverse = right_wrist_direction_matrix^(-1);
    
    % extract positions
    LASI_current = extractMarkerTrajectories(marker_current, marker_labels, 'LASI')';
    LPSI_current = extractMarkerTrajectories(marker_current, marker_labels, 'LPSI')';
    LKNE_current = extractMarkerTrajectories(marker_current, marker_labels, 'LKNE')';
    LANK_current = extractMarkerTrajectories(marker_current, marker_labels, 'LANK')';
    LTOE_current = extractMarkerTrajectories(marker_current, marker_labels, 'LTOE')';
    LTOEL_current = extractMarkerTrajectories(marker_current, marker_labels, 'LTOEL')';
    LSHO_current = extractMarkerTrajectories(marker_current, marker_labels, 'LSHO')';
    LFHD_current = extractMarkerTrajectories(marker_current, marker_labels, 'LFHD')';
    LBHD_current = extractMarkerTrajectories(marker_current, marker_labels, 'LBHD')';
    LELB_current = extractMarkerTrajectories(marker_current, marker_labels, 'LELB')';
    LWRB_current = extractMarkerTrajectories(marker_current, marker_labels, 'LWRB')';
    LFIN_current = extractMarkerTrajectories(marker_current, marker_labels, 'LFIN')';
    RASI_current = extractMarkerTrajectories(marker_current, marker_labels, 'RASI')';
    RPSI_current = extractMarkerTrajectories(marker_current, marker_labels, 'RPSI')';
    RKNE_current = extractMarkerTrajectories(marker_current, marker_labels, 'RKNE')';
    RANK_current = extractMarkerTrajectories(marker_current, marker_labels, 'RANK')';
    RTOE_current = extractMarkerTrajectories(marker_current, marker_labels, 'RTOE')';
    RTOEL_current = extractMarkerTrajectories(marker_current, marker_labels, 'RTOEL')';
    RSHO_current = extractMarkerTrajectories(marker_current, marker_labels, 'RSHO')';
    RFHD_current = extractMarkerTrajectories(marker_current, marker_labels, 'RFHD')';
    RBHD_current = extractMarkerTrajectories(marker_current, marker_labels, 'RBHD')';
    RELB_current = extractMarkerTrajectories(marker_current, marker_labels, 'RELB')';
    RWRB_current = extractMarkerTrajectories(marker_current, marker_labels, 'RWRB')';
    RFIN_current = extractMarkerTrajectories(marker_current, marker_labels, 'RFIN')';
    left_hip_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LHIPCOR')';
    left_knee_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LKNEECOR')';
    left_ankle_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LANKLECOR')';
    left_toes_eef_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LTOESEEF')';
    left_shoulder_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LSHOULDERCOR')';
    left_elbow_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LELBOWCOR')';
    left_wrist_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LWRISTCOR')';
    left_hand_eef_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LHANDEEF')';
    right_hip_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'RHIPCOR')';
    right_knee_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'RKNEECOR')';
    right_ankle_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'RANKLECOR')';
    right_toes_eef_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'RTOESEEF')';
    right_shoulder_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'RSHOULDERCOR')';
    right_elbow_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'RELBOWCOR')';
    right_wrist_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'RWRISTCOR')';
    right_hand_eef_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'RHANDEEF')';
    lumbar_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LUMBARCOR')';
    cervix_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'CERVIXCOR')';
    head_center_current = mean([LFHD_current RFHD_current LBHD_current RBHD_current], 2);
    forehead_center_current = mean([LFHD_current RFHD_current], 2);

    LASI_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'LASI')';
    LPSI_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'LPSI')';
    LKNE_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'LKNE')';
    LANK_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'LANK')';
    LTOE_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'LTOE')';
    LTOEL_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'LTOEL')';
    LSHO_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'LSHO')';
    LFHD_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'LFHD')';
    LBHD_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'LBHD')';
    LELB_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'LELB')';
    LWRB_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'LWRB')';
    LFIN_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'LFIN')';
    RASI_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'RASI')';
    RPSI_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'RPSI')';
    RKNE_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'RKNE')';
    RANK_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'RANK')';
    RTOE_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'RTOE')';
    RTOEL_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'RTOEL')';
    RSHO_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'RSHO')';
    RFHD_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'RFHD')';
    RBHD_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'RBHD')';
    RELB_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'RELB')';
    RWRB_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'RWRB')';
    RFIN_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'RFIN')';
    left_hip_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LHIPCOR')';
    left_knee_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LKNEECOR')';
    left_ankle_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LANKLECOR')';
    left_toes_eef_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LTOESEEF')';
    left_shoulder_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LSHOULDERCOR')';
    left_elbow_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LELBOWCOR')';
    left_wrist_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LWRISTCOR')';
    left_hand_eef_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LHANDEEF')';
    right_hip_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RHIPCOR')';
    right_knee_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RKNEECOR')';
    right_ankle_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RANKLECOR')';
    right_toes_eef_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RTOESEEF')';
    right_shoulder_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RSHOULDERCOR')';
    right_elbow_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RELBOWCOR')';
    right_wrist_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RWRISTCOR')';
    right_hand_eef_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RHANDEEF')';
    lumbar_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LUMBARCOR')';
    cervix_cor_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'CERVIXCOR')';
    head_center_reference = mean([LFHD_reference RFHD_reference LBHD_reference RBHD_reference], 2);
    forehead_center_reference = mean([LFHD_reference RFHD_reference], 2);
    
    %% pelvis
    % pelvis angles
    pelvis_transformation_current_cell = calculateMcsToWcsTransformations_detailed(marker_current, marker_labels, {'PELVIS'});    
    pelvis_transformation_current = pelvis_transformation_current_cell{1};
    pelvis_translation_current = pelvis_transformation_current(1:3, 4);

    pelvis_transformation_reference_cell = calculateMcsToWcsTransformations_detailed(marker_reference, marker_labels, {'PELVIS'});
    pelvis_transformation_reference = pelvis_transformation_reference_cell{1};
    pelvis_translation_reference = pelvis_transformation_reference(1:3, 4);
    pelvis_translation_reference_to_current = pelvis_translation_current - pelvis_translation_reference;

    pelvis_joint_transformation = pelvis_transformation_current * pelvis_transformation_reference^(-1);
    pelvis_joint_rotation = pelvis_joint_transformation(1:3, 1:3);
    pelvis_rotation_angles = eulerAnglesFromRotationMatrixZXY(pelvis_joint_rotation); % these are the correct values
    
    joint_angles(1:3) = pelvis_translation_reference_to_current;
    joint_angles(4:6) = pelvis_rotation_angles;
    
    %% left leg
    % calculate left hip flexion-extension
    left_hip_cor_to_left_knee_cor_current_world = left_knee_cor_current - left_hip_cor_current;
    R_world_to_7 = left_hip_direction_matrix_inverse * pelvis_joint_rotation^(-1);
    left_hip_cor_to_left_knee_cor_current_7 = R_world_to_7 * left_hip_cor_to_left_knee_cor_current_world;
    theta_7 = atan2(left_hip_cor_to_left_knee_cor_current_7(2), -left_hip_cor_to_left_knee_cor_current_7(3));
    joint_angles(7) = theta_7;
    
    % left hip ab-adduction
    R_7 = expAxis(left_hip_direction_matrix(:, 1), joint_angles(7));
    R_world_to_8 = left_hip_direction_matrix_inverse * R_7^(-1) * pelvis_joint_rotation^(-1);
    left_hip_cor_to_left_knee_cor_current_8 = R_world_to_8 * left_hip_cor_to_left_knee_cor_current_world;
    theta_8 = atan2(-left_hip_cor_to_left_knee_cor_current_8(1), -left_hip_cor_to_left_knee_cor_current_8(3));
    joint_angles(8) = theta_8;
    
    % left hip internal/external rotation - from knee marker
    left_knee_cor_to_LKNE_current_world = LKNE_current - left_knee_cor_current;
    R_8 = expAxis(left_hip_direction_matrix(:, 2), joint_angles(8));
    R_world_to_9 = left_hip_direction_matrix_inverse * R_8^(-1) * R_7^(-1) * pelvis_joint_rotation^(-1);
    left_knee_cor_to_LKNE_current_9 = R_world_to_9 * left_knee_cor_to_LKNE_current_world;
    angle_now = atan2(left_knee_cor_to_LKNE_current_9(2), -left_knee_cor_to_LKNE_current_9(1));
    left_knee_cor_to_LKNE_reference_world = LKNE_reference - left_knee_cor_reference;
    left_knee_cor_to_LKNE_reference_hip = left_hip_direction_matrix_inverse * left_knee_cor_to_LKNE_reference_world;
    angle_reference = atan2(left_knee_cor_to_LKNE_reference_hip(2), -left_knee_cor_to_LKNE_reference_hip(1));
    joint_angles(9) = angle_now - angle_reference;
    
    % left knee flexion
    left_knee_cor_to_left_ankle_cor_current_world = left_ankle_cor_current - left_knee_cor_current;
    R_9 = expAxis(-left_hip_direction_matrix(:, 3), joint_angles(9));
    left_hip_joint_rotation = R_7 * R_8 * R_9;
    R_world_to_10 = left_knee_direction_matrix_inverse * left_hip_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    left_knee_cor_to_left_ankle_cor_current_10 = R_world_to_10 * left_knee_cor_to_left_ankle_cor_current_world;
    theta_10 = atan2(-left_knee_cor_to_left_ankle_cor_current_10(2), -left_knee_cor_to_left_ankle_cor_current_10(3));
    joint_angles(10) = theta_10;
    
    % left knee internal rotation - from ankle marker
    left_ankle_cor_to_LANK_current_world = LANK_current - left_ankle_cor_current;
    R_10 = expAxis(-left_knee_direction_matrix(:, 1), joint_angles(10));
    R_world_to_11 = left_knee_direction_matrix_inverse * R_10^(-1) * left_hip_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    left_ankle_cor_to_LANK_current_11 = R_world_to_11 * left_ankle_cor_to_LANK_current_world;
    angle_now = atan2(-left_ankle_cor_to_LANK_current_11(2), -left_ankle_cor_to_LANK_current_11(1));
    left_ankle_cor_to_LANK_reference_world = LANK_reference - left_ankle_cor_reference;
    left_ankle_cor_to_LANK_reference_knee = left_knee_direction_matrix_inverse * left_ankle_cor_to_LANK_reference_world;
    angle_reference = atan2(-left_ankle_cor_to_LANK_reference_knee(2), -left_ankle_cor_to_LANK_reference_knee(1));
    theta_11 = angle_now - angle_reference;
    joint_angles(11) = theta_11;
    
    % left ankle plantar-dorsiflexion
    left_ankle_cor_to_left_toes_eef_current_world = left_toes_eef_current - left_ankle_cor_current;
    R_11 = expAxis(left_knee_direction_matrix(:, 3), joint_angles(11));
    left_knee_joint_rotation = R_10 * R_11;
    R_world_to_12 = left_ankle_direction_matrix_inverse * left_knee_joint_rotation^(-1) * left_hip_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    left_ankle_cor_to_left_toes_eef_current_12 = R_world_to_12 * left_ankle_cor_to_left_toes_eef_current_world;
    theta_12 = atan2(left_ankle_cor_to_left_toes_eef_current_12(3), left_ankle_cor_to_left_toes_eef_current_12(2));
    joint_angles(12) = theta_12;
   
    % left ankle in-eversion - from toes markers
    left_toes_eef_to_LTOEL_current_world = LTOEL_current - left_toes_eef_current;
    R_12 = expAxis(left_ankle_direction_matrix(:, 1), joint_angles(12));
    R_world_to_13 = left_ankle_direction_matrix_inverse * R_12^(-1) * left_knee_joint_rotation^(-1) * left_hip_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    left_toes_eef_to_LTOEL_current_13 = R_world_to_13 * left_toes_eef_to_LTOEL_current_world;
    angle_now = atan2(left_toes_eef_to_LTOEL_current_13(3), -left_toes_eef_to_LTOEL_current_13(1));
    left_toes_eef_to_LTOEL_reference_world = LTOEL_reference - left_toes_eef_reference;
    left_toes_eef_to_LTOEL_reference_ankle = left_ankle_direction_matrix_inverse * left_toes_eef_to_LTOEL_reference_world;
    angle_reference = atan2(left_toes_eef_to_LTOEL_reference_ankle(3), -left_toes_eef_to_LTOEL_reference_ankle(1));
    theta_13 = angle_now - angle_reference;
    joint_angles(13) = theta_13;
    
    %% right leg
    
    % calculate right hip flexion-extension
    right_hip_cor_to_right_knee_cor_current_world = right_knee_cor_current - right_hip_cor_current;
    R_world_to_14 = right_hip_direction_matrix_inverse * pelvis_joint_rotation^(-1);
    right_hip_cor_to_right_knee_cor_current_14 = R_world_to_14 * right_hip_cor_to_right_knee_cor_current_world;
    theta_14 = atan2(right_hip_cor_to_right_knee_cor_current_14(2), -right_hip_cor_to_right_knee_cor_current_14(3));
    joint_angles(14) = theta_14;
    
    % right hip ab-adduction
    R_14 = expAxis(right_hip_direction_matrix(:, 1), joint_angles(14));
    R_world_to_15 = right_hip_direction_matrix_inverse * R_14^(-1) * pelvis_joint_rotation^(-1);
    right_hip_cor_to_right_knee_cor_current_15 = R_world_to_15 * right_hip_cor_to_right_knee_cor_current_world;
    theta_15 = atan2(right_hip_cor_to_right_knee_cor_current_15(1), -right_hip_cor_to_right_knee_cor_current_15(3));
    joint_angles(15) = theta_15;
    
    % right hip internal/external rotation - from knee marker
    right_knee_cor_to_RKNE_current_world = RKNE_current - right_knee_cor_current;
    R_15 = expAxis(-right_hip_direction_matrix(:, 2), joint_angles(15));
    R_world_to_16 = right_hip_direction_matrix_inverse * R_15^(-1) * R_14^(-1) * pelvis_joint_rotation^(-1);
    right_knee_cor_to_RKNE_current_16 = R_world_to_16 * right_knee_cor_to_RKNE_current_world;
    angle_now = atan2(right_knee_cor_to_RKNE_current_16(2), right_knee_cor_to_RKNE_current_16(1));
    right_knee_cor_to_RKNE_reference_world = RKNE_reference - right_knee_cor_reference;
    right_knee_cor_to_RKNE_reference_hip = right_hip_direction_matrix_inverse * right_knee_cor_to_RKNE_reference_world;
    angle_reference = atan2(right_knee_cor_to_RKNE_reference_hip(2), right_knee_cor_to_RKNE_reference_hip(1));
    joint_angles(16) = angle_now - angle_reference;
    
    % right knee flexion
    right_knee_cor_to_right_ankle_cor_current_world = right_ankle_cor_current - right_knee_cor_current;
    R_16 = expAxis(right_hip_direction_matrix(:, 3), joint_angles(16));
    right_hip_joint_rotation = R_14 * R_15 * R_16;
    R_world_to_17 = right_knee_direction_matrix_inverse * right_hip_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    right_knee_cor_to_right_ankle_cor_current_17 = R_world_to_17 * right_knee_cor_to_right_ankle_cor_current_world;
    theta_17 = atan2(-right_knee_cor_to_right_ankle_cor_current_17(2), -right_knee_cor_to_right_ankle_cor_current_17(3));
    joint_angles(17) = theta_17;
    
    % right knee internal rotation - from ankle marker
    right_ankle_cor_to_RANK_current_world = RANK_current - right_ankle_cor_current;
    R_17 = expAxis(-right_knee_direction_matrix(:, 1), joint_angles(17));
    R_world_to_18 = right_knee_direction_matrix_inverse * R_17^(-1) * right_hip_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    right_ankle_cor_to_RANK_current_18 = R_world_to_18 * right_ankle_cor_to_RANK_current_world;
    angle_now = atan2(right_ankle_cor_to_RANK_current_18(2), -right_ankle_cor_to_RANK_current_18(1));
    right_ankle_cor_to_RANK_reference_world = RANK_reference - right_ankle_cor_reference;
    right_ankle_cor_to_RANK_reference_knee = right_knee_direction_matrix_inverse * right_ankle_cor_to_RANK_reference_world;
    angle_reference = atan2(right_ankle_cor_to_RANK_reference_knee(2), -right_ankle_cor_to_RANK_reference_knee(1));
    theta_18 = angle_now - angle_reference;
    joint_angles(18) = theta_18;
    
    % right ankle plantar-dorsiflexion
    right_ankle_cor_to_right_toes_eef_current_world = right_toes_eef_current - right_ankle_cor_current;
    R_18 = expAxis(-right_knee_direction_matrix(:, 3), joint_angles(18));
    right_knee_joint_rotation = R_17 * R_18;
    R_world_to_19 = right_ankle_direction_matrix_inverse * right_knee_joint_rotation^(-1) * right_hip_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    right_ankle_cor_to_right_toes_eef_current_19 = R_world_to_19 * right_ankle_cor_to_right_toes_eef_current_world;
    theta_19 = atan2(right_ankle_cor_to_right_toes_eef_current_19(3), right_ankle_cor_to_right_toes_eef_current_19(2));
    joint_angles(19) = theta_19;
   
    % right ankle in-eversion - from toes markers
    right_toes_eef_to_RTOEL_current_world = RTOEL_current - right_toes_eef_current;
    R_19 = expAxis(right_ankle_direction_matrix(:, 1), joint_angles(19));
    R_world_to_20 = right_ankle_direction_matrix_inverse * R_19^(-1) * right_knee_joint_rotation^(-1) * right_hip_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    right_toes_eef_to_RTOEL_current_20 = R_world_to_20 * right_toes_eef_to_RTOEL_current_world;
    angle_now = atan2(right_toes_eef_to_RTOEL_current_20(3), right_toes_eef_to_RTOEL_current_20(1));
    right_toes_eef_to_RTOEL_reference_world = RTOEL_reference - right_toes_eef_reference;
    right_toes_eef_to_RTOEL_reference_ankle = right_ankle_direction_matrix_inverse * right_toes_eef_to_RTOEL_reference_world;
    angle_reference = atan2(right_toes_eef_to_RTOEL_reference_ankle(3), right_toes_eef_to_RTOEL_reference_ankle(1));
    theta_20 = angle_now - angle_reference;
    joint_angles(20) = theta_20;    
    
    %% lumbar
    
    % calculate lumbar flexion-extension
    lumbar_cor_to_cervix_cor_current_world = cervix_cor_current - lumbar_cor_current;
    R_world_to_21 = lumbar_direction_matrix_inverse * pelvis_joint_rotation^(-1);
    lumbar_cor_to_cervix_cor_current_21 = R_world_to_21 * lumbar_cor_to_cervix_cor_current_world;
    theta_21 = atan2(lumbar_cor_to_cervix_cor_current_21(2), lumbar_cor_to_cervix_cor_current_21(3));
    joint_angles(21) = theta_21;
    
    % left lumbar ab-adduction
    R_21 = expAxis(-lumbar_direction_matrix(:, 1), joint_angles(21));
    R_world_to_22 = lumbar_direction_matrix_inverse * R_21^(-1) * pelvis_joint_rotation^(-1);
    lumbar_cor_to_cervix_cor_current_22 = R_world_to_22 * lumbar_cor_to_cervix_cor_current_world;
    theta_22 = atan2(lumbar_cor_to_cervix_cor_current_22(1), lumbar_cor_to_cervix_cor_current_22(3));
    joint_angles(22) = theta_22;
    
    % left lumbar internal/external rotation - from cervix CoR to left shoulder marker
    cervix_cor_to_LSHO_current_world = LSHO_current - cervix_cor_current;
    R_22 = expAxis(lumbar_direction_matrix(:, 2), joint_angles(22));
    R_world_to_23 = lumbar_direction_matrix_inverse * R_22^(-1) * R_21^(-1) * pelvis_joint_rotation^(-1);
    cervix_cor_to_LSHO_current_23 = R_world_to_23 * cervix_cor_to_LSHO_current_world;
    angle_now = atan2(cervix_cor_to_LSHO_current_23(2), -cervix_cor_to_LSHO_current_23(1));
    cervix_cor_to_LSHO_reference_world = LSHO_reference - cervix_cor_reference;
    cervix_cor_to_LSHO_reference_lumbar = lumbar_direction_matrix_inverse * cervix_cor_to_LSHO_reference_world;
    angle_reference = atan2(cervix_cor_to_LSHO_reference_lumbar(2), -cervix_cor_to_LSHO_reference_lumbar(1));
    joint_angles(23) = angle_now - angle_reference;    
    
    R_23 = expAxis(-lumbar_direction_matrix(:, 3), joint_angles(23));
    lumbar_joint_rotation = R_21 * R_22 * R_23;
    
    %% cervix
    
    % calculate cervix flexion-extension
    cervix_cor_to_head_center_current_world = head_center_current - cervix_cor_current;
    R_world_to_24 = cervix_direction_matrix_inverse * lumbar_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    cervix_cor_to_head_center_current_24 = R_world_to_24 * cervix_cor_to_head_center_current_world;
    theta_24 = atan2(cervix_cor_to_head_center_current_24(2), cervix_cor_to_head_center_current_24(3));
    joint_angles(24) = theta_24;
    
    % left cervix ab-adduction
    R_24 = expAxis(-cervix_direction_matrix(:, 1), joint_angles(24));
    R_world_to_25 = cervix_direction_matrix_inverse * R_24^(-1) * lumbar_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    cervix_cor_to_head_center_current_25 = R_world_to_25 * cervix_cor_to_head_center_current_world;
    theta_25 = atan2(cervix_cor_to_head_center_current_25(1), cervix_cor_to_head_center_current_25(3));
    joint_angles(25) = theta_25;
    
    % left cervix internal/external rotation - from cervix CoR
    head_center_to_forehead_center_current_world = forehead_center_current - head_center_current;
    R_25 = expAxis(cervix_direction_matrix(:, 2), joint_angles(25));
    R_world_to_26 = cervix_direction_matrix_inverse * R_25^(-1) * R_24^(-1) * lumbar_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    head_center_to_forehead_center_current_26 = R_world_to_26 * head_center_to_forehead_center_current_world;
    angle_now = atan2(head_center_to_forehead_center_current_26(2), -head_center_to_forehead_center_current_26(1));
    head_center_to_forehead_center_reference_world = forehead_center_reference - head_center_reference;
    head_center_to_forehead_center_reference_cervix = cervix_direction_matrix_inverse * head_center_to_forehead_center_reference_world;
    angle_reference = atan2(head_center_to_forehead_center_reference_cervix(2), -head_center_to_forehead_center_reference_cervix(1));
    joint_angles(26) = angle_now - angle_reference;    
    
    %% left arm
    % calculate left shoulder flexion-extension
    left_shoulder_cor_to_left_elbow_cor_current_world = left_elbow_cor_current - left_shoulder_cor_current;
    R_world_to_27 = left_shoulder_direction_matrix_inverse * lumbar_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    left_shoulder_cor_to_left_elbow_cor_current_27 = R_world_to_27 * left_shoulder_cor_to_left_elbow_cor_current_world;
    theta_27 = atan2(left_shoulder_cor_to_left_elbow_cor_current_27(2), -left_shoulder_cor_to_left_elbow_cor_current_27(3));
    joint_angles(27) = theta_27;
    
    % left shoulder ab-adduction
    R_27 = expAxis(left_shoulder_direction_matrix(:, 1), joint_angles(27));
    R_world_to_28 = left_shoulder_direction_matrix_inverse * R_27^(-1) * lumbar_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    left_shoulder_cor_to_left_elbow_cor_current_28 = R_world_to_28 * left_shoulder_cor_to_left_elbow_cor_current_world;
    theta_28 = atan2(-left_shoulder_cor_to_left_elbow_cor_current_28(1), -left_shoulder_cor_to_left_elbow_cor_current_28(3));
    joint_angles(28) = theta_28;
    
    % left shoulder internal/external rotation - from elbow marker
    left_elbow_cor_to_LELB_current_world = LELB_current - left_elbow_cor_current;
    R_28 = expAxis(left_shoulder_direction_matrix(:, 2), joint_angles(28));
    R_world_to_29 = left_shoulder_direction_matrix_inverse * R_28^(-1) * R_27^(-1) * lumbar_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    left_elbow_cor_to_LELB_current_29 = R_world_to_29 * left_elbow_cor_to_LELB_current_world;
    angle_now = atan2(left_elbow_cor_to_LELB_current_29(2), -left_elbow_cor_to_LELB_current_29(1));
    left_elbow_cor_to_LELB_reference_world = LELB_reference - left_elbow_cor_reference;
    left_elbow_cor_to_LELB_reference_shoulder = left_shoulder_direction_matrix_inverse * left_elbow_cor_to_LELB_reference_world;
    angle_reference = atan2(left_elbow_cor_to_LELB_reference_shoulder(2), -left_elbow_cor_to_LELB_reference_shoulder(1));
    joint_angles(29) = angle_now - angle_reference;
    
    % left elbow flexion
    left_elbow_cor_to_left_wrist_cor_current_world = left_wrist_cor_current - left_elbow_cor_current;
    R_29 = expAxis(-left_shoulder_direction_matrix(:, 3), joint_angles(29));
    left_shoulder_joint_rotation = R_27 * R_28 * R_29;
    R_world_to_30 = left_elbow_direction_matrix_inverse * left_shoulder_joint_rotation^(-1) * lumbar_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    left_elbow_cor_to_left_wrist_cor_current_30 = R_world_to_30 * left_elbow_cor_to_left_wrist_cor_current_world;
    theta_30 = atan2(left_elbow_cor_to_left_wrist_cor_current_30(3), left_elbow_cor_to_left_wrist_cor_current_30(2));
    joint_angles(30) = theta_30;
    
    % left elbow internal rotation - from wrist marker
    left_wrist_cor_to_LWRB_current_world = LWRB_current - left_wrist_cor_current;
    R_30 = expAxis(left_elbow_direction_matrix(:, 1), joint_angles(30));
    R_world_to_31 = left_elbow_direction_matrix_inverse * R_30^(-1) * left_shoulder_joint_rotation^(-1) * lumbar_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    left_wrist_cor_to_LWRB_current_31 = R_world_to_31 * left_wrist_cor_to_LWRB_current_world;
    angle_now = atan2(-left_wrist_cor_to_LWRB_current_31(1), -left_wrist_cor_to_LWRB_current_31(3));
    left_wrist_cor_to_LWRB_reference_world = LWRB_reference - left_wrist_cor_reference;
    left_wrist_cor_to_LWRB_reference_elbow = left_elbow_direction_matrix_inverse * left_wrist_cor_to_LWRB_reference_world;
    angle_reference = atan2(-left_wrist_cor_to_LWRB_reference_elbow(1), -left_wrist_cor_to_LWRB_reference_elbow(3));
    theta_31 = angle_now - angle_reference;
    joint_angles(31) = theta_31;
    
    % left wrist flexion and radial-ulnar deviation angles - from finger marker
    left_wrist_cor_to_LFIN_current_world = LFIN_current - left_wrist_cor_current;
    R_31 = expAxis(left_elbow_direction_matrix(:, 2), joint_angles(31));
    left_elbow_joint_rotation = R_30 * R_31;
    R_world_to_32 = left_wrist_direction_matrix_inverse * left_elbow_joint_rotation^(-1) * left_shoulder_joint_rotation^(-1) * lumbar_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    left_wrist_cor_to_LFIN_current_32 = R_world_to_32 * left_wrist_cor_to_LFIN_current_world;
    left_wrist_cor_to_LFIN_reference_world = LFIN_reference - left_wrist_cor_reference;
    left_wrist_cor_to_LFIN_reference_wrist = left_wrist_direction_matrix_inverse * left_wrist_cor_to_LFIN_reference_world;
    angle_now = atan2(left_wrist_cor_to_LFIN_current_32(1), left_wrist_cor_to_LFIN_current_32(2));
    angle_reference = atan2(left_wrist_cor_to_LFIN_reference_wrist(1), left_wrist_cor_to_LFIN_reference_wrist(2));
    theta_32 = angle_now - angle_reference;
    joint_angles(32) = theta_32;
    
    
    %% right arm
    % calculate right shoulder flexion-extension
    right_shoulder_cor_to_right_elbow_cor_current_world = right_elbow_cor_current - right_shoulder_cor_current;
    R_world_to_33 = right_shoulder_direction_matrix_inverse * lumbar_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    right_shoulder_cor_to_right_elbow_cor_current_33 = R_world_to_33 * right_shoulder_cor_to_right_elbow_cor_current_world;
    theta_33 = atan2(right_shoulder_cor_to_right_elbow_cor_current_33(2), -right_shoulder_cor_to_right_elbow_cor_current_33(3));
    joint_angles(33) = theta_33;
    
    % right shoulder ab-adduction
    R_33 = expAxis(right_shoulder_direction_matrix(:, 1), joint_angles(33));
    R_world_to_34 = right_shoulder_direction_matrix_inverse * R_33^(-1) * lumbar_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    right_shoulder_cor_to_right_elbow_cor_current_34 = R_world_to_34 * right_shoulder_cor_to_right_elbow_cor_current_world;
    theta_34 = atan2(right_shoulder_cor_to_right_elbow_cor_current_34(1), -right_shoulder_cor_to_right_elbow_cor_current_34(3));
    joint_angles(34) = theta_34;
    
    % right shoulder internal/external rotation - from elbow marker
    right_elbow_cor_to_RELB_current_world = RELB_current - right_elbow_cor_current;
    R_34 = expAxis(-right_shoulder_direction_matrix(:, 2), joint_angles(34));
    R_world_to_35 = right_shoulder_direction_matrix_inverse * R_34^(-1) * R_33^(-1) * lumbar_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    right_elbow_cor_to_RELB_current_35 = R_world_to_35 * right_elbow_cor_to_RELB_current_world;
    angle_now = atan2(right_elbow_cor_to_RELB_current_35(2), right_elbow_cor_to_RELB_current_35(1));
    right_elbow_cor_to_RELB_reference_world = RELB_reference - right_elbow_cor_reference;
    right_elbow_cor_to_RELB_reference_shoulder = right_shoulder_direction_matrix_inverse * right_elbow_cor_to_RELB_reference_world;
    angle_reference = atan2(right_elbow_cor_to_RELB_reference_shoulder(2), right_elbow_cor_to_RELB_reference_shoulder(1));
    joint_angles(35) = angle_now - angle_reference;
    
    % right elbow flexion
    right_elbow_cor_to_right_wrist_cor_current_world = right_wrist_cor_current - right_elbow_cor_current;
    R_35 = expAxis(right_shoulder_direction_matrix(:, 3), joint_angles(35));
    right_shoulder_joint_rotation = R_33 * R_34 * R_35;
    R_world_to_36 = right_elbow_direction_matrix_inverse * right_shoulder_joint_rotation^(-1) * lumbar_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    right_elbow_cor_to_right_wrist_cor_current_36 = R_world_to_36 * right_elbow_cor_to_right_wrist_cor_current_world;
    theta_36 = atan2(right_elbow_cor_to_right_wrist_cor_current_36(3), right_elbow_cor_to_right_wrist_cor_current_36(2));
    joint_angles(36) = theta_36;
    
    % right elbow internal rotation - from wrist marker
    right_wrist_cor_to_RWRB_current_world = RWRB_current - right_wrist_cor_current;
    R_36 = expAxis(right_elbow_direction_matrix(:, 1), joint_angles(36));
    R_world_to_37 = right_elbow_direction_matrix_inverse * R_36^(-1) * right_shoulder_joint_rotation^(-1) * lumbar_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    right_wrist_cor_to_RWRB_current_37 = R_world_to_37 * right_wrist_cor_to_RWRB_current_world;
    angle_now = atan2(right_wrist_cor_to_RWRB_current_37(1), -right_wrist_cor_to_RWRB_current_37(3));
    right_wrist_cor_to_RWRB_reference_world = RWRB_reference - right_wrist_cor_reference;
    right_wrist_cor_to_RWRB_reference_elbow = right_elbow_direction_matrix_inverse * right_wrist_cor_to_RWRB_reference_world;
    angle_reference = atan2(right_wrist_cor_to_RWRB_reference_elbow(1), -right_wrist_cor_to_RWRB_reference_elbow(3));
    theta_37 = angle_now - angle_reference;
    joint_angles(37) = theta_37;

    % right wrist flexion and radial-ulnar deviation angles - from finger marker
    right_wrist_cor_to_RFIN_current_world = RFIN_current - right_wrist_cor_current;
    R_37 = expAxis(-right_elbow_direction_matrix(:, 2), joint_angles(37));
    right_elbow_joint_rotation = R_36 * R_37;
    R_world_to_38 = right_wrist_direction_matrix_inverse * right_elbow_joint_rotation^(-1) * right_shoulder_joint_rotation^(-1) * lumbar_joint_rotation^(-1) * pelvis_joint_rotation^(-1);
    right_wrist_cor_to_RFIN_current_38 = R_world_to_38 * right_wrist_cor_to_RFIN_current_world;
    right_wrist_cor_to_RFIN_reference_world = RFIN_reference - right_wrist_cor_reference;
    right_wrist_cor_to_RFIN_reference_wrist = right_wrist_direction_matrix_inverse * right_wrist_cor_to_RFIN_reference_world;
    angle_now = atan2(-right_wrist_cor_to_RFIN_current_38(1), right_wrist_cor_to_RFIN_current_38(2));
    angle_reference = atan2(-right_wrist_cor_to_RFIN_reference_wrist(1), right_wrist_cor_to_RFIN_reference_wrist(2));
    theta_38 = angle_now - angle_reference;
    joint_angles(38) = theta_38;    
    
    
%     R_37
%     exp37
    
    
    return
    
    
% disp('----------------------------------------------------------------------')    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %% obsolete
    if false    
        % calculate pelvis translation and rotation
        pelvis_transformation_current_cell = calculateMcsToWcsTransformations_detailed(marker_current, marker_headers, {'PELVIS'});    
        pelvis_transformation_current = pelvis_transformation_current_cell{1};
        pelvis_rotation_current = pelvis_transformation_current(1:3, 1:3);
        pelvis_translation_current = pelvis_transformation_current(1:3, 4);

        pelvis_transformation_reference_cell = calculateMcsToWcsTransformations_detailed(marker_reference, marker_headers, {'PELVIS'});
        pelvis_transformation_reference = pelvis_transformation_reference_cell{1};
        pelvis_rotation_reference = pelvis_transformation_reference(1:3, 1:3);
        pelvis_translation_reference = pelvis_transformation_reference(1:3, 4);

        pelvis_transformation_reference_to_current = pelvis_transformation_reference^(-1) * pelvis_transformation_current;
        pelvis_rotation_reference_to_current = pelvis_transformation_reference_to_current(1:3, 1:3);
        pelvis_translation_reference_to_current = pelvis_translation_current - pelvis_translation_reference;

        joint_angles(1:3) = pelvis_translation_reference_to_current;
    %     joint_angles(4:6) = eulerAnglesFromRotationMatrixZXY(pelvis_rotation_reference_to_current); % this gives the angles relative to the pelvis frame, but I need them in world frame

        pelvis_euler_angles_current = eulerAnglesFromRotationMatrixZXY(pelvis_rotation_current);
        pelvis_euler_angles_reference = eulerAnglesFromRotationMatrixZXY(pelvis_rotation_reference);
        pelvis_euler_angles = pelvis_euler_angles_current - pelvis_euler_angles_reference;
        joint_angles(4:6) = pelvis_euler_angles;

        R_4 = expAxis([0; 0; 1], joint_angles(4));
        R_5 = expAxis([1; 0; 0], joint_angles(5));
        R_6 = expAxis([0; 1; 0], joint_angles(6));
        R_pelvis = R_4 * R_5 * R_6;

        % check whether the hip CoRs are in the correct spot
    %     load subjectModel.mat;
    %     pelvis_transformation_reference_from_tree = kinematic_tree.endEffectorTransformations{10};
    %     kinematic_tree.jointAngles = joint_angles;
    %     kinematic_tree.updateConfiguration;
    %     pelvis_transformation_current_from_tree = kinematic_tree.endEffectorTransformations{10};

    %     pelvis_transformation_current
    %     pelvis_transformation_current_from_tree
    %     pelvis_transformation_reference
    %     pelvis_transformation_reference_from_tree

    %     left_hip_cor_current
    %     left_hip_cor_pelvis = pelvis_transformation_reference^(-1) * [left_hip_cor_reference; 1];
    %     left_hip_cor_current_check = pelvis_transformation_current * left_hip_cor_pelvis
    %     
    %     right_hip_cor_current
    %     right_hip_cor_pelvis = pelvis_transformation_reference^(-1) * [right_hip_cor_reference; 1];
    %     right_hip_cor_current_check = pelvis_transformation_current * right_hip_cor_pelvis






        % calculate left hip flexion-extension
        left_hip_cor_to_left_knee_cor_current_world = left_knee_cor_current - left_hip_cor_current;
        R_world_to_7 = left_hip_direction_matrix_inverse * R_pelvis^(-1);
        left_hip_cor_to_left_knee_cor_current_7 = R_world_to_7 * left_hip_cor_to_left_knee_cor_current_world;
        theta_7 = atan2(left_hip_cor_to_left_knee_cor_current_7(2), -left_hip_cor_to_left_knee_cor_current_7(3));
        joint_angles(7) = theta_7;

    %     % left hip ab-adduction
    %     R_7_tilde = expAxis([1; 0; 0], joint_angles(7));
    %     R_world_to_8 = R_7_tilde^(-1) * left_hip_direction_matrix_inverse * R_pelvis^(-1);
    %     left_hip_cor_to_left_knee_cor_current_8 = R_world_to_8 * left_hip_cor_to_left_knee_cor_current_world;
    %     theta_8 = atan2(-left_hip_cor_to_left_knee_cor_current_8(1), -left_hip_cor_to_left_knee_cor_current_8(3));
    %     joint_angles(8) = theta_8;
    %     
    %     % left hip internal/external rotation
    %     left_knee_cor_to_LKNE_current_world = LKNE_current - left_knee_cor_current;
    %     R_8 = expAxis([0; -1; 0], joint_angles(8));
    %     R_world_to_9 = R_8^(-1) * R_7_tilde^(-1) * left_hip_direction_matrix_inverse * R_pelvis^(-1)
    %     left_knee_cor_to_LKNE_current_9 = R_world_to_9 * left_knee_cor_to_LKNE_current_world;
    %     angle_now = atan2(left_knee_cor_to_LKNE_current_9(2), -left_knee_cor_to_LKNE_current_9(1));
    %     left_knee_cor_to_LKNE_reference_world = LKNE_reference - left_knee_cor_reference;
    %     left_knee_cor_to_LKNE_reference_hip = left_hip_direction_matrix_inverse * left_knee_cor_to_LKNE_reference_world;
    %     angle_reference = atan2(left_knee_cor_to_LKNE_reference_hip(2), -left_knee_cor_to_LKNE_reference_hip(1));
    %     joint_angles(9) = angle_now - angle_reference;

        % check
    %     left_hip_direction_matrix * R_7_tilde * left_hip_direction_matrix_inverse
    %     expAxis(left_hip_direction_matrix(:, 1), joint_angles(7))
    %     left_hip_direction_matrix * R_8_tilde * left_hip_direction_matrix_inverse
    %     expAxis(-left_hip_direction_matrix(:, 2), joint_angles(8))

        % left hip ab-adduction
        R_7 = expAxis(left_hip_direction_matrix(:, 1), joint_angles(7));
        R_world_to_8 = left_hip_direction_matrix_inverse * R_7^(-1) * R_pelvis^(-1);
        left_hip_cor_to_left_knee_cor_current_8 = R_world_to_8 * left_hip_cor_to_left_knee_cor_current_world;
        theta_8 = atan2(-left_hip_cor_to_left_knee_cor_current_8(1), -left_hip_cor_to_left_knee_cor_current_8(3));
        joint_angles(8) = theta_8;

        % left knee flexion
        left_knee_cor_to_left_ankle_cor_current_world = left_ankle_cor_current - left_knee_cor_current;
        R_8 = expAxis(left_hip_direction_matrix(:, 2), joint_angles(8));
        R_world_to_10 = left_knee_direction_matrix_inverse * R_8^(-1) * R_7^(-1) * R_pelvis^(-1);
        left_knee_cor_to_left_ankle_cor_current_10 = R_world_to_10 * left_knee_cor_to_left_ankle_cor_current_world;
        joint_angles(10) = subspace(left_knee_cor_to_left_ankle_cor_current_10, [0; 0; 1]);

        % left hip internal/external rotation - from knee marker
        left_knee_cor_to_LKNE_current_world = LKNE_current - left_knee_cor_current;
        R_8 = expAxis(left_hip_direction_matrix(:, 2), joint_angles(8));
        R_world_to_9 = left_hip_direction_matrix_inverse * R_8^(-1) * R_7^(-1) * R_pelvis^(-1);
        left_knee_cor_to_LKNE_current_9 = R_world_to_9 * left_knee_cor_to_LKNE_current_world;
        angle_now = atan2(left_knee_cor_to_LKNE_current_9(2), -left_knee_cor_to_LKNE_current_9(1));
        left_knee_cor_to_LKNE_reference_world = LKNE_reference - left_knee_cor_reference;
        left_knee_cor_to_LKNE_reference_hip = left_hip_direction_matrix_inverse * left_knee_cor_to_LKNE_reference_world;
        angle_reference = atan2(left_knee_cor_to_LKNE_reference_hip(2), -left_knee_cor_to_LKNE_reference_hip(1));
        theta_9_1 = angle_now - angle_reference;
        % I don't use this so far. Check later if the calculation based on the ankle cor makes problems for almost fully
        % extended knee joints. In that case, gradually switch to this version for small knee flexion angles.

        % left hip internal/external rotation - from ankle CoR
        R_10 = expAxis(-left_knee_direction_matrix(:, 1), joint_angles(10));
        R_world_to_9 = left_hip_direction_matrix_inverse * R_8^(-1) * R_7^(-1) * R_pelvis^(-1);
        left_knee_cor_to_left_ankle_cor_current_9 = R_world_to_9 * left_knee_cor_to_left_ankle_cor_current_world;
        angle_desired = atan2(-left_knee_cor_to_left_ankle_cor_current_9(1), -left_knee_cor_to_left_ankle_cor_current_9(2)); % angle of the actual knee-to-ankle vector in hip rotation frame
        left_knee_cor_to_left_ankle_cor_reference_world = left_ankle_cor_reference - left_knee_cor_reference;
        left_knee_cor_to_left_ankle_cor_moved_world = R_pelvis * R_7 * R_8 * R_10 * left_knee_cor_to_left_ankle_cor_reference_world;
        left_knee_cor_to_left_ankle_cor_moved_9 = R_world_to_9 * left_knee_cor_to_left_ankle_cor_moved_world;
        angle_moved = atan2(-left_knee_cor_to_left_ankle_cor_moved_9(1), -left_knee_cor_to_left_ankle_cor_moved_9(2)); % angle of the knee-to-ankle vector in hip rotation frame as moved by the joint angles calculated so far (1-8, 10)
        joint_angles(9) = angle_desired - angle_moved; % moving joint 9 by the difference aligns the two vectors











        % left ankle plantar/dorsiflexion
        left_ankle_cor_to_LTOE_current_world = LTOE_current - left_ankle_cor_current;
        R_9 = expAxis(-left_hip_direction_matrix(:, 3), joint_angles(9));
        R_left_hip = R_7 * R_8 * R_9;
        R_world_to_12 = left_ankle_direction_matrix_inverse * R_10^(-1) * R_left_hip^(-1) * R_pelvis^(-1);
        left_ankle_cor_to_LTOE_current_12 = R_world_to_12 * left_ankle_cor_to_LTOE_current_world;
        angle_desired = subspace(left_ankle_cor_to_LTOE_current_12, [0; 0; 1]);

    % late night comment: can't use subspace here, I'm interested in an angle in only one plane    


        left_ankle_cor_to_LTOE_reference_world = LTOE_reference - left_ankle_cor_reference;
        left_ankle_cor_to_LTOE_moved_world = R_pelvis * R_left_hip * R_10 * left_ankle_cor_to_LTOE_reference_world;
        left_ankle_cor_to_LTOE_moved_12 = R_world_to_12 * left_ankle_cor_to_LTOE_moved_world;
        angle_moved = subspace(left_ankle_cor_to_LTOE_moved_12, [0; 0; 1]); % angle of the ankle-to-LTOE vector in ankle rotation frame as moved by the joint angles calculated so far (1-10)


        theta_11 = angle_desired - angle_moved;





        load subjectModel.mat;
        kinematic_tree.jointAngles = joint_angles;
        kinematic_tree.updateConfiguration;
        ankle_joint_moved = kinematic_tree.jointTransformations{13}(1:3, 4);
        LTOE_moved = eye(3, 4) * kinematic_tree.productsOfExponentials{13} * kinematic_tree.markerReferencePositions{13}(:, 2);

        left_ankle_cor_to_LTOE_moved_world_from_tree = LTOE_moved - ankle_joint_moved;
        left_ankle_cor_to_LTOE_moved_world_from_tree_2 = kinematic_tree.productsOfExponentials{13}(1:3, 1:3) * left_ankle_cor_to_LTOE_reference_world;
        left_ankle_cor_to_LTOE_moved_world = R_pelvis * R_left_hip * R_10 * left_ankle_cor_to_LTOE_reference_world;

        return

    %     joint_angles(11) = theta_11;








    %     
    % left_hip_flexion_axis = left_hip_direction_matrix(:, 1);
    % left_hip_abduction_axis = left_hip_direction_matrix(:, 2);
    % left_hip_internal_rotation_axis = -left_hip_direction_matrix(:, 3);    
    % left_ankle_inversion_axis = left_ankle_direction_matrix(:, 2);
    % left_ankle_dorsiflexion_axis = left_ankle_direction_matrix(:, 1);

































        % I should project the ankle joint onto the plane of rotation of the knee joint


        % no. I can figure out the flexion first, because it doesn't depend upon the hip axis. Or does it?



        % left hip internal/external rotation and knee flexion - part of this is an alternative to the above
        left_knee_cor_to_left_ankle_cor_current_world = left_ankle_cor_current - left_knee_cor_current;
        R_8 = expAxis(-left_hip_direction_matrix(:, 2), joint_angles(8));
        R_world_to_10 = left_knee_direction_matrix_inverse * R_8^(-1) * R_7^(-1) * R_pelvis^(-1);
        left_knee_cor_to_left_ankle_cor_current_10 = R_world_to_10 * left_knee_cor_to_left_ankle_cor_current_world;
        azimuth_now = cart2sph(left_knee_cor_to_left_ankle_cor_current_9(1), -left_knee_cor_to_left_ankle_cor_current_9(2), -left_knee_cor_to_left_ankle_cor_current_9(3));

        left_knee_cor_to_left_ankle_cor_reference_world = left_ankle_cor_reference - left_knee_cor_reference;
        left_knee_cor_to_LKNE_reference_hip = left_hip_direction_matrix_inverse * left_knee_cor_to_left_ankle_cor_reference_world;
        azimuth_reference = cart2sph(left_knee_cor_to_LKNE_reference_hip(1), -left_knee_cor_to_LKNE_reference_hip(2), -left_knee_cor_to_LKNE_reference_hip(3));

        azimuth_diff = azimuth_now - azimuth_reference;
    %     elevation_diff = elevation_now - elevation_reference;

    %     joint_angles(9) = azimuth_diff;
    %     joint_angles(10) = elevation_diff;









        return

        left_knee_cor_to_left_ankle_cor_current_hip = - R_hip_current * pelvis_rotation_current^(-1) * left_knee_cor_to_left_ankle_cor_current_world;
        pelvis_rotation = pelvis_transformation_reference_to_current(1:3, 1:3);
        hip_rotation_phi = expAxis(left_hip_direction_matrix(:, 1), left_hip_rotation_phi_reference_to_current);
        hip_rotation_psi = expAxis(left_hip_direction_matrix(:, 2), left_hip_rotation_psi_reference_to_current);
        hip_rotation_gamma = expAxis(left_hip_direction_matrix(:, 3), left_hip_rotation_gamma_reference_to_current);
        hip_rotation = hip_rotation_phi * hip_rotation_psi * hip_rotation_gamma;

        R_knee_phi_current = expAxis(-left_knee_direction_matrix(:, 1), left_knee_rotation_phi_reference_to_current);
        knee_rotation = R_knee_phi_current;

        R_knee = pelvis_rotation * hip_rotation * knee_rotation;
        left_knee_cor_to_left_ankle_cor_current_world_not_from_tree = R_knee * left_knee_cor_to_left_ankle_cor_reference_world;

        % take the knee-to-ankle in world frame and transform it into hip frame
        left_knee_cor_to_left_ankle_cor_current_hip_not_from_tree = - R_hip_current * pelvis_rotation_current^(-1) * left_knee_cor_to_left_ankle_cor_current_world_not_from_tree;

        % now, by how much do I have to rotate around the z-axis of the hip frame, i.e. internal rotation, to get
        % left_knee_cor_to_left_ankle_cor_current_hip_from_tree onto left_knee_cor_to_left_ankle_cor_current_hip
        hip_internal_rotation_correction = subspace([[0; 0; 1] left_knee_cor_to_left_ankle_cor_current_hip], [[0; 0; 1] left_knee_cor_to_left_ankle_cor_current_hip_not_from_tree]);

        joint_angles(9) = joint_angles(9) + hip_internal_rotation_correction;




















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

        % left hip abduction-adduction, but rotate around the correct axis
        R_hip_phi = expAxis(left_hip_direction_matrix(:, 1), left_hip_rotation_phi_reference_to_current);

        % this is the hip flexion matrix. What do I do with it?
        left_hip_cor_to_left_knee_cor_current_hip_phi = R_pelvis * R_hip_phi * left_hip_cor_to_left_knee_cor_reference_world;



        left_hip_cor_to_left_knee_cor_current_phi = R_hip_phi_current * left_hip_cor_to_left_knee_cor_current_hip;
        left_hip_cor_to_left_knee_cor_reference_phi = R_hip_phi_reference * left_hip_cor_to_left_knee_cor_reference_hip;
        left_hip_rotation_psi_current = atan2(-left_hip_cor_to_left_knee_cor_current_phi(1), -left_hip_cor_to_left_knee_cor_current_phi(3));
        left_hip_rotation_psi_reference = atan2(-left_hip_cor_to_left_knee_cor_reference_phi(1), -left_hip_cor_to_left_knee_cor_reference_phi(3));
        left_hip_rotation_psi_reference_to_current = left_hip_rotation_psi_current - left_hip_rotation_psi_reference;
        joint_angles(8) = left_hip_rotation_psi_reference_to_current;


    return

        % left hip internal rotation
        R_hip_psi_current = expAxis([0; 1; 0], -left_hip_rotation_psi_current);
        R_hip_psi_reference = expAxis([0; 1; 0], -left_hip_rotation_psi_reference);
        left_knee_cor_to_LKNE_current_world = LKNE_current - left_knee_cor_current;
        left_knee_cor_to_LKNE_reference_world = LKNE_reference - left_knee_cor_reference;
        left_knee_cor_to_LKNE_current_phipsi = R_hip_psi_current * R_hip_phi_current * left_hip_direction_matrix_inverse * pelvis_rotation_current^(-1) * left_knee_cor_to_LKNE_current_world;
        left_knee_cor_to_LKNE_reference_phipsi = R_hip_psi_reference * R_hip_phi_reference * left_hip_direction_matrix_inverse * pelvis_rotation_reference^(-1) * left_knee_cor_to_LKNE_reference_world;
        left_hip_rotation_gamma_current = atan2(-left_knee_cor_to_LKNE_current_phipsi(2), -left_knee_cor_to_LKNE_current_phipsi(1));
        left_hip_rotation_gamma_reference = atan2(-left_knee_cor_to_LKNE_reference_phipsi(2), -left_knee_cor_to_LKNE_reference_phipsi(1));
        left_hip_rotation_gamma_reference_to_current = left_hip_rotation_gamma_current - left_hip_rotation_gamma_reference;
        joint_angles(9) = -left_hip_rotation_gamma_reference_to_current;

        % left knee flexion-extension
        R_hip_gamma_current = expAxis([0; 0; 1], -left_hip_rotation_gamma_current);
        R_hip_gamma_reference = expAxis([0; 0; 1], -left_hip_rotation_gamma_reference);
        R_hip_current = R_hip_gamma_current * R_hip_psi_current * R_hip_phi_current * left_hip_direction_matrix_inverse;
        R_hip_reference = R_hip_gamma_reference * R_hip_psi_reference * R_hip_phi_reference * left_hip_direction_matrix_inverse;

        left_knee_cor_to_left_ankle_cor_current_world = left_ankle_cor_current - left_knee_cor_current;
        left_knee_cor_to_left_ankle_cor_reference_world = left_ankle_cor_reference - left_knee_cor_reference;
        left_knee_cor_to_left_ankle_cor_current_knee = left_knee_direction_matrix_inverse * R_hip_current * pelvis_rotation_current^(-1) * left_knee_cor_to_left_ankle_cor_current_world;
        left_knee_cor_to_left_ankle_cor_reference_knee = left_knee_direction_matrix_inverse * R_hip_reference * pelvis_rotation_reference^(-1) * left_knee_cor_to_left_ankle_cor_reference_world;
        left_knee_rotation_phi_current = atan2(-left_knee_cor_to_left_ankle_cor_current_knee(2), -left_knee_cor_to_left_ankle_cor_current_knee(3));
        left_knee_rotation_phi_reference = atan2(-left_knee_cor_to_left_ankle_cor_reference_knee(2), -left_knee_cor_to_left_ankle_cor_reference_knee(3));
        left_knee_rotation_phi_reference_to_current = left_knee_rotation_phi_current - left_knee_rotation_phi_reference;
        joint_angles(10) = left_knee_rotation_phi_reference_to_current;



        % change the hip internal rotation to correct for the error at the ankle CoR
        left_knee_cor_to_left_ankle_cor_current_hip = - R_hip_current * pelvis_rotation_current^(-1) * left_knee_cor_to_left_ankle_cor_current_world;
        pelvis_rotation = pelvis_transformation_reference_to_current(1:3, 1:3);
        hip_rotation_phi = expAxis(left_hip_direction_matrix(:, 1), left_hip_rotation_phi_reference_to_current);
        hip_rotation_psi = expAxis(left_hip_direction_matrix(:, 2), left_hip_rotation_psi_reference_to_current);
        hip_rotation_gamma = expAxis(left_hip_direction_matrix(:, 3), left_hip_rotation_gamma_reference_to_current);
        hip_rotation = hip_rotation_phi * hip_rotation_psi * hip_rotation_gamma;

        R_knee_phi_current = expAxis(-left_knee_direction_matrix(:, 1), left_knee_rotation_phi_reference_to_current);
        knee_rotation = R_knee_phi_current;

        R_knee = pelvis_rotation * hip_rotation * knee_rotation;
        left_knee_cor_to_left_ankle_cor_current_world_not_from_tree = R_knee * left_knee_cor_to_left_ankle_cor_reference_world;

        % take the knee-to-ankle in world frame and transform it into hip frame
        left_knee_cor_to_left_ankle_cor_current_hip_not_from_tree = - R_hip_current * pelvis_rotation_current^(-1) * left_knee_cor_to_left_ankle_cor_current_world_not_from_tree;

        % now, by how much do I have to rotate around the z-axis of the hip frame, i.e. internal rotation, to get
        % left_knee_cor_to_left_ankle_cor_current_hip_from_tree onto left_knee_cor_to_left_ankle_cor_current_hip
        hip_internal_rotation_correction = subspace([[0; 0; 1] left_knee_cor_to_left_ankle_cor_current_hip], [[0; 0; 1] left_knee_cor_to_left_ankle_cor_current_hip_not_from_tree]);

        joint_angles(9) = joint_angles(9) + hip_internal_rotation_correction;



        return

























        % experimental: correct internal rotation angle to bring the ankle of the kinematic chain CoR closer to the previously estimated one
        R_knee_phi_current = expAxis([-1; 0; 0], -left_knee_rotation_phi_current);
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

        pelvis_rotation = pelvis_transformation_reference_to_current(1:3, 1:3);

        hip_rotation_phi = expAxis(left_hip_direction_matrix(:, 1), left_hip_rotation_phi_reference_to_current);
        hip_rotation_psi = expAxis(left_hip_direction_matrix(:, 2), left_hip_rotation_psi_reference_to_current);
        hip_rotation_gamma = expAxis(left_hip_direction_matrix(:, 3), left_hip_rotation_gamma_reference_to_current);
        hip_rotation = hip_rotation_phi * hip_rotation_psi * hip_rotation_gamma;

        R_knee_phi_current = expAxis(-left_knee_direction_matrix(:, 1), left_knee_rotation_phi_reference_to_current);
        knee_rotation = R_knee_phi_current;

        R_knee = pelvis_rotation * hip_rotation * knee_rotation;
        left_knee_cor_to_left_ankle_cor_current_world_not_from_tree = R_knee_from_tree * left_knee_cor_to_left_ankle_cor_reference_world;

        % take the knee-to-ankle in world frame and transform it into hip frame
        left_knee_cor_to_left_ankle_cor_current_hip_from_tree = - R_hip_current * pelvis_rotation_current^(-1) * left_knee_cor_to_left_ankle_cor_current_world_from_tree(1:3);
        left_knee_cor_to_left_ankle_cor_current_hip_not_from_tree = - R_hip_current * pelvis_rotation_current^(-1) * left_knee_cor_to_left_ankle_cor_current_world_not_from_tree;

        % now, by how much do I have to rotate around the z-axis of the hip frame, i.e. internal rotation, to get
        % left_knee_cor_to_left_ankle_cor_current_hip_from_tree onto left_knee_cor_to_left_ankle_cor_current_hip
        hip_internal_rotation_correction = subspace([[0; 0; 1] left_knee_cor_to_left_ankle_cor_current_hip], [[0; 0; 1] left_knee_cor_to_left_ankle_cor_current_hip_from_tree]);

    %     hip_internal_rotation_correction = atan2(left_knee_cor_to_left_ankle_cor_current_knee_phi(1), -left_knee_cor_to_left_ankle_cor_current_knee_phi(3));
        joint_angles(9) = joint_angles(9) + hip_internal_rotation_correction;


        % this works! now, figure out how to calculate the above without having to update the kinematic tree object
        % I did. For now I'll keep the kinematic tree version around for testing while I develop the code









        % left knee internal rotation
        R_hip_gamma_current = expAxis([0; 0; 1], -left_hip_rotation_gamma_current);
        R_hip_gamma_reference = expAxis([0; 0; 1], -left_hip_rotation_gamma_reference);
        R_hip_current = R_hip_gamma_current * R_hip_psi_current * R_hip_phi_current * left_hip_direction_matrix_inverse;
        R_hip_reference = R_hip_gamma_reference * R_hip_psi_reference * R_hip_phi_reference * left_hip_direction_matrix_inverse;
        R_knee_phi_current = expAxis([-1; 0; 0], -left_knee_rotation_phi_current);
        R_knee_phi_reference = expAxis([-1; 0; 0], -left_knee_rotation_phi_reference);
        left_ankle_cor_to_LTOE_current_world = LTOE_current - left_ankle_cor_current;
        left_ankle_cor_to_LTOE_reference_world = LTOE_reference - left_ankle_cor_reference;
        left_ankle_cor_to_LTOE_current_knee_phi = R_knee_phi_current * left_knee_direction_matrix_inverse * R_hip_current * pelvis_rotation_current^(-1) * left_ankle_cor_to_LTOE_current_world;
        left_ankle_cor_to_LTOE_reference_knee_phi = R_knee_phi_reference * left_knee_direction_matrix_inverse * R_hip_reference * pelvis_rotation_reference^(-1) * left_ankle_cor_to_LTOE_reference_world;
        left_knee_rotation_psi_current = atan2(-left_ankle_cor_to_LTOE_current_knee_phi(1), left_ankle_cor_to_LTOE_current_knee_phi(2));
        left_knee_rotation_psi_reference = atan2(-left_ankle_cor_to_LTOE_reference_knee_phi(1), left_ankle_cor_to_LTOE_reference_knee_phi(2));
        left_knee_rotation_psi_reference_to_current = left_knee_rotation_psi_current - left_knee_rotation_psi_reference;
        joint_angles(11) = left_knee_rotation_psi_reference_to_current;

        % left ankle dorsi/plantarflexion
        R_knee_psi_current = expAxis([0; 0; 1], -left_knee_rotation_psi_current);
        R_knee_psi_reference = expAxis([0; 0; 1], -left_knee_rotation_psi_reference);
        R_knee_current = R_knee_psi_current * R_knee_phi_current * left_knee_direction_matrix_inverse;
        R_knee_reference = R_knee_psi_reference * R_knee_phi_reference * left_knee_direction_matrix_inverse;
        left_ankle_cor_to_LTOE_current_knee = R_knee_current * left_knee_direction_matrix_inverse * R_hip_current * pelvis_rotation_current^(-1) * left_ankle_cor_to_LTOE_current_world;
        left_ankle_cor_to_LTOE_reference_knee = R_knee_reference * left_knee_direction_matrix_inverse * R_hip_reference * pelvis_rotation_reference^(-1) * left_ankle_cor_to_LTOE_reference_world;
        left_ankle_rotation_phi_current = atan2(left_ankle_cor_to_LTOE_current_knee(3), left_ankle_cor_to_LTOE_current_knee(2));
        left_ankle_rotation_phi_reference = atan2(left_ankle_cor_to_LTOE_reference_knee(3), left_ankle_cor_to_LTOE_reference_knee(2));
        left_ankle_rotation_phi_reference_to_current = left_ankle_rotation_phi_current - left_ankle_rotation_phi_reference;
        joint_angles(12) = left_ankle_rotation_phi_reference_to_current;

        % left ankle inversion/eversion
        R_ankle_phi_current = expAxis([1; 0; 0], -left_ankle_rotation_phi_current);
        R_ankle_phi_reference = expAxis([1; 0; 0], -left_ankle_rotation_phi_reference);
        LTOE_to_LTOEL_current_world = LTOEL_current - LTOE_current;
        LTOE_to_LTOEL_reference_world = LTOEL_reference - LTOE_reference;
        LTOE_to_LTOEL_current_ankle = R_ankle_phi_current * left_ankle_direction_matrix_inverse * R_knee_current * left_knee_direction_matrix_inverse * R_hip_current * pelvis_rotation_current^(-1) * LTOE_to_LTOEL_current_world;
        LTOE_to_LTOEL_reference_ankle = R_ankle_phi_reference * left_ankle_direction_matrix_inverse * R_knee_reference * left_knee_direction_matrix_inverse * R_hip_reference * pelvis_rotation_reference^(-1) * LTOE_to_LTOEL_reference_world;
        left_ankle_rotation_psi_current = -atan2(LTOE_to_LTOEL_current_ankle(3), -LTOE_to_LTOEL_current_ankle(1));
        left_ankle_rotation_psi_reference = -atan2(LTOE_to_LTOEL_reference_ankle(3), -LTOE_to_LTOEL_reference_ankle(1));
        left_ankle_rotation_psi_reference_to_current = left_ankle_rotation_psi_current - left_ankle_rotation_psi_reference;
        joint_angles(13) = left_ankle_rotation_psi_reference_to_current;

        % status: I'm not completely convinced that this is right, but there seems to be a small error in the knee, so let's
        % fix that one first. That might be because I'm not using the right rotation matrices somewhere along the line, i.e.
        % not correcting for joint directions


























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










    end
    
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

