function joint_center_positions = ...
    calculateJointCenterPositions ...
      ( ...
        marker_reference, ...
        marker_positions, ...
        marker_headers, ...
        joint_center_reference, ...
        joint_center_headers ...
      )
  
  
    segment_labels = ...
    { ...
      'HEAD', ...
      'TORSO', ...
      'LFOREARM', ...
      'RFOREARM', ...
      'LHAND', ...
      'RHAND', ...
      'PELVIS', ...
      'LFOOT', ...
      'RFOOT' ...
    };    

%     % calculate transformations for segments that are fully determined by markers
%     markers_by_segment = ...
%       {
%         'RFHD', 'LFHD', 'RBHD'; ...                     % head
%         'C7', 'CLAV', 'T10'; ...                        % torso
%         'LELB', 'LWRA', 'LWRB'; ...                     % left forearm
%         'RELB', 'RWRA', 'RWRB'; ...                     % right forearm
%         'LWRA', 'LWRB', 'LFIN'; ...                     % left hand
%         'RWRA', 'RWRB', 'RFIN'; ...                     % right hand
%         'RASI', 'LASI', 'RPSI'; ...                     % pelvis
%         'LANK', 'LHEE', 'LTOE'; ...                     % left foot
%         'RANK', 'RHEE', 'RTOE'; ...                     % right foot
%       };
%     transformations_reference = calculateMcsToWcsTransformations(marker_reference, marker_headers, markers_by_segment);
%     transformations_current = calculateMcsToWcsTransformations(marker_positions, marker_headers, markers_by_segment);
    
    transformations_reference = calculateMcsToWcsTransformations_detailed(marker_reference, marker_headers, segment_labels);
    transformations_current = calculateMcsToWcsTransformations_detailed(marker_positions, marker_headers, segment_labels);    
    
    
    % cervix
    T_reference_wcs_to_head_mcs = transformations_reference{strcmp(segment_labels, 'HEAD')}^(-1);
    cervix_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'CERVIXCOR')';
    cervix_cor_head_mcs = eye(3, 4) * T_reference_wcs_to_head_mcs * [cervix_cor_reference_wcs; 1];
    T_current_head_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'HEAD')};
    cervix_cor_current_wcs_from_head = eye(3, 4) * T_current_head_mcs_to_wcs * [cervix_cor_head_mcs; 1];
    
    T_reference_wcs_to_torso_mcs = transformations_reference{strcmp(segment_labels, 'TORSO')}^(-1);
    cervix_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'CERVIXCOR')';
    cervix_cor_torso_mcs = eye(3, 4) * T_reference_wcs_to_torso_mcs * [cervix_cor_reference_wcs; 1];
    T_current_torso_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'TORSO')};
    cervix_cor_current_wcs_from_torso = eye(3, 4) * T_current_torso_mcs_to_wcs * [cervix_cor_torso_mcs; 1];
  
    cervix_cor_current_wcs = mean([cervix_cor_current_wcs_from_head cervix_cor_current_wcs_from_torso], 2);
    
    % left shoulder
    T_reference_wcs_to_torso_mcs = transformations_reference{strcmp(segment_labels, 'TORSO')}^(-1);
    lshoulder_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LSHOULDERCOR')';
    lshoulder_cor_torso_mcs = eye(3, 4) * T_reference_wcs_to_torso_mcs * [lshoulder_cor_reference_wcs; 1];
    T_current_torso_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'TORSO')};
    lshoulder_cor_current_wcs = eye(3, 4) * T_current_torso_mcs_to_wcs * [lshoulder_cor_torso_mcs; 1];
    
    % right shoulder
    T_reference_wcs_to_torso_mcs = transformations_reference{strcmp(segment_labels, 'TORSO')}^(-1);
    rshoulder_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RSHOULDERCOR')';
    rshoulder_cor_torso_mcs = eye(3, 4) * T_reference_wcs_to_torso_mcs * [rshoulder_cor_reference_wcs; 1];
    T_current_torso_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'TORSO')};
    rshoulder_cor_current_wcs = eye(3, 4) * T_current_torso_mcs_to_wcs * [rshoulder_cor_torso_mcs; 1];
    
    % left wrist
    T_reference_wcs_to_hand_mcs = transformations_reference{strcmp(segment_labels, 'LHAND')}^(-1);
    lwrist_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LWRISTCOR')';
    lwrist_cor_hand_mcs = eye(3, 4) * T_reference_wcs_to_hand_mcs * [lwrist_cor_reference_wcs; 1];
    T_current_hand_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'LHAND')};
    lwrist_cor_current_wcs = eye(3, 4) * T_current_hand_mcs_to_wcs * [lwrist_cor_hand_mcs; 1];
    
    % right wrist
    T_reference_wcs_to_hand_mcs = transformations_reference{strcmp(segment_labels, 'RHAND')}^(-1);
    rwrist_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RWRISTCOR')';
    rwrist_cor_hand_mcs = eye(3, 4) * T_reference_wcs_to_hand_mcs * [rwrist_cor_reference_wcs; 1];
    T_current_hand_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'RHAND')};
    rwrist_cor_current_wcs = eye(3, 4) * T_current_hand_mcs_to_wcs * [rwrist_cor_hand_mcs; 1];
    
    % lumbar
    T_reference_wcs_to_pelvis_mcs = transformations_reference{strcmp(segment_labels, 'PELVIS')}^(-1);
    lumbar_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LUMBARCOR')';
    lumbar_cor_pelvis_mcs = eye(3, 4) * T_reference_wcs_to_pelvis_mcs * [lumbar_cor_reference_wcs; 1];
    T_current_pelvis_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'PELVIS')};
    lumbar_cor_current_wcs_from_pelvis = eye(3, 4) * T_current_pelvis_mcs_to_wcs * [lumbar_cor_pelvis_mcs; 1];
    
    lumbar_cor_current_wcs = lumbar_cor_current_wcs_from_pelvis;
    
%     T_reference_wcs_to_torso_mcs = transformations_reference{strcmp(segment_labels, 'TORSO')}^(-1);
%     lumbar_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LUMBARCOR')';
%     lumbar_cor_torso_mcs = eye(3, 4) * T_reference_wcs_to_torso_mcs * [lumbar_cor_reference_wcs; 1];
%     T_current_torso_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'TORSO')};
%     lumbar_cor_current_wcs_from_torso = eye(3, 4) * T_current_torso_mcs_to_wcs * [lumbar_cor_torso_mcs; 1];
%   
%     lumbar_cor_current_wcs = mean([lumbar_cor_current_wcs_from_pelvis lumbar_cor_current_wcs_from_torso], 2);
    
    % left hip
    T_reference_wcs_to_pelvis_mcs = transformations_reference{strcmp(segment_labels, 'PELVIS')}^(-1);
    lhip_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LHIPCOR')';
    lhip_cor_pelvis_mcs = eye(3, 4) * T_reference_wcs_to_pelvis_mcs * [lhip_cor_reference_wcs; 1];
    T_current_pelvis_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'PELVIS')};
    lhip_cor_current_wcs = eye(3, 4) * T_current_pelvis_mcs_to_wcs * [lhip_cor_pelvis_mcs; 1];
    
    % right hip
    T_reference_wcs_to_pelvis_mcs = transformations_reference{strcmp(segment_labels, 'PELVIS')}^(-1);
    rhip_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RHIPCOR')';
    rhip_cor_pelvis_mcs = eye(3, 4) * T_reference_wcs_to_pelvis_mcs * [rhip_cor_reference_wcs; 1];
    T_current_pelvis_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'PELVIS')};
    rhip_cor_current_wcs = eye(3, 4) * T_current_pelvis_mcs_to_wcs * [rhip_cor_pelvis_mcs; 1];
        
    % left ankle
    T_reference_wcs_to_foot_mcs = transformations_reference{strcmp(segment_labels, 'LFOOT')}^(-1);
    lankle_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LANKLECOR')';
    lankle_cor_foot_mcs = eye(3, 4) * T_reference_wcs_to_foot_mcs * [lankle_cor_reference_wcs; 1];
    T_current_foot_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'LFOOT')};
    lankle_cor_current_wcs = eye(3, 4) * T_current_foot_mcs_to_wcs * [lankle_cor_foot_mcs; 1];
    
    % right ankle
    T_reference_wcs_to_foot_mcs = transformations_reference{strcmp(segment_labels, 'RFOOT')}^(-1);
    rankle_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RANKLECOR')';
    rankle_cor_foot_mcs = eye(3, 4) * T_reference_wcs_to_foot_mcs * [rankle_cor_reference_wcs; 1];
    T_current_foot_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'RFOOT')};
    rankle_cor_current_wcs = eye(3, 4) * T_current_foot_mcs_to_wcs * [rankle_cor_foot_mcs; 1];
    
    % left elbow
    LELB_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LELB')';
    LSHOULDERCOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LSHOULDERCOR')';
    LELBOWCOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LELBOWCOR')';
    LWRISTCOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LWRISTCOR')';
    LELB_current = extractMarkerTrajectories(marker_positions, marker_headers, 'LELB')';
    r_marker = norm(LELB_reference - LELBOWCOR_reference);
    r_shoulder = norm(LELB_reference - LSHOULDERCOR_reference);
    r_wrist = norm(LELB_reference - LWRISTCOR_reference);
    lelbow_cor_current_wcs = findHingeJointCenter(lshoulder_cor_current_wcs, lwrist_cor_current_wcs, LELB_current, r_marker, r_shoulder, r_wrist, 'negative');
    
    % right elbow
    RELB_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RELB')';
    RSHOULDERCOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RSHOULDERCOR')';
    RELBOWCOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RELBOWCOR')';
    RWRISTCOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RWRISTCOR')';
    RELB_current = extractMarkerTrajectories(marker_positions, marker_headers, 'RELB')';
    r_marker = norm(RELB_reference - RELBOWCOR_reference);
    r_shoulder = norm(RELB_reference - RSHOULDERCOR_reference);
    r_wrist = norm(RELB_reference - RWRISTCOR_reference);
    relbow_cor_current_wcs = findHingeJointCenter(rshoulder_cor_current_wcs, rwrist_cor_current_wcs, RELB_current, r_marker, r_shoulder, r_wrist, 'positive');

    % left knee
    LANK_current = extractMarkerTrajectories(marker_positions, marker_headers, 'LANK')';
    lknee_axis_direction = normVector(lankle_cor_current_wcs - LANK_current);
    LKNE_current = extractMarkerTrajectories(marker_positions, marker_headers, 'LKNE')';
    LKNE_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LKNE')';
    LKNEECOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LKNEECOR')';
    lknee_cor_to_marker = norm(LKNE_reference - LKNEECOR_reference);
    lknee_cor_current_wcs = LKNE_current + lknee_cor_to_marker * lknee_axis_direction;
    
    % right knee
    RANK_current = extractMarkerTrajectories(marker_positions, marker_headers, 'RANK')';
    rknee_axis_direction = normVector(rankle_cor_current_wcs - RANK_current);
    RKNE_current = extractMarkerTrajectories(marker_positions, marker_headers, 'RKNE')';
    RKNE_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RKNE')';
    RKNEECOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RKNEECOR')';
    rknee_cor_to_marker = norm(RKNE_reference - RKNEECOR_reference);
    rknee_cor_current_wcs = RKNE_current + rknee_cor_to_marker * rknee_axis_direction;
    
    % combine results
    joint_center_positions = ...
      [ ...
        cervix_cor_current_wcs' ...
        lshoulder_cor_current_wcs' ...
        rshoulder_cor_current_wcs' ...
        lelbow_cor_current_wcs', ...
        relbow_cor_current_wcs', ...
        lwrist_cor_current_wcs', ...
        rwrist_cor_current_wcs', ...
        lumbar_cor_current_wcs', ...
        lhip_cor_current_wcs', ...
        rhip_cor_current_wcs', ...
        lknee_cor_current_wcs', ...
        rknee_cor_current_wcs', ...
        lankle_cor_current_wcs', ...
        rankle_cor_current_wcs' ...
      ];
    
    
end
        
        