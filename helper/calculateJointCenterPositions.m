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

    
    % calculate transformations for segments that are fully determined by markers
    transformations_reference = calculateMcsToWcsTransformations_new(marker_reference, marker_headers, segment_labels);
    transformations_current = calculateMcsToWcsTransformations_new(marker_positions, marker_headers, segment_labels);    
    
    
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

% replace this with a version that uses the forearm instead of the hand markers
%     % left wrist
%     T_reference_wcs_to_hand_mcs = transformations_reference{strcmp(segment_labels, 'LHAND')}^(-1);
%     lwrist_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LWRISTCOR')';
%     lwrist_cor_hand_mcs = eye(3, 4) * T_reference_wcs_to_hand_mcs * [lwrist_cor_reference_wcs; 1];
%     T_current_hand_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'LHAND')};
%     lwrist_cor_current_wcs = eye(3, 4) * T_current_hand_mcs_to_wcs * [lwrist_cor_hand_mcs; 1];
%     
%     % right wrist
%     T_reference_wcs_to_hand_mcs = transformations_reference{strcmp(segment_labels, 'RHAND')}^(-1);
%     rwrist_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RWRISTCOR')';
%     rwrist_cor_hand_mcs = eye(3, 4) * T_reference_wcs_to_hand_mcs * [rwrist_cor_reference_wcs; 1];
%     T_current_hand_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'RHAND')};
%     rwrist_cor_current_wcs = eye(3, 4) * T_current_hand_mcs_to_wcs * [rwrist_cor_hand_mcs; 1];

    % left wrist
    T_reference_wcs_to_forearm_mcs = transformations_reference{strcmp(segment_labels, 'LFOREARM')}^(-1); % debugging: this does not depend on the wrist angle
    lwrist_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LWRISTCOR')'; % debugging: this does not depend on the wrist angle
    lwrist_cor_forearm_mcs = eye(3, 4) * T_reference_wcs_to_forearm_mcs * [lwrist_cor_reference_wcs; 1]; % debugging: this does not depend on the wrist angle
    T_current_forearm_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'LFOREARM')};
    lwrist_cor_current_wcs = eye(3, 4) * T_current_forearm_mcs_to_wcs * [lwrist_cor_forearm_mcs; 1];
    
    % right wrist
    T_reference_wcs_to_forearm_mcs = transformations_reference{strcmp(segment_labels, 'RFOREARM')}^(-1);
    rwrist_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RWRISTCOR')';
    rwrist_cor_forearm_mcs = eye(3, 4) * T_reference_wcs_to_forearm_mcs * [rwrist_cor_reference_wcs; 1];
    T_current_forearm_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'RFOREARM')};
    rwrist_cor_current_wcs = eye(3, 4) * T_current_forearm_mcs_to_wcs * [rwrist_cor_forearm_mcs; 1]; 
    
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
    
    % left toes eef
    T_reference_wcs_to_foot_mcs = transformations_reference{strcmp(segment_labels, 'LFOOT')}^(-1);
    ltoes_eef_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LTOESEEF')';
    ltoes_eef_foot_mcs = eye(3, 4) * T_reference_wcs_to_foot_mcs * [ltoes_eef_reference_wcs; 1];
    T_current_foot_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'LFOOT')};
    ltoes_eef_current_wcs = eye(3, 4) * T_current_foot_mcs_to_wcs * [ltoes_eef_foot_mcs; 1];
    
    % right toes eef
    T_reference_wcs_to_foot_mcs = transformations_reference{strcmp(segment_labels, 'RFOOT')}^(-1);
    rtoes_eef_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RTOESEEF')';
    rtoes_eef_foot_mcs = eye(3, 4) * T_reference_wcs_to_foot_mcs * [rtoes_eef_reference_wcs; 1];
    T_current_foot_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'RFOOT')};
    rtoes_eef_current_wcs = eye(3, 4) * T_current_foot_mcs_to_wcs * [rtoes_eef_foot_mcs; 1];
    
    % left hand eef
    T_reference_wcs_to_hand_mcs = transformations_reference{strcmp(segment_labels, 'LHAND')}^(-1);
    lhand_eef_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LHANDEEF')';
    lhand_eef_hand_mcs = eye(3, 4) * T_reference_wcs_to_hand_mcs * [lhand_eef_reference_wcs; 1];
    T_current_hand_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'LHAND')};
    lhand_eef_current_wcs = eye(3, 4) * T_current_hand_mcs_to_wcs * [lhand_eef_hand_mcs; 1];
    
    % right hand eef
    T_reference_wcs_to_hand_mcs = transformations_reference{strcmp(segment_labels, 'RHAND')}^(-1);
    rhand_eef_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RHANDEEF')';
    rhand_eef_hand_mcs = eye(3, 4) * T_reference_wcs_to_hand_mcs * [rhand_eef_reference_wcs; 1];
    T_current_hand_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'RHAND')};
    rhand_eef_current_wcs = eye(3, 4) * T_current_hand_mcs_to_wcs * [rhand_eef_hand_mcs; 1];
    
    % left elbow
    T_reference_wcs_to_forearm_mcs = transformations_reference{strcmp(segment_labels, 'LFOREARM')}^(-1);
    lelbow_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LELBOWCOR')';
    lelbow_cor_forearm_mcs = eye(3, 4) * T_reference_wcs_to_forearm_mcs * [lelbow_cor_reference_wcs; 1];
    T_current_forearm_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'LFOREARM')};
    lelbow_cor_current_wcs = eye(3, 4) * T_current_forearm_mcs_to_wcs * [lelbow_cor_forearm_mcs; 1];
    
    % right elbow
    T_reference_wcs_to_forearm_mcs = transformations_reference{strcmp(segment_labels, 'RFOREARM')}^(-1);
    relbow_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RELBOWCOR')';
    relbow_cor_forearm_mcs = eye(3, 4) * T_reference_wcs_to_forearm_mcs * [relbow_cor_reference_wcs; 1];
    T_current_forearm_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'RFOREARM')};
    relbow_cor_current_wcs = eye(3, 4) * T_current_forearm_mcs_to_wcs * [relbow_cor_forearm_mcs; 1]; 
    
    
    % exploratory: calculate knees and elbows also based on markers
    markers_and_cor_reference = [marker_reference lhip_cor_reference_wcs' rhip_cor_reference_wcs'];
    markers_and_cor_positions = [marker_positions lhip_cor_current_wcs' rhip_cor_current_wcs'];
    markers_and_cor_labels = [marker_headers {'LHIPCOR', 'RHIPCOR'}];
    segment_labels = {'LTHIGH', 'RTHIGH'};
    
    transformations_reference = calculateMcsToWcsTransformations_new(markers_and_cor_reference, markers_and_cor_labels, segment_labels);
    transformations_current = calculateMcsToWcsTransformations_new(markers_and_cor_positions, markers_and_cor_labels, segment_labels);    
    
    % left knee
    T_reference_wcs_to_thigh_mcs = transformations_reference{strcmp(segment_labels, 'LTHIGH')}^(-1);
    lknee_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LKNEECOR')'; % this seems wrong
    lknee_cor_thigh_mcs = eye(3, 4) * T_reference_wcs_to_thigh_mcs * [lknee_cor_reference_wcs; 1];
    T_current_thigh_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'LTHIGH')};
    lknee_cor_current_wcs = eye(3, 4) * T_current_thigh_mcs_to_wcs * [lknee_cor_thigh_mcs; 1];
    
    % right knee
    T_reference_wcs_to_thigh_mcs = transformations_reference{strcmp(segment_labels, 'RTHIGH')}^(-1);
    rknee_cor_reference_wcs = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RKNEECOR')';
    rknee_cor_thigh_mcs = eye(3, 4) * T_reference_wcs_to_thigh_mcs * [rknee_cor_reference_wcs; 1];
    T_current_thigh_mcs_to_wcs = transformations_current{strcmp(segment_labels, 'RTHIGH')};
    rknee_cor_current_wcs = eye(3, 4) * T_current_thigh_mcs_to_wcs * [rknee_cor_thigh_mcs; 1];
    
    
    
    
    

    
    
    
    
    % TODO check whether these values for lknee_cor_current_wcs and rknee_cor_current_wcs are correct
    % they are
    
    
    % calculate joint centers for segments that are NOT fully determined by markers
    

% the following is being replaced with what's above, using the LTHI and RTHI markers
if false
    % left knee
    LANK_current = extractMarkerTrajectories(marker_positions, marker_headers, 'LANK')'; % left ankle marker
    lknee_axis_direction = normVector(lankle_cor_current_wcs - LANK_current); 
    % knee axis direction, defined to be the same as the ankle axis direction
    % but that isn't true, because I have an additional joint for knee internal rotation
    
    
    
    % that means I can't use the lower leg to calculate the knee joint center position, because I only have to markers on the segment, ankle and anklecor
    % I also can't use the thigh, for the same reason, only knee marker and hipcor
    
    % so what do I do? Two options
    
    % 1. use the circle-approach that I developed a while ago, which turned out to be problematic
    
    % 2. combine inverse kinematics with the calculation of the joint angles
    % it should be possible to determine the hip and knee joint angles without using the knee cor - shouldn't it?
    
    % 3. option
    % use the thigh and shank markers
    % define thigh by 
    
    
    LKNE_current = extractMarkerTrajectories(marker_positions, marker_headers, 'LKNE')';
    LKNE_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LKNE')';
    LKNEECOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LKNEECOR')';
    lknee_cor_to_marker = norm(LKNE_reference - LKNEECOR_reference);
    lknee_cor_current_wcs = LKNE_current + lknee_cor_to_marker * lknee_axis_direction;
    lknee_cor_current_wcs_from_markers = lknee_cor_current_wcs;
    
    % right knee
    RANK_current = extractMarkerTrajectories(marker_positions, marker_headers, 'RANK')';
    rknee_axis_direction = normVector(rankle_cor_current_wcs - RANK_current);
    RKNE_current = extractMarkerTrajectories(marker_positions, marker_headers, 'RKNE')';
    RKNE_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RKNE')';
    RKNEECOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RKNEECOR')';
    rknee_cor_to_marker = norm(RKNE_reference - RKNEECOR_reference);
    rknee_cor_current_wcs = RKNE_current + rknee_cor_to_marker * rknee_axis_direction;
    
    % what follows uses findHingeJointCenter, i.e. the circle method
    % left knee
    LKNE_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LKNE')';
    LHIPCOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LHIPCOR')';
    LKNEECOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LKNEECOR')';
    LANKLECOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LANKLECOR')';
    LKNE_current = extractMarkerTrajectories(marker_positions, marker_headers, 'LKNE')';
    r_marker = norm(LKNE_reference - LKNEECOR_reference);
    r_hip = norm(LKNE_reference - LHIPCOR_reference);
    r_ankle = norm(LKNE_reference - LANKLECOR_reference);
    
    lknee_cor_current_wcs = findHingeJointCenter(lhip_cor_current_wcs, lankle_cor_current_wcs, LKNE_current, r_marker, r_hip, r_ankle, 'positive');
    
    lknee_cor_current_wcs_from_findHingeJointCenter = lknee_cor_current_wcs
    
    % right knee
    RKNE_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RKNE')';
    RHIPCOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RHIPCOR')';
    RKNEECOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RKNEECOR')';
    RANKLECOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RANKLECOR')';
    RKNE_current = extractMarkerTrajectories(marker_positions, marker_headers, 'RKNE')';
    r_marker = norm(RKNE_reference - RKNEECOR_reference);
    r_hip = norm(RKNE_reference - RHIPCOR_reference);
    r_ankle = norm(RKNE_reference - RANKLECOR_reference);
    rknee_cor_current_wcs = findHingeJointCenter(rhip_cor_current_wcs, rankle_cor_current_wcs, RKNE_current, r_marker, r_hip, r_ankle, 'negative');    
    
    knee_to_marker_distance_left = norm(lknee_cor_current_wcs - LKNE_current);
    knee_to_marker_distance_right = norm(rknee_cor_current_wcs - RKNE_current);
    
    % left elbow
    LELB_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LELB')';
    LSHOULDERCOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LSHOULDERCOR')';
    LELBOWCOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LELBOWCOR')';
    LWRISTCOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'LWRISTCOR')';
    LELB_current = extractMarkerTrajectories(marker_positions, marker_headers, 'LELB')';
    r_marker = norm(LELB_reference - LELBOWCOR_reference);
    r_shoulder = norm(LELB_reference - LSHOULDERCOR_reference);
    r_wrist = norm(LELB_reference - LWRISTCOR_reference);
%     lelbow_cor_current_wcs = findHingeJointCenter(lshoulder_cor_current_wcs, lwrist_cor_current_wcs, LELB_current, r_marker, r_shoulder, r_wrist, 'negative');
    lelbow_cor_current_wcs = zeros(3, 1);
    
    % right elbow
    RELB_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RELB')';
    RSHOULDERCOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RSHOULDERCOR')';
    RELBOWCOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RELBOWCOR')';
    RWRISTCOR_reference = extractMarkerTrajectories(joint_center_reference, joint_center_headers, 'RWRISTCOR')';
    RELB_current = extractMarkerTrajectories(marker_positions, marker_headers, 'RELB')';
    r_marker = norm(RELB_reference - RELBOWCOR_reference);
    r_shoulder = norm(RELB_reference - RSHOULDERCOR_reference);
    r_wrist = norm(RELB_reference - RWRISTCOR_reference);
%     relbow_cor_current_wcs = findHingeJointCenter(rshoulder_cor_current_wcs, rwrist_cor_current_wcs, RELB_current, r_marker, r_shoulder, r_wrist, 'positive');
    relbow_cor_current_wcs = zeros(3, 1);
    
end
    
    

    
    
    
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
        rankle_cor_current_wcs', ...
        ltoes_eef_current_wcs', ...
        rtoes_eef_current_wcs', ...
        lhand_eef_current_wcs', ...
        rhand_eef_current_wcs' ...
      ];
    
    
end
        
        