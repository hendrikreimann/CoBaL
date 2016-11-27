
static_reference_type = 'casual';
% static_reference_type = 'anatomic';
% static_reference_type = 'motorcycle';
% static_reference_type = 'ski';

% hip_joint_center_estimation_method = 'SCoRE';
hip_joint_center_estimation_method = 'Tylkowski';

% knee_joint_axis_estimation_method = 'SARA';
knee_joint_axis_estimation_method = 'markers';

create_kinematic_tree                   = 1;
show_visualization                      = 1;

static_reference_trial_type = 'walking';
static_reference_file_index = 1;

left_hip_calibration_file_index = 2;
right_hip_calibration_file_index = 3;
left_knee_calibration_file_index = 4;
right_knee_calibration_file_index = 5;

load subjectInfo.mat;

% if width measurements are not available, use best guess
if knee_width == 0
    knee_width = 0.07;
end
if ankle_width == 0
    ankle_width = 0.07;
end
if elbow_width == 0
    elbow_width = 0.06;
end


%% create static reference

% load static reference file
load(['processed' filesep makeFileName(date, subject_id, static_reference_trial_type, static_reference_file_index, 'markerTrajectories')]);

% find first time step where all markers are available
i_time = 1;
while any(isnan(marker_trajectories(i_time, :)))
    i_time = i_time + 1;
end
marker_reference = marker_trajectories(i_time, :);

%% extract marker reference positions

% head
LFHD_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LFHD')';
RFHD_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RFHD')';
LBHD_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LBHD')';
RBHD_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RBHD')';

% torso
C7_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'C7')';
T10_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'T10')';
CLAV_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'CLAV')';
STRN_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'STRN')';
RBAK_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RBAK')';

% left arm
LSHO_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LSHO')';
LELB_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LELB')';
LWRA_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LWRA')';
LWRB_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LWRB')';
LFIN_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LFIN')';

% right arm
RSHO_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RSHO')';
RELB_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RELB')';
RWRA_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RWRA')';
RWRB_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RWRB')';
RFIN_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RFIN')';

% pelvis
LASI_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LASI')';
RASI_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RASI')';
LPSI_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LPSI')';
RPSI_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RPSI')';

% left leg
LTHI_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LTHI')';
LTHIA_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LTHIA')';
LKNE_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LKNE')';
LTIB_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LTIB')';
LTIBA_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LTIBA')';
LANK_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LANK')';
LHEE_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LHEE')';
LTOE_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LTOE')';
LTOEL_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'LTOEL')';

% left leg
RTHI_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RTHI')';
RTHIA_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RTHIA')';
RKNE_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RKNE')';
RTIB_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RTIB')';
RTIBA_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RTIBA')';
RANK_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RANK')';
RHEE_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RHEE')';
RTOE_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RTOE')';
RTOEL_reference = extractMarkerTrajectories(marker_reference, marker_headers, 'RTOEL')';

% groups
head_markers_reference = [LFHD_reference' RFHD_reference' LBHD_reference' RBHD_reference'];
trunk_markers_reference = [C7_reference' T10_reference' CLAV_reference' STRN_reference' RBAK_reference'];
pelvis_markers_reference = [LASI_reference' RASI_reference' LPSI_reference' RPSI_reference'];
left_thigh_markers_reference = [LTHI_reference' LTHIA_reference' LKNE_reference'];
left_shank_markers_reference = [LTIB_reference' LTIBA_reference' LANK_reference'];
left_foot_markers_reference = [LHEE_reference' LTOE_reference' LTOEL_reference'];
right_thigh_markers_reference = [RTHI_reference' RTHIA_reference' RKNE_reference'];
right_shank_markers_reference = [RTIB_reference' RTIBA_reference' RANK_reference'];
right_foot_markers_reference = [RHEE_reference' RTOE_reference' RTOEL_reference'];

%% calculate anatomical points
MASIS_reference = mean([LASI_reference RASI_reference], 2);
MPSIS_reference = mean([LPSI_reference RPSI_reference], 2);

% define directions
e_1 = [1; 0; 0];
e_2 = [0; 1; 0];
e_3 = [0; 0; 1];
anterior_direction = normVector(MASIS_reference - MPSIS_reference);
posterior_direction = - anterior_direction;
left_direction = normVector(LASI_reference - RASI_reference);
right_direction = - left_direction;
proximal_direction = cross(right_direction, anterior_direction);
distal_direction = - proximal_direction;

% calculate anatomical landmarks and apply marker and flesh offsets
centroid_to_skin_correction = 0.0152; % in meters
skin_to_bone_correction_ASIS = 0.01; % in meters
skin_to_bone_correction_PSIS = 0.01; % in meters
c7 = C7_reference + centroid_to_skin_correction * anterior_direction;
suprasternale = CLAV_reference + centroid_to_skin_correction * posterior_direction;
left_acromion = LSHO_reference + centroid_to_skin_correction * distal_direction;
left_lateral_humeral_epicondyle = LELB_reference + centroid_to_skin_correction * distal_direction;
left_inner_wrist = LWRA_reference + centroid_to_skin_correction * right_direction;
left_outer_wrist = LWRB_reference + centroid_to_skin_correction * left_direction;
left_hand = LFIN_reference + centroid_to_skin_correction * distal_direction;
right_inner_wrist = RWRA_reference + centroid_to_skin_correction * left_direction;
right_outer_wrist = RWRB_reference + centroid_to_skin_correction * right_direction;
right_hand = RFIN_reference + centroid_to_skin_correction * distal_direction;
right_acromion = RSHO_reference + centroid_to_skin_correction * distal_direction;
right_lateral_humeral_epicondyle = RELB_reference + centroid_to_skin_correction * distal_direction;
LASIS = LASI_reference + (centroid_to_skin_correction + skin_to_bone_correction_ASIS) * posterior_direction;
RASIS = RASI_reference + (centroid_to_skin_correction + skin_to_bone_correction_ASIS) * posterior_direction;
LPSIS = LPSI_reference + (centroid_to_skin_correction + skin_to_bone_correction_PSIS) * anterior_direction;
RPSIS = RPSI_reference + (centroid_to_skin_correction + skin_to_bone_correction_PSIS) * anterior_direction;
left_lateral_femoral_epicondyle = LKNE_reference + centroid_to_skin_correction * right_direction;
left_lateral_malleolus = LANK_reference + centroid_to_skin_correction * right_direction;
left_calcaneus = LHEE_reference + centroid_to_skin_correction * anterior_direction;
left_toe_mid = LTOE_reference + centroid_to_skin_correction * distal_direction;
right_lateral_femoral_epicondyle = RKNE_reference + centroid_to_skin_correction * left_direction;
right_lateral_malleolus = RANK_reference + centroid_to_skin_correction * left_direction;
right_calcaneus = RHEE_reference + centroid_to_skin_correction * anterior_direction;
right_toe_mid = RTOE_reference + centroid_to_skin_correction * distal_direction;
head_width = norm(LFHD_reference - RFHD_reference);
head_center_to_vertex = head_width; % this is an ad-hoc assumption that seems to work out well in some examples
head_vertex = mean([LFHD_reference RFHD_reference LBHD_reference RBHD_reference], 2) + head_center_to_vertex * proximal_direction;
sellion = mean([LFHD_reference RFHD_reference], 2);

% calculate some distances
inter_ASIS_distance = norm(LASIS - RASIS);

%% estimate hip joint centers
if strcmp(hip_joint_center_estimation_method, 'SCoRE')
    pelvis_center_reference = mean(reshape(pelvis_markers_reference, 3, size(pelvis_markers_reference, 2)/3), 2);

    % find left hip CoR
    left_hip_reference_file_name = ['processed' filesep makeFileName(date, subject_id, 'calibration', left_hip_calibration_file_index, 'markerTrajectories')];
    disp(['Left hip reference file name: ' left_hip_reference_file_name]);
    load(left_hip_reference_file_name);
    hip_reference = marker_trajectories;

    LASI_trajectory = extractMarkerTrajectories(hip_reference, marker_headers, 'LASI')';
    RASI_trajectory = extractMarkerTrajectories(hip_reference, marker_headers, 'RASI')';
    LPSI_trajectory = extractMarkerTrajectories(hip_reference, marker_headers, 'LPSI')';
    RPSI_trajectory = extractMarkerTrajectories(hip_reference, marker_headers, 'RPSI')';
    LTHI_trajectory = extractMarkerTrajectories(hip_reference, marker_headers, 'LTHI')';
    LTHIA_trajectory = extractMarkerTrajectories(hip_reference, marker_headers, 'LTHIA')';
    LKNE_trajectory = extractMarkerTrajectories(hip_reference, marker_headers, 'LKNE')';
    pelvis_markers_trajectory = [LASI_trajectory' RASI_trajectory' LPSI_trajectory' RPSI_trajectory'];
    left_thigh_markers_trajectory = [LTHI_trajectory' LTHIA_trajectory' LKNE_trajectory'];

%     pelvis_markers_trajectory = hip_reference(:, pelvis_markers_indices);
%     left_thigh_markers_trajectory = hip_reference(:, left_thigh_markers_indices);
    [left_hip_cor, left_hip_cor_error] = estimateJointKinematics ...
      ( ...
        pelvis_markers_reference, ...
        left_thigh_markers_reference, ...
        pelvis_markers_trajectory, ...
        left_thigh_markers_trajectory ...
      );

    % find right hip CoR
    right_hip_reference_file_name = ['processed' filesep makeFileName(date, subject_id, 'calibration', right_hip_calibration_file_index, 'markerTrajectories')];
    disp(['Right hip reference file name: ' right_hip_reference_file_name]);
    load(right_hip_reference_file_name);
    hip_reference = marker_trajectories;

    LASI_trajectory = extractMarkerTrajectories(hip_reference, marker_headers, 'LASI')';
    RASI_trajectory = extractMarkerTrajectories(hip_reference, marker_headers, 'RASI')';
    LPSI_trajectory = extractMarkerTrajectories(hip_reference, marker_headers, 'LPSI')';
    RPSI_trajectory = extractMarkerTrajectories(hip_reference, marker_headers, 'RPSI')';
    RTHI_trajectory = extractMarkerTrajectories(hip_reference, marker_headers, 'RTHI')';
    RTHIA_trajectory = extractMarkerTrajectories(hip_reference, marker_headers, 'RTHIA')';
    RKNE_trajectory = extractMarkerTrajectories(hip_reference, marker_headers, 'RKNE')';
    pelvis_markers_trajectory = [LASI_trajectory' RASI_trajectory' LPSI_trajectory' RPSI_trajectory'];
    right_thigh_markers_trajectory = [RTHI_trajectory' RTHIA_trajectory' RKNE_trajectory'];
%     pelvis_markers_trajectory = hip_reference(:, pelvis_markers_indices);
%     right_thigh_markers_trajectory = hip_reference(:, right_thigh_markers_indices);
    [right_hip_cor, right_hip_cor_error] = estimateJointKinematics ...
      ( ...
        pelvis_markers_reference, ...
        right_thigh_markers_reference, ...
        pelvis_markers_trajectory, ...
        right_thigh_markers_trajectory ...
      );


elseif strcmp(hip_joint_center_estimation_method, 'Tylkowski')
    
%     hjc_correction_factor_lateral = 0.11; % Tylkowski, after Bell et al., 1990
%     hjc_correction_factor_distal = 0.12; % Tylkowski, after Bell et al., 1990
%     hjc_correction_factor_frontal = 0.21; % Tylkowski, after Bell et al., 1990
%     % estimate hip CoRs
%     left_hip_cor = LASIS ...
%                     + hjc_correction_factor_lateral * inter_ASIS_distance * right_direction ...
%                     + hjc_correction_factor_distal * inter_ASIS_distance * distal_direction ...
%                     + hjc_correction_factor_frontal * inter_ASIS_distance * posterior_direction;
%     right_hip_cor = RASIS ...
%                     + hjc_correction_factor_lateral * inter_ASIS_distance * left_direction ...
%                     + hjc_correction_factor_distal * inter_ASIS_distance * distal_direction ...
%                     + hjc_correction_factor_frontal * inter_ASIS_distance * posterior_direction;    
                
% not sure where I got these values, but looking at the Bell 1990 paper again, the following seems to make more sense

    MASIS = mean([LASIS RASIS], 2);
    MPSIS = mean([LPSIS RPSIS], 2);
    z_pelvis = normVector(LASIS - RASIS); % medial-lateral, positive is left
    x_pelvis = normVector(MPSIS - MASIS); % anterior-posterior, positive is backwards
    y_pelvis = normVector(cross(z_pelvis, x_pelvis)); % proximal-distal, positive is up
    
    hjc_correction_factor_lateral = 0.14; % Tylkowski, after Bell et al., 1990
    hjc_correction_factor_distal = 0.30; % Tylkowski, after Bell et al., 1990
    hjc_correction_factor_anterior = 0.19; % Tylkowski, after Bell et al., 1990

    % estimate hip CoRs
    left_hip_cor = LASIS ...
                    - hjc_correction_factor_lateral * inter_ASIS_distance * z_pelvis ...
                    - hjc_correction_factor_distal * inter_ASIS_distance * y_pelvis ...
                    + hjc_correction_factor_anterior * inter_ASIS_distance * x_pelvis;
    right_hip_cor = RASIS ...
                    + hjc_correction_factor_lateral * inter_ASIS_distance * z_pelvis ...
                    - hjc_correction_factor_distal * inter_ASIS_distance * y_pelvis ...
                    + hjc_correction_factor_anterior * inter_ASIS_distance * x_pelvis;    

else
    error('hip joint center estimation method not recognized. Options are "SCoRE" or "Tylkowski".');
end

%% estimate knee joint centers and axes
if strcmp(knee_joint_axis_estimation_method, 'SARA')
    % find left knee CoR
    left_knee_reference_file_name = ['processed' filesep makeFileName(date, subject_id, 'calibration', left_knee_calibration_file_index, 'markerTrajectories')];
    disp(['Left knee reference file name: ' left_knee_reference_file_name]);
    load(left_knee_reference_file_name);
    knee_reference = marker_trajectories;

    LTHI_trajectory = extractMarkerTrajectories(knee_reference, marker_headers, 'LTHI')';
    LTHIA_trajectory = extractMarkerTrajectories(knee_reference, marker_headers, 'LTHIA')';
    LKNE_trajectory = extractMarkerTrajectories(knee_reference, marker_headers, 'LKNE')';
    LTIB_trajectory = extractMarkerTrajectories(knee_reference, marker_headers, 'LTIB')';
    LTIBA_trajectory = extractMarkerTrajectories(knee_reference, marker_headers, 'LTIBA')';
    LANK_trajectory = extractMarkerTrajectories(knee_reference, marker_headers, 'LANK')';
    left_thigh_markers_trajectory = [LTHI_trajectory' LTHIA_trajectory' LKNE_trajectory'];
    left_shank_markers_trajectory = [LTIB_trajectory' LTIBA_trajectory' LANK_trajectory'];
    
    
%     left_thigh_markers_trajectory = knee_reference(:, left_thigh_markers_indices);
%     left_shank_markers_trajectory = knee_reference(:, left_shank_markers_indices);
    [left_knee_point, left_knee_cor_error, left_knee_aor] = estimateJointKinematics ...
      ( ...
        left_thigh_markers_reference, ...
        left_shank_markers_reference, ...
        left_thigh_markers_trajectory, ...
        left_shank_markers_trajectory, ...
        1 ...
      );

    % find right knee CoR
    right_knee_reference_file_name = ['processed' filesep makeFileName(date, subject_id, 'calibration', right_knee_calibration_file_index, 'markerTrajectories')];
    disp(['Right knee reference file name: ' right_knee_reference_file_name]);
    load(right_knee_reference_file_name);
    knee_reference = marker_trajectories;

    RTHI_trajectory = extractMarkerTrajectories(knee_reference, marker_headers, 'RTHI')';
    RTHIA_trajectory = extractMarkerTrajectories(knee_reference, marker_headers, 'RTHIA')';
    RKNE_trajectory = extractMarkerTrajectories(knee_reference, marker_headers, 'RKNE')';
    RTIB_trajectory = extractMarkerTrajectories(knee_reference, marker_headers, 'RTIB')';
    RTIBA_trajectory = extractMarkerTrajectories(knee_reference, marker_headers, 'RTIBA')';
    RANK_trajectory = extractMarkerTrajectories(knee_reference, marker_headers, 'RANK')';
    right_thigh_markers_trajectory = [RTHI_trajectory' RTHIA_trajectory' RKNE_trajectory'];
    right_shank_markers_trajectory = [RTIB_trajectory' RTIBA_trajectory' RANK_trajectory'];
    
%     right_thigh_markers_trajectory = knee_reference(:, right_thigh_markers_indices);
%     right_shank_markers_trajectory = knee_reference(:, right_shank_markers_indices);
    [right_knee_point, right_knee_cor_error, right_knee_aor] = estimateJointKinematics ...
      ( ...
        right_thigh_markers_reference, ...
        right_shank_markers_reference, ...
        right_thigh_markers_trajectory, ...
        right_shank_markers_trajectory, ...
        1 ...
      );
elseif strcmp(knee_joint_axis_estimation_method, 'markers')
    % assume that the knee axis of rotation is the vector between the knee markers
    left_knee_aor = normVector(left_lateral_femoral_epicondyle - right_lateral_femoral_epicondyle);
    right_knee_aor = normVector(left_lateral_femoral_epicondyle - right_lateral_femoral_epicondyle);
    
else
    error('hip joint center estimation method not recognized. Options are "SARA" or "markers".');
end
% correct directions
if dot(right_knee_aor, left_direction) < 0
    right_knee_aor = -right_knee_aor;
end
if dot(left_knee_aor, left_direction) < 0
    left_knee_aor = -left_knee_aor;
end

% estimate knee CoRs
kjc_correction_factor = 0.5;
left_knee_cor = left_lateral_femoral_epicondyle - kjc_correction_factor*knee_width*left_knee_aor;
right_knee_cor = right_lateral_femoral_epicondyle + kjc_correction_factor*knee_width*right_knee_aor;

%% calculate other joint centers and axes
% define correction factors for joint centers
ajc_correction_factor = 0.5;
ejc_correction_factor = 0.5;
if strcmp(gender, 'male')
    ljc_correction_factor_frontal = 0.264;
    ljc_correction_factor_vertical = 0.126;
    cjc_correction_angle = 8; % in degrees
    cjc_correction_factor = 0.55;
    sjc_correction_angle = 11;
    sjc_correction_factor = 0.43;
elseif strcmp(gender, 'female')
    ljc_correction_factor_frontal = 0.289;
    ljc_correction_factor_vertical = 0.172;
    cjc_correction_angle = 14; % in degrees
    cjc_correction_factor = 0.53;
    sjc_correction_angle = 5;
    sjc_correction_factor = 0.53;
else
    error('Gender must be either male or female');
end

% calculate some distances
thorax_width = norm(suprasternale - c7);

% estimate ankle CoRs
left_ankle_cor = left_lateral_malleolus - ajc_correction_factor * ankle_width * left_knee_aor;
right_ankle_cor = right_lateral_malleolus + ajc_correction_factor * ankle_width * right_knee_aor;

% estimate lumbar joint center
pelvis_acs_origin = mean([LASIS RASIS], 2);
lumbar_cor = pelvis_acs_origin + ljc_correction_factor_frontal*posterior_direction*inter_ASIS_distance + ljc_correction_factor_vertical*proximal_direction*inter_ASIS_distance;

% estimate cervical joint center
c7_to_suprasternale_direction = normVector(suprasternale - c7);
c7_to_cervix_cor = expAxis(right_direction, cjc_correction_angle*pi/180) * c7_to_suprasternale_direction;
cervix_cor = c7 + cjc_correction_factor*thorax_width*c7_to_cervix_cor;

% estimate shoulder joint centers
acromion_to_shoulder_cor = expAxis(left_direction, sjc_correction_angle*pi/180) * c7_to_suprasternale_direction;
left_shoulder_cor = left_acromion + sjc_correction_factor*thorax_width*acromion_to_shoulder_cor;
right_shoulder_cor = right_acromion + sjc_correction_factor*thorax_width*acromion_to_shoulder_cor;

% estimate wrist joint centers
left_wrist_cor = mean([left_outer_wrist, left_inner_wrist], 2);
right_wrist_cor = mean([right_outer_wrist, right_inner_wrist], 2);
left_wrist_flexion_axis = normVector(left_outer_wrist - left_inner_wrist);
right_wrist_flexion_axis = normVector(right_inner_wrist - right_outer_wrist);

% estimate elbow axes and joint centers
left_elbow_cor = left_lateral_humeral_epicondyle + ejc_correction_factor*elbow_width*right_direction;
left_elbow_axis = right_direction;
left_radioulnar_axis = normVector(left_wrist_cor - left_elbow_cor);

right_elbow_cor = right_lateral_humeral_epicondyle - ejc_correction_factor*elbow_width*right_direction;
right_elbow_axis = right_direction;
right_radioulnar_axis = - normVector(right_wrist_cor - right_elbow_cor);

left_wrist_inversion_axis = cross(left_wrist_flexion_axis, left_radioulnar_axis);
right_wrist_inversion_axis = cross(right_wrist_flexion_axis, right_radioulnar_axis);

% TODO: some of these assumptions are not valid for all types of reference configurations

%% define scaling factors
if strcmp(gender, 'male')
    % mass
    head_mass_scaling_factor      = 0.067;
    torso_mass_scaling_factor     = 0.333;
    arm_mass_scaling_factor       = 0.024;
    forearm_mass_scaling_factor   = 0.017;
    hand_mass_scaling_factor      = 0.006;
    pelvis_mass_scaling_factor    = 0.142;
    thigh_mass_scaling_factor     = 0.123;
    shank_mass_scaling_factor     = 0.048;
    foot_mass_scaling_factor      = 0.012;
    
    % CoM
    head_com_scaling_factor_x    = -0.062;  head_com_scaling_factor_y    =  0.555;  head_com_scaling_factor_z    =  0.001;
    torso_com_scaling_factor_x   = -0.036;  torso_com_scaling_factor_y   = -0.420;  torso_com_scaling_factor_z   = -0.002;
    arm_com_scaling_factor_x     =  0.017;  arm_com_scaling_factor_y     = -0.452;  arm_com_scaling_factor_z     = -0.026;
    forearm_com_scaling_factor_x =  0.010;  forearm_com_scaling_factor_y = -0.417;  forearm_com_scaling_factor_z =  0.014;
    hand_com_scaling_factor_x    =  0.082;  hand_com_scaling_factor_y    = -0.839;  hand_com_scaling_factor_z    =  0.074;
    pelvis_com_scaling_factor_x  =  0.028;  pelvis_com_scaling_factor_y  = -0.280;  pelvis_com_scaling_factor_z  = -0.006;
    thigh_com_scaling_factor_x   = -0.041;  thigh_com_scaling_factor_y   = -0.429;  thigh_com_scaling_factor_z   =  0.033;
    shank_com_scaling_factor_x   = -0.048;  shank_com_scaling_factor_y   = -0.410;  shank_com_scaling_factor_z   =  0.007;
    foot_com_scaling_factor_x    =  0.382;  foot_com_scaling_factor_y    = -0.151;  foot_com_scaling_factor_z    =  0.026;
    
    % rog
    head_rxx_scaling_factor    = 0.31;  head_ryy_scaling_factor    = 0.25;  head_rzz_scaling_factor    = 0.33;  head_rxy_scaling_factor    = 0.09*1i;	head_rxz_scaling_factor    = 0.02*1i;   head_ryz_scaling_factor    = 0.03;
    torso_rxx_scaling_factor   = 0.27;  torso_ryy_scaling_factor   = 0.25;  torso_rzz_scaling_factor   = 0.28;  torso_rxy_scaling_factor   = 0.18;      torso_rxz_scaling_factor   = 0.02;      torso_ryz_scaling_factor   = 0.04*1i;
    arm_rxx_scaling_factor     = 0.31;  arm_ryy_scaling_factor     = 0.14;  arm_rzz_scaling_factor     = 0.32;  arm_rxy_scaling_factor     = 0.06;      arm_rxz_scaling_factor     = 0.05;      arm_ryz_scaling_factor     = 0.02;
    forearm_rxx_scaling_factor = 0.28;  forearm_ryy_scaling_factor = 0.11;  forearm_rzz_scaling_factor = 0.27;  forearm_rxy_scaling_factor = 0.03;      forearm_rxz_scaling_factor = 0.02;      forearm_ryz_scaling_factor = 0.08*1i;
    hand_rxx_scaling_factor    = 0.61;  hand_ryy_scaling_factor    = 0.38;  hand_rzz_scaling_factor    = 0.56;  hand_rxy_scaling_factor    = 0.22;      hand_rxz_scaling_factor    = 0.15;      hand_ryz_scaling_factor    = 0.20*1i;
    pelvis_rxx_scaling_factor  = 1.01;  pelvis_ryy_scaling_factor  = 1.06;  pelvis_rzz_scaling_factor  = 0.95;  pelvis_rxy_scaling_factor  = 0.25*1i;   pelvis_rxz_scaling_factor  = 0.12*1i;   pelvis_ryz_scaling_factor  = 0.08*1i;
    thigh_rxx_scaling_factor   = 0.29;  thigh_ryy_scaling_factor   = 0.15;  thigh_rzz_scaling_factor   = 0.30;  thigh_rxy_scaling_factor   = 0.07;      thigh_rxz_scaling_factor   = 0.02*1i;   thigh_ryz_scaling_factor   = 0.07*1i;
    shank_rxx_scaling_factor   = 0.28;  shank_ryy_scaling_factor   = 0.10;  shank_rzz_scaling_factor   = 0.28;  shank_rxy_scaling_factor   = 0.04*1i;   shank_rxz_scaling_factor   = 0.02*1i;   shank_ryz_scaling_factor   = 0.05;
    foot_rxx_scaling_factor    = 0.17;  foot_ryy_scaling_factor    = 0.37;	foot_rzz_scaling_factor    = 0.36;  foot_rxy_scaling_factor    = 0.13;      foot_rxz_scaling_factor    = 0.08*1i;	foot_ryz_scaling_factor    = 0.00;
elseif strcmp(gender, 'female')
    % mass
    head_mass_scaling_factor    = 0.067;
    torso_mass_scaling_factor   = 0.304;
    arm_mass_scaling_factor     = 0.022;
    forearm_mass_scaling_factor = 0.013;
    hand_mass_scaling_factor    = 0.005;
    pelvis_mass_scaling_factor  = 0.146;
    thigh_mass_scaling_factor   = 0.146;
    shank_mass_scaling_factor   = 0.045;
    foot_mass_scaling_factor    = 0.010;
    
    % CoM
    head_com_scaling_factor_x    = -0.070;  head_com_scaling_factor_y    =  0.597;  head_com_scaling_factor_z    =  0.000;
    torso_com_scaling_factor_x   = -0.016;  torso_com_scaling_factor_y   = -0.436;  torso_com_scaling_factor_z   = -0.006;
    arm_com_scaling_factor_x     =  0.073;  arm_com_scaling_factor_y     = -0.454;  arm_com_scaling_factor_z     = -0.028;
    forearm_com_scaling_factor_x =  0.021;  forearm_com_scaling_factor_y = -0.411;  forearm_com_scaling_factor_z =  0.019;
    hand_com_scaling_factor_x    =  0.077;  hand_com_scaling_factor_y    = -0.768;  hand_com_scaling_factor_z    =  0.048;
    pelvis_com_scaling_factor_x  = -0.009;  pelvis_com_scaling_factor_y  = -0.232;  pelvis_com_scaling_factor_z  =  0.002;
    thigh_com_scaling_factor_x   = -0.077;  thigh_com_scaling_factor_y   = -0.377;  thigh_com_scaling_factor_z   =  0.009;
    shank_com_scaling_factor_x   = -0.049;  shank_com_scaling_factor_y   = -0.404;  shank_com_scaling_factor_z   =  0.031;
    foot_com_scaling_factor_x    =  0.270;  foot_com_scaling_factor_y    = -0.218;  foot_com_scaling_factor_z    =  0.039;

    % rog
    head_rxx_scaling_factor  = 0.32;         head_ryy_scaling_factor = 0.27;         head_rzz_scaling_factor = 0.34;
    head_rxy_scaling_factor  = 0.06*1i;      head_rxz_scaling_factor = 0.01;         head_ryz_scaling_factor = 0.01*1i;
    torso_rxx_scaling_factor = 0.29;        torso_ryy_scaling_factor = 0.27;        torso_rzz_scaling_factor = 0.29;    torso_rxy_scaling_factor = 0.22;        torso_rxz_scaling_factor = 0.05;        torso_ryz_scaling_factor = 0.05*1i;
    % ACHTUNG: used the values for males here because the values for females fails to visualize as an ellipsoid.
    % Explore later!
    torso_rxx_scaling_factor = 0.27;        torso_ryy_scaling_factor = 0.25;        torso_rzz_scaling_factor = 0.28;    torso_rxy_scaling_factor = 0.18;        torso_rxz_scaling_factor = 0.02;        torso_ryz_scaling_factor = 0.04*1i;

    arm_rxx_scaling_factor     = 0.33;  arm_ryy_scaling_factor     = 0.17;  arm_rzz_scaling_factor     = 0.33;  arm_rxy_scaling_factor     = 0.03;      arm_rxz_scaling_factor     = 0.05*1i;   arm_ryz_scaling_factor     = 0.14;
    forearm_rxx_scaling_factor = 0.26;  forearm_ryy_scaling_factor = 0.14;  forearm_rzz_scaling_factor = 0.25;  forearm_rxy_scaling_factor = 0.10;      forearm_rxz_scaling_factor = 0.04;      forearm_ryz_scaling_factor = 0.14*1i;
    hand_rxx_scaling_factor    = 0.63;  hand_ryy_scaling_factor    = 0.43;  hand_rzz_scaling_factor    = 0.56;  hand_rxy_scaling_factor    = 0.29;      hand_rxz_scaling_factor    = 0.23;      hand_ryz_scaling_factor    = 0.28*1i;
    
    pelvis_rxx_scaling_factor = 0.91;       pelvis_ryy_scaling_factor = 1.00;       pelvis_rzz_scaling_factor = 0.79;   pelvis_rxy_scaling_factor = 0.34*1i;    pelvis_rxz_scaling_factor = 0.01*1i;    pelvis_ryz_scaling_factor = 0.01*1i;
    thigh_rxx_scaling_factor  = 0.31;        thigh_ryy_scaling_factor = 0.19;        thigh_rzz_scaling_factor = 0.32;    thigh_rxy_scaling_factor = 0.07;        thigh_rxz_scaling_factor = 0.02*1i;     thigh_ryz_scaling_factor = 0.07*1i;
    shank_rxx_scaling_factor  = 0.28;        shank_ryy_scaling_factor = 0.10;        shank_rzz_scaling_factor = 0.28;    shank_rxy_scaling_factor = 0.02;        shank_rxz_scaling_factor = 0.01;        shank_ryz_scaling_factor = 0.06;
    foot_rxx_scaling_factor   = 0.17;         foot_ryy_scaling_factor = 0.36;         foot_rzz_scaling_factor = 0.35;     foot_rxy_scaling_factor = 0.10*1i;      foot_rxz_scaling_factor = 0.06;         foot_ryz_scaling_factor = 0.04*1i;
end

%% set up the segment coordinate systems (SCS)

% pelvis
pelvis_scs_z = normVector(RASI_reference - LASI_reference);
pelvis_scs_y = normVector(cross(RASI_reference - MPSIS_reference, LASI_reference - MPSIS_reference));
pelvis_scs_x = cross(pelvis_scs_y, pelvis_scs_z);

% left thigh
left_thigh_scs_y = normVector(left_hip_cor - left_knee_cor);
left_thigh_scs_x = normVector(cross(left_knee_aor, left_thigh_scs_y));
left_thigh_scs_z = cross(left_thigh_scs_x, left_thigh_scs_y);

% left leg
left_leg_scs_y = normVector(left_knee_cor - left_ankle_cor);
left_leg_scs_x = normVector(cross(left_knee_aor, left_leg_scs_y));
left_leg_scs_z = cross(left_leg_scs_x, left_leg_scs_y);

% left foot
left_foot_scs_x = normVector(left_toe_mid - left_calcaneus);
left_foot_scs_y = normVector(cross(left_toe_mid-left_calcaneus, left_direction));
left_foot_scs_z = cross(left_foot_scs_x, left_foot_scs_y);

 % right thigh
right_thigh_scs_y = normVector(right_hip_cor - right_knee_cor);
right_thigh_scs_x = normVector(cross(right_knee_aor, right_thigh_scs_y));
right_thigh_scs_z = cross(right_thigh_scs_x, right_thigh_scs_y);

% right leg
right_leg_scs_y = normVector(right_knee_cor - right_ankle_cor);
right_leg_scs_x = normVector(cross(right_knee_aor, right_leg_scs_y));
right_leg_scs_z = cross(right_leg_scs_x, right_leg_scs_y);

% right foot
right_foot_scs_x = normVector(right_toe_mid - right_calcaneus);
right_foot_scs_y = normVector(cross(right_toe_mid-right_calcaneus, left_direction));
right_foot_scs_z = cross(right_foot_scs_x, right_foot_scs_y);

% torso
torso_scs_y = normVector(cervix_cor - lumbar_cor);
torso_scs_z = normVector(cross(cervix_cor - suprasternale, lumbar_cor - suprasternale));
torso_scs_x = cross(torso_scs_y, torso_scs_z);

% head and neck
head_scs_y = normVector(head_vertex - cervix_cor);
head_scs_z = normVector(cross(sellion - cervix_cor, head_vertex - cervix_cor));
head_scs_x = cross(head_scs_y, head_scs_z);

% left upper arm
left_arm_scs_y = normVector(left_shoulder_cor - left_elbow_cor);
left_arm_scs_x = normVector(cross(left_arm_scs_y, left_elbow_axis));
left_arm_scs_z = cross(left_arm_scs_y, left_arm_scs_z);

% left lower arm
left_forearm_scs_y = normVector(left_elbow_cor - left_wrist_cor);
left_forearm_scs_x = cross(left_forearm_scs_y, left_wrist_flexion_axis);
left_forearm_scs_z = normVector(cross(left_forearm_scs_x, left_forearm_scs_y));

% left hand
left_hand_scs_y = normVector(left_inner_wrist - left_hand);
left_hand_scs_x = cross(left_hand_scs_y, left_wrist_flexion_axis);
left_hand_scs_z = normVector(cross(left_hand_scs_x, left_hand_scs_y));

% right upper arm
right_arm_scs_y = normVector(right_shoulder_cor - right_elbow_cor);
right_arm_scs_x = normVector(cross(right_arm_scs_y, right_elbow_axis));
right_arm_scs_z = cross(right_arm_scs_y, right_arm_scs_z);

% right lower arm
right_forearm_scs_y = normVector(right_elbow_cor - right_wrist_cor);
right_forearm_scs_x = cross(right_forearm_scs_y, right_wrist_flexion_axis);
right_forearm_scs_z = normVector(cross(right_forearm_scs_x, right_forearm_scs_y));

% right hand
right_hand_scs_y = normVector(right_inner_wrist - right_hand);
right_hand_scs_x = cross(right_hand_scs_y, right_wrist_flexion_axis);
right_hand_scs_z = normVector(cross(right_hand_scs_x, right_hand_scs_y));

%% calculate segment mass and CoM

% calculate segment masses and correct for mean rounding errors
pelvis_segment_mass = pelvis_mass_scaling_factor    * weight * 1.002^(-1);
thigh_segment_mass = thigh_mass_scaling_factor      * weight * 1.002^(-1);
leg_segment_mass = shank_mass_scaling_factor        * weight * 1.002^(-1);
foot_segment_mass = foot_mass_scaling_factor        * weight * 1.002^(-1);
head_segment_mass = head_mass_scaling_factor        * weight * 1.002^(-1);
torso_segment_mass = torso_mass_scaling_factor      * weight * 1.002^(-1);
arm_segment_mass = arm_mass_scaling_factor          * weight * 1.002^(-1);
forearm_segment_mass = forearm_mass_scaling_factor  * weight * 1.002^(-1);
hand_segment_mass = hand_mass_scaling_factor        * weight * 1.002^(-1);

% pelvis
pelvis_segment_length = norm(mean([right_hip_cor left_hip_cor], 2) - lumbar_cor);
pelvis_com = lumbar_cor ...
                   + pelvis_com_scaling_factor_x * pelvis_segment_length * pelvis_scs_x ...
                   + pelvis_com_scaling_factor_y * pelvis_segment_length * pelvis_scs_y ...
                   + pelvis_com_scaling_factor_z * pelvis_segment_length * pelvis_scs_z;

% left thigh
left_thigh_segment_length = norm(left_hip_cor - left_knee_cor);
left_thigh_com = left_hip_cor ...
                   + thigh_com_scaling_factor_x * left_thigh_segment_length * left_thigh_scs_x ...
                   + thigh_com_scaling_factor_y * left_thigh_segment_length * left_thigh_scs_y ...
                   - thigh_com_scaling_factor_z * left_thigh_segment_length * left_thigh_scs_z;

% left leg
left_leg_segment_length = norm(left_knee_cor - left_ankle_cor);
left_leg_com = left_knee_cor ...
                   + shank_com_scaling_factor_x * left_leg_segment_length * left_leg_scs_x ...
                   + shank_com_scaling_factor_y * left_leg_segment_length * left_leg_scs_y ...
                   - shank_com_scaling_factor_z * left_leg_segment_length * left_leg_scs_z;

% left foot
left_foot_segment_length = norm(left_toe_mid - left_ankle_cor);
left_foot_com = left_ankle_cor ...
                   + foot_com_scaling_factor_x * left_foot_segment_length * left_foot_scs_x ...
                   + foot_com_scaling_factor_y * left_foot_segment_length * left_foot_scs_y ...
                   - foot_com_scaling_factor_z * left_foot_segment_length * left_foot_scs_z;

% right thigh
right_thigh_segment_length = norm(right_hip_cor - right_knee_cor);
right_thigh_com = right_hip_cor ...
                   + thigh_com_scaling_factor_x * right_thigh_segment_length * right_thigh_scs_x ...
                   + thigh_com_scaling_factor_y * right_thigh_segment_length * right_thigh_scs_y ...
                   + thigh_com_scaling_factor_z * right_thigh_segment_length * right_thigh_scs_z;

% right leg
right_leg_segment_length = norm(right_knee_cor - right_ankle_cor);
right_leg_com = right_knee_cor ...
                   + shank_com_scaling_factor_x * right_leg_segment_length * right_leg_scs_x ...
                   + shank_com_scaling_factor_y * right_leg_segment_length * right_leg_scs_y ...
                   + shank_com_scaling_factor_z * right_leg_segment_length * right_leg_scs_z;

% right foot
right_foot_segment_length = norm(right_toe_mid - right_ankle_cor);
right_foot_com = right_ankle_cor ...
                   + foot_com_scaling_factor_x * right_foot_segment_length * right_foot_scs_x ...
                   + foot_com_scaling_factor_y * right_foot_segment_length * right_foot_scs_y ...
                   + foot_com_scaling_factor_z * right_foot_segment_length * right_foot_scs_z;

% torso               
torso_segment_length = norm(lumbar_cor - cervix_cor);
torso_com = cervix_cor ...
                   + torso_com_scaling_factor_x * torso_segment_length * torso_scs_x ...
                   + torso_com_scaling_factor_y * torso_segment_length * torso_scs_y ...
                   + torso_com_scaling_factor_z * torso_segment_length * torso_scs_z;

% head and neck
head_segment_length = norm(cervix_cor - head_vertex);
head_com = cervix_cor ...
                   + head_com_scaling_factor_x * head_segment_length * head_scs_x ...
                   + head_com_scaling_factor_y * head_segment_length * head_scs_y ...
                   + head_com_scaling_factor_z * head_segment_length * head_scs_z;
               
% left upper arm
left_arm_segment_length = norm(left_shoulder_cor - left_elbow_cor);
left_arm_com = left_shoulder_cor ...
                   + arm_com_scaling_factor_x * left_arm_segment_length * left_arm_scs_x ...
                   + arm_com_scaling_factor_y * left_arm_segment_length * left_arm_scs_y ...
                   - arm_com_scaling_factor_z * left_arm_segment_length * left_arm_scs_z;

% left forearm
left_forearm_segment_length = norm(left_elbow_cor - left_wrist_cor);
left_forearm_com = left_elbow_cor ...
                   + forearm_com_scaling_factor_x * left_forearm_segment_length * left_forearm_scs_x ...
                   + forearm_com_scaling_factor_y * left_forearm_segment_length * left_forearm_scs_y ...
                   - forearm_com_scaling_factor_z * left_forearm_segment_length * left_forearm_scs_z;

% left hand
left_hand_segment_length = norm(left_inner_wrist - left_hand);
left_hand_com = left_wrist_cor ...
                   + hand_com_scaling_factor_x * left_hand_segment_length * left_hand_scs_x ...
                   + hand_com_scaling_factor_y * left_hand_segment_length * left_hand_scs_y ...
                   - hand_com_scaling_factor_z * left_hand_segment_length * left_hand_scs_z;

% right upper arm
right_arm_segment_length = norm(right_shoulder_cor - right_elbow_cor);
right_arm_com = right_shoulder_cor ...
                   + arm_com_scaling_factor_x * right_arm_segment_length * right_arm_scs_x ...
                   + arm_com_scaling_factor_y * right_arm_segment_length * right_arm_scs_y ...
                   + arm_com_scaling_factor_z * right_arm_segment_length * right_arm_scs_z;

% right forearm
right_forearm_segment_length = norm(right_elbow_cor - right_wrist_cor);
right_forearm_com = right_elbow_cor ...
                   + forearm_com_scaling_factor_x * right_forearm_segment_length * right_forearm_scs_x ...
                   + forearm_com_scaling_factor_y * right_forearm_segment_length * right_forearm_scs_y ...
                   + forearm_com_scaling_factor_z * right_forearm_segment_length * right_forearm_scs_z;

% right hand
right_hand_segment_length = norm(right_inner_wrist - right_hand);
right_hand_com = right_wrist_cor ...
                   + hand_com_scaling_factor_x * right_hand_segment_length * right_hand_scs_x ...
                   + hand_com_scaling_factor_y * right_hand_segment_length * right_hand_scs_y ...
                   + hand_com_scaling_factor_z * right_hand_segment_length * right_hand_scs_z;

%% calculate inertia tensors

% pelvis
pelvis_I_xx = (pelvis_rxx_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
pelvis_I_yy = (pelvis_ryy_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
pelvis_I_zz = (pelvis_rzz_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
pelvis_I_xy = (pelvis_rxy_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
pelvis_I_xz = (pelvis_rxz_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
pelvis_I_yz = (pelvis_ryz_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
pelvis_inertia_tensor = [pelvis_I_xx pelvis_I_xy pelvis_I_xz; pelvis_I_xy pelvis_I_yy pelvis_I_yz; pelvis_I_xz pelvis_I_yz pelvis_I_zz];

% right thigh
right_thigh_I_xx = (thigh_rxx_scaling_factor*right_thigh_segment_length)^2 * thigh_segment_mass;
right_thigh_I_yy = (thigh_ryy_scaling_factor*right_thigh_segment_length)^2 * thigh_segment_mass;
right_thigh_I_zz = (thigh_rzz_scaling_factor*right_thigh_segment_length)^2 * thigh_segment_mass;
right_thigh_I_xy = (thigh_rxy_scaling_factor*right_thigh_segment_length)^2 * thigh_segment_mass;
right_thigh_I_xz = (thigh_rxz_scaling_factor*right_thigh_segment_length)^2 * thigh_segment_mass;
right_thigh_I_yz = (thigh_ryz_scaling_factor*right_thigh_segment_length)^2 * thigh_segment_mass;
right_thigh_inertia_tensor = [right_thigh_I_xx right_thigh_I_xy right_thigh_I_xz; right_thigh_I_xy right_thigh_I_yy right_thigh_I_yz; right_thigh_I_xz right_thigh_I_yz right_thigh_I_zz];

% right leg
right_leg_I_xx = (shank_rxx_scaling_factor*right_leg_segment_length)^2 * leg_segment_mass;
right_leg_I_yy = (shank_ryy_scaling_factor*right_leg_segment_length)^2 * leg_segment_mass;
right_leg_I_zz = (shank_rzz_scaling_factor*right_leg_segment_length)^2 * leg_segment_mass;
right_leg_I_xy = (shank_rxy_scaling_factor*right_leg_segment_length)^2 * leg_segment_mass;
right_leg_I_xz = (shank_rxz_scaling_factor*right_leg_segment_length)^2 * leg_segment_mass;
right_leg_I_yz = (shank_ryz_scaling_factor*right_leg_segment_length)^2 * leg_segment_mass;
right_leg_inertia_tensor = [right_leg_I_xx right_leg_I_xy right_leg_I_xz; right_leg_I_xy right_leg_I_yy right_leg_I_yz; right_leg_I_xz right_leg_I_yz right_leg_I_zz];

% right foot
right_foot_I_xx = (foot_rxx_scaling_factor*right_foot_segment_length)^2 * foot_segment_mass;
right_foot_I_yy = (foot_ryy_scaling_factor*right_foot_segment_length)^2 * foot_segment_mass;
right_foot_I_zz = (foot_rzz_scaling_factor*right_foot_segment_length)^2 * foot_segment_mass;
right_foot_I_xy = (foot_rxy_scaling_factor*right_foot_segment_length)^2 * foot_segment_mass;
right_foot_I_xz = (foot_rxz_scaling_factor*right_foot_segment_length)^2 * foot_segment_mass;
right_foot_I_yz = (foot_ryz_scaling_factor*right_foot_segment_length)^2 * foot_segment_mass;
right_foot_inertia_tensor = [right_foot_I_xx right_foot_I_xy right_foot_I_xz; right_foot_I_xy right_foot_I_yy right_foot_I_yz; right_foot_I_xz right_foot_I_yz right_foot_I_zz];

% left thigh
left_thigh_I_xx = (thigh_rxx_scaling_factor*left_thigh_segment_length)^2 * thigh_segment_mass;
left_thigh_I_yy = (thigh_ryy_scaling_factor*left_thigh_segment_length)^2 * thigh_segment_mass;
left_thigh_I_zz = (thigh_rzz_scaling_factor*left_thigh_segment_length)^2 * thigh_segment_mass;
left_thigh_I_xy = (thigh_rxy_scaling_factor*left_thigh_segment_length)^2 * thigh_segment_mass;
left_thigh_I_xz = -(thigh_rxz_scaling_factor*left_thigh_segment_length)^2 * thigh_segment_mass;
left_thigh_I_yz = -(thigh_ryz_scaling_factor*left_thigh_segment_length)^2 * thigh_segment_mass;
left_thigh_inertia_tensor = [left_thigh_I_xx left_thigh_I_xy left_thigh_I_xz; left_thigh_I_xy left_thigh_I_yy left_thigh_I_yz; left_thigh_I_xz left_thigh_I_yz left_thigh_I_zz];

% left leg
left_leg_I_xx = (shank_rxx_scaling_factor*left_leg_segment_length)^2 * leg_segment_mass;
left_leg_I_yy = (shank_ryy_scaling_factor*left_leg_segment_length)^2 * leg_segment_mass;
left_leg_I_zz = (shank_rzz_scaling_factor*left_leg_segment_length)^2 * leg_segment_mass;
left_leg_I_xy = (shank_rxy_scaling_factor*left_leg_segment_length)^2 * leg_segment_mass;
left_leg_I_xz = -(shank_rxz_scaling_factor*left_leg_segment_length)^2 * leg_segment_mass;
left_leg_I_yz = -(shank_ryz_scaling_factor*left_leg_segment_length)^2 * leg_segment_mass;
left_leg_inertia_tensor = [left_leg_I_xx left_leg_I_xy left_leg_I_xz; left_leg_I_xy left_leg_I_yy left_leg_I_yz; left_leg_I_xz left_leg_I_yz left_leg_I_zz];

% left foot
left_foot_I_xx = (foot_rxx_scaling_factor*left_foot_segment_length)^2 * foot_segment_mass;
left_foot_I_yy = (foot_ryy_scaling_factor*left_foot_segment_length)^2 * foot_segment_mass;
left_foot_I_zz = (foot_rzz_scaling_factor*left_foot_segment_length)^2 * foot_segment_mass;
left_foot_I_xy = (foot_rxy_scaling_factor*left_foot_segment_length)^2 * foot_segment_mass;
left_foot_I_xz = -(foot_rxz_scaling_factor*left_foot_segment_length)^2 * foot_segment_mass;
left_foot_I_yz = -(foot_ryz_scaling_factor*left_foot_segment_length)^2 * foot_segment_mass;
left_foot_inertia_tensor = [left_foot_I_xx left_foot_I_xy left_foot_I_xz; left_foot_I_xy left_foot_I_yy left_foot_I_yz; left_foot_I_xz left_foot_I_yz left_foot_I_zz];

% torso
torso_I_xx = (torso_rxx_scaling_factor*torso_segment_length)^2 * torso_segment_mass;
torso_I_yy = (torso_ryy_scaling_factor*torso_segment_length)^2 * torso_segment_mass;
torso_I_zz = (torso_rzz_scaling_factor*torso_segment_length)^2 * torso_segment_mass;
torso_I_xy = (torso_rxy_scaling_factor*torso_segment_length)^2 * torso_segment_mass;
torso_I_xz = (torso_rxz_scaling_factor*torso_segment_length)^2 * torso_segment_mass;
torso_I_yz = (torso_ryz_scaling_factor*torso_segment_length)^2 * torso_segment_mass;
torso_inertia_tensor = [torso_I_xx torso_I_xy torso_I_xz; torso_I_xy torso_I_yy torso_I_yz; torso_I_xz torso_I_yz torso_I_zz];

% head and neck
head_I_xx = (head_rxx_scaling_factor*head_segment_length)^2 * head_segment_mass;
head_I_yy = (head_ryy_scaling_factor*head_segment_length)^2 * head_segment_mass;
head_I_zz = (head_rzz_scaling_factor*head_segment_length)^2 * head_segment_mass;
head_I_xy = (head_rxy_scaling_factor*head_segment_length)^2 * head_segment_mass;
head_I_xz = (head_rxz_scaling_factor*head_segment_length)^2 * head_segment_mass;
head_I_yz = (head_ryz_scaling_factor*head_segment_length)^2 * head_segment_mass;
head_inertia_tensor = [head_I_xx head_I_xy head_I_xz; head_I_xy head_I_yy head_I_yz; head_I_xz head_I_yz head_I_zz];

% right arm
right_arm_I_xx = (arm_rxx_scaling_factor*right_arm_segment_length)^2 * arm_segment_mass;
right_arm_I_yy = (arm_ryy_scaling_factor*right_arm_segment_length)^2 * arm_segment_mass;
right_arm_I_zz = (arm_rzz_scaling_factor*right_arm_segment_length)^2 * arm_segment_mass;
right_arm_I_xy = (arm_rxy_scaling_factor*right_arm_segment_length)^2 * arm_segment_mass;
right_arm_I_xz = (arm_rxz_scaling_factor*right_arm_segment_length)^2 * arm_segment_mass;
right_arm_I_yz = (arm_ryz_scaling_factor*right_arm_segment_length)^2 * arm_segment_mass;
right_arm_inertia_tensor = [right_arm_I_xx right_arm_I_xy right_arm_I_xz; right_arm_I_xy right_arm_I_yy right_arm_I_yz; right_arm_I_xz right_arm_I_yz right_arm_I_zz];

% right forearm
right_forearm_I_xx = (forearm_rxx_scaling_factor*right_forearm_segment_length)^2 * forearm_segment_mass;
right_forearm_I_yy = (forearm_ryy_scaling_factor*right_forearm_segment_length)^2 * forearm_segment_mass;
right_forearm_I_zz = (forearm_rzz_scaling_factor*right_forearm_segment_length)^2 * forearm_segment_mass;
right_forearm_I_xy = (forearm_rxy_scaling_factor*right_forearm_segment_length)^2 * forearm_segment_mass;
right_forearm_I_xz = (forearm_rxz_scaling_factor*right_forearm_segment_length)^2 * forearm_segment_mass;
right_forearm_I_yz = (forearm_ryz_scaling_factor*right_forearm_segment_length)^2 * forearm_segment_mass;
right_forearm_inertia_tensor = [right_forearm_I_xx right_forearm_I_xy right_forearm_I_xz; right_forearm_I_xy right_forearm_I_yy right_forearm_I_yz; right_forearm_I_xz right_forearm_I_yz right_forearm_I_zz];

% right hand
right_hand_I_xx = (hand_rxx_scaling_factor*right_hand_segment_length)^2 * hand_segment_mass;
right_hand_I_yy = (hand_ryy_scaling_factor*right_hand_segment_length)^2 * hand_segment_mass;
right_hand_I_zz = (hand_rzz_scaling_factor*right_hand_segment_length)^2 * hand_segment_mass;
right_hand_I_xy = (hand_rxy_scaling_factor*right_hand_segment_length)^2 * hand_segment_mass;
right_hand_I_xz = (hand_rxz_scaling_factor*right_hand_segment_length)^2 * hand_segment_mass;
right_hand_I_yz = (hand_ryz_scaling_factor*right_hand_segment_length)^2 * hand_segment_mass;
right_hand_inertia_tensor = [right_hand_I_xx right_hand_I_xy right_hand_I_xz; right_hand_I_xy right_hand_I_yy right_hand_I_yz; right_hand_I_xz right_hand_I_yz right_hand_I_zz];

% left arm
left_arm_I_xx = (arm_rxx_scaling_factor*left_arm_segment_length)^2 * arm_segment_mass;
left_arm_I_yy = (arm_ryy_scaling_factor*left_arm_segment_length)^2 * arm_segment_mass;
left_arm_I_zz = (arm_rzz_scaling_factor*left_arm_segment_length)^2 * arm_segment_mass;
left_arm_I_xy = (arm_rxy_scaling_factor*left_arm_segment_length)^2 * arm_segment_mass;
left_arm_I_xz = -(arm_rxz_scaling_factor*left_arm_segment_length)^2 * arm_segment_mass;
left_arm_I_yz = -(arm_ryz_scaling_factor*left_arm_segment_length)^2 * arm_segment_mass;
left_arm_inertia_tensor = [left_arm_I_xx left_arm_I_xy left_arm_I_xz; left_arm_I_xy left_arm_I_yy left_arm_I_yz; left_arm_I_xz left_arm_I_yz left_arm_I_zz];

% left forearm
left_forearm_I_xx = (forearm_rxx_scaling_factor*left_forearm_segment_length)^2 * forearm_segment_mass;
left_forearm_I_yy = (forearm_ryy_scaling_factor*left_forearm_segment_length)^2 * forearm_segment_mass;
left_forearm_I_zz = (forearm_rzz_scaling_factor*left_forearm_segment_length)^2 * forearm_segment_mass;
left_forearm_I_xy = (forearm_rxy_scaling_factor*left_forearm_segment_length)^2 * forearm_segment_mass;
left_forearm_I_xz = -(forearm_rxz_scaling_factor*left_forearm_segment_length)^2 * forearm_segment_mass;
left_forearm_I_yz = -(forearm_ryz_scaling_factor*left_forearm_segment_length)^2 * forearm_segment_mass;
left_forearm_inertia_tensor = [left_forearm_I_xx left_forearm_I_xy left_forearm_I_xz; left_forearm_I_xy left_forearm_I_yy left_forearm_I_yz; left_forearm_I_xz left_forearm_I_yz left_forearm_I_zz];

% left hand
left_hand_I_xx = (hand_rxx_scaling_factor*left_hand_segment_length)^2 * hand_segment_mass;
left_hand_I_yy = (hand_ryy_scaling_factor*left_hand_segment_length)^2 * hand_segment_mass;
left_hand_I_zz = (hand_rzz_scaling_factor*left_hand_segment_length)^2 * hand_segment_mass;
left_hand_I_xy = (hand_rxy_scaling_factor*left_hand_segment_length)^2 * hand_segment_mass;
left_hand_I_xz = -(hand_rxz_scaling_factor*left_hand_segment_length)^2 * hand_segment_mass;
left_hand_I_yz = -(hand_ryz_scaling_factor*left_hand_segment_length)^2 * hand_segment_mass;
left_hand_inertia_tensor = [left_hand_I_xx left_hand_I_xy left_hand_I_xz; left_hand_I_xy left_hand_I_yy left_hand_I_yz; left_hand_I_xz left_hand_I_yz left_hand_I_zz];

%% assemble kinematic tree
if create_kinematic_tree
    % adjust end-effector positions for foot extension in z-direction
    right_heel = right_ankle_cor; right_heel(3) = right_toe_mid(3);
    left_heel = left_ankle_cor; left_heel(3) = left_toe_mid(3);


    joint_positions = ...
    { ...
        pelvis_com, pelvis_com, pelvis_com, pelvis_com, pelvis_com, pelvis_com, ...
        left_hip_cor, left_hip_cor, left_hip_cor, left_knee_cor, left_ankle_cor, left_ankle_cor, ...
        right_hip_cor, right_hip_cor, right_hip_cor, right_knee_cor, right_ankle_cor, right_ankle_cor, ...
        lumbar_cor, lumbar_cor, lumbar_cor, cervix_cor, cervix_cor, cervix_cor, ...
        left_shoulder_cor, left_shoulder_cor, left_shoulder_cor, left_elbow_cor, left_elbow_cor, left_wrist_cor, left_wrist_cor ...
        right_shoulder_cor, right_shoulder_cor, right_shoulder_cor, right_elbow_cor, right_elbow_cor, right_wrist_cor, right_wrist_cor ...
    };

    joint_axes = ...
    { ...
        e_1, e_2, e_3, e_1, e_2, e_3, ...       % pelvis free body DoFs
        right_direction, anterior_direction, proximal_direction, left_knee_aor, left_direction, posterior_direction, ...   % left leg
        right_direction, posterior_direction, distal_direction, right_knee_aor, left_direction, anterior_direction,...   % right leg
        left_direction, posterior_direction, proximal_direction, left_direction, posterior_direction, proximal_direction, ...       % trunk and neck
        anterior_direction, right_direction, distal_direction, left_elbow_axis, left_radioulnar_axis, left_wrist_flexion_axis, left_wrist_inversion_axis, ... % left arm
        posterior_direction, right_direction, proximal_direction, right_elbow_axis, right_radioulnar_axis, right_wrist_flexion_axis, right_wrist_inversion_axis, ... % right arm
    };


    joint_types = ...
      [ ...
        2 2 2 1 1 1 ... % virtual dofs
        1 1 1 1 1 1 ... % left leg
        1 1 1 1 1 1 ... % right leg
        1 1 1 1 1 1 ... % l5 and neck
        1 1 1 1 1 1 1 ... % left arm
        1 1 1 1 1 1 1 ... % right arm
      ]; % 1 = rotational, 2 = prismatic

    % link setup
    link_com_positions =  ...
    { ...
        pelvis_com, pelvis_com, pelvis_com, pelvis_com, pelvis_com, pelvis_com, ...
        left_thigh_com, left_thigh_com, left_thigh_com, left_leg_com, left_foot_com, left_foot_com, ...
        right_thigh_com, right_thigh_com, right_thigh_com, right_leg_com, right_foot_com, right_foot_com, ...
        torso_com, torso_com, torso_com, head_com, head_com, head_com ...
        left_arm_com, left_arm_com, left_arm_com, left_forearm_com, left_forearm_com, left_hand_com, left_hand_com, ...
        right_arm_com, right_arm_com, right_arm_com, right_forearm_com, right_forearm_com, right_hand_com, right_hand_com ...
    };

    link_orientations = ...
      {
        eye(3); eye(3); eye(3); eye(3); eye(3); [pelvis_scs_x pelvis_scs_y pelvis_scs_z]; ...
        eye(3); eye(3); [left_thigh_scs_x left_thigh_scs_y left_thigh_scs_z]; [left_leg_scs_x left_leg_scs_y left_leg_scs_z]; eye(3); [left_foot_scs_x left_foot_scs_y left_foot_scs_z]; ...
        eye(3); eye(3); [right_thigh_scs_x right_thigh_scs_y right_thigh_scs_z]; [right_leg_scs_x right_leg_scs_y right_leg_scs_z]; eye(3); [right_foot_scs_x right_foot_scs_y right_foot_scs_z]; ...
        eye(3); eye(3); [torso_scs_x torso_scs_y torso_scs_z]; eye(3); eye(3); [head_scs_x head_scs_y head_scs_z]; ...
        eye(3); eye(3); [left_arm_scs_x left_arm_scs_y left_arm_scs_z]; eye(3); [left_forearm_scs_x left_forearm_scs_y left_forearm_scs_z]; eye(3); [left_hand_scs_x left_hand_scs_y left_hand_scs_z]; ...
        eye(3); eye(3); [right_arm_scs_x right_arm_scs_y right_arm_scs_z]; eye(3); [right_forearm_scs_x right_forearm_scs_y right_forearm_scs_z]; eye(3); [right_hand_scs_x right_hand_scs_y right_hand_scs_z]; ...
      };

    generalized_inertia_matrix_pelvis           = [pelvis_segment_mass*eye(3) zeros(3); zeros(3) pelvis_inertia_tensor];
    generalized_inertia_matrix_left_thigh       = [thigh_segment_mass*eye(3) zeros(3); zeros(3) left_thigh_inertia_tensor];
    generalized_inertia_matrix_left_shank       = [leg_segment_mass*eye(3) zeros(3); zeros(3) left_leg_inertia_tensor];
    generalized_inertia_matrix_left_foot        = [foot_segment_mass*eye(3) zeros(3); zeros(3) left_foot_inertia_tensor];
    generalized_inertia_matrix_right_thigh      = [thigh_segment_mass*eye(3) zeros(3); zeros(3) right_thigh_inertia_tensor];
    generalized_inertia_matrix_right_shank      = [leg_segment_mass*eye(3) zeros(3); zeros(3) right_leg_inertia_tensor];
    generalized_inertia_matrix_right_foot       = [foot_segment_mass*eye(3) zeros(3); zeros(3) right_foot_inertia_tensor];
    generalized_inertia_matrix_torso            = [torso_segment_mass*eye(3) zeros(3); zeros(3) torso_inertia_tensor];
    generalized_inertia_matrix_head             = [head_segment_mass*eye(3) zeros(3); zeros(3) head_inertia_tensor];
    generalized_inertia_matrix_left_arm         = [arm_segment_mass*eye(3) zeros(3); zeros(3) left_arm_inertia_tensor];
    generalized_inertia_matrix_left_forearm     = [forearm_segment_mass*eye(3) zeros(3); zeros(3) left_forearm_inertia_tensor];
    generalized_inertia_matrix_left_hand        = [hand_segment_mass*eye(3) zeros(3); zeros(3) left_hand_inertia_tensor];
    generalized_inertia_matrix_right_arm        = [arm_segment_mass*eye(3) zeros(3); zeros(3) right_arm_inertia_tensor];
    generalized_inertia_matrix_right_forearm    = [forearm_segment_mass*eye(3) zeros(3); zeros(3) right_forearm_inertia_tensor];
    generalized_inertia_matrix_right_hand       = [hand_segment_mass*eye(3) zeros(3); zeros(3) right_hand_inertia_tensor];

    generalized_link_inertia_matrices = ...
      { ...
        zeros(6, 6), zeros(6, 6), zeros(6, 6), zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_pelvis, ...
        zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_left_thigh, generalized_inertia_matrix_left_shank, zeros(6, 6), generalized_inertia_matrix_left_foot, ...
        zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_right_thigh, generalized_inertia_matrix_right_shank, zeros(6, 6), generalized_inertia_matrix_right_foot, ...
        zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_torso, zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_head ...
        zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_left_arm, zeros(6, 6), generalized_inertia_matrix_left_forearm, zeros(6, 6), generalized_inertia_matrix_left_hand, ...
        zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_right_arm, zeros(6, 6), generalized_inertia_matrix_right_forearm, zeros(6, 6), generalized_inertia_matrix_right_hand, ...
      };

    end_effector_transformations = ...
      { ...
        [left_foot_scs_x left_foot_scs_y left_foot_scs_z left_heel; 0 0 0 1], ...
        [left_foot_scs_x left_foot_scs_y left_foot_scs_z left_toe_mid; 0 0 0 1], ...
        [left_foot_scs_x left_foot_scs_y left_foot_scs_z left_ankle_cor; 0 0 0 1], ...
        [right_foot_scs_x right_foot_scs_y right_foot_scs_z right_heel; 0 0 0 1], ...
        [right_foot_scs_x right_foot_scs_y right_foot_scs_z right_toe_mid; 0 0 0 1], ...
        [right_foot_scs_x right_foot_scs_y right_foot_scs_z right_ankle_cor; 0 0 0 1], ...
        [eye(3), head_vertex; 0 0 0 1], ...
        [eye(3), left_wrist_cor + [0; 0.05; 0]; 0 0 0 1], ...
        [eye(3), right_wrist_cor + [0; 0.05; 0]; 0 0 0 1], ...
      };

    branch_matrix = ...
      [ ...
        1 1 1 1 1 1    1 1 1 1 1 1   0 0 0 0 0 0   0 0 0   0 0 0   0 0 0 0 0 0 0   0 0 0 0 0 0 0; ... % left leg until heel
        1 1 1 1 1 1    1 1 1 1 1 1   0 0 0 0 0 0   0 0 0   0 0 0   0 0 0 0 0 0 0   0 0 0 0 0 0 0; ... % left leg until toes
        1 1 1 1 1 1    1 1 1 1 1 1   0 0 0 0 0 0   0 0 0   0 0 0   0 0 0 0 0 0 0   0 0 0 0 0 0 0; ... % left leg until ankle
        1 1 1 1 1 1    0 0 0 0 0 0   1 1 1 1 1 1   0 0 0   0 0 0   0 0 0 0 0 0 0   0 0 0 0 0 0 0; ... % right leg until heel
        1 1 1 1 1 1    0 0 0 0 0 0   1 1 1 1 1 1   0 0 0   0 0 0   0 0 0 0 0 0 0   0 0 0 0 0 0 0; ... % right leg until toes
        1 1 1 1 1 1    0 0 0 0 0 0   1 1 1 1 1 1   0 0 0   0 0 0   0 0 0 0 0 0 0   0 0 0 0 0 0 0; ... % right leg until ankle
        1 1 1 1 1 1    0 0 0 0 0 0   0 0 0 0 0 0   1 1 1   1 1 1   0 0 0 0 0 0 0   0 0 0 0 0 0 0; ... % head
        1 1 1 1 1 1    0 0 0 0 0 0   0 0 0 0 0 0   1 1 1   0 0 0   1 1 1 1 1 1 1   0 0 0 0 0 0 0; ... % left arm
        1 1 1 1 1 1    0 0 0 0 0 0   0 0 0 0 0 0   1 1 1   0 0 0   0 0 0 0 0 0 0   1 1 1 1 1 1 1; ... % right arm
      ]; % each row is a branch, listing the joints that move the end-effector of that branch

    kinematic_tree = GeneralKinematicTree ...
    ( ...
      joint_positions, ...
      joint_axes, ...
      joint_types, ...
      branch_matrix, ...
      end_effector_transformations, ...
      link_com_positions, ...
      link_orientations, ...
      generalized_link_inertia_matrices ...
    );

    % adjust the main axes of the toe frames
    kinematic_tree.referenceJointTransformations{12}(1:3, 1:3) = [left_foot_scs_x left_foot_scs_y left_foot_scs_z];
    kinematic_tree.referenceJointTransformations{18}(1:3, 1:3) = [right_foot_scs_x right_foot_scs_y right_foot_scs_z];

    % update
    kinematic_tree.updateInternals();

    kinematic_tree.jointLabels = ...
      { ...
        'pelvis, x-translation', 'pelvis, y-translation', 'pelvis, z-translation', 'pelvis, x-rotation', 'pelvis, y-rotation', 'pelvis, z-rotation', ...
        'left hip flexion/extension', 'left hip ab/adduction', 'left hip external/internal rotation', 'left knee flexion/extension', 'left ankle plantar/dorsiflexion', 'left ankle inversion/eversion', ...
        'right hip flexion/extension', 'right hip ab/adduction', 'right hip external/internal rotation', 'right knee flexion/extension', 'right ankle plantar/dorsiflexion', 'right ankle inversion/eversion', ...
        'lumbar joint - forward/backward bending', 'lumbar joint - sideways bending (left/right)', 'lumbar joint - internal rotation (left/right)', ...
        'cervical joint - forward/backward bending', 'cervical joint - sideways bending (left/right)', 'cervical joint - internal rotation (left/right)', ...
        'left shoulder ab/adduction', 'left shoulder flexion/extension', 'left shoulder in/external rotation', 'left elbow flexion/extension', 'left pronation/supination', 'left wrist flexion/extension', 'left wrist radial/ulnar deviation', ...
        'right shoulder ab/adduction', 'right shoulder flexion/extension', 'right shoulder in/external rotation', 'right elbow flexion/extension', 'right pronation/supination', 'right wrist flexion/extension', 'right radial/ulnar deviation', ...
      };
    kinematic_tree.markerLabels = marker_headers;

    % define markers
    red = [1 0 0];
    green = [0 1 0];
    blue = [0 0 1];
    yellow = [1 0.9 0];

    markerSegments = zeros(1, length(marker_headers));
    markerSegments(strcmp(marker_headers, 'LFHD')) = 24;
    markerSegments(strcmp(marker_headers, 'RFHD')) = 24;
    markerSegments(strcmp(marker_headers, 'LBHD')) = 24;
    markerSegments(strcmp(marker_headers, 'RBHD')) = 24;

    markerSegments(strcmp(marker_headers, 'C7')) = 21;
    markerSegments(strcmp(marker_headers, 'T10')) = 21;
    markerSegments(strcmp(marker_headers, 'CLAV')) = 21;
    markerSegments(strcmp(marker_headers, 'STRN')) = 21;
    markerSegments(strcmp(marker_headers, 'RBAK')) = 21;

    markerSegments(strcmp(marker_headers, 'LSHO')) = 21;
    markerSegments(strcmp(marker_headers, 'LUPA')) = 27;
    markerSegments(strcmp(marker_headers, 'LELB')) = 27;
    markerSegments(strcmp(marker_headers, 'LFRM')) = 29;
    markerSegments(strcmp(marker_headers, 'LWRA')) = 29;
    markerSegments(strcmp(marker_headers, 'LWRB')) = 29;
    markerSegments(strcmp(marker_headers, 'LFIN')) = 31;

    markerSegments(strcmp(marker_headers, 'RSHO')) = 21;
    markerSegments(strcmp(marker_headers, 'RUPA')) = 34;
    markerSegments(strcmp(marker_headers, 'RELB')) = 34;
    markerSegments(strcmp(marker_headers, 'RFRM')) = 36;
    markerSegments(strcmp(marker_headers, 'RWRA')) = 36;
    markerSegments(strcmp(marker_headers, 'RWRB')) = 36;
    markerSegments(strcmp(marker_headers, 'RFIN')) = 38;

    markerSegments(strcmp(marker_headers, 'LASI')) = 6;
    markerSegments(strcmp(marker_headers, 'RASI')) = 6;
    markerSegments(strcmp(marker_headers, 'LPSI')) = 6;
    markerSegments(strcmp(marker_headers, 'RPSI')) = 6;

    markerSegments(strcmp(marker_headers, 'LTHI')) = 9;
    markerSegments(strcmp(marker_headers, 'LTHIA')) = 9;
    markerSegments(strcmp(marker_headers, 'LKNE')) = 9;
    markerSegments(strcmp(marker_headers, 'LTIB')) = 10;
    markerSegments(strcmp(marker_headers, 'LTIBA')) = 10;
    markerSegments(strcmp(marker_headers, 'LANK')) = 10;
    markerSegments(strcmp(marker_headers, 'LHEE')) = 12;
    markerSegments(strcmp(marker_headers, 'LTOE')) = 12;
    markerSegments(strcmp(marker_headers, 'LTOEL')) = 12;

    markerSegments(strcmp(marker_headers, 'RTHI')) = 15;
    markerSegments(strcmp(marker_headers, 'RTHIA')) = 15;
    markerSegments(strcmp(marker_headers, 'RKNE')) = 15;
    markerSegments(strcmp(marker_headers, 'RTIB')) = 16;
    markerSegments(strcmp(marker_headers, 'RTIBA')) = 16;
    markerSegments(strcmp(marker_headers, 'RANK')) = 16;
    markerSegments(strcmp(marker_headers, 'RHEE')) = 18;
    markerSegments(strcmp(marker_headers, 'RTOE')) = 18;
    markerSegments(strcmp(marker_headers, 'RTOEL')) = 18;

    marker_color_list = cell(1, length(marker_headers));
    marker_color_list{strcmp(marker_headers, 'LFHD')} = red;
    marker_color_list{strcmp(marker_headers, 'RFHD')} = green;
    marker_color_list{strcmp(marker_headers, 'LBHD')} = red;
    marker_color_list{strcmp(marker_headers, 'RBHD')} = green;

    marker_color_list{strcmp(marker_headers, 'C7')} = blue;
    marker_color_list{strcmp(marker_headers, 'T10')} = blue;
    marker_color_list{strcmp(marker_headers, 'CLAV')} = blue;
    marker_color_list{strcmp(marker_headers, 'STRN')} = blue;
    marker_color_list{strcmp(marker_headers, 'RBAK')} = blue;

    marker_color_list{strcmp(marker_headers, 'LSHO')} = red;
    marker_color_list{strcmp(marker_headers, 'LUPA')} = red;
    marker_color_list{strcmp(marker_headers, 'LELB')} = red;
    marker_color_list{strcmp(marker_headers, 'LFRA')} = red;
%     marker_color_list{strcmp(marker_headers, 'LFRM')} = red;
    marker_color_list{strcmp(marker_headers, 'LWRA')} = red;
    marker_color_list{strcmp(marker_headers, 'LWRB')} = red;
    marker_color_list{strcmp(marker_headers, 'LFIN')} = red;

    marker_color_list{strcmp(marker_headers, 'RSHO')} = green;
    marker_color_list{strcmp(marker_headers, 'RUPA')} = green;
    marker_color_list{strcmp(marker_headers, 'RELB')} = green;
%     marker_color_list{strcmp(marker_headers, 'RFRM')} = green;
    marker_color_list{strcmp(marker_headers, 'RFRA')} = green;
    marker_color_list{strcmp(marker_headers, 'RWRA')} = green;
    marker_color_list{strcmp(marker_headers, 'RWRB')} = green;
    marker_color_list{strcmp(marker_headers, 'RFIN')} = green;

    marker_color_list{strcmp(marker_headers, 'LASI')} = red;
    marker_color_list{strcmp(marker_headers, 'RASI')} = green;
    marker_color_list{strcmp(marker_headers, 'LPSI')} = red;
    marker_color_list{strcmp(marker_headers, 'RPSI')} = green;

    marker_color_list{strcmp(marker_headers, 'LTHI')} = red;
%     marker_color_list{strcmp(marker_headers, 'LTHIA')} = red;
    marker_color_list{strcmp(marker_headers, 'LKNE')} = red;
    marker_color_list{strcmp(marker_headers, 'LTIB')} = red;
%     marker_color_list{strcmp(marker_headers, 'LTIBA')} = red;
    marker_color_list{strcmp(marker_headers, 'LANK')} = red;
    marker_color_list{strcmp(marker_headers, 'LHEE')} = red;
    marker_color_list{strcmp(marker_headers, 'LTOE')} = red;
%     marker_color_list{strcmp(marker_headers, 'LTOEL')} = red;

    marker_color_list{strcmp(marker_headers, 'RTHI')} = green;
%     marker_color_list{strcmp(marker_headers, 'RTHIA')} = green;
    marker_color_list{strcmp(marker_headers, 'RKNE')} = green;
    marker_color_list{strcmp(marker_headers, 'RTIB')} = green;
%     marker_color_list{strcmp(marker_headers, 'RTIBA')} = green;
    marker_color_list{strcmp(marker_headers, 'RANK')} = green;
    marker_color_list{strcmp(marker_headers, 'RHEE')} = green;
    marker_color_list{strcmp(marker_headers, 'RTOE')} = green;
%     marker_color_list{strcmp(marker_headers, 'RTOEL')} = green;

    number_of_markers = length(markerSegments);
    for i_marker = 1 : number_of_markers
        kinematic_tree.addMarker(markerSegments(i_marker), marker_reference((i_marker-1)*3+1 : (i_marker-1)*3+3)', marker_color_list{i_marker});
    end

    % kinematic_tree.addMarkerConnectionLineBody(1:3);
    % kinematic_tree.addMarkerConnectionLineBody(4:7);
    % kinematic_tree.addMarkerConnectionLineBody(8:10);
    % kinematic_tree.addMarkerConnectionLineBody(11:13);
    % kinematic_tree.addMarkerConnectionLineBody(14:16);
    % kinematic_tree.addMarkerConnectionLineBody(17:19);
    % kinematic_tree.addMarkerConnectionLineBody(21:23);
    % kinematic_tree.addMarkerConnectionLineBody(24:26);
    % kinematic_tree.addMarkerConnectionLineBody(27:29);
    % kinematic_tree.addMarkerConnectionLineBody(29:30);

    kinematic_tree.updateInternals();
end

%% set up marker coordinate systems (MCS)

% define joint center references
joint_center_reference = ...
  [ ...
    cervix_cor', ...
    left_shoulder_cor', ...
    right_shoulder_cor', ...
    left_elbow_cor', ...
    right_elbow_cor', ...
    left_wrist_cor', ...
    right_wrist_cor', ...
    lumbar_cor', ...
    left_hip_cor', ...
    right_hip_cor', ...
    left_knee_cor', ...
    right_knee_cor', ...
    left_ankle_cor', ...
    right_ankle_cor' ...
  ];

joint_center_headers = ...
  { ...
    'CERVIXCOR', ...
    'LSHOULDERCOR', ...
    'RSHOULDERCOR', ...
    'LELBOWCOR', ...
    'RELBOWCOR', ...
    'LWRISTCOR', ...
    'RWRISTCOR', ...
    'LUMBARCOR', ...
    'LHIPCOR', ...
    'RHIPCOR', ...
    'LKNEECOR', ...
    'RKNEECOR', ...
    'LANKLECOR', ...
    'RANKLECOR' ...
  };

% specify markers that define segments
segment_labels = ...
  { ...
    'HEAD', ...
    'TORSO', ...
    'LUPPERARM', ...
    'RUPPERARM', ...
    'LFOREARM', ...
    'RFOREARM', ...
    'LHAND', ...
    'RHAND', ...
    'PELVIS', ...
    'LTHIGH', ...
    'RTHIGH', ...
    'LSHANK', ...
    'RSHANK', ...
    'LFOOT', ...
    'RFOOT' ...
  };

markers_by_segment = ...
  {
    'RFHD', 'LFHD', 'RBHD'; ...                     % head
    'C7', 'CLAV', 'T10'; ...                        % torso
    'LSHOULDERCOR', 'LELBOWCOR', 'LELB'; ...        % left upper arm
    'RSHOULDERCOR', 'RELBOWCOR', 'RELB'; ...        % right upper arm
    'LELB', 'LWRA', 'LWRB'; ...                     % left forearm
    'RELB', 'RWRA', 'RWRB'; ...                     % right forearm
    'LWRA', 'LWRB', 'LFIN'; ...                     % left hand
    'RWRA', 'RWRB', 'RFIN'; ...                     % right hand
    'RASI', 'LASI', 'RPSI'; ...                     % pelvis
    'LHIPCOR', 'LKNEECOR', 'LKNE'; ...              % left thigh
    'RHIPCOR', 'RKNEECOR', 'RKNE'; ...              % right thigh
    'LKNEECOR', 'LKNE', 'LANK'; ...                 % left shank
    'RKNEECOR', 'RKNE', 'RANK'; ...                 % right shank
    'LANK', 'LHEE', 'LTOE'; ...                     % left foot
    'RANK', 'RHEE', 'RTOE'; ...                     % right foot
  };
number_of_segments = size(markers_by_segment, 1);

mcs_to_wcs_transformations = calculateMcsToWcsTransformations([marker_reference joint_center_reference], [marker_headers joint_center_headers], markers_by_segment);

% assemble segment masses
segment_masses = ...
  [ ...
    head_segment_mass, ...
    torso_segment_mass, ...
    arm_segment_mass, ...
    arm_segment_mass, ...
    forearm_segment_mass, ...
    forearm_segment_mass, ...
    hand_segment_mass, ...
    hand_segment_mass, ...
    pelvis_segment_mass, ...
    thigh_segment_mass, ...
    thigh_segment_mass, ...
    leg_segment_mass, ...
    leg_segment_mass, ...
    foot_segment_mass, ...
    foot_segment_mass ...
  ];


% calculate segment CoMs in marker coordinates
segment_coms_wcs = ...
  { ...
    head_com; ...
    torso_com; ...
    left_arm_com; ...
    right_arm_com; ...
    left_forearm_com; ...
    right_forearm_com; ...
    left_hand_com; ...
    right_hand_com; ...
    pelvis_com; ...
    left_thigh_com; ...
    right_thigh_com; ...
    left_leg_com; ...
    right_leg_com; ...
    left_foot_com; ...
    right_foot_com; ...
  };
segment_coms_mcs = cell(size(segment_coms_wcs));
for i_segment = 1 : number_of_segments
    segment_com_wcs = [segment_coms_wcs{i_segment}; 1];
    T_wcs_to_mcs = mcs_to_wcs_transformations{i_segment}^(-1);
    segment_coms_mcs{i_segment} = eye(3, 4) * T_wcs_to_mcs * segment_com_wcs;
end
save ...
  ( ...
    'subjectModel', ...
    'marker_headers', ...
    'marker_reference', ...
    'joint_center_headers', ...
    'joint_center_reference', ...
    'segment_labels', ...
    'segment_masses', ...
    'markers_by_segment', ...
    'segment_coms_mcs' ...
  );

%% show visualization
if show_visualization

    % show stick figure of geometric model
    hip_center = (right_hip_cor + left_hip_cor) * 0.5;
    scene_bound = repmat(hip_center, 1, 2) + 2*[-0.5 0.5; -0.5 0.5; -1 1];
    stick_figure = KinematicTreeController(kinematic_tree, scene_bound, 'ellipsoid');
    
    % show segment CoMs
    for i_segment = 1 : number_of_segments
        segment_com_mcs = segment_coms_mcs{i_segment};
        T_mcs_to_wcs = mcs_to_wcs_transformations{i_segment};
        segment_com_wcs = T_mcs_to_wcs * [segment_com_mcs; 1];
        plot3(stick_figure.sceneAxes, segment_com_wcs(1), segment_com_wcs(2), segment_com_wcs(3), 'ko', 'linewidth', 2, 'markersize', 10);
    end
    
    % show joint centers
    for i_marker = 1 : (length(joint_center_reference) / 3);
        plot3 ...
          ( ...
            stick_figure.sceneAxes, ...
            joint_center_reference((i_marker-1)*3+1), joint_center_reference((i_marker-1)*3+2), joint_center_reference((i_marker-1)*3+3), ...
            'd', ...
            'color', [1 1 1]*0.5, ...
            'linewidth', 2, ...
            'markersize', 10 ...
          );
    end
    
end






























































