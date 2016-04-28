function kinematic_tree = walkingModel(markerReference, jointAxes, bodyMass, gender)

% TODO: the link inertia matrices should be set up so that the principal
% axes coincide with the actual segment axes, not the world axes, i.e. the
% z-axis of the link frame should point in the same direction as line
% between the two adjacent joints

% TODO: the position of the markers should be part of the constructor

% define measurements that I don't have right now
knee_width = 0.10;
ankle_width = 0.08;
elbow_width = 0.07;
head_thickness = 0.2;

%% define marker numbers and indices
LFHD_marker = 1;
RFHD_marker = 2;
LBHD_marker = 3;
RBHD_marker = 4;
C7_marker = 5;
CLAV_marker = 7;
LSHO_marker = 10;
LELB_marker = 12;
LWRA_marker = 14;
LWRB_marker = 15;
LFIN_marker = 16;
RSHO_marker = 17;
RELB_marker = 19;
RWRA_marker = 21;
RWRB_marker = 22;
RFIN_marker = 23;
LASI_marker = 24;
RASI_marker = 25;
LPSI_marker = 26;
RPSI_marker = 27;
LKNE_marker = 30;
LANK_marker = 33;
LHEE_marker = 34;
LTOE_marker = 35;
RKNE_marker = 38;
RANK_marker = 41;
RHEE_marker = 42;
RTOE_marker = 43;

LFHD_markers_indices = reshape([(LFHD_marker - 1) * 3 + 1; (LFHD_marker - 1) * 3 + 2; (LFHD_marker - 1) * 3 + 3], 1, length(LFHD_marker)*3);
RFHD_markers_indices = reshape([(RFHD_marker - 1) * 3 + 1; (RFHD_marker - 1) * 3 + 2; (RFHD_marker - 1) * 3 + 3], 1, length(RFHD_marker)*3);
LBHD_markers_indices = reshape([(LBHD_marker - 1) * 3 + 1; (LBHD_marker - 1) * 3 + 2; (LBHD_marker - 1) * 3 + 3], 1, length(LBHD_marker)*3);
RBHD_markers_indices = reshape([(RBHD_marker - 1) * 3 + 1; (RBHD_marker - 1) * 3 + 2; (RBHD_marker - 1) * 3 + 3], 1, length(RBHD_marker)*3);
C7_markers_indices = reshape([(C7_marker - 1) * 3 + 1; (C7_marker - 1) * 3 + 2; (C7_marker - 1) * 3 + 3], 1, length(C7_marker)*3);
CLAV_markers_indices = reshape([(CLAV_marker - 1) * 3 + 1; (CLAV_marker - 1) * 3 + 2; (CLAV_marker - 1) * 3 + 3], 1, length(CLAV_marker)*3);
LSHO_markers_indices = reshape([(LSHO_marker - 1) * 3 + 1; (LSHO_marker - 1) * 3 + 2; (LSHO_marker - 1) * 3 + 3], 1, length(LSHO_marker)*3);
LELB_markers_indices = reshape([(LELB_marker - 1) * 3 + 1; (LELB_marker - 1) * 3 + 2; (LELB_marker - 1) * 3 + 3], 1, length(LELB_marker)*3);
LWRA_markers_indices = reshape([(LWRA_marker - 1) * 3 + 1; (LWRA_marker - 1) * 3 + 2; (LWRA_marker - 1) * 3 + 3], 1, length(LWRA_marker)*3);
LWRB_markers_indices = reshape([(LWRB_marker - 1) * 3 + 1; (LWRB_marker - 1) * 3 + 2; (LWRB_marker - 1) * 3 + 3], 1, length(LWRB_marker)*3);
LFIN_markers_indices = reshape([(LFIN_marker - 1) * 3 + 1; (LFIN_marker - 1) * 3 + 2; (LFIN_marker - 1) * 3 + 3], 1, length(LFIN_marker)*3);
RSHO_markers_indices = reshape([(RSHO_marker - 1) * 3 + 1; (RSHO_marker - 1) * 3 + 2; (RSHO_marker - 1) * 3 + 3], 1, length(RSHO_marker)*3);
RELB_markers_indices = reshape([(RELB_marker - 1) * 3 + 1; (RELB_marker - 1) * 3 + 2; (RELB_marker - 1) * 3 + 3], 1, length(RELB_marker)*3);
RWRA_markers_indices = reshape([(RWRA_marker - 1) * 3 + 1; (RWRA_marker - 1) * 3 + 2; (RWRA_marker - 1) * 3 + 3], 1, length(RWRA_marker)*3);
RWRB_markers_indices = reshape([(RWRB_marker - 1) * 3 + 1; (RWRB_marker - 1) * 3 + 2; (RWRB_marker - 1) * 3 + 3], 1, length(RWRB_marker)*3);
RFIN_markers_indices = reshape([(RFIN_marker - 1) * 3 + 1; (RFIN_marker - 1) * 3 + 2; (RFIN_marker - 1) * 3 + 3], 1, length(RFIN_marker)*3);
LASI_markers_indices = reshape([(LASI_marker - 1) * 3 + 1; (LASI_marker - 1) * 3 + 2; (LASI_marker - 1) * 3 + 3], 1, length(LASI_marker)*3);
RASI_markers_indices = reshape([(RASI_marker - 1) * 3 + 1; (RASI_marker - 1) * 3 + 2; (RASI_marker - 1) * 3 + 3], 1, length(RASI_marker)*3);
LPSI_markers_indices = reshape([(LPSI_marker - 1) * 3 + 1; (LPSI_marker - 1) * 3 + 2; (LPSI_marker - 1) * 3 + 3], 1, length(LPSI_marker)*3);
RPSI_markers_indices = reshape([(RPSI_marker - 1) * 3 + 1; (RPSI_marker - 1) * 3 + 2; (RPSI_marker - 1) * 3 + 3], 1, length(RPSI_marker)*3);
LKNE_markers_indices = reshape([(LKNE_marker - 1) * 3 + 1; (LKNE_marker - 1) * 3 + 2; (LKNE_marker - 1) * 3 + 3], 1, length(LKNE_marker)*3);
LANK_markers_indices = reshape([(LANK_marker - 1) * 3 + 1; (LANK_marker - 1) * 3 + 2; (LANK_marker - 1) * 3 + 3], 1, length(LANK_marker)*3);
LHEE_markers_indices = reshape([(LHEE_marker - 1) * 3 + 1; (LHEE_marker - 1) * 3 + 2; (LHEE_marker - 1) * 3 + 3], 1, length(LHEE_marker)*3);
LTOE_markers_indices = reshape([(LTOE_marker - 1) * 3 + 1; (LTOE_marker - 1) * 3 + 2; (LTOE_marker - 1) * 3 + 3], 1, length(LTOE_marker)*3);
RKNE_markers_indices = reshape([(RKNE_marker - 1) * 3 + 1; (RKNE_marker - 1) * 3 + 2; (RKNE_marker - 1) * 3 + 3], 1, length(RKNE_marker)*3);
RANK_markers_indices = reshape([(RANK_marker - 1) * 3 + 1; (RANK_marker - 1) * 3 + 2; (RANK_marker - 1) * 3 + 3], 1, length(RANK_marker)*3);
RHEE_markers_indices = reshape([(RHEE_marker - 1) * 3 + 1; (RHEE_marker - 1) * 3 + 2; (RHEE_marker - 1) * 3 + 3], 1, length(RHEE_marker)*3);
RTOE_markers_indices = reshape([(RTOE_marker - 1) * 3 + 1; (RTOE_marker - 1) * 3 + 2; (RTOE_marker - 1) * 3 + 3], 1, length(RTOE_marker)*3);

%% extract marker reference positions
LFHD_reference = markerReference(LFHD_markers_indices)';
RFHD_reference = markerReference(RFHD_markers_indices)';
LBHD_reference = markerReference(LBHD_markers_indices)';
RBHD_reference = markerReference(RBHD_markers_indices)';
C7_reference = markerReference(C7_markers_indices)';
CLAV_reference = markerReference(CLAV_markers_indices)';
LSHO_reference = markerReference(LSHO_markers_indices)';
LELB_reference = markerReference(LELB_markers_indices)';
LWRA_reference = markerReference(LWRA_markers_indices)';
LWRB_reference = markerReference(LWRB_markers_indices)';
LFIN_reference = markerReference(LFIN_markers_indices)';
RSHO_reference = markerReference(RSHO_markers_indices)';
RELB_reference = markerReference(RELB_markers_indices)';
RWRA_reference = markerReference(RWRA_markers_indices)';
RWRB_reference = markerReference(RWRB_markers_indices)';
RFIN_reference = markerReference(RFIN_markers_indices)';
LASI_reference = markerReference(LASI_markers_indices)';
RASI_reference = markerReference(RASI_markers_indices)';
LPSI_reference = markerReference(LPSI_markers_indices)';
RPSI_reference = markerReference(RPSI_markers_indices)';
LKNE_reference = markerReference(LKNE_markers_indices)';
LANK_reference = markerReference(LANK_markers_indices)';
LHEE_reference = markerReference(LHEE_markers_indices)';
LTOE_reference = markerReference(LTOE_markers_indices)';
RKNE_reference = markerReference(RKNE_markers_indices)';
RANK_reference = markerReference(RANK_markers_indices)';
RHEE_reference = markerReference(RHEE_markers_indices)';
RTOE_reference = markerReference(RTOE_markers_indices)';

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
skin_to_correction_ASIS = 0.01; % in meters
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
LASIS = LASI_reference + (centroid_to_skin_correction + skin_to_correction_ASIS) * posterior_direction;
RASIS = RASI_reference + (centroid_to_skin_correction + skin_to_correction_ASIS) * posterior_direction;
left_lateral_femoral_epicondyle = LKNE_reference + centroid_to_skin_correction * right_direction;
left_lateral_malleolus = LANK_reference + centroid_to_skin_correction * right_direction;
left_calcaneus = LHEE_reference + centroid_to_skin_correction * anterior_direction;
left_toe_mid = LTOE_reference + centroid_to_skin_correction * distal_direction;
right_lateral_femoral_epicondyle = RKNE_reference + centroid_to_skin_correction * left_direction;
right_lateral_malleolus = RANK_reference + centroid_to_skin_correction * left_direction;
right_calcaneus = RHEE_reference + centroid_to_skin_correction * anterior_direction;
right_toe_mid = RTOE_reference + centroid_to_skin_correction * distal_direction;
head_vertex = mean([LFHD_reference RFHD_reference LBHD_reference RBHD_reference], 2) + head_thickness/2 * proximal_direction;
sellion = mean([LFHD_reference RFHD_reference], 2);

%% calculate joint centers and axes
% define correction factors for joint centers
ajc_correction_factor = 0.5;
kjc_correction_factor = 0.5;
hjc_correction_factor_lateral = 0.11; % Tylkowski, after Bell et al., 1990
hjc_correction_factor_distal = 0.12; % Tylkowski, after Bell et al., 1990
hjc_correction_factor_frontal = 0.21; % Tylkowski, after Bell et al., 1990
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
inter_ASIS_distance = norm(LASIS - RASIS);
thorax_width = norm(suprasternale - c7);

% estimate hip CoRs
left_hip_cor = LASIS ...
                + hjc_correction_factor_lateral * inter_ASIS_distance * right_direction ...
                + hjc_correction_factor_distal * inter_ASIS_distance * distal_direction ...
                + hjc_correction_factor_frontal * inter_ASIS_distance * posterior_direction;
right_hip_cor = RASIS ...
                + hjc_correction_factor_lateral * inter_ASIS_distance * left_direction ...
                + hjc_correction_factor_distal * inter_ASIS_distance * distal_direction ...
                + hjc_correction_factor_frontal * inter_ASIS_distance * posterior_direction;    

% check that knee axes are pointing in the right direction (leftward)
left_knee_aor = normVector(jointAxes{1});
right_knee_aor = normVector(jointAxes{2});
if dot(right_knee_aor, left_direction) < 0
    right_knee_aor = -right_knee_aor;
end
if dot(left_knee_aor, left_direction) < 0
    left_knee_aor = -left_knee_aor;
end
            
% estimate knee CoRs
left_knee_cor = left_lateral_femoral_epicondyle - kjc_correction_factor*knee_width*left_knee_aor;
right_knee_cor = right_lateral_femoral_epicondyle + kjc_correction_factor*knee_width*right_knee_aor;

% estimate ankle CoRs
left_ankle_cor = left_lateral_malleolus - ajc_correction_factor*ankle_width*left_knee_aor;
right_ankle_cor = right_lateral_malleolus + ajc_correction_factor*ankle_width*right_knee_aor;

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
left_wrist_inversion_axis = distal_direction;
right_wrist_inversion_axis = proximal_direction;

% estimate elbow axes and joint centers
left_elbow_axis = distal_direction;
right_elbow_axis = proximal_direction;
left_elbow_cor = left_lateral_humeral_epicondyle + ejc_correction_factor*elbow_width*left_elbow_axis;
right_elbow_cor = right_lateral_humeral_epicondyle - ejc_correction_factor*elbow_width*right_elbow_axis;
left_radioulnar_axis = normVector(left_wrist_cor - left_elbow_cor);
right_radioulnar_axis = - normVector(right_wrist_cor - right_elbow_cor);

%% define scaling factors
if strcmp(gender, 'male')
    % mass
    head_mass_scaling_factor =      0.067;
    torso_mass_scaling_factor =     0.333;
    arm_mass_scaling_factor =       0.024;
    forearm_mass_scaling_factor =   0.017;
    hand_mass_scaling_factor =      0.006;
    pelvis_mass_scaling_factor =    0.142;
    thigh_mass_scaling_factor =     0.123;
    shank_mass_scaling_factor =     0.048;
    foot_mass_scaling_factor =      0.012;
    
    % CoM
    head_com_scaling_factor_x =    -0.062;  head_com_scaling_factor_y =     0.555;  head_com_scaling_factor_z =     0.001;
    torso_com_scaling_factor_x =   -0.036;  torso_com_scaling_factor_y =   -0.420;  torso_com_scaling_factor_z =   -0.002;
    arm_com_scaling_factor_x =      0.017;  arm_com_scaling_factor_y =     -0.452;  arm_com_scaling_factor_z =     -0.026;
    forearm_com_scaling_factor_x =  0.010;  forearm_com_scaling_factor_y = -0.417;  forearm_com_scaling_factor_z =  0.014;
    hand_com_scaling_factor_x =     0.082;  hand_com_scaling_factor_y =    -0.839;  hand_com_scaling_factor_z =     0.074;
    pelvis_com_scaling_factor_x =   0.028;  pelvis_com_scaling_factor_y =  -0.280;  pelvis_com_scaling_factor_z =  -0.006;
    thigh_com_scaling_factor_x =   -0.041;  thigh_com_scaling_factor_y =   -0.429;  thigh_com_scaling_factor_z =    0.033;
    shank_com_scaling_factor_x =   -0.048;  shank_com_scaling_factor_y =   -0.410;  shank_com_scaling_factor_z =    0.007;
    foot_com_scaling_factor_x =     0.382;  foot_com_scaling_factor_y =    -0.151;  foot_com_scaling_factor_z =     0.026;
    
    % rog
    head_rxx_scaling_factor =    0.31;  head_ryy_scaling_factor =    0.25;  head_rzz_scaling_factor =    0.33;  head_rxy_scaling_factor =    0.09*1i;	head_rxz_scaling_factor =    0.02*1i;   head_ryz_scaling_factor =    0.03;
    torso_rxx_scaling_factor =   0.27;  torso_ryy_scaling_factor =   0.25;  torso_rzz_scaling_factor =   0.28;  torso_rxy_scaling_factor =   0.18;      torso_rxz_scaling_factor =   0.02;      torso_ryz_scaling_factor =   0.04*1i;
    arm_rxx_scaling_factor =     0.31;  arm_ryy_scaling_factor =     0.14;  arm_rzz_scaling_factor =     0.32;  arm_rxy_scaling_factor =     0.06;      arm_rxz_scaling_factor =     0.05;      arm_ryz_scaling_factor =     0.02;
    forearm_rxx_scaling_factor = 0.28;  forearm_ryy_scaling_factor = 0.11;  forearm_rzz_scaling_factor = 0.27;  forearm_rxy_scaling_factor = 0.03;      forearm_rxz_scaling_factor = 0.02;      forearm_ryz_scaling_factor = 0.08*1i;
    hand_rxx_scaling_factor =    0.61;  hand_ryy_scaling_factor =    0.38;  hand_rzz_scaling_factor =    0.56;  hand_rxy_scaling_factor =    0.22;      hand_rxz_scaling_factor =    0.15;      hand_ryz_scaling_factor =    0.20*1i;
    pelvis_rxx_scaling_factor =  1.01;  pelvis_ryy_scaling_factor =  1.06;  pelvis_rzz_scaling_factor =  0.95;  pelvis_rxy_scaling_factor =  0.25*1i;   pelvis_rxz_scaling_factor =  0.12*1i;   pelvis_ryz_scaling_factor =  0.08*1i;
    thigh_rxx_scaling_factor =   0.29;  thigh_ryy_scaling_factor =   0.15;  thigh_rzz_scaling_factor =   0.30;  thigh_rxy_scaling_factor =   0.07;      thigh_rxz_scaling_factor =   0.02*1i;   thigh_ryz_scaling_factor =   0.07*1i;
    shank_rxx_scaling_factor =   0.28;  shank_ryy_scaling_factor =   0.10;  shank_rzz_scaling_factor =   0.28;  shank_rxy_scaling_factor =   0.04*1i;   shank_rxz_scaling_factor =   0.02*1i;   shank_ryz_scaling_factor =   0.05;
    foot_rxx_scaling_factor =    0.17;  foot_ryy_scaling_factor =    0.37;	foot_rzz_scaling_factor =    0.36;  foot_rxy_scaling_factor =    0.13;      foot_rxz_scaling_factor =    0.08*1i;	foot_ryz_scaling_factor =    0.00;
elseif strcmp(gender, 'female')
    % mass
    head_mass_scaling_factor = 0.067;
    torso_mass_scaling_factor = 0.304;
    arm_mass_scaling_factor = 0.022;
    forearm_mass_scaling_factor = 0.013;
    hand_mass_scaling_factor = 0.005;
    pelvis_mass_scaling_factor = 0.146;
    thigh_mass_scaling_factor = 0.146;
    shank_mass_scaling_factor = 0.045;
    foot_mass_scaling_factor = 0.010;
    
    % CoM
    head_com_scaling_factor_x = -0.070;     head_com_scaling_factor_y = 0.597;      head_com_scaling_factor_z = 0.000;
    torso_com_scaling_factor_x = -0.016;    torso_com_scaling_factor_y = -0.436;    torso_com_scaling_factor_z = -0.006;
    pelvis_com_scaling_factor_x = -0.009;   pelvis_com_scaling_factor_y = -0.232;   pelvis_com_scaling_factor_z = 0.002;
    thigh_com_scaling_factor_x = -0.077;    thigh_com_scaling_factor_y = -0.377;    thigh_com_scaling_factor_z =  0.009;
    shank_com_scaling_factor_x = -0.049;    shank_com_scaling_factor_y = -0.404;    shank_com_scaling_factor_z =  0.031;
    foot_com_scaling_factor_x =  0.270;     foot_com_scaling_factor_y = -0.218;     foot_com_scaling_factor_z =  0.039;
    head_rxx_scaling_factor = 0.32;         head_ryy_scaling_factor = 0.27;         head_rzz_scaling_factor = 0.34;
    head_rxy_scaling_factor = 0.06*1i;      head_rxz_scaling_factor = 0.01;         head_ryz_scaling_factor = 0.01*1i;

    torso_rxx_scaling_factor = 0.29;        torso_ryy_scaling_factor = 0.27;        torso_rzz_scaling_factor = 0.29;    torso_rxy_scaling_factor = 0.22;        torso_rxz_scaling_factor = 0.05;        torso_ryz_scaling_factor = 0.05*1i;
    % ACHTUNG: used the values for males here because the values for females fails to visualize as an ellipsoid.
    % Explore later!
    torso_rxx_scaling_factor = 0.27;        torso_ryy_scaling_factor = 0.25;        torso_rzz_scaling_factor = 0.28;    torso_rxy_scaling_factor = 0.18;        torso_rxz_scaling_factor = 0.02;        torso_ryz_scaling_factor = 0.04*1i;

    pelvis_rxx_scaling_factor = 0.91;       pelvis_ryy_scaling_factor = 1.00;       pelvis_rzz_scaling_factor = 0.79;   pelvis_rxy_scaling_factor = 0.34*1i;    pelvis_rxz_scaling_factor = 0.01*1i;    pelvis_ryz_scaling_factor = 0.01*1i;
    thigh_rxx_scaling_factor = 0.31;        thigh_ryy_scaling_factor = 0.19;        thigh_rzz_scaling_factor = 0.32;    thigh_rxy_scaling_factor = 0.07;        thigh_rxz_scaling_factor = 0.02*1i;     thigh_ryz_scaling_factor = 0.07*1i;
    shank_rxx_scaling_factor = 0.28;        shank_ryy_scaling_factor = 0.10;        shank_rzz_scaling_factor = 0.28;    shank_rxy_scaling_factor = 0.02;        shank_rxz_scaling_factor = 0.01;        shank_ryz_scaling_factor = 0.06;
    foot_rxx_scaling_factor = 0.17;         foot_ryy_scaling_factor = 0.36;         foot_rzz_scaling_factor = 0.35;     foot_rxy_scaling_factor = 0.10*1i;      foot_rxz_scaling_factor = 0.06;         foot_ryz_scaling_factor = 0.04*1i;
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
left_arm_scs_z = normVector(cross(left_arm_scs_y, left_elbow_axis));
left_arm_scs_x = cross(left_arm_scs_y, left_arm_scs_z);

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
right_arm_scs_z = normVector(cross(right_arm_scs_y, right_elbow_axis));
right_arm_scs_x = cross(right_arm_scs_y, right_arm_scs_z);

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
pelvis_segment_mass = pelvis_mass_scaling_factor    * bodyMass * 1.002^(-1);
thigh_segment_mass = thigh_mass_scaling_factor      * bodyMass * 1.002^(-1);
leg_segment_mass = shank_mass_scaling_factor        * bodyMass * 1.002^(-1);
foot_segment_mass = foot_mass_scaling_factor        * bodyMass * 1.002^(-1);
head_segment_mass = head_mass_scaling_factor        * bodyMass * 1.002^(-1);
torso_segment_mass = torso_mass_scaling_factor      * bodyMass * 1.002^(-1);
arm_segment_mass = arm_mass_scaling_factor          * bodyMass * 1.002^(-1);
forearm_segment_mass = forearm_mass_scaling_factor  * bodyMass * 1.002^(-1);
hand_segment_mass = hand_mass_scaling_factor        * bodyMass * 1.002^(-1);

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

%% assemble

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
    anterior_direction, distal_direction, left_direction, left_elbow_axis, left_radioulnar_axis, left_wrist_flexion_axis, left_wrist_inversion_axis, ... % left arm
    posterior_direction, proximal_direction, left_direction, right_elbow_axis, right_radioulnar_axis, right_wrist_flexion_axis, right_wrist_inversion_axis, ... % right arm
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


% markers
red = [1 0 0];
green = [0 1 0];
blue = [0 0 1];
yellow = [1 0.9 0];

markerSegments = ...
  [ ...
    24 24 24 24 ...             % head
    21 21 21 21 21 ...          % trunk
    21 ...                      % left shoulder
    27 27 ...                   % left upper arm
    29 29 29 ...                % left lower arm
    31 ...                      % left hand
    21 ...                      % right shoulder
    34 34 ...                   % right upper arm
    36 36 36 ...                % right lower arm
    38 ...                      % right hand
    6 6 6 6 ...                 % pelvis
    9 9 9 ...                   % left thigh
    10 10 10 ...                % left shank
    12 12 ...                   % left foot
    15 15 15 ...                % right thigh
    16 16 16 ...                % right shank
    18 18 ...                   % right foot
  ];

marker_color_list = ...
  { ...
    blue blue blue blue ...                         % head
    blue blue blue blue blue ...                    % trunk
    red red red red red red red ...                 % left arm (fixed to trunk for now)
    green green green green green green green ...   % right arm (fixed to trunk for now)
    blue blue blue blue ...                         % pelvis
    red red red ...                                 % left thigh
    red red red ...                                 % left shank
    red red ...                                     % left foot
    green green green  ...                          % right thigh
    green green green  ...                          % right shank
    green green ...                                 % right foot
  };

number_of_markers = length(markerSegments);
for i_marker = 1 : number_of_markers
    kinematic_tree.addMarker(markerSegments(i_marker), markerReference((i_marker-1)*3+1 : (i_marker-1)*3+3)', marker_color_list{i_marker});
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
    

    
    
    
