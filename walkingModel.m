function kinematic_tree = walkingModel(bodyMass, bodyHeight, jointPositions, jointAxes, markerReference, gender)

% TODO: the link inertia matrices should be set up so that the principal
% axes coincide with the actual segment axes, not the world axes, i.e. the
% z-axis of the link frame should point in the same direction as line
% between the two adjacent joints

% TODO: the position of the markers should be part of the constructor


if nargin < 4
    % set up joint positions from anthropometric estimates according to Winter, 1990
    ankle_height = 0.039*bodyHeight;
    knee_height = 0.285*bodyHeight;
    hip_height = 0.530*bodyHeight;
    hip_width = 0.191*bodyHeight;
    hip_excursion = hip_width * 0.5;
    pelvis_height = 0.630*bodyHeight;
    neck_height = 0.870*bodyHeight;
    shoulder_height = 0.818*bodyHeight;
    head_height = 0.936*bodyHeight;
    foot_length = 0.09*bodyHeight;
    toes_length = 0.02*bodyHeight;
    shank_length = knee_height - ankle_height;
    thigh_length = hip_height - knee_height;
    pelvis_length = pelvis_height - hip_height;
    gtgh_length = shoulder_height - hip_height; % greater trochanter to glenohumeral joint
    torso_length = gtgh_length;
    head_length = bodyHeight - neck_height;
    metatarsal_height = 0.01*bodyHeight;
    
    left_met_mean = [foot_length; hip_excursion; 0];
    left_metatarsal_cor = left_met_mean;
    left_toetip_position = [foot_length + toes_length; hip_excursion; 0];
    left_ankle_cor = [0; hip_excursion; ankle_height];
    left_knee_cor = [0; hip_excursion; knee_height];
    left_hip_cor = [0; hip_excursion; hip_height];
    right_hip_cor = [0; -hip_excursion; hip_height];
    right_knee_cor = [0; -hip_excursion; knee_height];
    right_ankle_cor = [0; -hip_excursion; ankle_height];
    right_met_mean = [foot_length; -hip_excursion; 0];
    right_metatarsal_cor = right_met_mean;
    right_toetip_position = [foot_length + toes_length; -hip_excursion; 0];
    hip_center_position = [0; 0; hip_height];
    l5_position = [0; 0; pelvis_height];
    neck_position = [0; 0; neck_height];
    head_vertex_position = [0; 0; bodyHeight];
    
    lumbar_cor = l5_position;
    cervix_cor = neck_position;
    
    shoulder_center_position = [0; 0; shoulder_height];
    head_center_position = [0; 0; head_height];
    
    foot_com_from_distal = 0.5*ankle_height;
    foot_com_offset_ap = 0.05;
    toes_com_from_distal = 0.5*metatarsal_height;
    toes_com_offset_ap = 0.05;
    shank_com_height = ankle_height + 0.567*shank_length;
    thigh_com_height = knee_height + 0.567*thigh_length;
    pelvis_com_height = hip_height + 0.105*gtgh_length;
    trunk_com_height = hip_height + 0.63*gtgh_length;
    head_com_height = head_height;
    
    right_thigh_com = [0; -hip_excursion; thigh_com_height];
    right_leg_com = [0; -hip_excursion; shank_com_height];
    right_foot_com = [foot_com_offset_ap; -hip_excursion; foot_com_from_distal];
    right_toes_com = [toes_com_offset_ap; -hip_excursion; toes_com_from_distal];
    left_thigh_com = [0; hip_excursion; thigh_com_height];
    left_leg_com = [0; hip_excursion; shank_com_height];
    left_foot_com = [foot_com_offset_ap; hip_excursion; foot_com_from_distal];
    left_toes_com = [toes_com_offset_ap; hip_excursion; toes_com_from_distal];
    pelvis_com = [0; 0; pelvis_com_height];
    torso_com = [0; 0; trunk_com_height];
    head_com = [0; 0; head_com_height];
    trunk_center_position = torso_com;
    
    foot_rog = 0.475 * foot_length;
    foot_segment_mass = 0.0145*bodyMass;
    foot_inertia_y = foot_segment_mass * foot_rog^2;
    foot_inertia_x = foot_inertia_y * 0.1;
    foot_inertia_z = foot_inertia_y;
    right_foot_inertia_tensor = [foot_inertia_x 0 0; 0 foot_inertia_y 0; 0 0 foot_inertia_z];
    left_foot_inertia_tensor = [foot_inertia_x 0 0; 0 foot_inertia_y 0; 0 0 foot_inertia_z];
    right_foot_scs_x = [1; 0; 0];
    right_foot_scs_y = [0; 1; 0];
    right_foot_scs_z = [0; 0; 1];
    left_foot_scs_x = [1; 0; 0];
    left_foot_scs_y = [0; 1; 0];
    left_foot_scs_z = [0; 0; 1];
    
    shank_rog = 0.302 * shank_length;
    leg_segment_mass = 0.0465*bodyMass;
    leg_inertia_x = leg_segment_mass * shank_rog^2;
    leg_inertia_y = leg_inertia_x;
    leg_inertia_z = leg_inertia_x * 0.1;
    right_leg_inertia_tensor = [leg_inertia_x 0 0; 0 leg_inertia_y 0; 0 0 leg_inertia_z];
    left_leg_inertia_tensor = [leg_inertia_x 0 0; 0 leg_inertia_y 0; 0 0 leg_inertia_z];
    right_leg_scs_x = [1; 0; 0];
    right_leg_scs_y = [0; 1; 0];
    right_leg_scs_z = [0; 0; 1];
    left_leg_scs_x = [1; 0; 0];
    left_leg_scs_y = [0; 1; 0];
    left_leg_scs_z = [0; 0; 1];
    
    thigh_rog = 0.323 * thigh_length;
    thigh_segment_mass = 0.1*bodyMass;
    thigh_inertia_x = thigh_segment_mass * thigh_rog^2;
    thigh_inertia_y = thigh_inertia_x;
    thigh_inertia_z = thigh_inertia_x * 0.1;
    right_thigh_inertia_tensor = [thigh_inertia_x 0 0; 0 thigh_inertia_y 0; 0 0 thigh_inertia_z];
    left_thigh_inertia_tensor = [thigh_inertia_x 0 0; 0 thigh_inertia_y 0; 0 0 thigh_inertia_z];
    right_thigh_scs_x = [1; 0; 0];
    right_thigh_scs_y = [0; 1; 0];
    right_thigh_scs_z = [0; 0; 1];
    left_thigh_scs_x = [1; 0; 0];
    left_thigh_scs_y = [0; 1; 0];
    left_thigh_scs_z = [0; 0; 1];
    
    pelvis_rog = 0.1 * pelvis_length;
    pelvis_segment_mass = 0.142*bodyMass;
    pelvis_inertia_x = pelvis_segment_mass * pelvis_rog^2;
    pelvis_inertia_y = pelvis_inertia_x;
    pelvis_inertia_z = pelvis_inertia_x * 0.1;
    pelvis_inertia_tensor = [pelvis_inertia_x 0 0; 0 pelvis_inertia_y 0; 0 0 pelvis_inertia_z];
    pelvis_scs_x = [1; 0; 0];
    pelvis_scs_y = [0; 1; 0];
    pelvis_scs_z = [0; 0; 1];
    
    torso_rog = 0.3 * torso_length;
    torso_segment_mass = 0.455*bodyMass;
    torso_inertia_x = torso_segment_mass * torso_rog^2;
    torso_inertia_y = torso_inertia_x * 0.9;
    torso_inertia_z = torso_inertia_x * 0.7;
    torso_inertia_tensor = [torso_inertia_x 0 0; 0 torso_inertia_y 0; 0 0 torso_inertia_z];
    torso_scs_x = [1; 0; 0];
    torso_scs_y = [0; 1; 0];
    torso_scs_z = [0; 0; 1];
    
    head_rog = 0.3 * head_length;
    head_segment_mass = 0.081*bodyMass;
    head_inertia_x = head_segment_mass * head_rog^2;
    head_inertia_y = head_inertia_x * 0.6;
    head_inertia_z = head_inertia_x * 0.8;
    head_inertia_tensor = [head_inertia_x 0 0; 0 head_inertia_y 0; 0 0 head_inertia_z];
    head_scs_x = [1; 0; 0];
    head_scs_y = [0; 1; 0];
    head_scs_z = [0; 0; 1];
    
    right_knee_aor = [0 1 0]';
    left_knee_aor = [0 1 0]';
    right_metatarsal_aor = [0 1 0]';
    left_metatarsal_aor = [0 1 0]';

%     right_toe = right_met_mean - [0; 0; 0.03]; 
    right_palm = right_met_mean - [0; 0; 0.03]; 
    right_heel = right_ankle_cor; right_heel(3) = right_palm(3);
    left_palm = left_met_mean - [0; 0; 0.03];
    left_heel = left_ankle_cor; left_heel(3) = left_palm(3);
    
else
    right_hip_cor = jointPositions{1};
    right_knee_cor = jointPositions{2};
    right_ankle_cor = jointPositions{3};
    left_hip_cor = jointPositions{4};
    left_knee_cor = jointPositions{5};
    left_ankle_cor = jointPositions{6};
    
    if strcmp(gender, 'male')
        head_mass_scaling_factor = 0.067;
        torso_mass_scaling_factor = 0.333;
        arm_mass_scaling_factor = 0.024;
        forearm_mass_scaling_factor = 0.017;
        hand_mass_scaling_factor = 0.006;
        pelvis_mass_scaling_factor = 0.142;
        thigh_mass_scaling_factor = 0.123;
        shank_mass_scaling_factor = 0.048;
        foot_mass_scaling_factor = 0.012;
        head_com_scaling_factor_x = -0.062;
        head_com_scaling_factor_y =  0.555;
        head_com_scaling_factor_z =  0.001;
        torso_com_scaling_factor_x = -0.036;
        torso_com_scaling_factor_y = -0.420;
        torso_com_scaling_factor_z = -0.002;
        pelvis_com_scaling_factor_x =  0.028;
        pelvis_com_scaling_factor_y = -0.280;
        pelvis_com_scaling_factor_z = -0.006;
        thigh_com_scaling_factor_x = -0.041;
        thigh_com_scaling_factor_y = -0.429;
        thigh_com_scaling_factor_z =  0.033;
        shank_com_scaling_factor_x = -0.048;
        shank_com_scaling_factor_y = -0.410;
        shank_com_scaling_factor_z =  0.007;
        foot_com_scaling_factor_x =  0.382;
        foot_com_scaling_factor_y = -0.151;
        foot_com_scaling_factor_z =  0.026;
        head_rxx_scaling_factor = 0.31;
        head_ryy_scaling_factor = 0.25;
        head_rzz_scaling_factor = 0.33;
        head_rxy_scaling_factor = 0.09*1i;
        head_rxz_scaling_factor = 0.02*1i;
        head_ryz_scaling_factor = 0.03;
        torso_rxx_scaling_factor = 0.27;
        torso_ryy_scaling_factor = 0.25;
        torso_rzz_scaling_factor = 0.28;
        torso_rxy_scaling_factor = 0.18;
        torso_rxz_scaling_factor = 0.02;
        torso_ryz_scaling_factor = 0.04*1i;
        pelvis_rxx_scaling_factor = 1.01;
        pelvis_ryy_scaling_factor = 1.06;
        pelvis_rzz_scaling_factor = 0.95;
        pelvis_rxy_scaling_factor = 0.25*1i;
        pelvis_rxz_scaling_factor = 0.12*1i;
        pelvis_ryz_scaling_factor = 0.08*1i;
        thigh_rxx_scaling_factor = 0.29;
        thigh_ryy_scaling_factor = 0.15;
        thigh_rzz_scaling_factor = 0.30;
        thigh_rxy_scaling_factor = 0.07;
        thigh_rxz_scaling_factor = 0.02*1i;
        thigh_ryz_scaling_factor = 0.07*1i;
        shank_rxx_scaling_factor = 0.28;
        shank_ryy_scaling_factor = 0.10;
        shank_rzz_scaling_factor = 0.28;
        shank_rxy_scaling_factor = 0.04*1i;
        shank_rxz_scaling_factor = 0.02*1i;
        shank_ryz_scaling_factor = 0.05;
        foot_rxx_scaling_factor = 0.17;
        foot_ryy_scaling_factor = 0.37;
        foot_rzz_scaling_factor = 0.36;
        foot_rxy_scaling_factor = 0.13;
        foot_rxz_scaling_factor = 0.08*1i;
        foot_ryz_scaling_factor = 0.00;
    elseif strcmp(gender, 'female')
        head_mass_scaling_factor = 0.067;
        torso_mass_scaling_factor = 0.304;
        arm_mass_scaling_factor = 0.022;
        forearm_mass_scaling_factor = 0.013;
        hand_mass_scaling_factor = 0.005;
        pelvis_mass_scaling_factor = 0.146;
        thigh_mass_scaling_factor = 0.146;
        shank_mass_scaling_factor = 0.045;
        foot_mass_scaling_factor = 0.010;
        head_com_scaling_factor_x = -0.070;
        head_com_scaling_factor_y = 0.597;
        head_com_scaling_factor_z = 0.000;
        torso_com_scaling_factor_x = -0.016;
        torso_com_scaling_factor_y = -0.436;
        torso_com_scaling_factor_z = -0.006;
        pelvis_com_scaling_factor_x = -0.009;
        pelvis_com_scaling_factor_y = -0.232;
        pelvis_com_scaling_factor_z = 0.002;
        thigh_com_scaling_factor_x = -0.077;
        thigh_com_scaling_factor_y = -0.377;
        thigh_com_scaling_factor_z =  0.009;
        shank_com_scaling_factor_x = -0.049;
        shank_com_scaling_factor_y = -0.404;
        shank_com_scaling_factor_z =  0.031;
        foot_com_scaling_factor_x =  0.270;
        foot_com_scaling_factor_y = -0.218;
        foot_com_scaling_factor_z =  0.039;
        head_rxx_scaling_factor = 0.32;
        head_ryy_scaling_factor = 0.27;
        head_rzz_scaling_factor = 0.34;
        head_rxy_scaling_factor = 0.06*1i;
        head_rxz_scaling_factor = 0.01;
        head_ryz_scaling_factor = 0.01*1i;
        
        torso_rxx_scaling_factor = 0.29;
        torso_ryy_scaling_factor = 0.27;
        torso_rzz_scaling_factor = 0.29;
        torso_rxy_scaling_factor = 0.22;
        torso_rxz_scaling_factor = 0.05;
        torso_ryz_scaling_factor = 0.05*1i;
        % ACHTUNG: used the values for males here because the values for females fails to visualize as an ellipsoid.
        % Explore later!
        torso_rxx_scaling_factor = 0.27;
        torso_ryy_scaling_factor = 0.25;
        torso_rzz_scaling_factor = 0.28;
        torso_rxy_scaling_factor = 0.18;
        torso_rxz_scaling_factor = 0.02;
        torso_ryz_scaling_factor = 0.04*1i;
        
        pelvis_rxx_scaling_factor = 0.91;
        pelvis_ryy_scaling_factor = 1.00;
        pelvis_rzz_scaling_factor = 0.79;
        pelvis_rxy_scaling_factor = 0.34*1i;
        pelvis_rxz_scaling_factor = 0.01*1i;
        pelvis_ryz_scaling_factor = 0.01*1i;
        thigh_rxx_scaling_factor = 0.31;
        thigh_ryy_scaling_factor = 0.19;
        thigh_rzz_scaling_factor = 0.32;
        thigh_rxy_scaling_factor = 0.07;
        thigh_rxz_scaling_factor = 0.02*1i;
        thigh_ryz_scaling_factor = 0.07*1i;
        shank_rxx_scaling_factor = 0.28;
        shank_ryy_scaling_factor = 0.10;
        shank_rzz_scaling_factor = 0.28;
        shank_rxy_scaling_factor = 0.02;
        shank_rxz_scaling_factor = 0.01;
        shank_ryz_scaling_factor = 0.06;
        foot_rxx_scaling_factor = 0.17;
        foot_ryy_scaling_factor = 0.36;
        foot_rzz_scaling_factor = 0.35;
        foot_rxy_scaling_factor = 0.10*1i;
        foot_rxz_scaling_factor = 0.06;
        foot_ryz_scaling_factor = 0.04*1i;
    end
    
    %% set up the segment coordinate systems (SCS)

    % define marker indices
    left_temple_marker_index = 1;
    right_temple_marker_index = 2;
    LBHD_marker_index = 3;
    RBHD_marker_index = 4;
    c7_marker_index = 5;
    suprasternale_marker_index = 7;
    left_ASIS_marker_index = 24;
    right_ASIS_marker_index = 25;
    left_PSIS_marker_index = 26;
    right_PSIS_marker_index = 27;
    left_calcaneous_marker_index = 34;
    left_toetip_marker_index = 35;
    right_calcaneous_marker_index = 42;
    right_toetip_marker_index = 43;

    % determine anatomic directions
    MPSIS = mean([markerReference((right_PSIS_marker_index-1)*3+1 : (right_PSIS_marker_index-1)*3+3)' markerReference((left_PSIS_marker_index-1)*3+1 : (left_PSIS_marker_index-1)*3+3)'], 2);
    MASIS = mean([markerReference((right_ASIS_marker_index-1)*3+1 : (right_ASIS_marker_index-1)*3+3)' markerReference((left_ASIS_marker_index-1)*3+1 : (left_ASIS_marker_index-1)*3+3)'], 2);
    MPSIS_to_MASIS = MASIS - MPSIS;
    RASIS = markerReference((right_ASIS_marker_index-1)*3+1 : (right_ASIS_marker_index-1)*3+3)';
    LASIS = markerReference((left_ASIS_marker_index-1)*3+1 : (left_ASIS_marker_index-1)*3+3)';
    
    anterior_direction = normVector(MASIS - MPSIS);
    posterior_direction = - anterior_direction;
    left_direction = normVector(LASIS - RASIS);
    right_direction = - left_direction;
    up_direction = cross(right_direction, anterior_direction);
    down_direction = - up_direction;
    
    centroid_to_skin_correction = 0.0152; % in meters
    skin_to_bone_correction = 0.01; % in meters
    centroid_to_bone_correction = centroid_to_skin_correction + skin_to_bone_correction;
    
    sellion_marker_offset           = centroid_to_skin_correction * posterior_direction;
    head_vertex_marker_offset       = centroid_to_skin_correction * down_direction;
    right_temple_marker_offset      = centroid_to_skin_correction * left_direction;
    left_temple_marker_offset       = centroid_to_skin_correction * right_direction;
    suprasternale_marker_offset     = centroid_to_skin_correction * posterior_direction;
    c7_marker_offset                = centroid_to_skin_correction * anterior_direction;
    MPSIS_marker_offset             = centroid_to_skin_correction * anterior_direction;
    right_ASIS_marker_offset        = centroid_to_skin_correction * posterior_direction;
    left_ASIS_marker_offset         = centroid_to_skin_correction * posterior_direction;
    right_calcaneous_marker_offset  = centroid_to_skin_correction * anterior_direction;
    right_met_one_marker_offset     = centroid_to_skin_correction * down_direction;
    right_met_five_marker_offset    = centroid_to_skin_correction * down_direction;
    right_toetip_marker_offset      = centroid_to_skin_correction * down_direction;
    left_calcaneous_marker_offset   = centroid_to_skin_correction * anterior_direction;
    left_met_one_marker_offset      = centroid_to_skin_correction * down_direction;
    left_met_five_marker_offset     = centroid_to_skin_correction * down_direction;
    left_toetip_marker_offset       = centroid_to_skin_correction * down_direction;
    
%     sellion_position            = markerReference((sellion_marker_index-1)*3+1 : (sellion_marker_index-1)*3+3)'                     + sellion_marker_offset;
%     head_vertex_position        = markerReference((head_vertex_marker_index-1)*3+1 : (head_vertex_marker_index-1)*3+3)'             + head_vertex_marker_offset;
%     right_temple_position       = markerReference((right_temple_marker_index-1)*3+1 : (right_temple_marker_index-1)*3+3)'           + right_temple_marker_offset;
%     left_temple_position        = markerReference((left_temple_marker_index-1)*3+1 : (left_temple_marker_index-1)*3+3)'             + left_temple_marker_offset;
%     suprasternale_position      = markerReference((suprasternale_marker_index-1)*3+1 : (suprasternale_marker_index-1)*3+3)'         + suprasternale_marker_offset;
%     c7_position                 = markerReference((c7_marker_index-1)*3+1 : (c7_marker_index-1)*3+3)'                               + c7_marker_offset;
%     MPSIS_position              = markerReference((MPSIS_marker_index-1)*3+1 : (MPSIS_marker_index-1)*3+3)'                         + MPSIS_marker_offset;
%     right_ASIS_position         = markerReference((right_ASIS_marker_index-1)*3+1 : (right_ASIS_marker_index-1)*3+3)'               + right_ASIS_marker_offset;
%     left_ASIS_position          = markerReference((left_ASIS_marker_index-1)*3+1 : (left_ASIS_marker_index-1)*3+3)'                 + left_ASIS_marker_offset;
%     right_calcaneous_position   = markerReference((right_calcaneous_marker_index-1)*3+1 : (right_calcaneous_marker_index-1)*3+3)'   + right_calcaneous_marker_offset;
%     right_met_one_position      = markerReference((right_met_one_marker_index-1)*3+1 : (right_met_one_marker_index-1)*3+3)'         + right_met_one_marker_offset;
%     right_met_five_position     = markerReference((right_met_five_marker_index-1)*3+1 : (right_met_five_marker_index-1)*3+3)'       + right_met_five_marker_offset;
%     right_toetip_position       = markerReference((right_toetip_marker_index-1)*3+1 : (right_toetip_marker_index-1)*3+3)'           + right_toetip_marker_offset;
%     left_calcaneous_position    = markerReference((left_calcaneous_marker_index-1)*3+1 : (left_calcaneous_marker_index-1)*3+3)'     + left_calcaneous_marker_offset;
%     left_met_one_position       = markerReference((left_met_one_marker_index-1)*3+1 : (left_met_one_marker_index-1)*3+3)'           + left_met_one_marker_offset;
%     left_met_five_position      = markerReference((left_met_five_marker_index-1)*3+1 : (left_met_five_marker_index-1)*3+3)'         + left_met_five_marker_offset;
%     left_toetip_position        = markerReference((left_toetip_marker_index-1)*3+1 : (left_toetip_marker_index-1)*3+3)'             + left_toetip_marker_offset;

    % rough and dirty version
    right_temple_position       = markerReference((right_temple_marker_index-1)*3+1 : (right_temple_marker_index-1)*3+3)';
    left_temple_position        = markerReference((left_temple_marker_index-1)*3+1 : (left_temple_marker_index-1)*3+3)';
    sellion_position            = mean([right_temple_position, left_temple_position], 2);
    BHD_position                = mean([markerReference((LBHD_marker_index-1)*3+1 : (LBHD_marker_index-1)*3+3)' markerReference((RBHD_marker_index-1)*3+1 : (RBHD_marker_index-1)*3+3)'], 2);
    head_vertex_position        = mean([sellion_position, BHD_position], 2) + [0; 0; 0.05];
    suprasternale_position      = markerReference((suprasternale_marker_index-1)*3+1 : (suprasternale_marker_index-1)*3+3)'         + suprasternale_marker_offset;
    c7_position                 = markerReference((c7_marker_index-1)*3+1 : (c7_marker_index-1)*3+3)'                               + c7_marker_offset;
    MPSIS_position              = MPSIS;
    right_ASIS_position         = markerReference((right_ASIS_marker_index-1)*3+1 : (right_ASIS_marker_index-1)*3+3)'               + right_ASIS_marker_offset;
    left_ASIS_position          = markerReference((left_ASIS_marker_index-1)*3+1 : (left_ASIS_marker_index-1)*3+3)'                 + left_ASIS_marker_offset;
    right_calcaneous_position   = markerReference((right_calcaneous_marker_index-1)*3+1 : (right_calcaneous_marker_index-1)*3+3)'   + right_calcaneous_marker_offset;
    right_toetip_position       = markerReference((right_toetip_marker_index-1)*3+1 : (right_toetip_marker_index-1)*3+3)'           + right_toetip_marker_offset;
    left_calcaneous_position    = markerReference((left_calcaneous_marker_index-1)*3+1 : (left_calcaneous_marker_index-1)*3+3)'     + left_calcaneous_marker_offset;
    left_toetip_position        = markerReference((left_toetip_marker_index-1)*3+1 : (left_toetip_marker_index-1)*3+3)'             + left_toetip_marker_offset;
    
    right_met_one_position      = right_toetip_position + 0.05*left_direction;
    right_met_five_position     = right_toetip_position + 0.05*right_direction;
    left_met_one_position      = left_toetip_position + 0.05*right_direction;
    left_met_five_position     = left_toetip_position + 0.05*left_direction;
    
    
    % define segment masses and correct for mean rounding errors
    pelvis_segment_mass = pelvis_mass_scaling_factor*bodyMass * 1.002^(-1);
    foot_segment_mass = foot_mass_scaling_factor*bodyMass * 1.002^(-1);
    leg_segment_mass = shank_mass_scaling_factor*bodyMass * 1.002^(-1);
    thigh_segment_mass = thigh_mass_scaling_factor*bodyMass * 1.002^(-1);
    torso_segment_mass = (torso_mass_scaling_factor + 2*arm_mass_scaling_factor + 2*forearm_mass_scaling_factor + 2*hand_mass_scaling_factor)*bodyMass * 1.002^(-1); % adding the arm mass to the torso here
    head_segment_mass = head_mass_scaling_factor*bodyMass * 1.002^(-1);
    
    % find the lumbar joint center
    pelvis_width = norm(left_ASIS_position - right_ASIS_position);
    pelvis_acs_origin = mean([right_ASIS_position left_ASIS_position], 2);
    
    lumbar_cor = pelvis_acs_origin + 0.264*posterior_direction*pelvis_width + 0.126*up_direction*pelvis_width;
 
    % find the cervical joint center
    c7_to_sup_direction = normVector(suprasternale_position - c7_position);
    c7_to_cervix_cor = expAxis(posterior_direction, 8*pi/180) * c7_to_sup_direction;
    thorax_width = norm(suprasternale_position - c7_position);
    cervix_cor = c7_position + 0.55*thorax_width*c7_to_cervix_cor;
 
    % pelvis
    pelvis_center = mean([MPSIS_position right_ASIS_position left_ASIS_position], 2);
    pelvis_scs_z = normVector(right_ASIS_position - left_ASIS_position);
    pelvis_scs_y = normVector(cross(right_ASIS_position - MPSIS_position, left_ASIS_position - MPSIS_position));
    pelvis_scs_x = cross(pelvis_scs_y, pelvis_scs_z);
    pelvis_segment_length = norm(mean([right_hip_cor left_hip_cor], 2) - lumbar_cor);
%     pelvis_segment_length = norm(mean([right_hip_cor left_hip_cor], 2) - lumbar_cor) - 0.03; % adjusted because inertia seems unrealistically large
    pelvis_com = lumbar_cor ...
                       + pelvis_com_scaling_factor_x * pelvis_segment_length * pelvis_scs_x ...
                       + pelvis_com_scaling_factor_y * pelvis_segment_length * pelvis_scs_y ...
                       + pelvis_com_scaling_factor_z * pelvis_segment_length * pelvis_scs_z;
    pelvis_I_xx = (pelvis_rxx_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
    pelvis_I_yy = (pelvis_ryy_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
    pelvis_I_zz = (pelvis_rzz_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
    pelvis_I_xy = (pelvis_rxy_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
    pelvis_I_xz = (pelvis_rxz_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
    pelvis_I_yz = (pelvis_ryz_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
    pelvis_inertia_tensor = [pelvis_I_xx pelvis_I_xy pelvis_I_xz; pelvis_I_xy pelvis_I_yy pelvis_I_yz; pelvis_I_xz pelvis_I_yz pelvis_I_zz];
        
    % right foot
    right_met_mean = mean([right_met_five_position right_met_one_position], 2);
    right_foot_scs_x = normVector(right_met_mean - right_calcaneous_position);
    right_foot_scs_y = normVector(cross(right_met_five_position-right_calcaneous_position, right_met_one_position-right_calcaneous_position));
    right_foot_scs_z = cross(right_foot_scs_x, right_foot_scs_y);
    right_foot_segment_length = norm(right_met_mean - right_ankle_cor);
    right_foot_com = right_ankle_cor ...
                       + foot_com_scaling_factor_x * right_foot_segment_length * right_foot_scs_x ...
                       + foot_com_scaling_factor_y * right_foot_segment_length * right_foot_scs_y ...
                       + foot_com_scaling_factor_z * right_foot_segment_length * right_foot_scs_z;
    right_foot_I_xx = (foot_rxx_scaling_factor*right_foot_segment_length)^2 * foot_segment_mass;
    right_foot_I_yy = (foot_ryy_scaling_factor*right_foot_segment_length)^2 * foot_segment_mass;
    right_foot_I_zz = (foot_rzz_scaling_factor*right_foot_segment_length)^2 * foot_segment_mass;
    right_foot_I_xy = (foot_rxy_scaling_factor*right_foot_segment_length)^2 * foot_segment_mass;
    right_foot_I_xz = (foot_rxz_scaling_factor*right_foot_segment_length)^2 * foot_segment_mass;
    right_foot_I_yz = (foot_ryz_scaling_factor*right_foot_segment_length)^2 * foot_segment_mass;
    right_foot_inertia_tensor = [right_foot_I_xx right_foot_I_xy right_foot_I_xz; right_foot_I_xy right_foot_I_yy right_foot_I_yz; right_foot_I_xz right_foot_I_yz right_foot_I_zz];
    
    % right leg
    right_knee_aor = jointAxes{1};
    if right_knee_aor(2) < 0
        right_knee_aor = - right_knee_aor;
    end
    right_leg_scs_y = normVector(right_knee_cor - right_ankle_cor);
    right_leg_scs_x = normVector(cross(right_knee_aor, right_leg_scs_y));
    right_leg_scs_z = cross(right_leg_scs_x, right_leg_scs_y);
    right_leg_segment_length = norm(right_knee_cor - right_ankle_cor);
    right_leg_com = right_knee_cor ...
                       + shank_com_scaling_factor_x * right_leg_segment_length * right_leg_scs_x ...
                       + shank_com_scaling_factor_y * right_leg_segment_length * right_leg_scs_y ...
                       + shank_com_scaling_factor_z * right_leg_segment_length * right_leg_scs_z;
    right_leg_I_xx = (shank_rxx_scaling_factor*right_leg_segment_length)^2 * leg_segment_mass;
    right_leg_I_yy = (shank_ryy_scaling_factor*right_leg_segment_length)^2 * leg_segment_mass;
    right_leg_I_zz = (shank_rzz_scaling_factor*right_leg_segment_length)^2 * leg_segment_mass;
    right_leg_I_xy = (shank_rxy_scaling_factor*right_leg_segment_length)^2 * leg_segment_mass;
    right_leg_I_xz = (shank_rxz_scaling_factor*right_leg_segment_length)^2 * leg_segment_mass;
    right_leg_I_yz = (shank_ryz_scaling_factor*right_leg_segment_length)^2 * leg_segment_mass;
    right_leg_inertia_tensor = [right_leg_I_xx right_leg_I_xy right_leg_I_xz; right_leg_I_xy right_leg_I_yy right_leg_I_yz; right_leg_I_xz right_leg_I_yz right_leg_I_zz];
    
     % right thigh
    right_thigh_scs_y = normVector(right_hip_cor - right_knee_cor);
    right_thigh_scs_x = normVector(cross(right_knee_aor, right_thigh_scs_y));
    right_thigh_scs_z = cross(right_thigh_scs_x, right_thigh_scs_y);
    right_thigh_segment_length = norm(right_hip_cor - right_knee_cor);
    right_thigh_com = right_hip_cor ...
                       + thigh_com_scaling_factor_x * right_thigh_segment_length * right_thigh_scs_x ...
                       + thigh_com_scaling_factor_y * right_thigh_segment_length * right_thigh_scs_y ...
                       + thigh_com_scaling_factor_z * right_thigh_segment_length * right_thigh_scs_z;
    right_thigh_I_xx = (thigh_rxx_scaling_factor*right_thigh_segment_length)^2 * thigh_segment_mass;
    right_thigh_I_yy = (thigh_ryy_scaling_factor*right_thigh_segment_length)^2 * thigh_segment_mass;
    right_thigh_I_zz = (thigh_rzz_scaling_factor*right_thigh_segment_length)^2 * thigh_segment_mass;
    right_thigh_I_xy = (thigh_rxy_scaling_factor*right_thigh_segment_length)^2 * thigh_segment_mass;
    right_thigh_I_xz = (thigh_rxz_scaling_factor*right_thigh_segment_length)^2 * thigh_segment_mass;
    right_thigh_I_yz = (thigh_ryz_scaling_factor*right_thigh_segment_length)^2 * thigh_segment_mass;
    right_thigh_inertia_tensor = [right_thigh_I_xx right_thigh_I_xy right_thigh_I_xz; right_thigh_I_xy right_thigh_I_yy right_thigh_I_yz; right_thigh_I_xz right_thigh_I_yz right_thigh_I_zz];
    
    % left foot
    left_met_mean = mean([left_met_five_position left_met_one_position], 2);
    left_foot_scs_x = normVector(left_met_mean - left_calcaneous_position);
    left_foot_scs_y = normVector(cross(left_met_one_position-left_calcaneous_position, left_met_five_position-left_calcaneous_position));
    left_foot_scs_z = cross(left_foot_scs_x, left_foot_scs_y);
    left_foot_segment_length = norm(left_met_mean - left_ankle_cor);
    left_foot_com = left_ankle_cor ...
                       + foot_com_scaling_factor_x * left_foot_segment_length * left_foot_scs_x ...
                       + foot_com_scaling_factor_y * left_foot_segment_length * left_foot_scs_y ...
                       - foot_com_scaling_factor_z * left_foot_segment_length * left_foot_scs_z;
    left_foot_I_xx = (foot_rxx_scaling_factor*left_foot_segment_length)^2 * foot_segment_mass;
    left_foot_I_yy = (foot_ryy_scaling_factor*left_foot_segment_length)^2 * foot_segment_mass;
    left_foot_I_zz = (foot_rzz_scaling_factor*left_foot_segment_length)^2 * foot_segment_mass;
    left_foot_I_xy = (foot_rxy_scaling_factor*left_foot_segment_length)^2 * foot_segment_mass;
    left_foot_I_xz = -(foot_rxz_scaling_factor*left_foot_segment_length)^2 * foot_segment_mass;
    left_foot_I_yz = -(foot_ryz_scaling_factor*left_foot_segment_length)^2 * foot_segment_mass;
    left_foot_inertia_tensor = [left_foot_I_xx left_foot_I_xy left_foot_I_xz; left_foot_I_xy left_foot_I_yy left_foot_I_yz; left_foot_I_xz left_foot_I_yz left_foot_I_zz];
    
    % left leg
    left_knee_aor = jointAxes{1};
    if left_knee_aor(2) < 0
        left_knee_aor = - left_knee_aor;
    end
    left_leg_scs_y = normVector(left_knee_cor - left_ankle_cor);
    left_leg_scs_x = normVector(cross(left_knee_aor, left_leg_scs_y));
    left_leg_scs_z = cross(left_leg_scs_x, left_leg_scs_y);
    left_leg_segment_length = norm(left_knee_cor - left_ankle_cor);
    left_leg_com = left_knee_cor ...
                       + shank_com_scaling_factor_x * left_leg_segment_length * left_leg_scs_x ...
                       + shank_com_scaling_factor_y * left_leg_segment_length * left_leg_scs_y ...
                       - shank_com_scaling_factor_z * left_leg_segment_length * left_leg_scs_z;
    left_leg_I_xx = (shank_rxx_scaling_factor*left_leg_segment_length)^2 * leg_segment_mass;
    left_leg_I_yy = (shank_ryy_scaling_factor*left_leg_segment_length)^2 * leg_segment_mass;
    left_leg_I_zz = (shank_rzz_scaling_factor*left_leg_segment_length)^2 * leg_segment_mass;
    left_leg_I_xy = (shank_rxy_scaling_factor*left_leg_segment_length)^2 * leg_segment_mass;
    left_leg_I_xz = -(shank_rxz_scaling_factor*left_leg_segment_length)^2 * leg_segment_mass;
    left_leg_I_yz = -(shank_ryz_scaling_factor*left_leg_segment_length)^2 * leg_segment_mass;
    left_leg_inertia_tensor = [left_leg_I_xx left_leg_I_xy left_leg_I_xz; left_leg_I_xy left_leg_I_yy left_leg_I_yz; left_leg_I_xz left_leg_I_yz left_leg_I_zz];
    
    % left thigh
    left_thigh_scs_y = normVector(left_hip_cor - left_knee_cor);
    left_thigh_scs_x = normVector(cross(left_knee_aor, left_thigh_scs_y));
    left_thigh_scs_z = cross(left_thigh_scs_x, left_thigh_scs_y);
    left_thigh_segment_length = norm(left_hip_cor - left_knee_cor);
    left_thigh_com = left_hip_cor ...
                       + thigh_com_scaling_factor_x * left_thigh_segment_length * left_thigh_scs_x ...
                       + thigh_com_scaling_factor_y * left_thigh_segment_length * left_thigh_scs_y ...
                       - thigh_com_scaling_factor_z * left_thigh_segment_length * left_thigh_scs_z;
    left_thigh_I_xx = (thigh_rxx_scaling_factor*left_thigh_segment_length)^2 * thigh_segment_mass;
    left_thigh_I_yy = (thigh_ryy_scaling_factor*left_thigh_segment_length)^2 * thigh_segment_mass;
    left_thigh_I_zz = (thigh_rzz_scaling_factor*left_thigh_segment_length)^2 * thigh_segment_mass;
    left_thigh_I_xy = (thigh_rxy_scaling_factor*left_thigh_segment_length)^2 * thigh_segment_mass;
    left_thigh_I_xz = -(thigh_rxz_scaling_factor*left_thigh_segment_length)^2 * thigh_segment_mass;
    left_thigh_I_yz = -(thigh_ryz_scaling_factor*left_thigh_segment_length)^2 * thigh_segment_mass;
    left_thigh_inertia_tensor = [left_thigh_I_xx left_thigh_I_xy left_thigh_I_xz; left_thigh_I_xy left_thigh_I_yy left_thigh_I_yz; left_thigh_I_xz left_thigh_I_yz left_thigh_I_zz];
   

    % torso
    torso_scs_y = normVector(cervix_cor - lumbar_cor);
    torso_scs_z = normVector(cross(cervix_cor - suprasternale_position, lumbar_cor - suprasternale_position));
    torso_scs_x = cross(torso_scs_y, torso_scs_z);
    torso_segment_length = norm(lumbar_cor - cervix_cor);
    torso_com = cervix_cor ...
                       + torso_com_scaling_factor_x * torso_segment_length * torso_scs_x ...
                       + torso_com_scaling_factor_y * torso_segment_length * torso_scs_y ...
                       + torso_com_scaling_factor_z * torso_segment_length * torso_scs_z;
    torso_I_xx = (torso_rxx_scaling_factor*torso_segment_length)^2 * torso_segment_mass;
    torso_I_yy = (torso_ryy_scaling_factor*torso_segment_length)^2 * torso_segment_mass;
    torso_I_zz = (torso_rzz_scaling_factor*torso_segment_length)^2 * torso_segment_mass;
    torso_I_xy = (torso_rxy_scaling_factor*torso_segment_length)^2 * torso_segment_mass;
    torso_I_xz = (torso_rxz_scaling_factor*torso_segment_length)^2 * torso_segment_mass;
    torso_I_yz = (torso_ryz_scaling_factor*torso_segment_length)^2 * torso_segment_mass;
    torso_inertia_tensor = [torso_I_xx torso_I_xy torso_I_xz; torso_I_xy torso_I_yy torso_I_yz; torso_I_xz torso_I_yz torso_I_zz];
    
    % head and neck
    head_scs_y = normVector(head_vertex_position - cervix_cor);
    head_scs_z = normVector(cross(sellion_position - cervix_cor, head_vertex_position - cervix_cor));
    head_scs_x = cross(head_scs_y, head_scs_z);
    head_segment_length = norm(cervix_cor - head_vertex_position);
%     head_segment_length = norm(cervix_cor - head_vertex) - 0.03; % adjusted because inertia seems unrealistically large
    head_com = cervix_cor ...
                       + head_com_scaling_factor_x * head_segment_length * head_scs_x ...
                       + head_com_scaling_factor_y * head_segment_length * head_scs_y ...
                       + head_com_scaling_factor_z * head_segment_length * head_scs_z;
    head_I_xx = (head_rxx_scaling_factor*head_segment_length)^2 * head_segment_mass;
    head_I_yy = (head_ryy_scaling_factor*head_segment_length)^2 * head_segment_mass;
    head_I_zz = (head_rzz_scaling_factor*head_segment_length)^2 * head_segment_mass;
    head_I_xy = (head_rxy_scaling_factor*head_segment_length)^2 * head_segment_mass;
    head_I_xz = (head_rxz_scaling_factor*head_segment_length)^2 * head_segment_mass;
    head_I_yz = (head_ryz_scaling_factor*head_segment_length)^2 * head_segment_mass;
    head_inertia_tensor = [head_I_xx head_I_xy head_I_xz; head_I_xy head_I_yy head_I_yz; head_I_xz head_I_yz head_I_zz];
    
    
    % check that knee axes are pointing in the right direction
    right_knee_aor = jointAxes{1};
    left_knee_aor = jointAxes{2};
    e_1 = [1; 0; 0];
    e_2 = [0; 1; 0];
    e_3 = [0; 0; 1];
%     left_direction_vector = left_direction * e_1;
%     right_direction_vector = - left_direction_vector;
%     anterior_direction_vector = anterior_direction * e_2;
%     posterior_direction_vector = - anterior_direction_vector;
    if dot(right_knee_aor, left_direction) < 0
        right_knee_aor = -right_knee_aor;
    end
    if dot(left_knee_aor, left_direction) < 0
        left_knee_aor = -left_knee_aor;
    end
    
    
    
    % adjust end-effector positions for foot extension in z-direction
    right_heel = right_ankle_cor; right_heel(3) = right_toetip_position(3);
    left_heel = left_ankle_cor; left_heel(3) = left_toetip_position(3);
    
end

joint_positions = ...
{ ...
    pelvis_com, pelvis_com, pelvis_com, pelvis_com, pelvis_com, pelvis_com, ...
    left_hip_cor, left_hip_cor, left_hip_cor, left_knee_cor, left_ankle_cor, left_ankle_cor, ...
    right_hip_cor, right_hip_cor, right_hip_cor, right_knee_cor, right_ankle_cor, right_ankle_cor, ...
    lumbar_cor, lumbar_cor, lumbar_cor, cervix_cor, cervix_cor, cervix_cor ...
};

joint_axes = ...
{ ...
    e_1, e_2, e_3, e_1, e_2, e_3, ...       % pelvis free body DoFs
    right_direction, anterior_direction, up_direction, left_knee_aor, left_direction, posterior_direction, ...   % left leg
    right_direction, posterior_direction, down_direction, right_knee_aor, left_direction, anterior_direction,...   % right leg
    left_direction, posterior_direction, up_direction, left_direction, posterior_direction, up_direction, ...       % trunk and neck
};

joint_types = ...
  [ ...
    2 2 2 1 1 1 ... % virtual dofs
    1 1 1 1 1 1 ... % left leg
    1 1 1 1 1 1 ... % right leg
    1 1 1 1 1 1 ... % l5 and neck
  ]; % 1 = rotational, 2 = prismatic

% link setup
link_com_positions =  ...
{ ...
    pelvis_com, pelvis_com, pelvis_com, pelvis_com, pelvis_com, pelvis_com, ...
    left_thigh_com, left_thigh_com, left_thigh_com, left_leg_com, left_foot_com, left_foot_com, ...
    right_thigh_com, right_thigh_com, right_thigh_com, right_leg_com, right_foot_com, right_foot_com, ...
    torso_com, torso_com, torso_com, head_com, head_com, head_com ...
};

link_orientations = ...
  {
    eye(3); eye(3); eye(3); eye(3); eye(3); [pelvis_scs_x pelvis_scs_y pelvis_scs_z]; ...
    eye(3); eye(3); [left_thigh_scs_x left_thigh_scs_y left_thigh_scs_z]; [left_leg_scs_x left_leg_scs_y left_leg_scs_z]; eye(3); [left_foot_scs_x left_foot_scs_y left_foot_scs_z]; ...
    eye(3); eye(3); [right_thigh_scs_x right_thigh_scs_y right_thigh_scs_z]; [right_leg_scs_x right_leg_scs_y right_leg_scs_z]; eye(3); [right_foot_scs_x right_foot_scs_y right_foot_scs_z]; ...
    eye(3); eye(3); [torso_scs_x torso_scs_y torso_scs_z]; eye(3); eye(3); [head_scs_x head_scs_y head_scs_z]; ...
  };

generalized_inertia_matrix_pelvis = ...
  [ ...
    pelvis_segment_mass*eye(3) zeros(3); ...
    zeros(3) pelvis_inertia_tensor ...
  ];
generalized_inertia_matrix_left_thigh = ...
  [ ...
    thigh_segment_mass*eye(3) zeros(3); ...
    zeros(3) left_thigh_inertia_tensor ...
  ];
generalized_inertia_matrix_left_shank = ...
  [ ...
    leg_segment_mass*eye(3) zeros(3); ...
    zeros(3) left_leg_inertia_tensor ...
  ];
generalized_inertia_matrix_left_foot = ...
  [ ...
    foot_segment_mass*eye(3) zeros(3); ...
    zeros(3) left_foot_inertia_tensor ...
  ];
generalized_inertia_matrix_right_thigh = ...
  [ ...
    thigh_segment_mass*eye(3) zeros(3); ...
    zeros(3) right_thigh_inertia_tensor ...
  ];
generalized_inertia_matrix_right_shank = ...
  [ ...
    leg_segment_mass*eye(3) zeros(3); ...
    zeros(3) right_leg_inertia_tensor ...
  ];
generalized_inertia_matrix_right_foot = ...
  [ ...
    foot_segment_mass*eye(3) zeros(3); ...
    zeros(3) right_foot_inertia_tensor ...
  ];
generalized_inertia_matrix_torso = ...
  [ ...
    torso_segment_mass*eye(3) zeros(3); ...
    zeros(3) torso_inertia_tensor ...
  ];
generalized_inertia_matrix_head = ...
  [ ...
    head_segment_mass*eye(3) zeros(3); ...
    zeros(3) head_inertia_tensor ...
  ];

generalized_link_inertia_matrices = ...
  { ...
    zeros(6, 6), zeros(6, 6), zeros(6, 6), zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_pelvis, ...
    zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_left_thigh, generalized_inertia_matrix_left_shank, zeros(6, 6), generalized_inertia_matrix_left_foot, ...
    zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_right_thigh, generalized_inertia_matrix_right_shank, zeros(6, 6), generalized_inertia_matrix_right_foot, ...
    zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_torso, zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_head ...
  };
    
end_effector_transformations = ...
  { ...
    [eye(3), left_heel; 0 0 0 1], ...
    [eye(3), left_toetip_position; 0 0 0 1], ...
    [left_foot_scs_x left_foot_scs_y left_foot_scs_z left_ankle_cor; 0 0 0 1], ...
    [eye(3), right_heel; 0 0 0 1], ...
    [eye(3), right_toetip_position; 0 0 0 1], ...
    [right_foot_scs_x right_foot_scs_y right_foot_scs_z right_ankle_cor; 0 0 0 1], ...
    [eye(3), head_vertex_position; 0 0 0 1], ...
  };

branch_matrix = ...
  [ ...
    1 1 1 1 1 1    1 1 1 1 1 1   0 0 0 0 0 0   0 0 0 0 0 0; ... % left leg until heel
    1 1 1 1 1 1    1 1 1 1 1 1   0 0 0 0 0 0   0 0 0 0 0 0; ... % left leg until toes
    1 1 1 1 1 1    1 1 1 1 1 1   0 0 0 0 0 0   0 0 0 0 0 0; ... % left leg until ankle
    1 1 1 1 1 1    0 0 0 0 0 0   1 1 1 1 1 1   0 0 0 0 0 0; ... % right leg until heel
    1 1 1 1 1 1    0 0 0 0 0 0   1 1 1 1 1 1   0 0 0 0 0 0; ... % right leg until toes
    1 1 1 1 1 1    0 0 0 0 0 0   1 1 1 1 1 1   0 0 0 0 0 0; ... % right leg until ankle
    1 1 1 1 1 1    0 0 0 0 0 0   0 0 0 0 0 0   1 1 1 1 1 1; ... % head
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
    'cervical joint - forward/backward bending', 'cervical joint - sideways bending (left/right)', 'cervical joint - internal rotation (left/right)' ...
  };


% markers
red = [1 0 0];
green = [0 1 0];
blue = [0 0 1];
yellow = [1 0.9 0];
if nargin < 5
    
    kinematic_tree.addMarker(24, head_center_position + [0.1; -0.08; 0.0], red);    % head right
    kinematic_tree.addMarker(24, head_center_position + [0.1; 0.08; 0.0], red);    % head left
    
    kinematic_tree.addMarker(21, trunk_center_position + [0.1; 0; 0.1], green);    % manubrium
    kinematic_tree.addMarker(21, trunk_center_position + [-0.1; 0; 0.1], green);    % C7
    kinematic_tree.addMarker(21, trunk_center_position + [0; -0.2; 0.15], green);    % right shoulder
    kinematic_tree.addMarker(21, trunk_center_position + [0; 0.2; 0.15], green);    % left shoulder
    
    kinematic_tree.addMarker(6, lumbar_cor + [0.1; 0; 0.0], blue);    % pubis
    kinematic_tree.addMarker(6, right_hip_cor + [0.05; -0.02; 0.05], blue); % right pelvis
    kinematic_tree.addMarker(6, left_hip_cor + [0.05; 0.02; 0.05], blue); % left pelvis

    kinematic_tree.addMarker(9, (right_knee_cor+right_hip_cor)*0.5 + [0.01; 0.05; 0.05], red);    % right thigh
%     kinematic_tree.addMarker(9, (right_knee_cor+right_hip_cor)*0.5 + [0.01;0.05; -0.05], green);    % right thigh

%     kinematic_tree.addMarker(10, (right_knee_cor+right_ankle_cor)*0.5 + [0.01; 0.05; 0.05], red);    % right shank
    kinematic_tree.addMarker(10, (right_knee_cor+right_ankle_cor)*0.5 + [0.01;0.05; -0.05], green);    % right shank

    kinematic_tree.addMarker(12, (right_met_mean+right_ankle_cor)*0.5 + [0.05; 0.05; 0.01], blue);    % right foot
    kinematic_tree.addMarker(12, (right_met_mean+right_ankle_cor)*0.5 + [-0.05; 0.05; 0.01], yellow);    % right foot
    
    kinematic_tree.addMarker(15, (left_knee_cor+left_hip_cor)*0.5 + [0.02; 0.05; 0.05], red);      % left thigh
%     kinematic_tree.addMarker(15, (left_knee_cor+left_hip_cor)*0.5 + [0.02; 0.05; -0.05], green);      % left thigh
    
%     kinematic_tree.addMarker(16, (left_knee_cor+left_ankle_cor)*0.5 + [0.02; 0.05; 0.05], red);    % left shank
    kinematic_tree.addMarker(16, (left_knee_cor+left_ankle_cor)*0.5 + [0.02; 0.05; -0.05], green);    % left shank
    
    kinematic_tree.addMarker(18, (left_met_mean+left_ankle_cor)*0.5 + [0.05; 0.05; 0.03], blue);    % left foot
    kinematic_tree.addMarker(18, (left_met_mean+left_ankle_cor)*0.5 + [-0.05; 0.05; 0.03], yellow);    % left foot
    
    
else
    markerSegments = ...
      [ ...
        24 24 24 24 ...             % head
        21 21 21 21 21 ...          % trunk
        21 21 21 21 21 21 21 ...    % left arm (fixed to trunk for now)
        21 21 21 21 21 21 21 ...    % right arm (fixed to trunk for now)
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
    

    
    
    
