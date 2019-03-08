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

% this function creates a model of the subject

% input: 
% calibration files
% anthropometric data in subjectInfo.mat
% parameters in studySettings.txt
% - static reference type: ski, casual, motorcycle, anatomic (currently not used, will implement later)
% - hip_joint_center_estimation_method: Tylkowski, SCoRE
% - knee_joint_axis_estimation_method: markers, SARA
%
% output:
% subjectModel.mat

function stanceModel_4DoF(varargin)
    load('subjectInfo.mat');
    
    % parse arguments
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;
    
    subject_settings = SettingsCustodian('subjectSettings.txt');
    

    load('subjectInfo.mat', 'date', 'subject_id');
    weight = subject_settings.get('weight');

    %% create static reference

    % load static reference file
    load(['processed' filesep makeFileName(date, subject_id, subject_settings.get('static_reference_trial_type'), subject_settings.get('static_reference_trial_number'), 'markerTrajectories')]);

    % find first time step where all markers are available
    i_time = 1;
    while any(isnan(marker_trajectories(i_time, :)))
        i_time = i_time + 1;
    end
    marker_reference = marker_trajectories(i_time, :);

    %% extract marker reference positions
    
    % head
    eye_L_reference = extractMarkerData(marker_reference, marker_labels, 'LEYE')';
    eye_R_reference = extractMarkerData(marker_reference, marker_labels, 'REYE')';
    ear_L_reference = extractMarkerData(marker_reference, marker_labels, 'LMAS')';
    ear_R_reference = extractMarkerData(marker_reference, marker_labels, 'RMAS')';

    % left arm
    shoulder_L_reference = extractMarkerData(marker_reference, marker_labels, 'LACR')';
    shoulder_R_reference = extractMarkerData(marker_reference, marker_labels, 'RACR')';

    % pelvis
    lumbar_joint_reference = extractMarkerData(marker_reference, marker_labels, 'L5S1')';

    % left leg
    hip_L_reference = extractMarkerData(marker_reference, marker_labels, 'LGTR')';
    knee_L_reference = extractMarkerData(marker_reference, marker_labels, 'LCND')';
    ankle_L_reference = extractMarkerData(marker_reference, marker_labels, 'LMAL')';
    toe_L_reference = extractMarkerData(marker_reference, marker_labels, 'L5MT')';

    % right leg
    hip_R_reference = extractMarkerData(marker_reference, marker_labels, 'RGTR')';
    knee_R_reference = extractMarkerData(marker_reference, marker_labels, 'RCND')';
    ankle_R_reference = extractMarkerData(marker_reference, marker_labels, 'RMAL')';
    toe_R_reference = extractMarkerData(marker_reference, marker_labels, 'R5MT')';

%     inter_ASIS_distance = norm(LASI_reference - RASI_reference);
    
    % define directions
    e_1 = [1; 0; 0];
    e_2 = [0; 1; 0];
    e_3 = [0; 0; 1];
    right_direction = e_2;
    anterior_direction = -e_1;
    proximal_direction = e_3;

    left_direction = - right_direction;
    posterior_direction = - anterior_direction;
    distal_direction = - proximal_direction;
    up_direction = proximal_direction;
    down_direction = distal_direction;

    %% estimate joint centers
    ankle_cor = mean([ankle_L_reference ankle_R_reference], 2);
    knee_cor = mean([knee_L_reference knee_R_reference], 2);
    hip_cor = mean([hip_L_reference hip_R_reference], 2);
    
    toe_mid = mean([toe_L_reference toe_R_reference], 2);
    shoulders_mid = mean([shoulder_L_reference shoulder_R_reference], 2);
    ears_mid = mean([ear_L_reference ear_R_reference], 2);

    lumbar_cor = lumbar_joint_reference;
    neck_cor = shoulders_mid;

    ankle_dorsiflexion_axis = left_direction;
    knee_flexion_axis = left_direction;
    hip_flexion_axis = left_direction;
    neck_flexion_axis = left_direction;

    %% define scaling factors
    % according to R. Dumas , L. Cheze, J.-P. Verriest: "Adjustments to McConville et al. and Young et al. body
    % segment inertial parameters", Journal of Biomechanics 40 (2007) 543?553
    % mass
    head_mass_scaling_factor      = 0.067;
    trunk_mass_scaling_factor     = 0.333;
    arm_mass_scaling_factor       = 0.024;
    forearm_mass_scaling_factor   = 0.017;
    hand_mass_scaling_factor      = 0.006;
    pelvis_mass_scaling_factor    = 0.142;
    thigh_mass_scaling_factor     = 0.123;
    shank_mass_scaling_factor     = 0.048;
    foot_mass_scaling_factor      = 0.012;

    % CoM
    head_com_scaling_factor_x    = -0.062;  head_com_scaling_factor_y    =  0.555;  head_com_scaling_factor_z    =  0.001;
    trunk_com_scaling_factor_x   = -0.036;  trunk_com_scaling_factor_y   = -0.420;  trunk_com_scaling_factor_z   = -0.002;
    arm_com_scaling_factor_x     =  0.017;  arm_com_scaling_factor_y     = -0.452;  arm_com_scaling_factor_z     = -0.026;
    forearm_com_scaling_factor_x =  0.010;  forearm_com_scaling_factor_y = -0.417;  forearm_com_scaling_factor_z =  0.014;
    hand_com_scaling_factor_x    =  0.082;  hand_com_scaling_factor_y    = -0.839;  hand_com_scaling_factor_z    =  0.074;
    pelvis_com_scaling_factor_x  =  0.028;  pelvis_com_scaling_factor_y  = -0.280;  pelvis_com_scaling_factor_z  = -0.006;
    thigh_com_scaling_factor_x   = -0.041;  thigh_com_scaling_factor_y   = -0.429;  thigh_com_scaling_factor_z   =  0.033;
    shank_com_scaling_factor_x   = -0.048;  shank_com_scaling_factor_y   = -0.410;  shank_com_scaling_factor_z   =  0.007;
    foot_com_scaling_factor_x    =  0.382;  foot_com_scaling_factor_y    = -0.151;  foot_com_scaling_factor_z    =  0.026;

    % rog
    head_rxx_scaling_factor    = 0.31;  head_ryy_scaling_factor    = 0.25;  head_rzz_scaling_factor    = 0.33;  head_rxy_scaling_factor    = 0.09*1i;	head_rxz_scaling_factor    = 0.02*1i;   head_ryz_scaling_factor    = 0.03;
    trunk_rxx_scaling_factor   = 0.27;  trunk_ryy_scaling_factor   = 0.25;  trunk_rzz_scaling_factor   = 0.28;  trunk_rxy_scaling_factor   = 0.18;      trunk_rxz_scaling_factor   = 0.02;      trunk_ryz_scaling_factor   = 0.04*1i;
    arm_rxx_scaling_factor     = 0.31;  arm_ryy_scaling_factor     = 0.14;  arm_rzz_scaling_factor     = 0.32;  arm_rxy_scaling_factor     = 0.06;      arm_rxz_scaling_factor     = 0.05;      arm_ryz_scaling_factor     = 0.02;
    forearm_rxx_scaling_factor = 0.28;  forearm_ryy_scaling_factor = 0.11;  forearm_rzz_scaling_factor = 0.27;  forearm_rxy_scaling_factor = 0.03;      forearm_rxz_scaling_factor = 0.02;      forearm_ryz_scaling_factor = 0.08*1i;
    hand_rxx_scaling_factor    = 0.61;  hand_ryy_scaling_factor    = 0.38;  hand_rzz_scaling_factor    = 0.56;  hand_rxy_scaling_factor    = 0.22;      hand_rxz_scaling_factor    = 0.15;      hand_ryz_scaling_factor    = 0.20*1i;
    pelvis_rxx_scaling_factor  = 1.01;  pelvis_ryy_scaling_factor  = 1.06;  pelvis_rzz_scaling_factor  = 0.95;  pelvis_rxy_scaling_factor  = 0.25*1i;   pelvis_rxz_scaling_factor  = 0.12*1i;   pelvis_ryz_scaling_factor  = 0.08*1i;
    thigh_rxx_scaling_factor   = 0.29;  thigh_ryy_scaling_factor   = 0.15;  thigh_rzz_scaling_factor   = 0.30;  thigh_rxy_scaling_factor   = 0.07;      thigh_rxz_scaling_factor   = 0.02*1i;   thigh_ryz_scaling_factor   = 0.07*1i;
    shank_rxx_scaling_factor   = 0.28;  shank_ryy_scaling_factor   = 0.10;  shank_rzz_scaling_factor   = 0.28;  shank_rxy_scaling_factor   = 0.04*1i;   shank_rxz_scaling_factor   = 0.02*1i;   shank_ryz_scaling_factor   = 0.05;
    foot_rxx_scaling_factor    = 0.17;  foot_ryy_scaling_factor    = 0.37;	foot_rzz_scaling_factor    = 0.36;  foot_rxy_scaling_factor    = 0.13;      foot_rxz_scaling_factor    = 0.08*1i;	foot_ryz_scaling_factor    = 0.00;


    %% set up the segment coordinate systems (SCS)

    % pelvis
    pelvis_scs_z = right_direction;
    pelvis_scs_y = up_direction;
    pelvis_scs_x = cross(pelvis_scs_y, pelvis_scs_z);

    % thighs
    thigh_scs_y = normVector(hip_cor - knee_cor);
    thigh_scs_x = normVector(cross(knee_flexion_axis, thigh_scs_y));
    thigh_scs_z = cross(thigh_scs_x, thigh_scs_y);

    % left leg
    leg_scs_y = normVector(knee_cor - ankle_cor);
    leg_scs_x = normVector(cross(knee_flexion_axis, leg_scs_y));
    leg_scs_z = cross(leg_scs_x, leg_scs_y);

    % left foot
    foot_scs_x = normVector(toe_mid - ankle_cor);
    foot_scs_y = normVector(cross(toe_mid-ankle_cor, left_direction));
    foot_scs_z = cross(foot_scs_x, foot_scs_y);

    % arms and trunk
    trunk_scs_y = normVector(shoulders_mid - lumbar_cor);
    trunk_scs_z = right_direction;
    trunk_scs_x = cross(trunk_scs_y, trunk_scs_z);

    % head
    head_scs_y = normVector(ears_mid - shoulders_mid);
    head_scs_z = right_direction;
    head_scs_x = cross(head_scs_y, head_scs_z);

    %% calculate segment mass and CoM

    % calculate segment masses and correct for mean rounding errors
    pelvis_segment_mass = pelvis_mass_scaling_factor    * weight;
    thigh_segment_mass = thigh_mass_scaling_factor      * weight * 2;
    leg_segment_mass = shank_mass_scaling_factor        * weight * 2;
    foot_segment_mass = foot_mass_scaling_factor        * weight * 2;
    head_segment_mass = head_mass_scaling_factor        * weight;
    trunk_segment_mass_single = trunk_mass_scaling_factor      * weight;
    arm_segment_mass = arm_mass_scaling_factor          * weight;
    forearm_segment_mass = forearm_mass_scaling_factor  * weight;
    hand_segment_mass = hand_mass_scaling_factor        * weight;
    trunk_segment_mass = head_segment_mass + trunk_segment_mass_single + 2*arm_segment_mass + 2*forearm_segment_mass + 2*hand_segment_mass + pelvis_segment_mass;


    % pelvis
    pelvis_segment_length = norm(hip_cor - lumbar_cor);
    pelvis_com = lumbar_cor ...
                       + pelvis_com_scaling_factor_x * pelvis_segment_length * pelvis_scs_x ...
                       + pelvis_com_scaling_factor_y * pelvis_segment_length * pelvis_scs_y ...
                       + pelvis_com_scaling_factor_z * pelvis_segment_length * pelvis_scs_z;

    % thigh
    thigh_segment_length = norm(hip_cor - knee_cor);
    thigh_com = hip_cor ...
                       + thigh_com_scaling_factor_x * thigh_segment_length * thigh_scs_x ...
                       + thigh_com_scaling_factor_y * thigh_segment_length * thigh_scs_y ...
                       - 0 * thigh_segment_length * thigh_scs_z;

    % leg
    leg_segment_length = norm(knee_cor - ankle_cor);
    leg_com = knee_cor ...
                       + shank_com_scaling_factor_x * leg_segment_length * leg_scs_x ...
                       + shank_com_scaling_factor_y * leg_segment_length * leg_scs_y ...
                       - 0 * leg_segment_length * leg_scs_z;

    % foot
    foot_segment_length = norm(toe_mid - ankle_cor);
    foot_com = ankle_cor ...
                       + foot_com_scaling_factor_x * foot_segment_length * foot_scs_x ...
                       + foot_com_scaling_factor_y * foot_segment_length * foot_scs_y ...
                       - 0 * foot_segment_length * foot_scs_z;

    % trunk
    trunk_segment_length = norm(lumbar_cor - shoulders_mid);
    trunk_com = shoulders_mid ...
                       + trunk_com_scaling_factor_x * trunk_segment_length * trunk_scs_x ...
                       + trunk_com_scaling_factor_y * trunk_segment_length * trunk_scs_y ...
                       + trunk_com_scaling_factor_z * trunk_segment_length * trunk_scs_z;
    % head
    head_segment_length = norm(shoulders_mid - ears_mid);
    head_com = shoulders_mid ...
                       + head_com_scaling_factor_x * head_segment_length * head_scs_x ...
                       + head_com_scaling_factor_y * head_segment_length * head_scs_y ...
                       + head_com_scaling_factor_z * head_segment_length * head_scs_z;


    %% calculate inertia tensors

    % pelvis
    pelvis_I_xx = (pelvis_rxx_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
    pelvis_I_yy = (pelvis_ryy_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
    pelvis_I_zz = (pelvis_rzz_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
    pelvis_I_xy = (pelvis_rxy_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
    pelvis_I_xz = (pelvis_rxz_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
    pelvis_I_yz = (pelvis_ryz_scaling_factor*pelvis_segment_length)^2 * pelvis_segment_mass;
    pelvis_inertia_tensor = [pelvis_I_xx pelvis_I_xy pelvis_I_xz; pelvis_I_xy pelvis_I_yy pelvis_I_yz; pelvis_I_xz pelvis_I_yz pelvis_I_zz];


    % thigh
    thigh_I_xx = (thigh_rxx_scaling_factor*thigh_segment_length)^2 * thigh_segment_mass;
    thigh_I_yy = (thigh_ryy_scaling_factor*thigh_segment_length)^2 * thigh_segment_mass;
    thigh_I_zz = (thigh_rzz_scaling_factor*thigh_segment_length)^2 * thigh_segment_mass;
%     thigh_I_xy = (thigh_rxy_scaling_factor*thigh_segment_length)^2 * thigh_segment_mass;
%     thigh_I_xz = -(thigh_rxz_scaling_factor*thigh_segment_length)^2 * thigh_segment_mass;
%     thigh_I_yz = -(thigh_ryz_scaling_factor*thigh_segment_length)^2 * thigh_segment_mass;
    thigh_I_xy = 0;
    thigh_I_xz = 0;
    thigh_I_yz = 0;
    thigh_inertia_tensor = [thigh_I_xx thigh_I_xy thigh_I_xz; thigh_I_xy thigh_I_yy thigh_I_yz; thigh_I_xz thigh_I_yz thigh_I_zz];

    % leg
    leg_I_xx = (shank_rxx_scaling_factor*leg_segment_length)^2 * leg_segment_mass;
    leg_I_yy = (shank_ryy_scaling_factor*leg_segment_length)^2 * leg_segment_mass;
    leg_I_zz = (shank_rzz_scaling_factor*leg_segment_length)^2 * leg_segment_mass;
%     leg_I_xy = (shank_rxy_scaling_factor*leg_segment_length)^2 * leg_segment_mass;
%     leg_I_xz = -(shank_rxz_scaling_factor*leg_segment_length)^2 * leg_segment_mass;
%     leg_I_yz = -(shank_ryz_scaling_factor*leg_segment_length)^2 * leg_segment_mass;
    leg_I_xy = 0;
    leg_I_xz = 0;
    leg_I_yz = 0;
    leg_inertia_tensor = [leg_I_xx leg_I_xy leg_I_xz; leg_I_xy leg_I_yy leg_I_yz; leg_I_xz leg_I_yz leg_I_zz];

    % foot
    foot_I_xx = (foot_rxx_scaling_factor*foot_segment_length)^2 * foot_segment_mass;
    foot_I_yy = (foot_ryy_scaling_factor*foot_segment_length)^2 * foot_segment_mass;
    foot_I_zz = (foot_rzz_scaling_factor*foot_segment_length)^2 * foot_segment_mass;
    foot_I_xy = (foot_rxy_scaling_factor*foot_segment_length)^2 * foot_segment_mass;
    foot_I_xz = -(foot_rxz_scaling_factor*foot_segment_length)^2 * foot_segment_mass;
    foot_I_yz = -(foot_ryz_scaling_factor*foot_segment_length)^2 * foot_segment_mass;
    foot_I_xy = 0;
    foot_I_xz = 0;
    foot_I_yz = 0;
    foot_inertia_tensor = [foot_I_xx foot_I_xy foot_I_xz; foot_I_xy foot_I_yy foot_I_yz; foot_I_xz foot_I_yz foot_I_zz];

    % trunk
    trunk_I_xx = (trunk_rxx_scaling_factor*trunk_segment_length)^2 * trunk_segment_mass;
    trunk_I_yy = (trunk_ryy_scaling_factor*trunk_segment_length)^2 * trunk_segment_mass;
    trunk_I_zz = (trunk_rzz_scaling_factor*trunk_segment_length)^2 * trunk_segment_mass;
    trunk_I_xy = (trunk_rxy_scaling_factor*trunk_segment_length)^2 * trunk_segment_mass;
    trunk_I_xz = (trunk_rxz_scaling_factor*trunk_segment_length)^2 * trunk_segment_mass;
    trunk_I_yz = (trunk_ryz_scaling_factor*trunk_segment_length)^2 * trunk_segment_mass;
    trunk_I_xy = 0;
    trunk_I_xz = 0;
    trunk_I_yz = 0;
    trunk_inertia_tensor = [trunk_I_xx trunk_I_xy trunk_I_xz; trunk_I_xy trunk_I_yy trunk_I_yz; trunk_I_xz trunk_I_yz trunk_I_zz];

    head_I_xx = (head_rxx_scaling_factor*head_segment_length)^2 * head_segment_mass;
    head_I_yy = (head_ryy_scaling_factor*head_segment_length)^2 * head_segment_mass;
    head_I_zz = (head_rzz_scaling_factor*head_segment_length)^2 * head_segment_mass;
    head_I_xy = (head_rxy_scaling_factor*head_segment_length)^2 * head_segment_mass;
    head_I_xz = -(head_rxz_scaling_factor*head_segment_length)^2 * head_segment_mass;
    head_I_yz = -(head_ryz_scaling_factor*head_segment_length)^2 * head_segment_mass;
    head_I_xy = 0;
    head_I_xz = 0;
    head_I_yz = 0;
    head_inertia_tensor = [head_I_xx head_I_xy head_I_xz; head_I_xy head_I_yy head_I_yz; head_I_xz head_I_yz head_I_zz];
    
    %% set up marker coordinate systems (MCS)

    % define joint center references
    joint_center_reference = ...
      [ ...
        ankle_cor', ...
        knee_cor', ...
        hip_cor', ...
        neck_cor' ...
      ];

    joint_center_labels_single = ...
      { ...
        'ANKLECOR', ...
        'KNEECOR', ...
        'HIPCOR', ...
        'NECKCOR' ...
      };

    joint_center_labels = cell(3, length(joint_center_labels_single));
    for i_joint = 1 : length(joint_center_labels)
        joint_center_labels{1, i_joint} = [joint_center_labels_single{i_joint} '_x'];
        joint_center_labels{2, i_joint} = [joint_center_labels_single{i_joint} '_y'];
        joint_center_labels{3, i_joint} = [joint_center_labels_single{i_joint} '_z'];
    end
    joint_center_labels = reshape(joint_center_labels, 1, length(joint_center_labels_single)*3);
    
    number_of_joint_centers = length(joint_center_labels);
    joint_center_directions = cell(2, number_of_joint_centers);
    [joint_center_directions{1, 1 : 3 : number_of_joint_centers}] = deal(marker_directions{1, 1});
    [joint_center_directions{2, 1 : 3 : number_of_joint_centers}] = deal(marker_directions{2, 1});
    [joint_center_directions{1, 2 : 3 : number_of_joint_centers}] = deal(marker_directions{1, 2});
    [joint_center_directions{2, 2 : 3 : number_of_joint_centers}] = deal(marker_directions{2, 2});
    [joint_center_directions{1, 3 : 3 : number_of_joint_centers}] = deal(marker_directions{1, 3});
    [joint_center_directions{2, 3 : 3 : number_of_joint_centers}] = deal(marker_directions{2, 3});
    
    % specify markers that define segments
    segment_labels = ...
      { ...
        'LEGS', ...
        'THIGHS', ...
        'TRUNK', ...
        'HEAD' ...
      };

    number_of_segments = length(segment_labels);

    % assemble segment masses
    segment_masses = ...
      [ ...
        leg_segment_mass, ...
        thigh_segment_mass, ...
        trunk_segment_mass, ...
        head_segment_mass ...
      ];


    % calculate segment CoMs in marker coordinates
    segment_coms_wcs = ...
      { ...
        thigh_com; ...
        leg_com; ...
        trunk_com; ...
        head_com; ...
      };
    segment_coms_mcs = cell(size(segment_coms_wcs));

    %% assemble kinematic tree
    joint_positions = ...
    { ...
        ankle_cor, ...
        knee_cor, ...
        hip_cor, ...
        neck_cor ...
    };

    joint_axes = ...
    { ...
        ankle_dorsiflexion_axis, ...
        knee_flexion_axis, ...
        hip_flexion_axis, ...
        neck_flexion_axis ...
    };

    joint_types = [1 1 1 1];

    % link setup
    link_com_positions =  ...
    { ...
        leg_com, ...
        thigh_com, ...
        trunk_com, ...
        head_com ...
    };

    link_orientations = ...
    {
        [leg_scs_x leg_scs_y leg_scs_z]; ...
        [thigh_scs_x thigh_scs_y thigh_scs_z]; ...
        [trunk_scs_x trunk_scs_y trunk_scs_z]; ...
        [head_scs_x head_scs_y head_scs_z]; ...
    };

    generalized_inertia_matrix_shank        = [leg_segment_mass*eye(3) zeros(3); zeros(3) leg_inertia_tensor];
    generalized_inertia_matrix_thigh        = [thigh_segment_mass*eye(3) zeros(3); zeros(3) thigh_inertia_tensor];
    generalized_inertia_matrix_trunk        = [trunk_segment_mass*eye(3) zeros(3); zeros(3) trunk_inertia_tensor];
    generalized_inertia_matrix_head         = [head_segment_mass*eye(3) zeros(3); zeros(3) head_inertia_tensor];

    generalized_link_inertia_matrices = ...
      { ...
        generalized_inertia_matrix_shank, ...
        generalized_inertia_matrix_thigh, ...
        generalized_inertia_matrix_trunk, ...
        generalized_inertia_matrix_head ...
      };

    end_effector_transformations = ...
      { ...
        [eye(3) ears_mid; 0 0 0 1], ...
      };

    branch_matrix = [1 1 1 1]; % each row is a branch, listing the joints that move the end-effector of that branch

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

    % update
    kinematic_tree.updateInternals();

    
    
    joint_labels = ...
      { ...
        'ankle_joint', ...
        'knee_joint', ...
        'hip_joint', ...
        'neck_joint', ...
      };
    kinematic_tree.jointLabels = joint_labels;
  
    kinematic_tree.endEffectorLabels = ...
      { ...
        'head' ...
      };
    kinematic_tree.markerLabels = marker_labels;
    joint_directions = ...
      { ...
        'forward', 'backward'; ...
        'forward', 'backward'; ...
        'forward', 'backward'; ...
        'forward', 'backward'; ...
      }';

    % add segment labels and joint centers
    kinematic_tree.addSegmentLabel('SHANKS', 1);
    kinematic_tree.addSegmentLabel('THIGHS', 2);
    kinematic_tree.addSegmentLabel('TRUNK', 3);
    kinematic_tree.addSegmentLabel('HEAD', 4);
    
    % define markers
    marker_segment_list = createMarkerSegmentList(marker_labels, subject_settings);
    marker_color_list = createMarkerColorList(marker_labels, subject_settings);
    number_of_markers = length(marker_segment_list);
    for i_marker = 1 : number_of_markers
        kinematic_tree.addMarker(marker_segment_list(i_marker), marker_reference((i_marker-1)*3+1 : (i_marker-1)*3+3)', marker_color_list{i_marker});
    end

    kinematic_tree.updateInternals();

    %% save
    save ...
      ( ...
        'subjectModel', ...
        'kinematic_tree', ...
        'joint_labels', ...
        'joint_directions', ...
        'marker_labels', ...
        'marker_reference', ...
        'joint_center_labels', ...
        'joint_center_directions', ...
        'joint_center_reference', ...
        'segment_labels', ...
        'segment_masses', ...
        'segment_coms_wcs' ...
      );
  
    %% show visualization

    if visualize
        % show stick figure of geometric model
        scene_bound = repmat(hip_cor, 1, 2) + 2*[-0.5 0.5; -0.5 0.5; -1 1];
        stick_figure = KinematicTreeController(kinematic_tree, scene_bound, 'ellipsoid');

        % show joint centers
        for i_marker = 1 : (length(joint_center_reference) / 3)
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

end

function marker_segments = createMarkerSegmentList(marker_labels, subject_settings)
    marker_to_segment_map = subject_settings.get('marker_to_segment_map');
    number_of_markers = length(marker_labels) / 3;
    marker_segments = zeros(1, number_of_markers);
    for i_label = 1 : number_of_markers
        this_marker_label = marker_labels{i_label*3};
        this_marker_id = this_marker_label(1:end-2);
        marker_segments(i_label) = str2num(marker_to_segment_map{strcmp(marker_to_segment_map(:, 1), this_marker_id), 2}); %#ok<ST2NM>
    end
end

function marker_color_list = createMarkerColorList(marker_labels, subject_settings)
    marker_to_color_map = subject_settings.get('marker_to_color_map');
    number_of_markers = length(marker_labels) / 3;

    marker_color_list = cell(1, number_of_markers);
    for i_marker = 1 : number_of_markers
        this_marker_label = marker_labels{i_marker*3};
        this_marker_id = this_marker_label(1:end-2);
        this_marker_color = ...
          [
            str2num(marker_to_color_map{strcmp(marker_to_color_map(:, 1), this_marker_id), 2}) ...
            str2num(marker_to_color_map{strcmp(marker_to_color_map(:, 1), this_marker_id), 3}) ...
            str2num(marker_to_color_map{strcmp(marker_to_color_map(:, 1), this_marker_id), 4}) ...
          ]; %#ok<ST2NM>
        marker_color_list{i_marker} = this_marker_color;
    end
end













































