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

function createModel(varargin)
    load('subjectInfo.mat');
    
    % parse arguments
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;
    
%     subject_settings = loadSettingsFile('subjectSettings.txt');
    subject_settings = SettingsCustodian('subjectSettings.txt');

    % TODO: change stuff depending upon the static reference type

    load('subjectInfo.mat', 'date', 'subject_id', 'knee_width', 'ankle_width', 'elbow_width');

    % if width measurements are not available, use best guess
    if knee_width == 0
        knee_width = 0.1;
        disp('Knee width not specified, using 0.1m as a proxy. Please enter measurement in subjects.csv.');
    end
    if ankle_width == 0
        ankle_width = 0.07;
        disp('Ankle width not specified, using 0.07m as a proxy. Please enter measurement in subjects.csv.');
    end%
    if elbow_width == 0
        elbow_width = 0.07;
        disp('Elbow width not specified, using 0.07m as a proxy. Please enter measurement in subjects.csv.');
    end


    %% create static reference

    % load static reference file
    load(['processed' filesep makeFileName(date, subject_id, subject_settings.get('static_reference_trial_type'), subject_settings.get('static_reference_trial_number'), 'markerTrajectories')]);

    % find first time step where all markers are available
    i_time = 1;
    
    % For TU data
%     while any(isnan(marker_trajectories(i_time, :)))
%         i_time = i_time + 1;
%     end
%     marker_reference = marker_trajectories(i_time, :);

    % For UD data
    while any(isnan(marker_trajectories(i_time, :)))
        i_time = i_time + 1;
    end
    marker_reference = marker_trajectories(i_time, :);

    
    %% extract marker reference positions
    
    % head
    LFHD_reference = extractMarkerData(marker_reference, marker_labels, 'LFHD')';
    RFHD_reference = extractMarkerData(marker_reference, marker_labels, 'RFHD')';
    LBHD_reference = extractMarkerData(marker_reference, marker_labels, 'LBHD')';
    RBHD_reference = extractMarkerData(marker_reference, marker_labels, 'RBHD')';

    % torso
    C7_reference = extractMarkerData(marker_reference, marker_labels, 'C7')';
    T10_reference = extractMarkerData(marker_reference, marker_labels, 'T10')';
    CLAV_reference = extractMarkerData(marker_reference, marker_labels, 'CLAV')';
    STRN_reference = extractMarkerData(marker_reference, marker_labels, 'STRN')';
    RBAK_reference = extractMarkerData(marker_reference, marker_labels, 'RBAK')';

    % left arm
    LSHO_reference = extractMarkerData(marker_reference, marker_labels, 'LSHO')';
    LELB_reference = extractMarkerData(marker_reference, marker_labels, 'LELB')';
    LWRA_reference = extractMarkerData(marker_reference, marker_labels, 'LWRA')';
    LWRB_reference = extractMarkerData(marker_reference, marker_labels, 'LWRB')';
    LFIN_reference = extractMarkerData(marker_reference, marker_labels, 'LFIN')';

    % right arm
    RSHO_reference = extractMarkerData(marker_reference, marker_labels, 'RSHO')';
    RELB_reference = extractMarkerData(marker_reference, marker_labels, 'RELB')';
    RWRA_reference = extractMarkerData(marker_reference, marker_labels, 'RWRA')';
    RWRB_reference = extractMarkerData(marker_reference, marker_labels, 'RWRB')';
    RFIN_reference = extractMarkerData(marker_reference, marker_labels, 'RFIN')';

    % pelvis
    LASI_reference = extractMarkerData(marker_reference, marker_labels, 'LASI')';
    RASI_reference = extractMarkerData(marker_reference, marker_labels, 'RASI')';
    LPSI_reference = extractMarkerData(marker_reference, marker_labels, 'LPSI')';
    RPSI_reference = extractMarkerData(marker_reference, marker_labels, 'RPSI')';

    % left leg
    LTHI_reference = extractMarkerData(marker_reference, marker_labels, 'LTHI')';
    LTHIA_reference = extractMarkerData(marker_reference, marker_labels, 'LTHIA')';
    LKNE_reference = extractMarkerData(marker_reference, marker_labels, 'LKNE')';
    LTIB_reference = extractMarkerData(marker_reference, marker_labels, 'LTIB')';
    LTIBA_reference = extractMarkerData(marker_reference, marker_labels, 'LTIBA')';
    LANK_reference = extractMarkerData(marker_reference, marker_labels, 'LANK')';
    LHEE_reference = extractMarkerData(marker_reference, marker_labels, 'LHEE')';
    LTOE_reference = extractMarkerData(marker_reference, marker_labels, 'LTOE')';
    LTOEL_reference = extractMarkerData(marker_reference, marker_labels, 'LTOEL')';

    % left leg
    RTHI_reference = extractMarkerData(marker_reference, marker_labels, 'RTHI')';
    RTHIA_reference = extractMarkerData(marker_reference, marker_labels, 'RTHIA')';
    RKNE_reference = extractMarkerData(marker_reference, marker_labels, 'RKNE')';
    RTIB_reference = extractMarkerData(marker_reference, marker_labels, 'RTIB')';
    RTIBA_reference = extractMarkerData(marker_reference, marker_labels, 'RTIBA')';
    RANK_reference = extractMarkerData(marker_reference, marker_labels, 'RANK')';
    RHEE_reference = extractMarkerData(marker_reference, marker_labels, 'RHEE')';
    RTOE_reference = extractMarkerData(marker_reference, marker_labels, 'RTOE')';
    RTOEL_reference = extractMarkerData(marker_reference, marker_labels, 'RTOEL')';

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
    right_direction_prime = normVector(RASI_reference - LASI_reference);
    anterior_direction_prime = normVector(MASIS_reference - MPSIS_reference);
    proximal_direction_prime = cross(right_direction_prime, anterior_direction_prime);

    body_direction_matrix = orthogonalizeBasis([right_direction_prime anterior_direction_prime proximal_direction_prime]);
    right_direction = body_direction_matrix(:, 1);
    anterior_direction = body_direction_matrix(:, 2);
    proximal_direction = body_direction_matrix(:, 3);

    left_direction = - right_direction;
    posterior_direction = - anterior_direction;
    distal_direction = - proximal_direction;
    
    up_direction = proximal_direction;
    down_direction = distal_direction;

    % calculate anatomical landmarks and apply marker and flesh offsets
    % TODO: some of these corrections depend upon the type of reference (ski, motorcycle, casual, anatomical)
    centroid_to_skin_correction = 0.0152; % in meters
    skin_to_bone_correction_ASIS = 0.01; % in meters
    skin_to_bone_correction_PSIS = 0.01; % in meters
    c7 = C7_reference + centroid_to_skin_correction * anterior_direction;
    suprasternale = CLAV_reference + centroid_to_skin_correction * posterior_direction;
    left_acromion = LSHO_reference + centroid_to_skin_correction * distal_direction;
    left_lateral_humeral_epicondyle = LELB_reference + centroid_to_skin_correction * distal_direction;
    
    if strcmp(subject_settings.get('static_reference_posture'), 'ski')
        left_inner_wrist = LWRA_reference + centroid_to_skin_correction * down_direction;
        left_outer_wrist = LWRB_reference + centroid_to_skin_correction * up_direction;
        left_hand = LFIN_reference + centroid_to_skin_correction * right_direction;
        right_inner_wrist = RWRA_reference + centroid_to_skin_correction * down_direction;
        right_outer_wrist = RWRB_reference + centroid_to_skin_correction * up_direction;
        right_hand = RFIN_reference + centroid_to_skin_correction * left_direction;
    end
    if strcmp(subject_settings.get('static_reference_posture'), 'casual')
        left_wrist_marker_direction = normVector(LWRB_reference - LWRA_reference);
        left_inner_wrist = LWRA_reference + centroid_to_skin_correction * left_wrist_marker_direction;
        left_outer_wrist = LWRB_reference - centroid_to_skin_correction * left_wrist_marker_direction;
        left_hand = LFIN_reference;% + centroid_to_skin_correction * right_direction;
        right_wrist_marker_direction = normVector(RWRB_reference - RWRA_reference);
        right_inner_wrist = RWRA_reference + centroid_to_skin_correction * right_wrist_marker_direction;
        right_outer_wrist = RWRB_reference - centroid_to_skin_correction * right_wrist_marker_direction;
        right_hand = RFIN_reference;% + centroid_to_skin_correction * left_direction;
    end
    if strcmp(subject_settings.get('static_reference_posture'), 'motorcycle')
        % TODO: not tested yet
        left_inner_wrist = LWRA_reference + centroid_to_skin_correction * left_direction;
        left_outer_wrist = LWRB_reference + centroid_to_skin_correction * right_direction;
        left_hand = LFIN_reference + centroid_to_skin_correction * down_direction;
        right_inner_wrist = RWRA_reference + centroid_to_skin_correction * right_direction;
        right_outer_wrist = RWRB_reference + centroid_to_skin_correction * left_direction;
        right_hand = RFIN_reference + centroid_to_skin_correction * down_direction;
    end
    
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
    head_center_to_vertex = head_width / 2; % this is an ad-hoc assumption that seems to work out well in some examples
    head_center = mean([LFHD_reference RFHD_reference LBHD_reference RBHD_reference], 2);
    head_vertex = head_center + head_center_to_vertex * proximal_direction;
    sellion = mean([LFHD_reference RFHD_reference], 2);

    % calculate some distances
    inter_ASIS_distance = norm(LASIS - RASIS);

    %% estimate hip joint centers
    if strcmp(subject_settings.get('hip_joint_center_estimation_method'), 'SCoRE')
        pelvis_center_reference = mean(reshape(pelvis_markers_reference, 3, size(pelvis_markers_reference, 2)/3), 2);

        % find left hip CoR
        left_hip_reference_file_name = ['processed' filesep makeFileName(date, subject_id, 'calibration', left_hip_calibration_file_index, 'markerTrajectories')];
        disp(['Left hip reference file name: ' left_hip_reference_file_name]);
        load(left_hip_reference_file_name);
        hip_reference = marker_trajectories;

        LASI_trajectory = extractMarkerData(hip_reference, marker_labels, 'LASI')';
        RASI_trajectory = extractMarkerData(hip_reference, marker_labels, 'RASI')';
        LPSI_trajectory = extractMarkerData(hip_reference, marker_labels, 'LPSI')';
        RPSI_trajectory = extractMarkerData(hip_reference, marker_labels, 'RPSI')';
        LTHI_trajectory = extractMarkerData(hip_reference, marker_labels, 'LTHI')';
        LTHIA_trajectory = extractMarkerData(hip_reference, marker_labels, 'LTHIA')';
        LKNE_trajectory = extractMarkerData(hip_reference, marker_labels, 'LKNE')';
        pelvis_markers_trajectory = [LASI_trajectory' RASI_trajectory' LPSI_trajectory' RPSI_trajectory'];
        left_thigh_markers_trajectory = [LTHI_trajectory' LTHIA_trajectory' LKNE_trajectory'];

        left_hip_cor = estimateJointKinematics ...
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

        LASI_trajectory = extractMarkerData(hip_reference, marker_labels, 'LASI')';
        RASI_trajectory = extractMarkerData(hip_reference, marker_labels, 'RASI')';
        LPSI_trajectory = extractMarkerData(hip_reference, marker_labels, 'LPSI')';
        RPSI_trajectory = extractMarkerData(hip_reference, marker_labels, 'RPSI')';
        RTHI_trajectory = extractMarkerData(hip_reference, marker_labels, 'RTHI')';
        RTHIA_trajectory = extractMarkerData(hip_reference, marker_labels, 'RTHIA')';
        RKNE_trajectory = extractMarkerData(hip_reference, marker_labels, 'RKNE')';
        pelvis_markers_trajectory = [LASI_trajectory' RASI_trajectory' LPSI_trajectory' RPSI_trajectory'];
        right_thigh_markers_trajectory = [RTHI_trajectory' RTHIA_trajectory' RKNE_trajectory'];
        right_hip_cor = estimateJointKinematics ...
          ( ...
            pelvis_markers_reference, ...
            right_thigh_markers_reference, ...
            pelvis_markers_trajectory, ...
            right_thigh_markers_trajectory ...
          );


    elseif strcmp(subject_settings.get('hip_joint_center_estimation_method'), 'Tylkowski')

        MASIS = mean([LASIS RASIS], 2);
        MPSIS = mean([LPSIS RPSIS], 2);
        z_pelvis = normVector(LASIS - RASIS); % medial-lateral, positive is left
        x_pelvis = normVector(MPSIS - MASIS); % anterior-posterior, positive is backwards
        y_pelvis = normVector(cross(z_pelvis, x_pelvis)); % proximal-distal, positive is up

        hjc_correction_factor_lateral = 0.14; % Tylkowski, after Bell et al., 1990
        hjc_correction_factor_distal = 0.30; % Tylkowski, after Bell et al., 1990
        hjc_correction_factor_anterior = 0.19; % Tylkowski, after Bell et al., 1990

        % XXX playing around with this factor to see if it improves the problem with the kinematic chain leg markers being
        % higher than the measured markers
%         hjc_correction_factor_anterior = 0.09; % Tylkowski, after Bell et al., 1990
        % TODO: figure out a systematic way to do this, and also do it for some other joints

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
    if strcmp(subject_settings.get('knee_joint_axis_estimation_method'), 'SARA')
        % find left knee CoR
        left_knee_reference_file_name = ['processed' filesep makeFileName(date, subject_id, 'calibration', subject_settings.get('left_knee_calibration_file_index'), 'markerTrajectories')];
        disp(['Left knee reference file name: ' left_knee_reference_file_name]);
        load(left_knee_reference_file_name);
        knee_reference = marker_trajectories;

        LTHI_trajectory = extractMarkerData(knee_reference, marker_labels, 'LTHI')';
        LTHIA_trajectory = extractMarkerData(knee_reference, marker_labels, 'LTHIA')';
        LKNE_trajectory = extractMarkerData(knee_reference, marker_labels, 'LKNE')';
        LTIB_trajectory = extractMarkerData(knee_reference, marker_labels, 'LTIB')';
        LTIBA_trajectory = extractMarkerData(knee_reference, marker_labels, 'LTIBA')';
        LANK_trajectory = extractMarkerData(knee_reference, marker_labels, 'LANK')';
        left_thigh_markers_trajectory = [LTHI_trajectory' LTHIA_trajectory' LKNE_trajectory'];
        left_shank_markers_trajectory = [LTIB_trajectory' LTIBA_trajectory' LANK_trajectory'];


        [~, ~, left_knee_flexion_axis] = estimateJointKinematics ...
          ( ...
            left_thigh_markers_reference, ...
            left_shank_markers_reference, ...
            left_thigh_markers_trajectory, ...
            left_shank_markers_trajectory, ...
            1 ...
          );

        % find right knee CoR
        right_knee_reference_file_name = ['processed' filesep makeFileName(date, subject_id, 'calibration', subject_settings.get('right_knee_calibration_file_index'), 'markerTrajectories')];
        disp(['Right knee reference file name: ' right_knee_reference_file_name]);
        load(right_knee_reference_file_name);
        knee_reference = marker_trajectories;

        RTHI_trajectory = extractMarkerData(knee_reference, marker_labels, 'RTHI')';
        RTHIA_trajectory = extractMarkerData(knee_reference, marker_labels, 'RTHIA')';
        RKNE_trajectory = extractMarkerData(knee_reference, marker_labels, 'RKNE')';
        RTIB_trajectory = extractMarkerData(knee_reference, marker_labels, 'RTIB')';
        RTIBA_trajectory = extractMarkerData(knee_reference, marker_labels, 'RTIBA')';
        RANK_trajectory = extractMarkerData(knee_reference, marker_labels, 'RANK')';
        right_thigh_markers_trajectory = [RTHI_trajectory' RTHIA_trajectory' RKNE_trajectory'];
        right_shank_markers_trajectory = [RTIB_trajectory' RTIBA_trajectory' RANK_trajectory'];

        [~, ~, right_knee_flexion_axis] = estimateJointKinematics ...
          ( ...
            right_thigh_markers_reference, ...
            right_shank_markers_reference, ...
            right_thigh_markers_trajectory, ...
            right_shank_markers_trajectory, ...
            1 ...
          );
    elseif strcmp(subject_settings.get('knee_joint_axis_estimation_method'), 'markers')
        % assume that the knee axis of rotation is the vector between the knee markers
        left_knee_flexion_axis = normVector(left_lateral_femoral_epicondyle - right_lateral_femoral_epicondyle);
        right_knee_flexion_axis = normVector(left_lateral_femoral_epicondyle - right_lateral_femoral_epicondyle);

    else
        error('hip joint center estimation method not recognized. Options are "SARA" or "markers".');
    end
    % correct directions
    if dot(right_knee_flexion_axis, left_direction) < 0
        right_knee_flexion_axis = -right_knee_flexion_axis;
    end
    if dot(left_knee_flexion_axis, left_direction) < 0
        left_knee_flexion_axis = -left_knee_flexion_axis;
    end

    % estimate knee CoRs
    kjc_correction_factor = 0.5;
    left_knee_cor = left_lateral_femoral_epicondyle - kjc_correction_factor*knee_width*left_knee_flexion_axis;
    right_knee_cor = right_lateral_femoral_epicondyle + kjc_correction_factor*knee_width*right_knee_flexion_axis;

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
    left_ankle_cor = left_lateral_malleolus - ajc_correction_factor * ankle_width * left_knee_flexion_axis;
    right_ankle_cor = right_lateral_malleolus + ajc_correction_factor * ankle_width * right_knee_flexion_axis;

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
    if strcmp(subject_settings.get('static_reference_posture'), 'ski') || strcmp(subject_settings.get('static_reference_posture'), 'casual')
        left_elbow_cor = left_lateral_humeral_epicondyle + ejc_correction_factor*elbow_width*right_direction;
        left_elbow_axis = right_direction;
        right_elbow_cor = right_lateral_humeral_epicondyle - ejc_correction_factor*elbow_width*right_direction;
        right_elbow_axis = right_direction;
    end
    if strcmp(subject_settings.get('static_reference_posture'), 'motorcycle')
        % TODO: not tested yet
        left_elbow_cor = left_lateral_humeral_epicondyle + ejc_correction_factor*elbow_width*distal_direction;
        left_elbow_axis = distal_direction;
        right_elbow_cor = right_lateral_humeral_epicondyle + ejc_correction_factor*elbow_width*distal_direction;
        right_elbow_axis = proximal_direction;
    end

    % define hip directions
    left_hip_internal_rotation_axis_prime = left_knee_cor - left_hip_cor;
    left_hip_direction_matrix_prime = orthogonalizeBasis([left_hip_internal_rotation_axis_prime, left_direction, cross(left_direction, left_hip_internal_rotation_axis_prime)]);
    left_hip_direction_matrix = [-left_hip_direction_matrix_prime(:, 2), -left_hip_direction_matrix_prime(:, 3), -left_hip_direction_matrix_prime(:, 1)];
    left_hip_flexion_axis = left_hip_direction_matrix(:, 1);
    left_hip_abduction_axis = left_hip_direction_matrix(:, 2);
    left_hip_internal_rotation_axis = -left_hip_direction_matrix(:, 3);

    right_hip_internal_rotation_axis_prime = right_knee_cor - right_hip_cor;
    right_hip_direction_matrix_prime = orthogonalizeBasis([right_hip_internal_rotation_axis_prime, left_direction, cross(left_direction, right_hip_internal_rotation_axis_prime)]);
    right_hip_direction_matrix = [-right_hip_direction_matrix_prime(:, 2), -right_hip_direction_matrix_prime(:, 3), -right_hip_direction_matrix_prime(:, 1)];
    right_hip_flexion_axis = right_hip_direction_matrix(:, 1);
    right_hip_abduction_axis = -right_hip_direction_matrix(:, 2);
    right_hip_internal_rotation_axis = right_hip_direction_matrix(:, 3);

    % define knee directions
    left_knee_internal_rotation_axis_prime = left_ankle_cor - left_knee_cor;
    left_knee_direction_matrix_prime = orthogonalizeBasis([left_knee_internal_rotation_axis_prime, left_knee_flexion_axis, cross(left_knee_flexion_axis, left_knee_internal_rotation_axis_prime)]);
    left_knee_direction_matrix = [-left_knee_direction_matrix_prime(:, 2), -left_knee_direction_matrix_prime(:, 3), -left_knee_direction_matrix_prime(:, 1)];
    left_knee_flexion_axis = -left_knee_direction_matrix(:, 1);
    left_knee_internal_rotation_axis = left_knee_direction_matrix(:, 3);

    right_knee_internal_rotation_axis_prime = right_ankle_cor - right_knee_cor;
    right_knee_direction_matrix_prime = orthogonalizeBasis([right_knee_internal_rotation_axis_prime, right_knee_flexion_axis, cross(right_knee_flexion_axis, right_knee_internal_rotation_axis_prime)]);
    right_knee_direction_matrix = [-right_knee_direction_matrix_prime(:, 2), -right_knee_direction_matrix_prime(:, 3), -right_knee_direction_matrix_prime(:, 1)];
    right_knee_flexion_axis = -right_knee_direction_matrix(:, 1);
    right_knee_internal_rotation_axis = -right_knee_direction_matrix(:, 3);

    % define ankle axes
    left_ankle_inversion_axis_prime = left_toe_mid - left_ankle_cor;
    left_ankle_direction_matrix_prime = orthogonalizeBasis([left_ankle_inversion_axis_prime, right_direction, cross(left_ankle_inversion_axis_prime, right_direction)]);
    left_ankle_direction_matrix = [left_ankle_direction_matrix_prime(:, 2), left_ankle_direction_matrix_prime(:, 1), -left_ankle_direction_matrix_prime(:, 3)];
    left_ankle_inversion_axis = left_ankle_direction_matrix(:, 2);
    left_ankle_dorsiflexion_axis = left_ankle_direction_matrix(:, 1);

    % define ankle axes
    right_ankle_inversion_axis_prime = right_toe_mid - right_ankle_cor;
    right_ankle_direction_matrix_prime = orthogonalizeBasis([right_ankle_inversion_axis_prime, right_direction, cross(right_ankle_inversion_axis_prime, right_direction)]);
    right_ankle_direction_matrix = [right_ankle_direction_matrix_prime(:, 2), right_ankle_direction_matrix_prime(:, 1), -right_ankle_direction_matrix_prime(:, 3)];
    right_ankle_inversion_axis = -right_ankle_direction_matrix(:, 2);
    right_ankle_dorsiflexion_axis = right_ankle_direction_matrix(:, 1);

    % define lumbar directions
    lumbar_internal_rotation_axis_prime = cervix_cor - lumbar_cor;
    lumbar_direction_matrix_prime = orthogonalizeBasis([lumbar_internal_rotation_axis_prime, left_direction, cross(left_direction, lumbar_internal_rotation_axis_prime)]);
    lumbar_direction_matrix = [-lumbar_direction_matrix_prime(:, 2), lumbar_direction_matrix_prime(:, 3), lumbar_direction_matrix_prime(:, 1)];
    lumbar_flexion_axis = -lumbar_direction_matrix(:, 1);
    lumbar_abduction_axis = lumbar_direction_matrix(:, 2);
    lumbar_internal_rotation_axis = -lumbar_direction_matrix(:, 3);

    % define cervix directions
    cervix_internal_rotation_axis_prime = head_center - cervix_cor;
    cervix_direction_matrix_prime = orthogonalizeBasis([cervix_internal_rotation_axis_prime, left_direction, cross(left_direction, cervix_internal_rotation_axis_prime)]);
    cervix_direction_matrix = [-cervix_direction_matrix_prime(:, 2), cervix_direction_matrix_prime(:, 3), cervix_direction_matrix_prime(:, 1)];
    cervix_flexion_axis = -cervix_direction_matrix(:, 1);
    cervix_abduction_axis = cervix_direction_matrix(:, 2);
    cervix_internal_rotation_axis = -cervix_direction_matrix(:, 3);

    % define shoulder directions
    left_shoulder_internal_rotation_axis_prime = left_elbow_cor - left_shoulder_cor;
    left_shoulder_direction_matrix_prime = orthogonalizeBasis([left_shoulder_internal_rotation_axis_prime, left_direction, cross(left_direction, left_shoulder_internal_rotation_axis_prime)]);
    left_shoulder_direction_matrix = [-left_shoulder_direction_matrix_prime(:, 2), -left_shoulder_direction_matrix_prime(:, 3), -left_shoulder_direction_matrix_prime(:, 1)];
    left_shoulder_flexion_axis = left_shoulder_direction_matrix(:, 1);
    left_shoulder_abduction_axis = left_shoulder_direction_matrix(:, 2);
    left_shoulder_internal_rotation_axis = -left_shoulder_direction_matrix(:, 3);
    
    right_shoulder_internal_rotation_axis_prime = right_elbow_cor - right_shoulder_cor;
    right_shoulder_direction_matrix_prime = orthogonalizeBasis([right_shoulder_internal_rotation_axis_prime, left_direction, cross(left_direction, right_shoulder_internal_rotation_axis_prime)]);
    right_shoulder_direction_matrix = [-right_shoulder_direction_matrix_prime(:, 2), -right_shoulder_direction_matrix_prime(:, 3), -right_shoulder_direction_matrix_prime(:, 1)];
    right_shoulder_flexion_axis = right_shoulder_direction_matrix(:, 1);
    right_shoulder_abduction_axis = -right_shoulder_direction_matrix(:, 2);
    right_shoulder_internal_rotation_axis = right_shoulder_direction_matrix(:, 3);
    
    % define elbow directions
    left_elbow_internal_rotation_axis_prime = left_wrist_cor - left_elbow_cor;
    left_elbow_direction_matrix_prime = orthogonalizeBasis([left_elbow_internal_rotation_axis_prime, left_elbow_axis, cross(left_elbow_axis, left_elbow_internal_rotation_axis_prime)]);
    left_elbow_direction_matrix = [left_elbow_direction_matrix_prime(:, 2), left_elbow_direction_matrix_prime(:, 1), left_elbow_direction_matrix_prime(:, 3)];
    left_elbow_flexion_axis = left_elbow_direction_matrix(:, 1);
    left_elbow_internal_rotation_axis = left_elbow_direction_matrix(:, 2);

    right_elbow_internal_rotation_axis_prime = right_wrist_cor - right_elbow_cor;
    right_elbow_direction_matrix_prime = orthogonalizeBasis([right_elbow_internal_rotation_axis_prime, right_elbow_axis, cross(right_elbow_axis, right_elbow_internal_rotation_axis_prime)]);
    right_elbow_direction_matrix = [right_elbow_direction_matrix_prime(:, 2), right_elbow_direction_matrix_prime(:, 1), right_elbow_direction_matrix_prime(:, 3)];
    right_elbow_flexion_axis = right_elbow_direction_matrix(:, 1);
    right_elbow_internal_rotation_axis = -right_elbow_direction_matrix(:, 2);
    
    % define wrist axes
    left_wrist_inversion_axis_prime = cross(left_wrist_flexion_axis, left_elbow_internal_rotation_axis);
    left_wrist_direction_matrix_prime = orthogonalizeBasis([left_wrist_flexion_axis, left_wrist_inversion_axis_prime, cross(left_wrist_flexion_axis, left_wrist_inversion_axis_prime)]);
    left_wrist_direction_matrix = [left_wrist_direction_matrix_prime(:, 2), -left_wrist_direction_matrix_prime(:, 3), -left_wrist_direction_matrix_prime(:, 1)];
%     left_wrist_inversion_axis = left_wrist_direction_matrix(:, 1);
    left_wrist_to_eef = norm(LFIN_reference - LWRA_reference);
    left_hand_mid = left_wrist_cor + left_wrist_direction_matrix(:, 2) * left_wrist_to_eef;
    
    right_wrist_inversion_axis_prime = cross(right_wrist_flexion_axis, -right_elbow_internal_rotation_axis);
    right_wrist_direction_matrix_prime = orthogonalizeBasis([right_wrist_flexion_axis, right_wrist_inversion_axis_prime, cross(right_wrist_flexion_axis, right_wrist_inversion_axis_prime)]);
    right_wrist_direction_matrix = [-right_wrist_direction_matrix_prime(:, 2), -right_wrist_direction_matrix_prime(:, 3), right_wrist_direction_matrix_prime(:, 1)];
%     right_wrist_inversion_axis = right_wrist_direction_matrix(:, 1);
    right_wrist_to_eef = norm(RFIN_reference - RWRA_reference);
    right_hand_mid = right_wrist_cor + right_wrist_direction_matrix(:, 2) * right_wrist_to_eef;
    
    % TODO: some of these assumptions are not valid for all types of reference configurations

    %% define scaling factors
    % according to R. Dumas , L. Cheze, J.-P. Verriest: "Adjustments to McConville et al. and Young et al. body
    % segment inertial parameters", Journal of Biomechanics 40 (2007) 543?553
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

%         arm_rxx_scaling_factor     = 0.33;  arm_ryy_scaling_factor     = 0.17;  arm_rzz_scaling_factor     = 0.33;  arm_rxy_scaling_factor     = 0.03;      arm_rxz_scaling_factor     = 0.05*1i;   arm_ryz_scaling_factor     = 0.14;
%         forearm_rxx_scaling_factor = 0.26;  forearm_ryy_scaling_factor = 0.14;  forearm_rzz_scaling_factor = 0.25;  forearm_rxy_scaling_factor = 0.10;      forearm_rxz_scaling_factor = 0.04;      forearm_ryz_scaling_factor = 0.14*1i;
        % ACHTUNG: used the values for males here because the values for females give rather weird looking results
        arm_rxx_scaling_factor     = 0.31;  arm_ryy_scaling_factor     = 0.14;  arm_rzz_scaling_factor     = 0.32;  arm_rxy_scaling_factor     = 0.06;      arm_rxz_scaling_factor     = 0.05;      arm_ryz_scaling_factor     = 0.02;
        forearm_rxx_scaling_factor = 0.28;  forearm_ryy_scaling_factor = 0.11;  forearm_rzz_scaling_factor = 0.27;  forearm_rxy_scaling_factor = 0.03;      forearm_rxz_scaling_factor = 0.02;      forearm_ryz_scaling_factor = 0.08*1i;
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
    left_thigh_scs_x = normVector(cross(left_knee_flexion_axis, left_thigh_scs_y));
    left_thigh_scs_z = cross(left_thigh_scs_x, left_thigh_scs_y);

    % left leg
    left_leg_scs_y = normVector(left_knee_cor - left_ankle_cor);
    left_leg_scs_x = normVector(cross(left_knee_flexion_axis, left_leg_scs_y));
    left_leg_scs_z = cross(left_leg_scs_x, left_leg_scs_y);

    % left foot
    left_foot_scs_x = normVector(left_toe_mid - left_calcaneus);
    left_foot_scs_y = normVector(cross(left_toe_mid-left_calcaneus, left_direction));
    left_foot_scs_z = cross(left_foot_scs_x, left_foot_scs_y);

     % right thigh
    right_thigh_scs_y = normVector(right_hip_cor - right_knee_cor);
    right_thigh_scs_x = normVector(cross(right_knee_flexion_axis, right_thigh_scs_y));
    right_thigh_scs_z = cross(right_thigh_scs_x, right_thigh_scs_y);

    % right leg
    right_leg_scs_y = normVector(right_knee_cor - right_ankle_cor);
    right_leg_scs_x = normVector(cross(right_knee_flexion_axis, right_leg_scs_y));
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
    left_arm_scs_z = cross(left_arm_scs_x, left_arm_scs_y);

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
    right_arm_scs_z = cross(right_arm_scs_x, right_arm_scs_y);

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
        right_ankle_cor', ...
        left_toe_mid', ...
        right_toe_mid', ...
        left_hand_mid', ...
        right_hand_mid' ...
      ];

    joint_center_labels_single = ...
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
        'LTOESEEF', ...
        'RTOESEEF', ...
        'LHANDEEF', ...
        'RHANDEEF' ...
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
    [joint_center_directions{1, 1 : 3 : number_of_joint_centers}] = deal('right');
    [joint_center_directions{2, 1 : 3 : number_of_joint_centers}] = deal('left');
    [joint_center_directions{1, 2 : 3 : number_of_joint_centers}] = deal('forward');
    [joint_center_directions{2, 2 : 3 : number_of_joint_centers}] = deal('backward');
    [joint_center_directions{1, 3 : 3 : number_of_joint_centers}] = deal('up');
    [joint_center_directions{2, 3 : 3 : number_of_joint_centers}] = deal('down');
    

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
    number_of_segments = length(segment_labels);

    mcs_to_wcs_transformations = calculateMcsToWcsTransformations_new([marker_reference joint_center_reference], [marker_labels joint_center_labels], segment_labels, subject_settings);
    pelvis_transformation_current = mcs_to_wcs_transformations{strcmp(segment_labels, 'PELVIS')};

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

    %% assemble kinematic tree

    % calculate pelvis base point as center between the two hip CoRs
%     pelvis_base_point = [left_hip_cor + right_hip_cor] / 2;
    pelvis_base_point = mean([LASI_reference RASI_reference LPSI_reference RPSI_reference], 2);

    % adjust end-effector positions for foot extension in z-direction
    right_heel = right_ankle_cor; right_heel(3) = right_toe_mid(3);
    left_heel = left_ankle_cor; left_heel(3) = left_toe_mid(3);

    joint_positions = ...
    { ...
        pelvis_base_point, pelvis_base_point, pelvis_base_point, pelvis_base_point, pelvis_base_point, pelvis_base_point, ...
        left_hip_cor, left_hip_cor, left_hip_cor, left_knee_cor, left_knee_cor, left_ankle_cor, left_ankle_cor, ...
        right_hip_cor, right_hip_cor, right_hip_cor, right_knee_cor, right_knee_cor, right_ankle_cor, right_ankle_cor, ...
        lumbar_cor, lumbar_cor, lumbar_cor, cervix_cor, cervix_cor, cervix_cor, ...
        left_shoulder_cor, left_shoulder_cor, left_shoulder_cor, left_elbow_cor, left_elbow_cor, left_wrist_cor ...
        right_shoulder_cor, right_shoulder_cor, right_shoulder_cor, right_elbow_cor, right_elbow_cor, right_wrist_cor ...
    };

    joint_axes = ...
    { ...
        e_1, e_2, e_3, e_3, e_1, e_2, ...       % pelvis free body DoFs
        left_hip_flexion_axis, left_hip_abduction_axis, left_hip_internal_rotation_axis, left_knee_flexion_axis, left_knee_internal_rotation_axis, left_ankle_dorsiflexion_axis, left_ankle_inversion_axis, ...   % left leg
        right_hip_flexion_axis, right_hip_abduction_axis, right_hip_internal_rotation_axis, right_knee_flexion_axis, right_knee_internal_rotation_axis, right_ankle_dorsiflexion_axis, right_ankle_inversion_axis,...   % right leg
        lumbar_flexion_axis, lumbar_abduction_axis, lumbar_internal_rotation_axis, cervix_flexion_axis, cervix_abduction_axis, cervix_internal_rotation_axis, ...       % trunk and neck
        left_shoulder_flexion_axis, left_shoulder_abduction_axis, left_shoulder_internal_rotation_axis, left_elbow_flexion_axis, left_elbow_internal_rotation_axis, left_wrist_flexion_axis, ... % left arm
        right_shoulder_flexion_axis, right_shoulder_abduction_axis, right_shoulder_internal_rotation_axis, right_elbow_flexion_axis, right_elbow_internal_rotation_axis, right_wrist_flexion_axis, ... % right arm
    };

    joint_types = ...
      [ ...
        2 2 2 1 1 1 ... % virtual dofs
        1 1 1 1 1 1 1 ... % left leg
        1 1 1 1 1 1 1 ... % right leg
        1 1 1 1 1 1 ... % l5 and neck
        1 1 1 1 1 1 ... % left arm
        1 1 1 1 1 1 ... % right arm
      ]; % 1 = rotational, 2 = prismatic

    % link setup
    link_com_positions =  ...
    { ...
        pelvis_com, pelvis_com, pelvis_com, pelvis_com, pelvis_com, pelvis_com, ...
        left_thigh_com, left_thigh_com, left_thigh_com, left_leg_com, left_leg_com, left_foot_com, left_foot_com, ...
        right_thigh_com, right_thigh_com, right_thigh_com, right_leg_com, right_leg_com, right_foot_com, right_foot_com, ...
        torso_com, torso_com, torso_com, head_com, head_com, head_com ...
        left_arm_com, left_arm_com, left_arm_com, left_forearm_com, left_forearm_com, left_hand_com, ...
        right_arm_com, right_arm_com, right_arm_com, right_forearm_com, right_forearm_com, right_hand_com ...
    };

    link_orientations = ...
      {
        eye(3); eye(3); eye(3); eye(3); eye(3); [pelvis_scs_x pelvis_scs_y pelvis_scs_z]; ...
        eye(3); eye(3); [left_thigh_scs_x left_thigh_scs_y left_thigh_scs_z]; eye(3); [left_leg_scs_x left_leg_scs_y left_leg_scs_z]; eye(3); [left_foot_scs_x left_foot_scs_y left_foot_scs_z]; ...
        eye(3); eye(3); [right_thigh_scs_x right_thigh_scs_y right_thigh_scs_z]; eye(3); [right_leg_scs_x right_leg_scs_y right_leg_scs_z]; eye(3); [right_foot_scs_x right_foot_scs_y right_foot_scs_z]; ...
        eye(3); eye(3); [torso_scs_x torso_scs_y torso_scs_z]; eye(3); eye(3); [head_scs_x head_scs_y head_scs_z]; ...
        eye(3); eye(3); [left_arm_scs_x left_arm_scs_y left_arm_scs_z]; eye(3); [left_forearm_scs_x left_forearm_scs_y left_forearm_scs_z]; [left_hand_scs_x left_hand_scs_y left_hand_scs_z]; ...
        eye(3); eye(3); [right_arm_scs_x right_arm_scs_y right_arm_scs_z]; eye(3); [right_forearm_scs_x right_forearm_scs_y right_forearm_scs_z]; [right_hand_scs_x right_hand_scs_y right_hand_scs_z]; ...
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
        zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_left_thigh, zeros(6, 6), generalized_inertia_matrix_left_shank, zeros(6, 6), generalized_inertia_matrix_left_foot, ...
        zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_right_thigh, zeros(6, 6), generalized_inertia_matrix_right_shank, zeros(6, 6), generalized_inertia_matrix_right_foot, ...
        zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_torso, zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_head ...
        zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_left_arm, zeros(6, 6), generalized_inertia_matrix_left_forearm, generalized_inertia_matrix_left_hand, ...
        zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_right_arm, zeros(6, 6), generalized_inertia_matrix_right_forearm, generalized_inertia_matrix_right_hand, ...
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
        [left_hand_scs_x left_hand_scs_y left_hand_scs_z, left_hand_mid; 0 0 0 1], ...
        [right_hand_scs_x right_hand_scs_y right_hand_scs_z, right_hand_mid; 0 0 0 1], ...
        pelvis_transformation_current, ...
      };

    branch_matrix = ...
      [ ...
        1 1 1 1 1 1    1 1 1 1 1 1 1   0 0 0 0 0 0 0   0 0 0   0 0 0   0 0 0 0 0 0   0 0 0 0 0 0; ... % left leg until heel
        1 1 1 1 1 1    1 1 1 1 1 1 1   0 0 0 0 0 0 0   0 0 0   0 0 0   0 0 0 0 0 0   0 0 0 0 0 0; ... % left leg until toes
        1 1 1 1 1 1    1 1 1 1 1 1 1   0 0 0 0 0 0 0   0 0 0   0 0 0   0 0 0 0 0 0   0 0 0 0 0 0; ... % left leg until ankle
        1 1 1 1 1 1    0 0 0 0 0 0 0   1 1 1 1 1 1 1   0 0 0   0 0 0   0 0 0 0 0 0   0 0 0 0 0 0; ... % right leg until heel
        1 1 1 1 1 1    0 0 0 0 0 0 0   1 1 1 1 1 1 1   0 0 0   0 0 0   0 0 0 0 0 0   0 0 0 0 0 0; ... % right leg until toes
        1 1 1 1 1 1    0 0 0 0 0 0 0   1 1 1 1 1 1 1   0 0 0   0 0 0   0 0 0 0 0 0   0 0 0 0 0 0; ... % right leg until ankle
        1 1 1 1 1 1    0 0 0 0 0 0 0   0 0 0 0 0 0 0   1 1 1   1 1 1   0 0 0 0 0 0   0 0 0 0 0 0; ... % head
        1 1 1 1 1 1    0 0 0 0 0 0 0   0 0 0 0 0 0 0   1 1 1   0 0 0   1 1 1 1 1 1   0 0 0 0 0 0; ... % left arm
        1 1 1 1 1 1    0 0 0 0 0 0 0   0 0 0 0 0 0 0   1 1 1   0 0 0   0 0 0 0 0 0   1 1 1 1 1 1; ... % right arm
        1 1 1 1 1 1    0 0 0 0 0 0 0   0 0 0 0 0 0 0   0 0 0   0 0 0   0 0 0 0 0 0   0 0 0 0 0 0; ... % pelvis
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
        'pelvis, x-translation', ...
        'pelvis, y-translation', ...
        'pelvis, z-translation', ...
        'pelvis, z-rotation', ...
        'pelvis, x-rotation', ...
        'pelvis, y-rotation', ...
        'left hip flexion/extension', ...
        'left hip ab/adduction', ...
        'left hip internal/external rotation', ...
        'left knee flexion/extension', ...
        'left knee external/internal rotation', ...
        'left ankle dorsi/plantarflexion', ...
        'left ankle eversion/inversion', ...
        'right hip flexion/extension', ...
        'right hip ab/adduction', ...
        'right hip internal/external rotation', ...
        'right knee flexion/extension', ...
        'right knee external/internal rotation', ...
        'right ankle dorsi/plantarflexion', ...
        'right ankle eversion/inversion', ...
        'lumbar joint - forward/backward bending', ...
        'lumbar joint - sideways bending (right/left)', ...
        'lumbar joint - internal rotation (right/left)', ...
        'cervical joint - forward/backward bending', ...
        'cervical joint - sideways bending (right/left)', ...
        'cervical joint - internal rotation (right/left)', ...
        'left shoulder flexion/extension', ...
        'left shoulder ab/adduction', ...
        'left shoulder in/external rotation', ...
        'left elbow flexion/extension', ...
        'left pronation/supination', ...
        'left wrist flexion/extension', ...
        'right shoulder flexion/extension', ...
        'right shoulder ab/adduction', ...
        'right shoulder in/external rotation', ...
        'right elbow flexion/extension', ...
        'right pronation/supination', ...
        'right wrist flexion/extension', ...
      };
    joint_directions = ...
      { ...
        'right', 'left'; ...
        'forward', 'backward'; ...
        'up', 'down'; ...
        'up', 'down'; ...
        'right', 'left'; ...
        'forward', 'backward'; ...
        'flexion', 'extension'; ...
        'abduction', 'adduction'; ...
        'internal', 'external'; ...
        'flexion', 'extension'; ...
        'external', 'internal'; ...
        'dorsi', 'plantar'; ...
        'eversion', 'inversion'; ...
        'flexion', 'extension'; ...
        'abduction', 'adduction'; ...
        'internal', 'external'; ...
        'flexion', 'extension'; ...
        'external', 'internal'; ...
        'dorsi', 'plantar'; ...
        'eversion', 'inversion'; ...
        'forward', 'backward'; ...
        'right', 'left'; ...
        'right', 'left'; ...
        'forward', 'backward'; ...
        'right', 'left'; ...
        'right', 'left'; ...
        'flexion', 'extension'; ...
        'abduction', 'adduction'; ...
        'internal', 'external'; ...
        'flexion', 'extension'; ...
        'pronation', 'supination'; ...
        'flexion', 'extension'; ...
        'flexion', 'extension'; ...
        'abduction', 'adduction'; ...
        'internal', 'external'; ...
        'flexion', 'extension'; ...
        'pronation', 'supination'; ...
        'flexion', 'extension'; ...
      }';

%     kinematic_tree.endEffectorLabels = ...
%       { ...
%         'left heel', ...
%         'left toes', ...
%         'left ankle', ...
%         'right heel', ...
%         'right toes', ...
%         'right ankle', ...
%         'head', ...
%         'left hand', ...
%         'right hand', ...
%         'pelvis' ...
%       };
    kinematic_tree.endEffectorLabels = ...
      { ...
        'left heel', ...
        'LTOESEEF', ...
        'left ankle', ...
        'right heel', ...
        'RTOESEEF', ...
        'right ankle', ...
        'head', ...
        'LHANDEEF', ...
        'RHANDEEF', ...
        'pelvis' ...
      };
    kinematic_tree.markerLabels = marker_labels;
%     kinematic_tree.addPointOfInterest('LTOESEEF', left_toe_mid, 26);
%     kinematic_tree.addPointOfInterest('RTOESEEF', right_toe_mid, 26);
%     kinematic_tree.addPointOfInterest('LHANDEEF', left_hand_mid, 26);
%     kinematic_tree.addPointOfInterest('RHANDEEF', right_hand_mid, 26);


    % add segment labels and joint centers
    kinematic_tree.addSegmentLabel('HEAD', 26);
    kinematic_tree.addSegmentLabel('TORSO', 23);
    kinematic_tree.addSegmentLabel('LUPPERARM', 29);
    kinematic_tree.addSegmentLabel('RUPPERARM', 35);
    kinematic_tree.addSegmentLabel('LFOREARM', 31);
    kinematic_tree.addSegmentLabel('RFOREARM', 37);
    kinematic_tree.addSegmentLabel('LHAND', 32);
    kinematic_tree.addSegmentLabel('RHAND', 38);
    kinematic_tree.addSegmentLabel('PELVIS', 6);
    kinematic_tree.addSegmentLabel('LTHIGH', 9);
    kinematic_tree.addSegmentLabel('RTHIGH', 16);
    kinematic_tree.addSegmentLabel('LSHANK', 11);
    kinematic_tree.addSegmentLabel('RSHANK', 18);
    kinematic_tree.addSegmentLabel('LFOOT', 13);
    kinematic_tree.addSegmentLabel('RFOOT', 20);
    
    kinematic_tree.addPointOfInterest('CERVIXCOR', cervix_cor, 26);
    kinematic_tree.addPointOfInterest('LSHOULDERCOR', left_shoulder_cor, 26);
    kinematic_tree.addPointOfInterest('RSHOULDERCOR', right_shoulder_cor, 26);
    kinematic_tree.addPointOfInterest('LELBOWCOR', left_elbow_cor, 26);
    kinematic_tree.addPointOfInterest('RELBOWCOR', right_elbow_cor, 26);
    kinematic_tree.addPointOfInterest('LWRISTCOR', left_wrist_cor, 26);
    kinematic_tree.addPointOfInterest('RWRISTCOR', right_wrist_cor, 26);
    kinematic_tree.addPointOfInterest('LUMBARCOR', lumbar_cor, 26);
    kinematic_tree.addPointOfInterest('LHIPCOR', left_hip_cor, 26);
    kinematic_tree.addPointOfInterest('RHIPCOR', right_hip_cor, 26);
    kinematic_tree.addPointOfInterest('LKNEECOR', left_knee_cor, 26);
    kinematic_tree.addPointOfInterest('RKNEECOR', right_knee_cor, 26);
    kinematic_tree.addPointOfInterest('LANKLECOR', left_ankle_cor, 26);
    kinematic_tree.addPointOfInterest('RANKLECOR', right_ankle_cor, 26);
    kinematic_tree.addPointOfInterest('LTOESEEF', left_toe_mid, 26);
    kinematic_tree.addPointOfInterest('RTOESEEF', right_toe_mid, 26);
    kinematic_tree.addPointOfInterest('LHANDEEF', left_hand_mid, 26);
    kinematic_tree.addPointOfInterest('RHANDEEF', right_hand_mid, 26);
    
    % define joint groups
    kinematic_tree.addJointGroup('pelvis', 1:6);
    kinematic_tree.addJointGroup('left leg', 7:13);
    kinematic_tree.addJointGroup('right leg', 14:20);
    kinematic_tree.addJointGroup('torso', 21:26);
    kinematic_tree.addJointGroup('left arm', 27:32);
    kinematic_tree.addJointGroup('right arm', 33:38);

    kinematic_tree.addJointGroup('left hip', 7:9);
    kinematic_tree.addJointGroup('left knee', 10:11);
    kinematic_tree.addJointGroup('left ankle', 12:13);
    kinematic_tree.addJointGroup('right hip', 14:16);
    kinematic_tree.addJointGroup('right knee', 17:18);
    kinematic_tree.addJointGroup('right ankle', 19:20);
    kinematic_tree.addJointGroup('lumbar', 21:23);
    kinematic_tree.addJointGroup('cervix', 24:26);
    kinematic_tree.addJointGroup('left shoulder', 27:29);
    kinematic_tree.addJointGroup('left elbow', 30:31);
    kinematic_tree.addJointGroup('left wrist', 32);
    kinematic_tree.addJointGroup('right shoulder', 33:35);
    kinematic_tree.addJointGroup('right elbow', 36:37);
    kinematic_tree.addJointGroup('right wrist', 38);

    
    kinematic_tree.addJointGroup('LHIPCOR', 7:9);
    kinematic_tree.addJointGroup('LKNEECOR', 10:11);
    kinematic_tree.addJointGroup('LANKLECOR', 12:13);
    kinematic_tree.addJointGroup('RHIPCOR', 14:16);
    kinematic_tree.addJointGroup('RKNEECOR', 17:18);
    kinematic_tree.addJointGroup('RANKLECOR', 19:20);
    kinematic_tree.addJointGroup('LUMBARCOR', 21:23);
    kinematic_tree.addJointGroup('CERVIXCOR', 24:26);
    kinematic_tree.addJointGroup('LSHOULDERCOR', 27:29);
    kinematic_tree.addJointGroup('LELBOWCOR', 30:31);
    kinematic_tree.addJointGroup('LWRISTCOR', 32);
    kinematic_tree.addJointGroup('RSHOULDERCOR', 33:35);
    kinematic_tree.addJointGroup('RELBOWCOR', 36:37);
    kinematic_tree.addJointGroup('RWRISTCOR', 38);

    % TODO: deal with cases where some of these don't exist or are misnamed
    
    % define markers
    marker_segment_list = createMarkerSegmentList(marker_labels);
    marker_color_list = createMarkerColorList(marker_labels);
    number_of_markers = length(marker_segment_list);
    for i_marker = 1 : number_of_markers
        kinematic_tree.addMarker(marker_segment_list(i_marker), marker_reference((i_marker-1)*3+1 : (i_marker-1)*3+3)', marker_color_list{i_marker});
    end

    kinematic_tree.updateInternals();

    %% save
    for i_segment = 1 : number_of_segments
        segment_com_wcs = [segment_coms_wcs{i_segment}; 1];
        T_wcs_to_mcs = mcs_to_wcs_transformations{i_segment}^(-1);
        segment_coms_mcs{i_segment} = eye(3, 4) * T_wcs_to_mcs * segment_com_wcs;
    end
    direction_matrices = ...
      { ...
        body_direction_matrix, ...
        left_hip_direction_matrix, ...
        right_hip_direction_matrix, ...
        left_knee_direction_matrix, ...
        right_knee_direction_matrix, ...
        left_ankle_direction_matrix, ...
        right_ankle_direction_matrix, ...
        lumbar_direction_matrix, ...
        cervix_direction_matrix, ...
        left_shoulder_direction_matrix, ...
        right_shoulder_direction_matrix, ...
        left_elbow_direction_matrix, ...
        right_elbow_direction_matrix, ...
        left_wrist_direction_matrix, ...
        right_wrist_direction_matrix, ...
      }; %#ok<NASGU>
    direction_matrix_labels = ...
      { ...
        'body', ...
        'left_hip', ...
        'right_hip', ...
        'left_knee', ...
        'right_knee', ...
        'left_ankle', ...
        'right_ankle', ...
        'lumbar', ...
        'cervix', ...
        'left_shoulder', ...
        'right_shoulder', ...
        'left_elbow', ...
        'right_elbow', ...
        'left_wrist', ...
        'right_wrist' ...
      }; %#ok<NASGU>
    
    save ...
      ( ...
        'subjectModel', ...
        'kinematic_tree', ...
        'joint_directions', ...
        'marker_labels', ...
        'marker_reference', ...
        'joint_center_labels', ...
        'joint_center_directions', ...
        'joint_center_reference', ...
        'segment_labels', ...
        'segment_masses', ...
        'segment_coms_wcs', ...
        'direction_matrices', ...
        'direction_matrix_labels' ...
      );
%         'segment_coms_mcs', ...
  
  % TODO: the segment_coms should be saved in world frame, not in marker frame. Then when I need it in marker frame, I
  % can calculate it. The way it currently is done, I have to keep track of which marker frame I have it represented in,
  % which is bad (mkay)

    %% show visualization

    if visualize
        % show stick figure of geometric model
        hip_center = (right_hip_cor + left_hip_cor) * 0.5;
        scene_bound = repmat(hip_center, 1, 2) + 2*[-0.5 0.5; -0.5 0.5; -1 1];
        stick_figure = KinematicTreeController(kinematic_tree, scene_bound, 'ellipsoid');
%         stick_figure = KinematicTreeController(kinematic_tree, scene_bound, 'none');

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

end


function marker_segments = createMarkerSegmentList(marker_labels)
    marker_to_segment_map = ...
      { ...
        'LFHD', 26; ...
        'RFHD', 26; ...
        'LBHD', 26; ...
        'RBHD', 26; ...
        'C7', 23; ...
        'T10', 23; ...
        'CLAV', 23; ...
        'STRN', 23; ...
        'RBAK', 23; ...
        'LSHO', 23; ...
        'LUPA', 29; ...
        'LELB', 29; ...
        'LFRM', 31; ...
        'LWRA', 31; ...
        'LWRB', 31; ...
        'LFIN', 32; ...
        'RSHO', 23; ...
        'RUPA', 35; ...
        'RELB', 35; ...
        'RFRM', 37; ...
        'RWRA', 37; ...
        'RWRB', 37; ...
        'RFIN', 38; ...
        'LASI', 6; ...
        'RASI', 6; ...
        'LPSI', 6; ...
        'RPSI', 6; ...
        'LTHI', 9; ...
        'LTHIA', 9; ...
        'LKNE', 9; ...
        'LTIB', 11; ...
        'LTIBA', 11; ...
        'LANK', 11; ...
        'LHEE', 13; ...
        'LTOE', 13; ...
        'LTOEL', 13; ...
        'RTHI', 16; ...
        'RTHIA', 16; ...
        'RKNE', 16; ...
        'RTIB', 18; ...
        'RTIBA', 18; ...
        'RANK', 18; ...
        'RHEE', 20; ...
        'RTOE', 20; ...
        'RTOEL', 20; ...
      }; % TODO: this should be loaded from the subject file or something
    number_of_markers = length(marker_labels) / 3;
    marker_segments = zeros(1, number_of_markers);
    for i_label = 1 : number_of_markers
        this_marker_label = marker_labels{i_label*3};
        this_marker_id = this_marker_label(1:end-2);
        marker_segments(i_label) = marker_to_segment_map{strcmp(marker_to_segment_map(:, 1), this_marker_id), 2};
    end
    
    
%     marker_segments(strcmp(marker_labels, 'LFHD')) = 26;
%     marker_segments(strcmp(marker_labels, 'RFHD')) = 26;
%     marker_segments(strcmp(marker_labels, 'LBHD')) = 26;
%     marker_segments(strcmp(marker_labels, 'RBHD')) = 26;
% 
%     marker_segments(strcmp(marker_labels, 'C7')) = 23;
%     marker_segments(strcmp(marker_labels, 'T10')) = 23;
%     marker_segments(strcmp(marker_labels, 'CLAV')) = 23;
%     marker_segments(strcmp(marker_labels, 'STRN')) = 23;
%     marker_segments(strcmp(marker_labels, 'RBAK')) = 23;
% 
%     marker_segments(strcmp(marker_labels, 'LSHO')) = 23;
%     marker_segments(strcmp(marker_labels, 'LUPA')) = 29;
%     marker_segments(strcmp(marker_labels, 'LELB')) = 29;
%     marker_segments(strcmp(marker_labels, 'LFRM')) = 31;
%     marker_segments(strcmp(marker_labels, 'LWRA')) = 31;
%     marker_segments(strcmp(marker_labels, 'LWRB')) = 31;
%     marker_segments(strcmp(marker_labels, 'LFIN')) = 32;
% 
%     marker_segments(strcmp(marker_labels, 'RSHO')) = 23;
%     marker_segments(strcmp(marker_labels, 'RUPA')) = 35;
%     marker_segments(strcmp(marker_labels, 'RELB')) = 35;
%     marker_segments(strcmp(marker_labels, 'RFRM')) = 37;
%     marker_segments(strcmp(marker_labels, 'RWRA')) = 37;
%     marker_segments(strcmp(marker_labels, 'RWRB')) = 37;
%     marker_segments(strcmp(marker_labels, 'RFIN')) = 38;
% 
%     marker_segments(strcmp(marker_labels, 'LASI')) = 6;
%     marker_segments(strcmp(marker_labels, 'RASI')) = 6;
%     marker_segments(strcmp(marker_labels, 'LPSI')) = 6;
%     marker_segments(strcmp(marker_labels, 'RPSI')) = 6;
% 
%     marker_segments(strcmp(marker_labels, 'LTHI')) = 9;
%     marker_segments(strcmp(marker_labels, 'LTHIA')) = 9;
%     marker_segments(strcmp(marker_labels, 'LKNE')) = 9;
%     marker_segments(strcmp(marker_labels, 'LTIB')) = 11;
%     marker_segments(strcmp(marker_labels, 'LTIBA')) = 11;
%     marker_segments(strcmp(marker_labels, 'LANK')) = 11;
%     marker_segments(strcmp(marker_labels, 'LHEE')) = 13;
%     marker_segments(strcmp(marker_labels, 'LTOE')) = 13;
%     marker_segments(strcmp(marker_labels, 'LTOEL')) = 13;
% 
%     marker_segments(strcmp(marker_labels, 'RTHI')) = 16;
%     marker_segments(strcmp(marker_labels, 'RTHIA')) = 16;
%     marker_segments(strcmp(marker_labels, 'RKNE')) = 16;
%     marker_segments(strcmp(marker_labels, 'RTIB')) = 18;
%     marker_segments(strcmp(marker_labels, 'RTIBA')) = 18;
%     marker_segments(strcmp(marker_labels, 'RANK')) = 18;
%     marker_segments(strcmp(marker_labels, 'RHEE')) = 20;
%     marker_segments(strcmp(marker_labels, 'RTOE')) = 20;
%     marker_segments(strcmp(marker_labels, 'RTOEL')) = 20;

end

function marker_color_list = createMarkerColorList(marker_labels)
    red = [1 0 0];
    green = [0 1 0];
    blue = [0 0 1];
    yellow = [1 0.9 0];
    
    marker_to_color_map = ...
      { ...
        'LFHD', red; ...
        'RFHD', green; ...
        'LBHD', red; ...
        'RBHD', green; ...
        'C7', blue; ...
        'T10', blue; ...
        'CLAV', blue; ...
        'STRN', blue; ...
        'RBAK', blue; ...
        'LSHO', red; ...
        'LUPA', red; ...
        'LELB', red; ...
        'LFRM', red; ...
        'LWRA', red; ...
        'LWRB', red; ...
        'LFIN', red; ...
        'RSHO', green; ...
        'RUPA', green; ...
        'RELB', green; ...
        'RFRM', green; ...
        'RWRA', green; ...
        'RWRB', green; ...
        'RFIN', green; ...
        'LASI', red; ...
        'RASI', green; ...
        'LPSI', red; ...
        'RPSI', green; ...
        'LTHI', red; ...
        'LTHIA', red; ...
        'LKNE', red; ...
        'LTIB', red; ...
        'LTIBA', red; ...
        'LANK', red; ...
        'LHEE', red; ...
        'LTOE', red; ...
        'LTOEL', red; ...
        'RTHI', green; ...
        'RTHIA', green; ...
        'RKNE', green; ...
        'RTIB', green; ...
        'RTIBA', green; ...
        'RANK', green; ...
        'RHEE', green; ...
        'RTOE', green; ...
        'RTOEL', green; ...
      }; % TODO: this should be loaded from the subject file or something
    number_of_markers = length(marker_labels) / 3;
    marker_color_array = zeros(number_of_markers, 3) * 0.3;
    for i_label = 1 : number_of_markers
        this_marker_label = marker_labels{i_label*3};
        this_marker_id = this_marker_label(1:end-2);
        marker_color_array(i_label, :) = marker_to_color_map{strcmp(marker_to_color_map(:, 1), this_marker_id), 2};
    end
  
%     if any(strcmp(marker_labels, 'LFHD')) marker_color_array(strcmp(marker_labels, 'LFHD'), :) = red; end
%     if any(strcmp(marker_labels, 'RFHD')) marker_color_array(strcmp(marker_labels, 'RFHD'), :) = green; end
%     if any(strcmp(marker_labels, 'LBHD')) marker_color_array(strcmp(marker_labels, 'LBHD'), :) = red; end
%     if any(strcmp(marker_labels, 'RBHD')) marker_color_array(strcmp(marker_labels, 'RBHD'), :) = green; end
% 
%     if any(strcmp(marker_labels, 'C7')) marker_color_array(strcmp(marker_labels, 'C7'), :) = blue; end
%     if any(strcmp(marker_labels, 'T10')) marker_color_array(strcmp(marker_labels, 'T10'), :) = blue; end
%     if any(strcmp(marker_labels, 'CLAV')) marker_color_array(strcmp(marker_labels, 'CLAV'), :) = blue; end
%     if any(strcmp(marker_labels, 'STRN')) marker_color_array(strcmp(marker_labels, 'STRN'), :) = blue; end
%     if any(strcmp(marker_labels, 'RBAK')) marker_color_array(strcmp(marker_labels, 'RBAK'), :) = blue; end
% 
%     if any(strcmp(marker_labels, 'LSHO')) marker_color_array(strcmp(marker_labels, 'LSHO'), :) = red; end
%     if any(strcmp(marker_labels, 'LUPA')) marker_color_array(strcmp(marker_labels, 'LUPA'), :) = red; end
%     if any(strcmp(marker_labels, 'LELB')) marker_color_array(strcmp(marker_labels, 'LELB'), :) = red; end
%     if any(strcmp(marker_labels, 'LFRA')) marker_color_array(strcmp(marker_labels, 'LFRA'), :) = red; end
%     if any(strcmp(marker_labels, 'LFRM')) marker_color_array(strcmp(marker_labels, 'LFRM'), :) = red; end
%     if any(strcmp(marker_labels, 'LWRA')) marker_color_array(strcmp(marker_labels, 'LWRA'), :) = red; end
%     if any(strcmp(marker_labels, 'LWRB')) marker_color_array(strcmp(marker_labels, 'LWRB'), :) = red; end
%     if any(strcmp(marker_labels, 'LFIN')) marker_color_array(strcmp(marker_labels, 'LFIN'), :) = red; end
% 
%     if any(strcmp(marker_labels, 'RSHO')) marker_color_array(strcmp(marker_labels, 'RSHO'), :) = green; end
%     if any(strcmp(marker_labels, 'RUPA')) marker_color_array(strcmp(marker_labels, 'RUPA'), :) = green; end
%     if any(strcmp(marker_labels, 'RELB')) marker_color_array(strcmp(marker_labels, 'RELB'), :) = green; end
%     if any(strcmp(marker_labels, 'RFRA')) marker_color_array(strcmp(marker_labels, 'RFRA'), :) = green; end
%     if any(strcmp(marker_labels, 'RFRM')) marker_color_array(strcmp(marker_labels, 'RFRM'), :) = green; end
%     if any(strcmp(marker_labels, 'RWRA')) marker_color_array(strcmp(marker_labels, 'RWRA'), :) = green; end
%     if any(strcmp(marker_labels, 'RWRB')) marker_color_array(strcmp(marker_labels, 'RWRB'), :) = green; end
%     if any(strcmp(marker_labels, 'RFIN')) marker_color_array(strcmp(marker_labels, 'RFIN'), :) = green; end
% 
%     if any(strcmp(marker_labels, 'LASI')) marker_color_array(strcmp(marker_labels, 'LASI'), :) = red; end
%     if any(strcmp(marker_labels, 'RASI')) marker_color_array(strcmp(marker_labels, 'RASI'), :) = green; end
%     if any(strcmp(marker_labels, 'LPSI')) marker_color_array(strcmp(marker_labels, 'LPSI'), :) = red; end
%     if any(strcmp(marker_labels, 'RPSI')) marker_color_array(strcmp(marker_labels, 'RPSI'), :) = green; end
% 
%     if any(strcmp(marker_labels, 'LTHI')) marker_color_array(strcmp(marker_labels, 'LTHI'), :) = red; end
%     if any(strcmp(marker_labels, 'LTHIA')) marker_color_array(strcmp(marker_labels, 'LTHIA'), :) = red; end
%     if any(strcmp(marker_labels, 'LKNE')) marker_color_array(strcmp(marker_labels, 'LKNE'), :) = red; end
%     if any(strcmp(marker_labels, 'LTIB')) marker_color_array(strcmp(marker_labels, 'LTIB'), :) = red; end
%     if any(strcmp(marker_labels, 'LTIBA')) marker_color_array(strcmp(marker_labels, 'LTIBA'), :) = red; end
%     if any(strcmp(marker_labels, 'LANK')) marker_color_array(strcmp(marker_labels, 'LANK'), :) = red; end
%     if any(strcmp(marker_labels, 'LHEE')) marker_color_array(strcmp(marker_labels, 'LHEE'), :) = red; end
%     if any(strcmp(marker_labels, 'LTOE')) marker_color_array(strcmp(marker_labels, 'LTOE'), :) = red; end
%     if any(strcmp(marker_labels, 'LTOEL')) marker_color_array(strcmp(marker_labels, 'LTOEL'), :) = red; end
% 
%     if any(strcmp(marker_labels, 'RTHI')) marker_color_array(strcmp(marker_labels, 'RTHI'), :) = green; end
%     if any(strcmp(marker_labels, 'RTHIA')) marker_color_array(strcmp(marker_labels, 'RTHIA'), :) = green; end
%     if any(strcmp(marker_labels, 'RKNE')) marker_color_array(strcmp(marker_labels, 'RKNE'), :) = green; end
%     if any(strcmp(marker_labels, 'RTIB')) marker_color_array(strcmp(marker_labels, 'RTIB'), :) = green; end
%     if any(strcmp(marker_labels, 'RTIBA')) marker_color_array(strcmp(marker_labels, 'RTIBA'), :) = green; end
%     if any(strcmp(marker_labels, 'RANK')) marker_color_array(strcmp(marker_labels, 'RANK'), :) = green; end
%     if any(strcmp(marker_labels, 'RHEE')) marker_color_array(strcmp(marker_labels, 'RHEE'), :) = green; end
%     if any(strcmp(marker_labels, 'RTOE')) marker_color_array(strcmp(marker_labels, 'RTOE'), :) = green; end
%     if any(strcmp(marker_labels, 'RTOEL')) marker_color_array(strcmp(marker_labels, 'RTOEL'), :) = green; end

    % export to list
    marker_color_list = cell(1, number_of_markers);
    for i_marker = 1 : number_of_markers
        marker_color_list{i_marker} = marker_color_array(i_marker, :);
    end

end















































