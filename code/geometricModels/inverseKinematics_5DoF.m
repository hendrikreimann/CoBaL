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

% this function calculates kinematic variables from the marker data

% input: 
% subjectInfo.mat
% subjectModel.mat
% markerTrajectories
%
% output:
% file kinematicTrajectories.mat, containing
% - joint_center_trajectories
% - com_trajectories
% - com_labels
% - joint_angle_trajectories


function inverseKinematics_stance(varargin)
    % parse arguments
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'use_parallel', false)
    parse(parser, varargin{:})
    use_parallel = parser.Results.use_parallel;
    
    load('subjectInfo.mat', 'date', 'subject_id');
    load('subjectModel.mat');
    study_settings_file = '';
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    subject_settings = SettingsCustodian('subjectSettings.txt');
    
    % extract references
    number_of_joint_angles = kinematic_tree.numberOfJoints;
    REAR_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'R_ear')';
    LEAR_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'L_ear')';
    LSHO_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'L_shoulder')';
    RSHO_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'R_shoulder')';
    LIC_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'L_ic')';
    RIC_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'R_ic')';
    LGT_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'L_gt')';
    LKNE_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'L_knee')';
    LANK_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'L_ankle')';
    RGT_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'R_gt')';
    RKNE_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'R_knee')';
    RANK_reference = extractMarkerTrajectories(marker_reference, marker_labels, 'R_ankle')';
    
    % calculate segment angles for reference
    ankle_reference = mean([LANK_reference RANK_reference], 2);
    knee_reference = mean([LKNE_reference RKNE_reference], 2);
    hip_reference = mean([LGT_reference RGT_reference], 2);
    lumbar_reference = mean([LIC_reference RIC_reference], 2);
    neck_reference = mean([LSHO_reference RSHO_reference], 2);
    head_reference = mean([LEAR_reference REAR_reference], 2);
    
    lower_leg_vector_reference = knee_reference - ankle_reference;
    thigh_vector_reference = hip_reference - knee_reference;
    pelvis_vector_reference = lumbar_reference - hip_reference;
    trunk_vector_reference = neck_reference - lumbar_reference;
    head_vector_reference = head_reference - neck_reference;
    
    lower_leg_segment_angle_reference = atan2(lower_leg_vector_reference(3), -lower_leg_vector_reference(1));
    thigh_segment_angle_reference = atan2(thigh_vector_reference(3), -thigh_vector_reference(1));
    pelvis_segment_angle_reference = atan2(pelvis_vector_reference(3), -pelvis_vector_reference(1));
    trunk_segment_angle_reference = atan2(trunk_vector_reference(3), -trunk_vector_reference(1));
    head_segment_angle_reference = atan2(head_vector_reference(3), -head_vector_reference(1));
    
    lower_leg_length = norm(lower_leg_vector_reference);
    thigh_length = norm(thigh_vector_reference);
    pelvis_length = norm(pelvis_vector_reference);
    trunk_length = norm(trunk_vector_reference);
    head_length = norm(head_vector_reference);
    
    if use_parallel
        % get or open pool of workers
        poolobject = gcp;
        number_of_labs = poolobject.NumWorkers;
    end
    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            condition = condition_list{i_condition};
            
            % give feedback
%             disp([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Condition ' condition ', Trial ' num2str(i_trial)])
%             fprintf([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Calculating kinematic trajectories... \n'])
%             disp([' - Condition ' condition ', Trial ' num2str(i_trial)])
%             fprintf([' - Calculating kinematic trajectories... \n'])
            
            % load data
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);
            number_of_time_steps = size(marker_trajectories, 1);
            time_steps_to_process = 1 : number_of_time_steps;
%             time_steps_to_process = 29999 : 30000;
            number_of_time_steps_to_process = length(time_steps_to_process);
            
            com_labels = [segment_labels 'BODY'];
            for i_label = 1 : length(com_labels)
                com_labels{i_label} = [com_labels{i_label} 'COM'];
            end
            
            % determine indices for optional markers
            optional_marker_indices = [];
            optional_marker_list = study_settings.get('optional_markers');
            for i_marker = 1 : length(optional_marker_list)
                marker = find(strcmp(marker_labels, optional_marker_list{i_marker}));
                marker_indices = reshape([(marker - 1) * 3 + 1; (marker - 1) * 3 + 2; (marker - 1) * 3 + 3], 1, length(marker)*3);
                optional_marker_indices = [optional_marker_indices marker_indices];
            end
            essential_marker_indicator = ~ismember(1 : size(marker_trajectories, 2), optional_marker_indices);
            
            % calculate
            joint_center_trajectories = zeros(number_of_time_steps, length(joint_center_headers)*3);
            joint_center_labels = joint_center_headers; % TODO: fix this
            com_trajectories = zeros(number_of_time_steps, length(com_labels)*3);
            joint_angle_trajectories = zeros(number_of_time_steps, number_of_joint_angles);
            
            % extract marker positions
            REAR_position = extractMarkerTrajectories(marker_trajectories, marker_labels, 'R_ear');
            LEAR_position = extractMarkerTrajectories(marker_trajectories, marker_labels, 'L_ear');
            LSHO_position = extractMarkerTrajectories(marker_trajectories, marker_labels, 'L_shoulder');
            RSHO_position = extractMarkerTrajectories(marker_trajectories, marker_labels, 'R_shoulder');
            LIC_position = extractMarkerTrajectories(marker_trajectories, marker_labels, 'L_ic');
            RIC_position = extractMarkerTrajectories(marker_trajectories, marker_labels, 'R_ic');
            LGT_position = extractMarkerTrajectories(marker_trajectories, marker_labels, 'L_gt');
            LKNE_position = extractMarkerTrajectories(marker_trajectories, marker_labels, 'L_knee');
            LANK_position = extractMarkerTrajectories(marker_trajectories, marker_labels, 'L_ankle');
            RGT_position = extractMarkerTrajectories(marker_trajectories, marker_labels, 'R_gt');
            RKNE_position = extractMarkerTrajectories(marker_trajectories, marker_labels, 'R_knee');
            RANK_position = extractMarkerTrajectories(marker_trajectories, marker_labels, 'R_ankle');

            % calculate segment angles for reference
            ankle_cor = (LANK_position + RANK_position) * 0.5;
            knee_cor = (LKNE_position + RKNE_position) * 0.5;
            hip_cor = (LGT_position + RGT_position) * 0.5;
            lumbar_cor = (LIC_position + RIC_position) * 0.5;
            shoulders_mid = (LSHO_position + RSHO_position) * 0.5;
            ears_mid = (LEAR_position + REAR_position) * 0.5;

            lower_leg_vector = knee_cor - ankle_cor;
            thigh_vector = hip_cor - knee_cor;
            pelvis_vector = lumbar_cor - hip_cor;
            trunk_vector = shoulders_mid - lumbar_cor;
            head_vector = ears_mid - shoulders_mid;

            lower_leg_segment_angle_position = atan2(lower_leg_vector(:, 3), -lower_leg_vector(:, 1));
            thigh_segment_angle_position = atan2(thigh_vector(:, 3), -thigh_vector(:, 1));
            pelvis_segment_angle_position = atan2(pelvis_vector(:, 3), -pelvis_vector(:, 1));
            trunk_segment_angle_position = atan2(trunk_vector(:, 3), -trunk_vector(:, 1));
            head_segment_angle_position = atan2(head_vector(:, 3), -head_vector(:, 1));

            lower_leg_segment_angle_change = lower_leg_segment_angle_position - lower_leg_segment_angle_reference;
            thigh_segment_angle_change = thigh_segment_angle_position - thigh_segment_angle_reference;
            pelvis_segment_angle_change = pelvis_segment_angle_position - pelvis_segment_angle_reference;
            trunk_segment_angle_change = trunk_segment_angle_position - trunk_segment_angle_reference;
            head_segment_angle_change = head_segment_angle_position - head_segment_angle_reference;

            ankle_joint_angle = lower_leg_segment_angle_change;
            knee_joint_angle = thigh_segment_angle_change - lower_leg_segment_angle_change;
            hip_joint_angle = pelvis_segment_angle_change - thigh_segment_angle_change;
            lumbar_joint_angle = trunk_segment_angle_change - pelvis_segment_angle_change;
            neck_joint_angle = head_segment_angle_change - trunk_segment_angle_change;
            
            joint_center_trajectories = [ankle_cor knee_cor hip_cor lumbar_cor shoulders_mid];
            joint_angle_trajectories = [ankle_joint_angle knee_joint_angle hip_joint_angle lumbar_joint_angle neck_joint_angle];
            

%             % calculate segment centers of mass
%             segment_coms_wcs_reference = segment_coms_wcs;
%             number_of_segments = length(segment_coms_wcs_reference);
%             mcs_to_wcs_transformations_reference = calculateMcsToWcsTransformations_detailed([marker_reference joint_center_reference], [marker_labels joint_center_headers], segment_labels);
%             mcs_to_wcs_transformations_current = calculateMcsToWcsTransformations_detailed([marker_current joint_center_current], [marker_labels joint_center_headers], segment_labels);
%             segment_coms_wcs_current = cell(number_of_segments, 1);
%             for i_segment = 1 : number_of_segments
%                 segment_com_wcs_reference = [segment_coms_wcs_reference{i_segment}; 1];
%                 T_mcs_to_wcs_reference = mcs_to_wcs_transformations_reference{i_segment};
%                 T_mcs_to_wcs_current = mcs_to_wcs_transformations_current{i_segment};
%                 segment_com_mcs = T_mcs_to_wcs_reference^(-1) * segment_com_wcs_reference;
%                 segment_coms_wcs_current{i_segment} = eye(3, 4) * T_mcs_to_wcs_current * segment_com_mcs;
%             end
% 
%             % calculate whole body center of mass
%             body_com = [0; 0; 0];
%             for i_segment = 1 : number_of_segments
%                 body_com = body_com + segment_masses(i_segment) * segment_coms_wcs_current{i_segment};
%             end
%             body_com = body_com * 1 / sum(segment_masses);
% 
%             % export centers of mass
%             for i_segment = 1 : number_of_segments
%                 com_trajectories(i_time, (i_segment-1)*3+1 : (i_segment-1)*3+3) = segment_coms_wcs_current{i_segment};
%             end
%             com_trajectories(i_time, number_of_segments*3+1 : number_of_segments*3+3) = body_com;
% 
%             % calculate joint angles
%             joint_angle_trajectories(i_time, :) = ...
%                 markerToAngles ...
%                   ( ...
%                     marker_reference, ...
%                     marker_current, ...
%                     marker_labels, ...
%                     joint_center_reference, ...
%                     joint_center_current, ...
%                     joint_center_headers, ...
%                     direction_matrices, ...
%                     direction_matrix_labels ...
%                   );
% 
%             % give progress feedback
%             display_step = 1;
%             if (i_time_step / display_step) == floor(i_time_step / display_step)
%                 disp([num2str(i_time_step) '(' num2str(length(time_steps_to_process)) ')']);
%             end                        
            
%             fprintf([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Calculating kinematic trajectories... done\n'])
%             fprintf([' - Calculating kinematic trajectories... done\n'])

            
            % save
            variables_to_save = struct;
            variables_to_save.joint_center_trajectories = joint_center_trajectories;
            variables_to_save.joint_center_labels = joint_center_labels;
            variables_to_save.joint_angle_trajectories = joint_angle_trajectories;
            variables_to_save.joint_labels = kinematic_tree.jointLabels;
            variables_to_save.com_trajectories = com_trajectories;
            variables_to_save.com_labels = com_labels;
            variables_to_save.time_mocap = time_mocap;
            variables_to_save.sampling_rate_mocap = sampling_rate_mocap;
            
            save_folder = 'processed';
            save_file_name = makeFileName(date, subject_id, condition, i_trial, 'kinematicTrajectories.mat');
            saveDataToFile([save_folder filesep save_file_name], variables_to_save);
            disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' save_folder filesep save_file_name]);

            addAvailableData('joint_center_trajectories', 'time_mocap', 'sampling_rate_mocap', 'joint_center_labels', save_folder, save_file_name);
            addAvailableData('com_trajectories', 'time_mocap', 'sampling_rate_mocap', 'com_labels', save_folder, save_file_name);
            addAvailableData('joint_angle_trajectories', 'time_mocap', 'sampling_rate_mocap', 'joint_labels', save_folder, save_file_name);
        end
    end
end