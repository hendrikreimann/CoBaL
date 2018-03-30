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


function inverseKinematics_3DoF(varargin)
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
    REAR_reference = extractMarkerData(marker_reference, marker_labels, 'R_ear')';
    LEAR_reference = extractMarkerData(marker_reference, marker_labels, 'L_ear')';
    LSHO_reference = extractMarkerData(marker_reference, marker_labels, 'L_shoulder')';
    RSHO_reference = extractMarkerData(marker_reference, marker_labels, 'R_shoulder')';
    LIC_reference = extractMarkerData(marker_reference, marker_labels, 'L_ic')';
    RIC_reference = extractMarkerData(marker_reference, marker_labels, 'R_ic')';
    LGT_reference = extractMarkerData(marker_reference, marker_labels, 'L_gt')';
    LKNE_reference = extractMarkerData(marker_reference, marker_labels, 'L_knee')';
    LANK_reference = extractMarkerData(marker_reference, marker_labels, 'L_ankle')';
    RGT_reference = extractMarkerData(marker_reference, marker_labels, 'R_gt')';
    RKNE_reference = extractMarkerData(marker_reference, marker_labels, 'R_knee')';
    RANK_reference = extractMarkerData(marker_reference, marker_labels, 'R_ankle')';
    
    % calculate segment angles for reference
    ankle_reference = mean([LANK_reference RANK_reference], 2);
    knee_reference = mean([LKNE_reference RKNE_reference], 2);
    hip_reference = mean([LGT_reference RGT_reference], 2);
    shoulder_reference = mean([LSHO_reference RSHO_reference], 2);
    
    lower_leg_vector_reference = knee_reference - ankle_reference;
    thigh_vector_reference = hip_reference - knee_reference;
    trunk_vector_reference = shoulder_reference - hip_reference;
    
    lower_leg_segment_angle_reference = atan2(lower_leg_vector_reference(3), -lower_leg_vector_reference(1));
    thigh_segment_angle_reference = atan2(thigh_vector_reference(3), -thigh_vector_reference(1));
    trunk_segment_angle_reference = atan2(trunk_vector_reference(3), -trunk_vector_reference(1));
    
    lower_leg_length = norm(lower_leg_vector_reference);
    thigh_length = norm(thigh_vector_reference);
    trunk_length = norm(trunk_vector_reference);
    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            condition = condition_list{i_condition};
            
            % give feedback
%             disp([' - Condition ' condition ', Trial ' num2str(i_trial)])
%             fprintf([' - Calculating kinematic trajectories... \n'])
            
            % load data
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);
            number_of_time_steps = size(marker_trajectories, 1);
            time_steps_to_process = 1 : number_of_time_steps;
            number_of_time_steps_to_process = length(time_steps_to_process);
            
            com_labels = [segment_labels 'BODY'];
            for i_label = 1 : length(com_labels)
                com_labels{i_label} = [com_labels{i_label} 'COM'];
            end
            
            % calculate
            joint_center_labels = joint_center_labels; % TODO: fix this
            com_trajectories = zeros(number_of_time_steps, length(com_labels)*3);
            
            % extract marker positions
            REAR_position = extractMarkerData(marker_trajectories, marker_labels, 'R_ear');
            LEAR_position = extractMarkerData(marker_trajectories, marker_labels, 'L_ear');
            LSHO_position = extractMarkerData(marker_trajectories, marker_labels, 'L_shoulder');
            RSHO_position = extractMarkerData(marker_trajectories, marker_labels, 'R_shoulder');
            LIC_position = extractMarkerData(marker_trajectories, marker_labels, 'L_ic');
            RIC_position = extractMarkerData(marker_trajectories, marker_labels, 'R_ic');
            LGT_position = extractMarkerData(marker_trajectories, marker_labels, 'L_gt');
            LKNE_position = extractMarkerData(marker_trajectories, marker_labels, 'L_knee');
            LANK_position = extractMarkerData(marker_trajectories, marker_labels, 'L_ankle');
            RGT_position = extractMarkerData(marker_trajectories, marker_labels, 'R_gt');
            RKNE_position = extractMarkerData(marker_trajectories, marker_labels, 'R_knee');
            RANK_position = extractMarkerData(marker_trajectories, marker_labels, 'R_ankle');

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
            trunk_segment_angle_position = atan2(trunk_vector(:, 3), -trunk_vector(:, 1));

            lower_leg_segment_angle_change = lower_leg_segment_angle_position - lower_leg_segment_angle_reference;
            thigh_segment_angle_change = thigh_segment_angle_position - thigh_segment_angle_reference;
            trunk_segment_angle_change = trunk_segment_angle_position - trunk_segment_angle_reference;

            ankle_joint_angle = lower_leg_segment_angle_change;
            knee_joint_angle = thigh_segment_angle_change - lower_leg_segment_angle_change;
            hip_joint_angle = trunk_segment_angle_change - thigh_segment_angle_change;
            
            joint_center_trajectories = [ankle_cor knee_cor hip_cor shoulders_mid];
            joint_angle_trajectories = [ankle_joint_angle knee_joint_angle hip_joint_angle];
            
            %% save
            variables_to_save = struct;
            
            variables_to_save.joint_center_trajectories = joint_center_trajectories;
            variables_to_save.joint_center_labels = joint_center_labels;
            variables_to_save.joint_center_directions = joint_center_directions;
            
            variables_to_save.joint_angle_trajectories = joint_angle_trajectories;
            variables_to_save.joint_labels = kinematic_tree.jointLabels;
            variables_to_save.joint_directions = joint_directions;
            
            variables_to_save.time_mocap = time_mocap;
            variables_to_save.sampling_rate_mocap = sampling_rate_mocap;
            
            save_folder = 'processed';
            save_file_name = makeFileName(date, subject_id, condition, i_trial, 'kinematicTrajectories.mat');
            saveDataToFile([save_folder filesep save_file_name], variables_to_save);

            addAvailableData_new ...
              ( ...
                'joint_center_trajectories', ...
                'time_mocap', ...
                'sampling_rate_mocap', ...
                '_joint_center_labels', ...
                '_joint_center_directions', ...
                save_folder, ...
                save_file_name ...
              );
            addAvailableData_new ...
              ( ...
                'joint_angle_trajectories', ...
                'time_mocap', ...
                'sampling_rate_mocap', ...
                '_joint_labels', ...
                '_joint_directions', ...
                save_folder, ...
                save_file_name ...
              );            
            
            disp(['Calculating kinematic variables, condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' save_folder filesep save_file_name]);
        end
    end
end