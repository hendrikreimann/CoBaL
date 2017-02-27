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


function calculateKinematicTrajectories(varargin)
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    
    load('subjectInfo.mat', 'date', 'subject_id');
    load('subjectModel.mat');
    
    number_of_joint_angles = kinematic_tree.numberOfJoints;
    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            % load data
            condition = condition_list{i_condition};
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);
            
            number_of_time_steps = size(marker_trajectories, 1);

            com_labels = [segment_labels 'BODY'];
            for i_label = 1 : length(com_labels)
                com_labels{i_label} = [com_labels{i_label} 'COM'];
            end
            
            % calculate
            joint_center_trajectories = zeros(number_of_time_steps, length(joint_center_headers)*3);
            joint_center_labels = joint_center_headers; % TODO: fix this
            com_trajectories = zeros(number_of_time_steps, length(com_labels)*3);
            joint_angle_trajectories = zeros(number_of_time_steps, number_of_joint_angles);
            for i_time = 1 : number_of_time_steps

                % calculate joint center positions
                marker_current = marker_trajectories(i_time, :);
                joint_center_current = ...
                    calculateJointCenterPositions ...
                      ( ...
                        marker_reference, ...
                        marker_current, ...
                        marker_labels, ...
                        joint_center_reference, ...
                        joint_center_headers ...
                      );
                joint_center_trajectories(i_time, :) = joint_center_current;
                
                % calculate segment centers of mass
                segment_coms_wcs_reference = segment_coms_wcs;
                number_of_segments = length(segment_coms_wcs_reference);
                mcs_to_wcs_transformations_reference = calculateMcsToWcsTransformations_detailed([marker_reference joint_center_reference], [marker_labels joint_center_headers], segment_labels);
                mcs_to_wcs_transformations_current = calculateMcsToWcsTransformations_detailed([marker_current joint_center_current], [marker_labels joint_center_headers], segment_labels);
                segment_coms_wcs_current = cell(number_of_segments, 1);
                for i_segment = 1 : number_of_segments
                    segment_com_wcs_reference = [segment_coms_wcs_reference{i_segment}; 1];
                    T_mcs_to_wcs_reference = mcs_to_wcs_transformations_reference{i_segment};
                    T_mcs_to_wcs_current = mcs_to_wcs_transformations_current{i_segment};
                    segment_com_mcs = T_mcs_to_wcs_reference^(-1) * segment_com_wcs_reference;
                    segment_coms_wcs_current{i_segment} = eye(3, 4) * T_mcs_to_wcs_current * segment_com_mcs;
                end
                
                % calculate whole body center of mass
                body_com = [0; 0; 0];
                for i_segment = 1 : number_of_segments
                    body_com = body_com + segment_masses(i_segment) * segment_coms_wcs_current{i_segment};
                end
                body_com = body_com * 1 / sum(segment_masses);
                
                % export centers of mass
                for i_segment = 1 : number_of_segments
                    com_trajectories(i_time, (i_segment-1)*3+1 : (i_segment-1)*3+3) = segment_coms_wcs_current{i_segment};
                end
                com_trajectories(i_time, number_of_segments*3+1 : number_of_segments*3+3) = body_com;
                
                % calculate joint angles
                joint_angle_trajectories(i_time, :) = ...
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
                      );
            end
            
            % save
            variables_to_save = struct;
            variables_to_save.joint_center_trajectories = joint_center_trajectories;
            variables_to_save.joint_center_labels = joint_center_labels;
            variables_to_save.joint_angle_trajectories = joint_angle_trajectories;
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
            addAvailableData('joint_angle_trajectories', 'time_mocap', 'sampling_rate_mocap', 'joint_angle_trajectories', save_folder, save_file_name);
        end
    end
end