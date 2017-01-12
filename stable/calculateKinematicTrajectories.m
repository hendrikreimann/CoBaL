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
            com_trajectories = zeros(number_of_time_steps, length(com_labels)*3);
            joint_angle_trajectories = zeros(number_of_time_steps, number_of_joint_angles);
            for i_time = 1 : number_of_time_steps
%             for i_time = 294 : 300
%             for i_time = 514 : 550
                % calculate joint center positions
                marker_current = marker_trajectories(i_time, :);
                joint_center_current = ...
                    calculateJointCenterPositions ...
                      ( ...
                        marker_reference, ...
                        marker_current, ...
                        marker_headers, ...
                        joint_center_reference, ...
                        joint_center_headers ...
                      );
                joint_center_trajectories(i_time, :) = joint_center_current;
                
                % calculate segment centers of mass
                number_of_segments = length(segment_coms_mcs);
%                 mcs_to_wcs_transformations = calculateMcsToWcsTransformations([marker_current joint_center_current], [marker_headers joint_center_headers], markers_by_segment);
                mcs_to_wcs_transformations = calculateMcsToWcsTransformations_detailed([marker_current joint_center_current], [marker_headers joint_center_headers], segment_labels);
                segment_coms_wcs = cell(number_of_segments, 1);
                for i_segment = 1 : number_of_segments
                    segment_com_mcs = [segment_coms_mcs{i_segment}; 1];
                    T_mcs_to_wcs = mcs_to_wcs_transformations{i_segment};
                    segment_coms_wcs{i_segment} = eye(3, 4) * T_mcs_to_wcs * segment_com_mcs;
                end
                
                % calculate whole body center of mass
                body_com = [0; 0; 0];
                for i_segment = 1 : number_of_segments
                    body_com = body_com + segment_masses(i_segment) * segment_coms_wcs{i_segment};
                end
                body_com = body_com * 1 / sum(segment_masses);
                
                % export centers of mass
                for i_segment = 1 : number_of_segments
                    com_trajectories(i_time, (i_segment-1)*3+1 : (i_segment-1)*3+3) = segment_coms_wcs{i_segment};
                end
                com_trajectories(i_time, number_of_segments*3+1 : number_of_segments*3+3) = body_com;
                
                % calculate joint angles
%                 joint_angle_trajectories(i_time, :) = ...
%                     markerToAngles ...
%                       ( ...
%                         marker_reference, ...
%                         marker_current, ...
%                         marker_headers, ...
%                         joint_center_reference, ...
%                         joint_center_current, ...
%                         joint_center_headers, ...
%                         body_direction_matrix, ...
%                         left_hip_direction_matrix, ...
%                         right_hip_direction_matrix, ...
%                         left_knee_direction_matrix, ...
%                         right_knee_direction_matrix, ...
%                         left_ankle_direction_matrix, ...
%                         right_ankle_direction_matrix ...
%                       );
            end
            
            com_file_name = ['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'kinematicTrajectories')];
            save ...
              ( ...
                com_file_name, ...
                'joint_center_trajectories', ...
                'com_trajectories', ...
                'com_labels', ...
                'joint_angle_trajectories' ...
              );
            disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' com_file_name]);

        end
    end
end