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
%             for i_time = 294 : 300
%             for i_time = 514 : 550

                % calculate ground truth for simulated data
                load('/Users/reimajbi/Box Sync/inverseKinematics/BRC/processed/00000000_XXX_simulation_001_markerTrajectories.mat')
                load('subjectModel.mat')
                knee_cor_reference = kinematic_tree.jointPositions{10};
                kinematic_tree.jointAngles = joint_angle_trajectories_simulated(i_time, :)';
                kinematic_tree.updateConfiguration;


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
                
                % check some stuff
%                 left_hip_cor_ground_truth = kinematic_tree.jointPositions{9}
%                 left_knee_cor_ground_truth = kinematic_tree.jointPositions{10}
%                 left_ankle_cor_ground_truth = kinematic_tree.jointPositions{13}
                left_toes_eef_ground_truth = kinematic_tree.endEffectorPositions{2};
                left_shoulder_cor_ground_truth = kinematic_tree.jointPositions{29};
                left_elbow_cor_ground_truth = kinematic_tree.jointPositions{31};
                left_wrist_cor_ground_truth = kinematic_tree.jointPositions{33};
                left_hand_eef_ground_truth = kinematic_tree.endEffectorPositions{8};
                right_hand_eef_ground_truth = kinematic_tree.endEffectorPositions{9};
                
                joint_angles_ground_truth = joint_angle_trajectories_simulated(i_time, :)';
                hip_flexion_ground_truth = joint_angles_ground_truth(7);
                knee_flexion_ground_truth = joint_angles_ground_truth(11);

%                 left_hip_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LHIPCOR')'
%                 left_knee_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LKNEECOR')'
%                 left_ankle_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LANKLECOR')'
%                 left_toes_eef_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LTOESEEF')'
%                 left_toes_eef_ground_truth
%                 left_shoulder_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LSHOULDERCOR')'
%                 left_shoulder_cor_ground_truth
%                 left_elbow_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LELBOWCOR')'
%                 left_elbow_cor_ground_truth
%                 left_wrist_cor_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LWRISTCOR')'
%                 left_wrist_cor_ground_truth
%                 left_hand_eef_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'LHANDEEF')'
%                 left_hand_eef_ground_truth

%                 right_hand_eef_current = extractMarkerTrajectories(joint_center_current, joint_center_headers, 'RHANDEEF')'
%                 right_hand_eef_ground_truth
                
                % TODO: this should be equal, but for some reason it's not
                % change the way that joint_centers_current is calculated, where exactly do the headers come from? How
                % is determined what this function returns?
                % this is not equal because the finger marker and the wrist markers are NOT on a rigid body, which is
                % what's implicitly assumed here. Use the finger marker directly instead
                
                % calculate segment centers of mass
                
                % TODO: this was using the segment_com_mcs as loaded from file. I now save the segment_com_wcs in reference, so I
                % have to transform that into segment_com_mcs first. The rationale for that is that I have to decide
                % only here which markers I use, and don't run into danger of confusing that.
                
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
%                 number_of_segments = length(segment_coms_mcs);
% %                 mcs_to_wcs_transformations = calculateMcsToWcsTransformations([marker_current joint_center_current], [marker_labels joint_center_headers], markers_by_segment);
%                 mcs_to_wcs_transformations = calculateMcsToWcsTransformations_detailed([marker_current joint_center_current], [marker_labels joint_center_headers], segment_labels);
%                 segment_coms_wcs = cell(number_of_segments, 1);
%                 for i_segment = 1 : number_of_segments
%                     segment_com_mcs = [segment_coms_mcs{i_segment}; 1];
%                     T_mcs_to_wcs = mcs_to_wcs_transformations{i_segment};
%                     segment_coms_wcs{i_segment} = eye(3, 4) * T_mcs_to_wcs * segment_com_mcs;
%                 end
                
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
%                         body_direction_matrix, ...
%                         left_hip_direction_matrix, ...
%                         right_hip_direction_matrix, ...
%                         left_knee_direction_matrix, ...
%                         right_knee_direction_matrix, ...
%                         left_ankle_direction_matrix, ...
%                         right_ankle_direction_matrix, ...
%                         lumbar_direction_matrix, ...
%                         cervix_direction_matrix ...
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                % calculate ground truth for simulated data
                load('/Users/reimajbi/Box Sync/inverseKinematics/BRC/processed/00000000_XXX_simulation_001_markerTrajectories.mat')
                load('subjectModel.mat')
                kinematic_tree.jointAngles = joint_angle_trajectories_simulated(i_time, :)';
                kinematic_tree.updateConfiguration;
                pelvis_transformation_from_simulation = kinematic_tree.productsOfExponentials{6};
                lthigh_transformation_from_simulation = kinematic_tree.productsOfExponentials{9};
                exp1 = kinematic_tree.twistExponentials{1};
                exp2 = kinematic_tree.twistExponentials{2};
                exp3 = kinematic_tree.twistExponentials{3};
                
                pelvis_transformation_current_cell = calculateMcsToWcsTransformations_detailed(marker_current, marker_labels, {'PELVIS'});    
                pelvis_transformation_current = pelvis_transformation_current_cell{1};
                pelvis_rotation_current = pelvis_transformation_current(1:3, 1:3);
                pelvis_translation_current = pelvis_transformation_current(1:3, 4);

                pelvis_transformation_reference_cell = calculateMcsToWcsTransformations_detailed(marker_reference, marker_labels, {'PELVIS'});
                pelvis_transformation_reference = pelvis_transformation_reference_cell{1};
                pelvis_rotation_reference = pelvis_transformation_reference(1:3, 1:3);
                pelvis_translation_reference = pelvis_transformation_reference(1:3, 4);

%                 pelvis_transformation_reference_to_current = pelvis_transformation_reference^(-1) * pelvis_transformation_current;
%                 pelvis_rotation_reference_to_current = pelvis_transformation_reference_to_current(1:3, 1:3)
                pelvis_translation_reference_to_current = pelvis_translation_current - pelvis_translation_reference;
                
                pelvis_joint_transformation = pelvis_transformation_current * pelvis_transformation_reference^(-1);
                pelvis_joint_rotation = pelvis_joint_transformation(1:3, 1:3); % this is equal to the rotation in kinematic_tree.productsOfExponentials{6}
                pelvis_rotation_angles = eulerAnglesFromRotationMatrixZXY(pelvis_joint_rotation); % these are the correct values

                joint_angles(1:3) = pelvis_translation_reference_to_current;

                pelvis_euler_angles_current = eulerAnglesFromRotationMatrixZXY(pelvis_rotation_current);
                pelvis_euler_angles_reference = eulerAnglesFromRotationMatrixZXY(pelvis_rotation_reference);
                pelvis_euler_angles = pelvis_euler_angles_current - pelvis_euler_angles_reference;
                joint_angles(4:6) = pelvis_euler_angles;

                R_4 = expAxis([0; 0; 1], joint_angles(4));
                R_5 = expAxis([1; 0; 0], joint_angles(5));
                R_6 = expAxis([0; 1; 0], joint_angles(6));
                R_pelvis = R_4 * R_5 * R_6;
                
                % thigh
                lthigh_transformation_current_cell = calculateMcsToWcsTransformations_detailed([marker_current joint_center_current], [marker_labels joint_center_labels], {'LTHIGH'});    
                lthigh_transformation_current = lthigh_transformation_current_cell{1};
                lthigh_rotation_current = lthigh_transformation_current(1:3, 1:3);
                lthigh_translation_current = lthigh_transformation_current(1:3, 4);

                lthigh_transformation_reference_cell = calculateMcsToWcsTransformations_detailed([marker_reference joint_center_reference], [marker_labels joint_center_labels], {'LTHIGH'});
                lthigh_transformation_reference = lthigh_transformation_reference_cell{1};
                lthigh_rotation_reference = lthigh_transformation_reference(1:3, 1:3);
                lthigh_translation_reference = lthigh_transformation_reference(1:3, 4);

%                 lthigh_transformation_reference_to_current = lthigh_transformation_reference^(-1) * lthigh_transformation_current;
%                 lthigh_rotation_reference_to_current = lthigh_transformation_reference_to_current(1:3, 1:3)
                lthigh_translation_reference_to_current = lthigh_translation_current - lthigh_translation_reference;
                
%                 lthigh_joint_transformation = lthigh_transformation_current * lthigh_transformation_reference^(-1);
                lthigh_joint_rotation = R_pelvis^(-1) * lthigh_rotation_current * lthigh_rotation_reference^(-1);
%                 lthigh_joint_rotation = lthigh_joint_transformation(1:3, 1:3) % this is equal to the rotation in kinematic_tree.productsOfExponentials{9} ? it's not
                lthigh_rotation_angles = eulerAnglesFromRotationMatrixZXY(lthigh_joint_rotation); % these should be the correct values
                
                
                % test equations
%                 T_cur = pelvis_transformation_current
%                 T_ref = pelvis_transformation_reference;
%                 exp1_6 = kinematic_tree.productsOfExponentials{6};
%                 T_cur_check = exp1_6 * T_ref
%                 
%                 T_cur = lthigh_transformation_current
%                 T_ref = lthigh_transformation_reference;
%                 exp1_9 = kinematic_tree.productsOfExponentials{9};
%                 T_cur_check = exp1_9 * T_ref
                
%                 R_cur = lthigh_rotation_current
%                 R_ref = lthigh_rotation_reference;
%                 exp1_9 = kinematic_tree.productsOfExponentials{9};
%                 R_cur_check = exp1_9(1:3, 1:3) * R_ref
                
                joint_angles_simulated_current = joint_angle_trajectories_simulated(i_time, :);
                
                
                
                
                
                
                
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
%             save ...
%               ( ...
%                 [save_folder filesep save_file_name], ...
%                 'joint_center_trajectories', ...
%                 'joint_center_labels', ...
%                 'joint_angle_trajectories', ...
%                 'com_trajectories', ...
%                 'com_labels', ...
%                 'time_mocap', ...
%                 'sampling_rate_mocap' ...
%               );
            disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' save_folder filesep save_file_name]);

            addAvailableData('joint_center_trajectories', 'time_mocap', 'sampling_rate_mocap', 'joint_center_labels', save_folder, save_file_name);
            addAvailableData('com_trajectories', 'time_mocap', 'sampling_rate_mocap', 'com_labels', save_folder, save_file_name);
            addAvailableData('joint_angle_trajectories', 'time_mocap', 'sampling_rate_mocap', 'joint_angle_trajectories', save_folder, save_file_name);
        end
    end
end