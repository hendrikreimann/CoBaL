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
    % parse arguments
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'use_parallel', false)
    parse(parser, varargin{:})
    use_parallel = parser.Results.use_parallel;
    
    load('subjectInfo.mat', 'date', 'subject_id');
    load('subjectModel.mat');
    
    number_of_joint_angles = kinematic_tree.numberOfJoints;
    
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
            disp([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Condition ' condition ', Trial ' num2str(i_trial)])
            fprintf([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Calculating kinematic trajectories... \n'])
            
            % load data
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);
            number_of_time_steps = size(marker_trajectories, 1); %#ok<NODEF>
            time_steps_to_process = 1 : number_of_time_steps;
            com_labels = [segment_labels 'BODY'];
            for i_label = 1 : length(com_labels)
                com_labels{i_label} = [com_labels{i_label} 'COM'];
            end
            
            % calculate
            joint_center_trajectories = zeros(number_of_time_steps, length(joint_center_headers)*3);
            joint_center_labels = joint_center_headers; % TODO: fix this
            com_trajectories = zeros(number_of_time_steps, length(com_labels)*3);
            joint_angle_trajectories = zeros(number_of_time_steps, number_of_joint_angles);
            
            if use_parallel
                % make variables accessible to workers by declaring them
                joint_center_trajectories_pool = zeros(number_of_time_steps, length(joint_center_headers)*3);
                com_trajectories_pool = zeros(number_of_time_steps, length(com_labels)*3);
                joint_angle_trajectories_pool = zeros(number_of_time_steps, number_of_joint_angles);
                marker_reference_pool = marker_reference;
                marker_labels_pool = marker_labels;
                joint_center_reference_pool = joint_center_reference;
                joint_center_labels_pool = joint_center_headers;
                segment_labels_pool = segment_labels;
                segment_coms_wcs_pool = segment_coms_wcs;
                segment_masses_pool = segment_masses;
                direction_matrices_pool = direction_matrices;
                direction_matrix_labels_pool = direction_matrix_labels;
                spmd
                    for i_time = time_steps_to_process(1)+labindex-1 : numlabs : time_steps_to_process(end)
                        % calculate joint center positions
                        marker_current = marker_trajectories(i_time, :);
                        joint_center_current = ...
                            calculateJointCenterPositions ...
                              ( ...
                                marker_reference_pool, ...
                                marker_current, ...
                                marker_labels_pool, ...
                                joint_center_reference_pool, ...
                                joint_center_labels_pool ...
                              );
                        joint_center_trajectories_pool(i_time, :) = joint_center_current;                    
                        
                        % calculate segment centers of mass
                        segment_coms_wcs_reference = segment_coms_wcs_pool;
                        number_of_segments = length(segment_coms_wcs_reference);
                        mcs_to_wcs_transformations_reference = calculateMcsToWcsTransformations_detailed([marker_reference_pool joint_center_reference_pool], [marker_labels_pool joint_center_labels_pool], segment_labels_pool);
                        mcs_to_wcs_transformations_current = calculateMcsToWcsTransformations_detailed([marker_current joint_center_current], [marker_labels_pool joint_center_labels_pool], segment_labels_pool);
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
                            body_com = body_com + segment_masses_pool(i_segment) * segment_coms_wcs_current{i_segment};
                        end
                        body_com = body_com * 1 / sum(segment_masses_pool);

                        % export centers of mass
                        for i_segment = 1 : number_of_segments
                            com_trajectories_pool(i_time, (i_segment-1)*3+1 : (i_segment-1)*3+3) = segment_coms_wcs_current{i_segment};
                        end
                        com_trajectories_pool(i_time, number_of_segments*3+1 : number_of_segments*3+3) = body_com;

                        % calculate joint angles
                        joint_angle_trajectories_pool(i_time, :) = ...
                            markerToAngles ...
                              ( ...
                                marker_reference_pool, ...
                                marker_current, ...
                                marker_labels_pool, ...
                                joint_center_reference_pool, ...
                                joint_center_current, ...
                                joint_center_labels_pool, ...
                                direction_matrices_pool, ...
                                direction_matrix_labels_pool ...
                              );                        
                    end
                end
                % reassemble
                for i_lab = 1 : number_of_labs
                    joint_center_trajectories_lab = joint_center_trajectories_pool{i_lab};
                    com_trajectories_lab = com_trajectories_pool{i_lab};
                    joint_angle_trajectories_lab = joint_angle_trajectories_pool{i_lab};

                    joint_center_trajectories(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
                        = joint_center_trajectories_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
                    com_trajectories(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
                        = com_trajectories_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
                    joint_angle_trajectories(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
                        = joint_angle_trajectories_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
                end               
            end
            if ~use_parallel
                for i_time = time_steps_to_process

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
            end
            
            fprintf([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Calculating kinematic trajectories... done\n'])

            
            % save
            variables_to_save = struct;
            variables_to_save.joint_center_trajectories = joint_center_trajectories;
            variables_to_save.joint_center_labels = joint_center_labels;
            variables_to_save.joint_angle_trajectories = joint_angle_trajectories;
            variables_to_save.joint_labels = ...
              { ...
                'pelvis, x-translation', 'pelvis, y-translation', 'pelvis, z-translation', 'pelvis, z-rotation', 'pelvis, x-rotation', 'pelvis, y-rotation', ...
                'left hip flexion/extension', 'left hip ab/adduction', 'left hip internal/external rotation', 'left knee flexion/extension', 'left knee external/internal rotation', 'left ankle dorsi/plantarflexion', 'left ankle inversion/eversion', ...
                'right hip flexion/extension', 'right hip ab/adduction', 'right hip internal/external rotation', 'right knee flexion/extension', 'right knee external/internal rotation', 'right ankle dorsi/plantarflexion', 'right ankle inversion/eversion', ...
                'lumbar joint - forward/backward bending', 'lumbar joint - sideways bending (right/left)', 'lumbar joint - internal rotation (right/left)', ...
                'cervical joint - forward/backward bending', 'cervical joint - sideways bending (right/left)', 'cervical joint - internal rotation (right/left)', ...
                'left shoulder flexion/extension', 'left shoulder ab/adduction', 'left shoulder in/external rotation', 'left elbow flexion/extension', 'left pronation/supination', 'left wrist flexion/extension', ...
                'right shoulder flexion/extension', 'right shoulder ab/adduction', 'right shoulder in/external rotation', 'right elbow flexion/extension', 'right pronation/supination', 'right wrist flexion/extension', ...
              };
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