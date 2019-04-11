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


function calculateComFromAngles(varargin)
    % parse arguments
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'use_parallel', false)
    parse(parser, varargin{:})
    use_parallel = parser.Results.use_parallel;
    
    load('subjectInfo.mat', 'date', 'subject_id');
    model_data = load('subjectModel.mat');
    kinematic_tree = model_data.kinematic_tree;
    segment_labels = model_data.segment_labels;
    eef_labels_single = kinematic_tree.endEffectorLabels;
    joint_center_labels = model_data.joint_center_labels;
    joint_center_directions = model_data.joint_center_directions;
    joint_directions = model_data.joint_directions;
    
    
    study_settings_file = '';
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    subject_settings = SettingsCustodian('subjectSettings.txt');
    
    
    number_of_joint_angles = kinematic_tree.numberOfJoints;
    
    if use_parallel & ~register_without_calculating
        % get or open pool of workers
        poolobject = gcp;
        number_of_labs = poolobject.NumWorkers;
    end
    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            %% prepare
            trial_type = condition_list{i_condition};
            
            % give feedback
%             disp([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Condition ' condition ', Trial ' num2str(i_trial)])
%             fprintf([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Calculating kinematic trajectories... \n'])
            
            % load data
            loaded_marker_data = load(['processed' filesep makeFileName(date, subject_id, trial_type, i_trial, 'markerTrajectories')]);
            marker_trajectories_trial = loaded_marker_data.marker_trajectories;
            marker_labels_trial = loaded_marker_data.marker_labels;
            time_mocap = loaded_marker_data.time_mocap;
            sampling_rate_mocap = loaded_marker_data.sampling_rate_mocap;
            loaded_kinematic_data = load(['processed' filesep makeFileName(date, subject_id, trial_type, i_trial, 'kinematicTrajectories')]);
            angle_trajectories_trial = loaded_kinematic_data.joint_angle_trajectories;
            joint_labels_trial = loaded_kinematic_data.joint_labels;
            
            number_of_time_steps = size(marker_trajectories_trial, 1); %#ok<NODEF>
            time_steps_to_process = 1 : number_of_time_steps;
%             time_steps_to_process = 1 : 10;
            number_of_time_steps_to_process = length(time_steps_to_process);
            
            com_labels_single = [segment_labels 'BODY'];
            for i_label = 1 : length(com_labels_single)
                com_labels_single{i_label} = [com_labels_single{i_label} 'COM'];
            end
            
            
            %% process
%             disp([' - Condition ' trial_type ', Trial ' num2str(i_trial)])
%             fprintf([' - Calculating kinematic trajectories... \n'])

            % calculate
            com_trajectories = zeros(number_of_time_steps, length(com_labels_single)*3) * NaN;
            eef_trajectories = zeros(number_of_time_steps, 3 * kinematic_tree.numberOfBranches);
            if use_parallel
%                 % make variables accessible to workers by declaring them
%                 joint_center_trajectories_pool = zeros(number_of_time_steps, length(joint_center_labels));
%                 com_trajectories_pool = zeros(number_of_time_steps, length(com_labels_single)*3);
%                 joint_angle_trajectories_pool = zeros(number_of_time_steps, number_of_joint_angles);
%                 marker_positions_reference_pool = marker_positions_reference;
%                 marker_labels_reference_pool = marker_labels_reference;
%                 marker_labels_trial_pool = marker_labels_trial;
%                 joint_center_positions_reference_pool = joint_center_positions_reference;
%                 joint_center_labels_pool = joint_center_labels;
%                 segment_labels_pool = segment_labels;
%                 segment_coms_wcs_pool = segment_coms_wcs;
%                 segment_masses_pool = segment_masses;
%                 direction_matrices_pool = direction_matrices;
%                 direction_matrix_labels_pool = direction_matrix_labels;
%                 subject_settings_pool = subject_settings;
%                 new_ignore_times_pool = new_ignore_times;
%                 time_mocap_pool = time_mocap;
%                 spmd
%                     for i_time_index = labindex : numlabs : number_of_time_steps_to_process
%                         i_time = time_steps_to_process(i_time_index);
%                         % check for missing markers
%                         marker_positions_current = marker_trajectories_trial(i_time, :);
%                         if any(any(isnan(marker_trajectories_trial(i_time, essential_marker_indicator))))
% 
%                             joint_center_trajectories_pool(i_time, :) = NaN;
%                             com_trajectories_pool(i_time, :) = NaN;
%                             joint_angle_trajectories_pool(i_time, :) = NaN;
%                             new_ignore_times_pool = [new_ignore_times_pool; time_mocap_pool(i_time)]; %#ok<AGROW>
%                         else
%                             % calculate joint center positions
%                             subject_settings_pool.verbose = false;
%                             joint_center_positions_current = ...
%                                 calculateJointCenterPositions ...
%                                   ( ...
%                                     marker_positions_reference_pool, ...
%                                     marker_positions_current, ...
%                                     marker_labels_reference_pool, ...
%                                     marker_labels_trial_pool, ...
%                                     joint_center_positions_reference_pool, ...
%                                     joint_center_labels_pool, ...
%                                     subject_settings_pool ...
%                                   );
%                             joint_center_trajectories_pool(i_time, :) = joint_center_positions_current;                    
% 
%                             % calculate segment centers of mass
%                             segment_coms_wcs_reference = segment_coms_wcs_pool;
%                             number_of_segments = length(segment_coms_wcs_reference);
%                             mcs_to_wcs_transformations_reference = ...
%                                 calculateMcsToWcsTransformations_detailed ...
%                                   ( ...
%                                     [marker_positions_reference_pool joint_center_positions_reference_pool], ...
%                                     [marker_labels_reference_pool joint_center_labels_pool], ...
%                                     segment_labels_pool ...
%                                   );
%                             mcs_to_wcs_transformations_current = ...
%                                 calculateMcsToWcsTransformations_detailed ...
%                                   ( ...
%                                     [marker_positions_current joint_center_positions_current], ...
%                                     [marker_labels_reference_pool joint_center_labels_pool], ...
%                                     segment_labels_pool ...
%                                   );
% 
%                             segment_coms_wcs_current = cell(number_of_segments, 1);
%                             for i_segment = 1 : number_of_segments
%                                 segment_com_wcs_reference = [segment_coms_wcs_reference{i_segment}; 1];
%                                 T_mcs_to_wcs_reference = mcs_to_wcs_transformations_reference{i_segment};
%                                 T_mcs_to_wcs_current = mcs_to_wcs_transformations_current{i_segment};
%                                 segment_com_mcs = T_mcs_to_wcs_reference^(-1) * segment_com_wcs_reference;
%                                 segment_coms_wcs_current{i_segment} = eye(3, 4) * T_mcs_to_wcs_current * segment_com_mcs;
%                             end
% 
%                             % calculate whole body center of mass
%                             body_com = [0; 0; 0];
%                             for i_segment = 1 : number_of_segments
%                                 body_com = body_com + segment_masses_pool(i_segment) * segment_coms_wcs_current{i_segment};
%                             end
%                             body_com = body_com * 1 / sum(segment_masses_pool);
% 
%                             % export centers of mass
%                             for i_segment = 1 : number_of_segments
%                                 com_trajectories_pool(i_time, (i_segment-1)*3+1 : (i_segment-1)*3+3) = segment_coms_wcs_current{i_segment};
%                             end
%                             com_trajectories_pool(i_time, number_of_segments*3+1 : number_of_segments*3+3) = body_com;
% 
%                             % calculate joint angles
%                             joint_angle_trajectories_pool(i_time, :) = ...
%                                 markerToAngles ...
%                                   ( ...
%                                     marker_positions_reference_pool, ...
%                                     marker_positions_current, ...
%                                     marker_labels_reference_pool, ...
%                                     marker_labels_trial_pool, ...
%                                     joint_center_positions_reference_pool, ...
%                                     joint_center_positions_current, ...
%                                     joint_center_labels_pool, ...
%                                     joint_center_labels_pool, ...
%                                     direction_matrices_pool, ...
%                                     direction_matrix_labels_pool ...
%                                   );    
%                         end
%                     end
%                 end
%                 % reassemble
%                 for i_lab = 1 : number_of_labs
%                     joint_center_trajectories_lab = joint_center_trajectories_pool{i_lab};
%                     com_trajectories_lab = com_trajectories_pool{i_lab};
%                     joint_angle_trajectories_lab = joint_angle_trajectories_pool{i_lab};
%                     new_ignore_times_lab = new_ignore_times_pool{i_lab};
%                     joint_center_trajectories(time_steps_to_process(i_lab : number_of_labs : number_of_time_steps_to_process), :) ...
%                         = joint_center_trajectories_lab(time_steps_to_process(i_lab : number_of_labs : number_of_time_steps_to_process), :);
%                     com_trajectories(time_steps_to_process(i_lab : number_of_labs : number_of_time_steps_to_process), :) ...
%                         = com_trajectories_lab(time_steps_to_process(i_lab : number_of_labs : number_of_time_steps_to_process), :);
%                     joint_angle_trajectories(time_steps_to_process(i_lab : number_of_labs : number_of_time_steps_to_process), :) ...
%                         = joint_angle_trajectories_lab(time_steps_to_process(i_lab : number_of_labs : number_of_time_steps_to_process), :);
%                     new_ignore_times = [new_ignore_times; new_ignore_times_lab];
%                 end               
            end
            if ~use_parallel
                for i_time_step = 1 : length(time_steps_to_process)
                    % update model to current time step
                    i_time = time_steps_to_process(i_time_step);
                    joint_angles = angle_trajectories_trial(i_time, :)';
                    kinematic_tree.jointAngles = joint_angles;
                    kinematic_tree.updateInternals();
                    
                    % calculate center of mass for segments and whole body
                    for i_joint = 1 : number_of_joint_angles
                        this_segment_com = kinematic_tree.linkTransformations{i_joint}(1:3, 4);
                        com_trajectories(i_time, i_joint*3 + [1 2 3] - 3) = this_segment_com;
                    end
                    whole_body_com = kinematic_tree.calculateCenterOfMassPosition();
                    com_trajectories(i_time, end-2:end) = whole_body_com;
                    
                    % store end-effector
                    for i_branch = 1 : kinematic_tree.numberOfBranches
                        this_eef = kinematic_tree.endEffectorPositions{i_branch};
                        eef_trajectories(i_time, i_branch*3 + [1 2 3] - 3) = this_eef;
                    end

                    % give progress feedback
%                     display_step = 100;
%                     if (i_time_step / display_step) == floor(i_time_step / display_step)
%                         disp([num2str(i_time_step) '(' num2str(length(time_steps_to_process)) ')']);
%                     end

                end                
            end
%             fprintf([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Calculating kinematic trajectories... done\n'])

            % triplicate com labels and add direction
            number_of_com_labels = length(com_labels_single);
            com_labels = cell(3, number_of_com_labels);
            for i_com = 1 : length(com_labels)
                com_labels{1, i_com} = [com_labels_single{i_com} '_x'];
                com_labels{2, i_com} = [com_labels_single{i_com} '_y'];
                com_labels{3, i_com} = [com_labels_single{i_com} '_z'];
            end
            com_labels = reshape(com_labels, 1, number_of_com_labels*3);
            
            number_of_com_trajectories = size(com_trajectories, 2);
            com_directions = cell(2, number_of_com_trajectories);
            [com_directions{1, 1 : 3 : number_of_com_trajectories}] = deal('right');
            [com_directions{2, 1 : 3 : number_of_com_trajectories}] = deal('left');
            [com_directions{1, 2 : 3 : number_of_com_trajectories}] = deal('forward');
            [com_directions{2, 2 : 3 : number_of_com_trajectories}] = deal('backward');
            [com_directions{1, 3 : 3 : number_of_com_trajectories}] = deal('up');
            [com_directions{2, 3 : 3 : number_of_com_trajectories}] = deal('down');
            
            % triplicate eef labels and add direction
            number_of_eef_labels = length(eef_labels_single);
            eef_labels = cell(3, number_of_eef_labels);
            for i_eef = 1 : length(eef_labels_single)
                eef_labels{1, i_eef} = [eef_labels_single{i_eef} '_x'];
                eef_labels{2, i_eef} = [eef_labels_single{i_eef} '_y'];
                eef_labels{3, i_eef} = [eef_labels_single{i_eef} '_z'];
            end
            eef_labels = reshape(eef_labels, 1, number_of_eef_labels*3);
            
            number_of_eef_trajectories = size(eef_trajectories, 2);
            eef_directions = cell(2, number_of_eef_trajectories);
            [eef_directions{1, 1 : 3 : number_of_eef_trajectories}] = deal('right');
            [eef_directions{2, 1 : 3 : number_of_eef_trajectories}] = deal('left');
            [eef_directions{1, 2 : 3 : number_of_eef_trajectories}] = deal('forward');
            [eef_directions{2, 2 : 3 : number_of_eef_trajectories}] = deal('backward');
            [eef_directions{1, 3 : 3 : number_of_eef_trajectories}] = deal('up');
            [eef_directions{2, 3 : 3 : number_of_eef_trajectories}] = deal('down');
            
            
            
            % TODO: directions and labels for the end-effector
            
            
            %% save com trajectories
            variables_to_save = struct;
            variables_to_save.com_trajectories = com_trajectories;
            variables_to_save.com_labels = com_labels;
            variables_to_save.com_directions = com_directions;
            variables_to_save.time_mocap = time_mocap;
            variables_to_save.sampling_rate_mocap = sampling_rate_mocap;
            
            save_folder = 'processed';
            save_file_name = makeFileName(date, subject_id, trial_type, i_trial, 'comTrajectories.mat');
            saveDataToFile([save_folder filesep save_file_name], variables_to_save);

            addAvailableData ...
              ( ...
                'com_trajectories', ...
                'time_mocap', ...
                'sampling_rate_mocap', ...
                '_com_labels', ...
                '_com_directions', ...
                save_folder, ...
                save_file_name ...
              );
          
            %% save eef trajectories
            variables_to_save = struct;
            variables_to_save.eef_trajectories = eef_trajectories;
            variables_to_save.eef_labels = eef_labels;
            variables_to_save.eef_directions = eef_directions;
            variables_to_save.time_mocap = time_mocap;
            variables_to_save.sampling_rate_mocap = sampling_rate_mocap;
            
            save_folder = 'processed';
            save_file_name = makeFileName(date, subject_id, trial_type, i_trial, 'eefTrajectories.mat');
            saveDataToFile([save_folder filesep save_file_name], variables_to_save);

            addAvailableData ...
              ( ...
                'eef_trajectories', ...
                'time_mocap', ...
                'sampling_rate_mocap', ...
                '_eef_labels', ...
                '_eef_directions', ...
                save_folder, ...
                save_file_name ...
              );
            
            % feedback
            disp([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ': type ' trial_type ', trial ' num2str(i_trial) ' completed, saved as ' save_folder filesep save_file_name]);
        end
    end
end