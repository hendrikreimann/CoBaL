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


function rigidBodyFill(marker_to_fill, marker_source_1, marker_source_2, marker_source_3, varargin)
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;
    
    
    
    
    % load settings
    study_settings_file = '';
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    
    load('subjectInfo.mat', 'date', 'subject_id');
    load('subjectModel.mat');

    % get reference positions
    marker_reference = kinematic_tree.exportMarkerPositions;
    marker_labels = kinematic_tree.markerLabels;
    marker_to_fill_reference = extractMarkerTrajectories(marker_reference, marker_labels, marker_to_fill)';
    marker_source_1_reference = extractMarkerTrajectories(marker_reference, marker_labels, marker_source_1)';
    marker_source_2_reference = extractMarkerTrajectories(marker_reference, marker_labels, marker_source_2)';
    marker_source_3_reference = extractMarkerTrajectories(marker_reference, marker_labels, marker_source_3)';
    
    % define coordinate frame based on positions
    p = mean([marker_source_1_reference marker_source_2_reference marker_source_3_reference], 2);
    u1 = marker_source_1_reference - p;
    u2 = marker_source_2_reference - p;
    u3 = cross(u1, u2);
    R = orthogonalizeBasis([u1 u2 u3]);

    % define transformation
    T_wm = [R p; [0 0 0 1]]; % this matrix transforms coordinates of a point from marker frame to world frame
    T_mw = T_wm^(-1);
    q_world = [marker_to_fill_reference; 1];
    q_marker = T_mw * q_world;
    
    %% optimize
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            %% load data
            condition = condition_list{i_condition};
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);
            number_of_time_steps = size(marker_trajectories, 1);
            
            marker_source_1_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, marker_source_1);
            marker_source_2_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, marker_source_2);
            marker_source_3_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, marker_source_3);
            
            marker_to_fill_trajectory = zeros(number_of_time_steps, 3);
            
            % process
            for i_time = 1 : number_of_time_steps
                % define coordinate frame based on positions
                p = mean([marker_source_1_trajectory(i_time, :)' marker_source_2_trajectory(i_time, :)' marker_source_3_trajectory(i_time, :)'], 2);
                u1 = marker_source_1_trajectory(i_time, :)' - p;
                u2 = marker_source_2_trajectory(i_time, :)' - p;
                u3 = cross(u1, u2);
                R = orthogonalizeBasis([u1 u2 u3]);

                % transform and store
                T_wm = [R p; [0 0 0 1]]; % this matrix transforms coordinates of a point from marker frame to world frame
                q_world = T_wm * q_marker;
                marker_to_fill_trajectory(i_time, :) = q_world(1:3);
            end
            
            if visualize
                marker_to_fill_trajectory_original = extractMarkerTrajectories(marker_trajectories, marker_labels, marker_to_fill);
                figure; hold on
                x_plot = plot(marker_to_fill_trajectory_original(:, 1), 'linewidth', 2);
                y_plot = plot(marker_to_fill_trajectory_original(:, 2), 'linewidth', 2);
                z_plot = plot(marker_to_fill_trajectory_original(:, 3), 'linewidth', 2);
                plot(marker_to_fill_trajectory(:, 1), 'linewidth', 1, 'color', get(x_plot, 'color'));
                plot(marker_to_fill_trajectory(:, 2), 'linewidth', 1, 'color', get(y_plot, 'color'));
                plot(marker_to_fill_trajectory(:, 3), 'linewidth', 1, 'color', get(z_plot, 'color'));
                
                difference = marker_to_fill_trajectory - marker_to_fill_trajectory_original;
                error = sum(difference.^2, 2).^0.5;
                figure; hold on
                plot(error, 'linewidth', 3)
            end
            
            % insert reconstructed trajectory back into array
            marker_number = find(strcmp(marker_labels, marker_to_fill));
            markers_indices = reshape([(marker_number - 1) * 3 + 1; (marker_number - 1) * 3 + 2; (marker_number - 1) * 3 + 3], 1, length(marker_number)*3);
            marker_trajectories(:, markers_indices) = marker_to_fill_trajectory;
            
            % save
            variables_to_save = struct;
            variables_to_save.marker_trajectories = marker_trajectories;
            
            save_folder = 'processed';
            save_file_name = makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories.mat');
            saveDataToFile([save_folder filesep save_file_name], variables_to_save);
            disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' save_folder filesep save_file_name]);
            
        end
    end
    


end