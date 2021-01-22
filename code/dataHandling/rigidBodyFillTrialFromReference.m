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


function rigidBodyFillTrialFromReference(marker_to_fill, marker_source_1, marker_source_2, marker_source_3, varargin)
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;
    
    % load settings
     subject_settings = SettingsCustodian('subjectSettings.txt');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');


    % get reference positions
%     marker_reference = load(['processed' filesep makeFileName(date, subject_id, 'calibration', '1', 'markerTrajectories')]);
     % load static reference file
    load(['processed' filesep makeFileName(collection_date, subject_id, subject_settings.get('static_reference_trial_type'), subject_settings.get('static_reference_trial_number'), 'markerTrajectories')]);

        % find first time step where all markers are available
    i_time = 1;
    while any(isnan(marker_trajectories(i_time, :)))
        i_time = i_time + 1;
    end
    marker_reference = marker_trajectories(i_time, :);
    
    % marker_vlabels = 
    marker_to_fill_reference = extractMarkerData(marker_reference, marker_labels, marker_to_fill)';
    marker_source_1_reference = extractMarkerData(marker_reference,marker_labels, marker_source_1)';
    marker_source_2_reference = extractMarkerData(marker_reference, marker_labels, marker_source_2)';
    marker_source_3_reference = extractMarkerData(marker_reference, marker_labels, marker_source_3)';

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
            load(['processed' filesep makeFileName(collection_date, subject_id, condition, i_trial, 'markerTrajectories')]);
            number_of_time_steps = size(marker_trajectories, 1);
            
            marker_source_1_trajectory = extractMarkerData(marker_trajectories, marker_labels, marker_source_1);
            marker_source_2_trajectory = extractMarkerData(marker_trajectories, marker_labels, marker_source_2);
            marker_source_3_trajectory = extractMarkerData(marker_trajectories, marker_labels, marker_source_3);
            
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
                marker_to_fill_trajectory_original = extractMarkerData(marker_trajectories, marker_labels, marker_to_fill);
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
            marker_indices = extractMarkerData(marker_trajectories, marker_labels, marker_to_fill, 'indices');
            if isempty(marker_indices)
                % marker was not present at all in this trial, so add it
                marker_labels = ...
                  [ ...
                    marker_labels, ...
                    [marker_to_fill '_x'], ...
                    [marker_to_fill '_y'], ...
                    [marker_to_fill '_z'] ...
                  ]; %#ok<AGROW>
                single_marker_directions = marker_directions(:, 1:3); % re-use this, assuming the direction is the same for all markers
                marker_directions = [marker_directions, single_marker_directions]; %#ok<AGROW>
                marker_indices = (1 : 3) + size(marker_trajectories, 2);
            end
            
            marker_trajectories(:, marker_indices) = marker_to_fill_trajectory; %#ok<AGROW>
            
            % save
            variables_to_save = struct;
            variables_to_save.marker_trajectories = marker_trajectories;
            variables_to_save.marker_labels = marker_labels;
            variables_to_save.marker_directions = marker_directions;
            
            save_folder = 'processed';
            save_file_name = makeFileName(collection_date, subject_id, condition, i_trial, 'markerTrajectories.mat');
            saveDataToFile([save_folder filesep save_file_name], variables_to_save);
            disp(['Rigid body filling marker ' marker_to_fill ', condition ' condition ', Trial ' num2str(i_trial)]);
            
        end
    end
end