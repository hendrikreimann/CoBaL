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


function [marker_trajectories, marker_labels, marker_directions] = rigidBodyFillMarkerFromReference ...
           ( ...
             marker_trajectories, ...
             marker_labels, ...
             marker_directions, ...
             marker_trajectories_reference, ...
             marker_labels_reference, ...
             marker_directions_reference, ...
             marker_to_fill, ...
             marker_source_1, ...
             marker_source_2, ...
             marker_source_3, ...
             visualize ...
           )
    
    % optional arguments
    if nargin < 8
        visualize = 0;
    end
       
    % extract
    marker_to_fill_reference = extractMarkerData(marker_trajectories_reference, marker_labels_reference, marker_to_fill)';
    marker_source_1_reference = extractMarkerData(marker_trajectories_reference,marker_labels_reference, marker_source_1)';
    marker_source_2_reference = extractMarkerData(marker_trajectories_reference, marker_labels_reference, marker_source_2)';
    marker_source_3_reference = extractMarkerData(marker_trajectories_reference, marker_labels_reference, marker_source_3)';

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
    
    % optimize
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
        title(['rigid body fill marker ' marker_to_fill ', based on ' marker_source_1 ', ' marker_source_2 ' and ' marker_source_3 ' - components'])
        x_plot = plot(marker_to_fill_trajectory_original(:, 1), 'linewidth', 4, 'DisplayName', 'original - x');
        y_plot = plot(marker_to_fill_trajectory_original(:, 2), 'linewidth', 4, 'DisplayName', 'original - y');
        z_plot = plot(marker_to_fill_trajectory_original(:, 3), 'linewidth', 4, 'DisplayName', 'original - z');
        plot(marker_to_fill_trajectory(:, 1), 'linewidth', 1, 'color', get(x_plot, 'color'), 'DisplayName', 'filled - x');
        plot(marker_to_fill_trajectory(:, 2), 'linewidth', 1, 'color', get(y_plot, 'color'), 'DisplayName', 'filled - y');
        plot(marker_to_fill_trajectory(:, 3), 'linewidth', 1, 'color', get(z_plot, 'color'), 'DisplayName', 'filled - z');
        legend('toggle')
        
        difference = marker_to_fill_trajectory - marker_to_fill_trajectory_original;
        error = sum(difference.^2, 2).^0.5;
        figure; axes; hold on
        title(['rigid body fill marker ' marker_to_fill ', based on ' marker_source_1 ', ' marker_source_2 ' and ' marker_source_3 ' - error'])
        plot(error, 'linewidth', 3)
        xlabel('time steps')
        ylabel('error')
        
        drawnow
    end

    % insert reconstructed trajectory back into array
    marker_indices = extractMarkerData(marker_trajectories, marker_labels, marker_to_fill, 'indices');
    if isempty(marker_indices)
        marker_indices_reference = extractMarkerData(marker_trajectories_reference, marker_labels_reference, marker_to_fill, 'indices');
        % marker was not present at all in this trial, so add it
        marker_labels = [marker_labels, marker_labels_reference(marker_indices_reference)];
        marker_directions = [marker_directions, marker_directions_reference(marker_indices_reference)];
        marker_indices = (1 : 3) + size(marker_trajectories, 2);
    end

    marker_trajectories(:, marker_indices) = marker_to_fill_trajectory;
end