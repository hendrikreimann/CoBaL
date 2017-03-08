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



% 
% load('/Users/reimajbi/Box Sync/inverseKinematics/BRC/processed/00000000_XXX_simulation_001_markerTrajectories.mat')
% load('/Users/reimajbi/Box Sync/inverseKinematics/BRC/processed/00000000_XXX_simulation_001_kinematicTrajectories.mat')
load('/Users/reimajbi/Box Sync/inverseKinematics/DCH/processed/20170110_DCH_baselineOG_001_markerTrajectories.mat')
load('/Users/reimajbi/Box Sync/inverseKinematics/DCH/processed/20170110_DCH_baselineOG_001_kinematicTrajectories.mat')

load('subjectInfo.mat', 'date', 'subject_id');
load('subjectModel.mat');

% joint_groups = ...
%   { ...
%     1:3 ...
%     4:6 ...
%   };
marker_groups = ...
  { ...
    1:4, ...
    [5:9 10 17], ...
    11:16, ...
    18:23, ...
    24:27, ...
    28:30, ...
    31:33, ...
    34:36, ...
    37:39, ...
    40:42, ...
    43:45, ...
  };

time_to_process = 1 : 500;

% reconstruct
number_of_time_steps = size(marker_trajectories, 1);
number_of_markers = size(marker_trajectories, 2)/3;
marker_positions_recorded = marker_trajectories;
marker_positions_reconstructed_calculated = zeros(size(marker_trajectories));
marker_positions_reconstructed_optimized = zeros(size(marker_trajectories));
marker_position_errors_calculated = zeros(number_of_time_steps, number_of_markers);
marker_position_errors_optimized = zeros(number_of_time_steps, number_of_markers);
for i_time = 1 : length(time_to_process)
% for i_time = 1 : 500
    kinematic_tree.jointAngles = joint_angle_trajectories_calculated(i_time, :)';
    kinematic_tree.updateConfiguration;
    marker_positions_reconstructed_calculated(i_time, :) = kinematic_tree.exportMarkerPositions;
    for i_marker = 1 : number_of_markers
        this_marker_position_reconstructed = marker_positions_reconstructed_calculated(i_time, 3*(i_marker-1) + [1 2 3]);
        this_marker_position_recorded = marker_positions_recorded(i_time, 3*(i_marker-1) + [1 2 3]);
        delta_position = this_marker_position_reconstructed - this_marker_position_recorded;
        error = norm(delta_position);
        marker_position_errors_calculated(i_time, i_marker) = error;
    end
    
    kinematic_tree.jointAngles = joint_angle_trajectories_optimized(i_time, :)';
    kinematic_tree.updateConfiguration;
    marker_positions_reconstructed_optimized(i_time, :) = kinematic_tree.exportMarkerPositions;
    for i_marker = 1 : number_of_markers
        this_marker_position_reconstructed = marker_positions_reconstructed_optimized(i_time, 3*(i_marker-1) + [1 2 3]);
        this_marker_position_recorded = marker_positions_recorded(i_time, 3*(i_marker-1) + [1 2 3]);
        delta_position = this_marker_position_reconstructed - this_marker_position_recorded;
        error = norm(delta_position);
        marker_position_errors_optimized(i_time, i_marker) = error;
    end
end

for i_group = 1 : length(marker_groups)
    figure;
    axes_x = subplot(4, 1, 1); hold on;
    axes_y = subplot(4, 1, 2); hold on;
    axes_z = subplot(4, 1, 3); hold on;
    axes_e = subplot(4, 1, 4); hold on;
    
    for i_marker = marker_groups{i_group}
        marker_positions_recorded_plot = plot(axes_x, time_to_process, marker_positions_recorded(time_to_process, 3*(i_marker-1) + 1), 'linewidth', 2, 'linestyle', '-', 'HandleVisibility', 'off');
        plot(axes_x, time_to_process, marker_positions_reconstructed_calculated(time_to_process, 3*(i_marker-1) + 1), 'color', marker_positions_recorded_plot.Color, 'linestyle', '-', 'linewidth', 1, 'displayname', num2str(i_marker));
        plot(axes_x, time_to_process, marker_positions_reconstructed_optimized(time_to_process, 3*(i_marker-1) + 1), 'color', marker_positions_recorded_plot.Color, 'linestyle', '--', 'linewidth', 1, 'displayname', num2str(i_marker));
        
        marker_positions_recorded_plot = plot(axes_y, time_to_process, marker_positions_recorded(time_to_process, 3*(i_marker-1) + 2), 'linewidth', 2, 'linestyle', '-', 'HandleVisibility', 'off');
        plot(axes_y, time_to_process, marker_positions_reconstructed_calculated(time_to_process, 3*(i_marker-1) + 2), 'color', marker_positions_recorded_plot.Color, 'linestyle', '-', 'linewidth', 1, 'displayname', num2str(i_marker));
        plot(axes_y, time_to_process, marker_positions_reconstructed_optimized(time_to_process, 3*(i_marker-1) + 2), 'color', marker_positions_recorded_plot.Color, 'linestyle', '--', 'linewidth', 1, 'displayname', num2str(i_marker));
        
        marker_positions_recorded_plot = plot(axes_z, time_to_process, marker_positions_recorded(time_to_process, 3*(i_marker-1) + 3), 'linewidth', 2, 'linestyle', '-', 'HandleVisibility', 'off');
        plot(axes_z, time_to_process, marker_positions_reconstructed_calculated(time_to_process, 3*(i_marker-1) + 3), 'color', marker_positions_recorded_plot.Color, 'linestyle', '-', 'linewidth', 1, 'displayname', num2str(i_marker));
        plot(axes_z, time_to_process, marker_positions_reconstructed_optimized(time_to_process, 3*(i_marker-1) + 3), 'color', marker_positions_recorded_plot.Color, 'linestyle', '--', 'linewidth', 1, 'displayname', num2str(i_marker));
        
        error_plot = plot(axes_e, time_to_process, marker_position_errors_calculated(time_to_process, i_marker), 'linewidth', 2, 'linestyle', '-', 'displayname', 'calculated');
        plot(axes_e, time_to_process, marker_position_errors_optimized(time_to_process, i_marker), 'color', error_plot.Color, 'linewidth', 2, 'linestyle', '--', 'displayname', 'optimized');
    end
%     legend('show')
end
