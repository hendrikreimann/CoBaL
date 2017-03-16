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




load('/Users/reimajbi/Box Sync/inverseKinematics/YJI/processed/20170307_YJI_walking_001_dynamicTrajectories.mat')
load('/Users/reimajbi/Box Sync/inverseKinematics/YJI/processed/20170307_YJI_walking_001_forceplateTrajectories.mat')


color_x_a = [1 0 0];
color_y_a = [0 1 0];
color_z_a = [0 0 1];
color_x_b = [0.7 0 0.3];
color_y_b = [0.3 0.7 0];
color_z_b = [0 0.3 0.7];

figure; total_grf_axes = axes; hold on; title('total ground reaction forces');
plot(time_forceplate, total_forceplate_wrench_world(:, 1), 'color', color_x_a);
plot(time_forceplate, total_forceplate_wrench_world(:, 2), 'color', color_y_a);
plot(time_forceplate, total_forceplate_wrench_world(:, 3), 'color', color_z_a);
plot(time_mocap, total_ground_reaction_wrench_trajectory(:, 1), 'color', color_x_b, 'linewidth', 2);
plot(time_mocap, total_ground_reaction_wrench_trajectory(:, 2), 'color', color_y_b, 'linewidth', 2);
plot(time_mocap, total_ground_reaction_wrench_trajectory(:, 3), 'color', color_z_b, 'linewidth', 2);

figure; total_grm_axes = axes; hold on; title('total ground reaction moments');
plot(time_forceplate, total_forceplate_wrench_world(:, 4), 'color', color_x_a);
plot(time_forceplate, total_forceplate_wrench_world(:, 5), 'color', color_y_a);
plot(time_forceplate, total_forceplate_wrench_world(:, 6), 'color', color_z_a);
plot(time_mocap, total_ground_reaction_wrench_trajectory(:, 4), 'color', color_x_b, 'linewidth', 2);
plot(time_mocap, total_ground_reaction_wrench_trajectory(:, 5), 'color', color_y_b, 'linewidth', 2);
plot(time_mocap, total_ground_reaction_wrench_trajectory(:, 6), 'color', color_z_b, 'linewidth', 2);

figure; left_grf_axes = axes; hold on; title('left ground reaction forces');
plot(time_forceplate, left_forceplate_wrench_world(:, 1), 'color', color_x_a);
plot(time_forceplate, left_forceplate_wrench_world(:, 2), 'color', color_y_a);
plot(time_forceplate, left_forceplate_wrench_world(:, 3), 'color', color_z_a);
plot(time_mocap, left_ground_reaction_wrench_trajectory(:, 1), 'color', color_x_b, 'linewidth', 2);
plot(time_mocap, left_ground_reaction_wrench_trajectory(:, 2), 'color', color_y_b, 'linewidth', 2);
plot(time_mocap, left_ground_reaction_wrench_trajectory(:, 3), 'color', color_z_b, 'linewidth', 2);

figure; left_grm_axes = axes; hold on; title('left ground reaction moments');
plot(time_forceplate, left_forceplate_wrench_world(:, 4), 'color', color_x_a);
plot(time_forceplate, left_forceplate_wrench_world(:, 5), 'color', color_y_a);
plot(time_forceplate, left_forceplate_wrench_world(:, 6), 'color', color_z_a);
plot(time_mocap, left_ground_reaction_wrench_trajectory(:, 4), 'color', color_x_b, 'linewidth', 2);
plot(time_mocap, left_ground_reaction_wrench_trajectory(:, 5), 'color', color_y_b, 'linewidth', 2);
plot(time_mocap, left_ground_reaction_wrench_trajectory(:, 6), 'color', color_z_b, 'linewidth', 2);

figure; right_grf_axes = axes; hold on; title('right ground reaction forces');
plot(time_forceplate, right_forceplate_wrench_world(:, 1), 'color', color_x_a);
plot(time_forceplate, right_forceplate_wrench_world(:, 2), 'color', color_y_a);
plot(time_forceplate, right_forceplate_wrench_world(:, 3), 'color', color_z_a);
plot(time_mocap, right_ground_reaction_wrench_trajectory(:, 1), 'color', color_x_b, 'linewidth', 2);
plot(time_mocap, right_ground_reaction_wrench_trajectory(:, 2), 'color', color_y_b, 'linewidth', 2);
plot(time_mocap, right_ground_reaction_wrench_trajectory(:, 3), 'color', color_z_b, 'linewidth', 2);

figure; right_grm_axes = axes; hold on; title('right ground reaction moments');
plot(time_forceplate, right_forceplate_wrench_world(:, 4), 'color', color_x_a);
plot(time_forceplate, right_forceplate_wrench_world(:, 5), 'color', color_y_a);
plot(time_forceplate, right_forceplate_wrench_world(:, 6), 'color', color_z_a);
plot(time_mocap, right_ground_reaction_wrench_trajectory(:, 4), 'color', color_x_b, 'linewidth', 2);
plot(time_mocap, right_ground_reaction_wrench_trajectory(:, 5), 'color', color_y_b, 'linewidth', 2);
plot(time_mocap, right_ground_reaction_wrench_trajectory(:, 6), 'color', color_z_b, 'linewidth', 2);

linkaxes([total_grf_axes total_grm_axes left_grf_axes left_grm_axes right_grf_axes right_grm_axes], 'x')
distFig('rows', 3, 'tight', true);

