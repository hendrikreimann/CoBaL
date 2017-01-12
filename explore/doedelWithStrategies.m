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

% figure; plot(left_ankle_acceleration_strategy(:, 1:6));
figure; plot(left_ankle_acceleration_strategy(:, 7:12));
% figure; plot(left_hip_acceleration_strategy(:, 1:6));
figure; plot(left_hip_acceleration_strategy(:, 7:12));
% figure; plot(right_ankle_acceleration_strategy(:, 1:6));
figure; plot(right_ankle_acceleration_strategy(:, 13:18));
% figure; plot(right_hip_acceleration_strategy(:, 1:6));
figure; plot(right_hip_acceleration_strategy(:, 13:18));
return
% A_left_trajectory_unfolded = zeros(256, 5, 38);
% J_c_left_trajectory_unfolded = zeros(256, 2, 38);
% for i_time = 1 : 256
%     A_left_trajectory_unfolded(i_time, :, :) = A_left_trajectory{i_time};
%     J_c_left_trajectory_unfolded(i_time, :, :) = J_c_left_trajectory{i_time};
% end

% for i_row = 1 : 5
%     figure; hold on
%     for i_joint = 1 : 38;
%         plot3(1:256, A_left_trajectory_unfolded(:, i_row, i_joint), i_joint*ones(256, 1));
%     end
% end

% for i_row = 1 : 2
%     figure; hold on
%     for i_joint = 1 : 38;
%         plot3(1:256, J_c_left_trajectory_unfolded(:, i_row, i_joint), i_joint*ones(256, 1));
%     end
% end

scene_limits = [-0.8 0.8; 56 58; -0.1 2];
stick_figure = showMarkerComparisonStickFigure(plant, joint_angles_mean_left_control, [], scene_limits);
stick_figure.showLinkMassEllipsoids = false;
stick_figure.update;

scene_limits = [-0.8 0.8; 50 52; -0.1 2];
stick_figure = showMarkerComparisonStickFigure(plant, joint_angles_mean_right_control, [], scene_limits);
stick_figure.showLinkMassEllipsoids = false;
stick_figure.update;




