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
% load('/Users/reimajbi/Box Sync/inverseKinematics/YJI/processed/20170307_YJI_walking_001_dynamicTrajectories.mat')


% joint_groups = ...
%   { ...
%     1:3 ...
%     4:6 ...
%   };
joint_groups = ...
  { ...
    1:6, ...
    7:13, ...
    14:20, ...
    21:26, ...
    27:32, ...
    33:38, ...
  };

% pelvis
for i_group = 1 : length(joint_groups)
    figure; axes; hold on
    for i_joint = joint_groups{i_group}
        joint_angle_simulated_plot = plot(time_mocap, joint_torque_trajectories_belt(:, i_joint), 'linewidth', 2, 'linestyle', '-', 'displayname', num2str(i_joint));
    end
    legend('show')
end
