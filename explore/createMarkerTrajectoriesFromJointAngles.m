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

% this script creates marker trajectories from known joint angles and the kinematic model

% load model
load('subjectModel.mat')

% time
time_step = 0.01;
total_time = 0.1;
time_mocap = time_step : time_step : total_time;

% define joint angle trajectories
%     0.1; 0.2; 0.3; 0.4; 0.5; 0.6; % pelvis free body DoFs
%     0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; ... % left leg
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; % pelvis free body DoFs
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.5; ... % left leg
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; % lumbar and cervical
%     0.14; 0.24; 0.34; 0.44; 0.54; 0.64; 0.74; ... % left arm
theta_init = ...
  [ ...
    0.0; 0.0; 0.0; 0.0; 0.0; 0.0; % pelvis free body DoFs
    0.11; 0.21; 0.31; 0.41; 0.51; 0.61; 0.71; ... % left leg
    0.12; 0.22; 0.32; 0.42; 0.52; 0.62; 0.72; ... % right leg
    0.13; 0.23; 0.33; 0.43; 0.53; 0.63; % lumbar and cervical
    0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; ... % left arm
    0; 0; 0; 0; 0; 0; 0; ... right arm
  ];

%     0.1; 0.2; 0.3; 0.5; 0.0; 0.0; 0.0; ... % left leg
theta_direction_1 = ...
  [ ...
    0.1; 0.2; 0.3; 0.4; 0.5; 0.6; % pelvis free body DoFs
    0.11; 0.21; 0.31; 0.41; 0.51; 0.61; 0.71; ... % left leg
    0.12; 0.22; 0.32; 0.42; 0.52; 0.62; 0.72; ... % right leg
    0.13; 0.23; 0.33; 0.43; 0.53; 0.63; % lumbar and cervical
    0.14; 0.24; 0.34; 0.44; 0.54; 0.64; 0; ... left arm
    0; 0; 0; 0; 0; 0; 0; ... right arm
  ];

% theta_init = rand(size(theta_init));
% theta_direction_1 = rand(size(theta_direction_1));
% theta_init = ...
%   [ ...
%     0.1; 0.2; 0.3; 0.4; 0.5; 0.6; % pelvis free body DoFs
%     .5; .1; .3; .15; .2; .4; 0.3; ... % left leg
%     .0; .0; .0; .0; .0; .0; 0; ... % right leg
%     0; 0; 0; 0; 0; 0; ... % lumbar and cervical
%     0; 0; 0; 0; 0; 0; 0; ... left arm
%     0; 0; 0; 0; 0; 0; 0; ... right arm
%   ];
% 
% theta_direction_1 = ...
%   [ ...
%     0.1; 0.2; 0.3; 0.4; 0.5; 0.6; % pelvis free body DoFs
%     .5; .1; .3; .15; .2; .4; 0.3; ... % left leg
%     .0; .0; .0; .0; .0; .0; 0; ... % right leg
%     0; 0; 0; 0; 0; 0; ... % lumbar and cervical
%     0; 0; 0; 0; 0; 0; 0; ... left arm
%     0; 0; 0; 0; 0; 0; 0; ... right arm
%   ];

amplitude_1 = 0.5;
frequency_1 = 0.5;
% amplitude_2 = -0.2;
% frequency_2 = 1.2;
% amplitude_3 = 0.3;
% frequency_3 = -2.2;

sinusoid_pos_1 = amplitude_1 * (-cos(time_mocap * 2 * pi * frequency_1 / total_time)+1);
% sinusoid_pos_2 = amplitude_2 * (-cos(time * 2 * pi * frequency_2 / total_time)+1);
% sinusoid_pos_3 = amplitude_3 * (-cos(time * 2 * pi * frequency_3 / total_time)+1);
% sinusoid_vel_1 = amplitude_1 * 2 * pi * frequency_1 / total_time * sin(time * 2 * pi * frequency_1 / total_time);
% sinusoid_vel_2 = amplitude_2 * 2 * pi * frequency_2 / total_time * sin(time * 2 * pi * frequency_2 / total_time);
% sinusoid_vel_3 = amplitude_3 * 2 * pi * frequency_3 / total_time * sin(time * 2 * pi * frequency_3 / total_time);
% sinusoid_acc_1 = amplitude_1 * (2 * pi * frequency_1 / total_time)^2 * cos(time * 2 * pi * frequency_1 / total_time);
% sinusoid_acc_2 = amplitude_2 * (2 * pi * frequency_2 / total_time)^2 * cos(time * 2 * pi * frequency_2 / total_time);
% sinusoid_acc_3 = amplitude_3 * (2 * pi * frequency_3 / total_time)^2 * cos(time * 2 * pi * frequency_3 / total_time);

joint_angle_trajectories = ...
  ( ...
    theta_init * ones(size(sinusoid_pos_1)) ...
    + theta_direction_1 * sinusoid_pos_1 ...
  )';
%     + theta_direction_2 * sinusoid_pos_2 ...
%     + theta_direction_3 * sinusoid_pos_3 ...
% joint_velocity_trajectories{2} = ...
%   ( ...
%     theta_direction_1 * sinusoid_vel_1 ...
%     + theta_direction_2 * sinusoid_vel_2 ...
%     + theta_direction_3 * sinusoid_vel_3 ...
%   )';
% joint_acceleration_trajectories{2} = ...
%   ( ...
%     theta_direction_1 * sinusoid_acc_1 ...
%     + theta_direction_2 * sinusoid_acc_2 ...
%     + theta_direction_3 * sinusoid_acc_3 ...
%   )';



% create marker trajectories
marker_trajectories = zeros(length(time_mocap), size(kinematic_tree.exportMarkerPositions, 2));
for i_time = 1 : length(time_mocap)
    theta = joint_angle_trajectories(i_time, :)';
    kinematic_tree.jointAngles = theta;
    kinematic_tree.updateConfiguration;
    
    marker_positions = kinematic_tree.exportMarkerPositions;
    marker_trajectories(i_time, :) = marker_positions;
end

joint_angle_trajectories_simulated = joint_angle_trajectories;

% save trajectories
save_folder = 'processed';
date = '00000000';
subject_id = 'XXX';
marker_labels = kinematic_tree.markerLabels;
sampling_rate_mocap = 1 / time_step;
save_file_name = makeFileName(date, subject_id, 'simulation', 1, 'markerTrajectories.mat');
save ...
  ( ...
    [save_folder filesep save_file_name], ...
    'marker_trajectories', ...
    'joint_angle_trajectories_simulated', ...
    'time_mocap', ...
    'sampling_rate_mocap', ...
    'marker_labels' ...
  );
addAvailableData('marker_trajectories', 'time_mocap', 'sampling_rate_mocap', 'marker_labels', save_folder, save_file_name);


clear







