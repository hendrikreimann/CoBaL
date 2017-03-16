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
total_time = 0.5;
time_mocap = time_step : time_step : total_time;
number_of_time_steps = length(time_mocap);

% define joint t trajectories
%     0.1; 0.2; 0.3; 0.4; 0.5; 0.6; % pelvis free body DoFs
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; % pelvis free body DoFs
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; ... % left leg
%     0.11; 0.21; 0.31; 0.41; 0.51; 0.61; 0.71; ... % left leg
%     0.12; 0.22; 0.32; 0.42; 0.52; 0.62; 0.72; ... % right leg
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; % pelvis free body DoFs
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.5; ... % left leg
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; % lumbar and cervical
%     0.13; 0.23; 0.33; 0.43; 0.53; 0.63; % lumbar and cervical
%     0.14; 0.24; 0.34; 0.44; 0.54; 0.64; 0.74; ... % left arm
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; ... left arm
% tau_init = ...
%   [ ...
%     0.1; 0.2; 0.3; 0.4; 0.5; 0.6; % pelvis free body DoFs
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; ... % left leg
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; ... % right leg
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; % lumbar and cervical
%     0.14; 0.24; 0.34; 0.44; 0.54; 0.64; ... left arm
%     0; 0; 0; 0; 0; 0; ... right arm
%   ];
tau_init = ...
  [ ...
    0.1; 0.2; 0.3; 0.4; 0.5; 0.6; % pelvis free body DoFs
    0.11; 0.21; 0.31; 0.41; 0.51; 0.61; 0.71; ... % left leg
    0.12; 0.22; 0.32; 0.42; 0.52; 0.62; 0.72; ... % right leg
    0.13; 0.23; 0.33; 0.43; 0.53; 0.63; % lumbar and cervical
    0.14; 0.24; 0.34; 0.44; 0.54; 0.64; ... % left arm
    0.15; 0.25; 0.35; 0.45; 0.55; 0.65; ... % right arm
  ];

%     0.1; 0.2; 0.3; 0.5; 0.0; 0.0; 0.0; ... % left leg
% tau_direction_1 = ...
%   [ ...
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; % pelvis free body DoFs
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; ... % left leg
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; ... % right leg
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; % lumbar and cervical
%     0.0; 0.0; 0.0; 0.0; 0.0; 0.0; ... left arm
%     0; 0; 0; 0; 0; 0; ... right arm
%   ];
tau_direction_1 = ...
  [ ...
    0.1; 0.2; 0.3; 0.4; 0.5; 0.6; % pelvis free body DoFs
    0.11; 0.21; 0.31; 0.41; 0.51; 0.61; 0.71; ... % left leg
    0.12; 0.22; 0.32; 0.42; 0.52; 0.62; 0.72; ... % right leg
    0.13; 0.23; 0.33; 0.43; 0.53; 0.63; % lumbar and cervical
    0.14; 0.24; 0.34; 0.44; 0.54; 0.64; ... % left arm
    0.15; 0.25; 0.35; 0.45; 0.55; 0.65; ... % right arm
  ];

% tau_init = rand(size(tau_init));
% tau_direction_1 = rand(size(tau_direction_1));
% tau_init = ...
%   [ ...
%     0.1; 0.2; 0.3; 0.4; 0.5; 0.6; % pelvis free body DoFs
%     .5; .1; .3; .15; .2; .4; 0.3; ... % left leg
%     .0; .0; .0; .0; .0; .0; 0; ... % right leg
%     0; 0; 0; 0; 0; 0; ... % lumbar and cervical
%     0; 0; 0; 0; 0; 0; 0; ... left arm
%     0; 0; 0; 0; 0; 0; 0; ... right arm
%   ];
% 
% tau_direction_1 = ...
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

torque_unconstrained_trajectories = ...
  ( ...
    tau_init * ones(size(sinusoid_pos_1)) ...
    + tau_direction_1 * sinusoid_pos_1 ...
  )';

right_foot_constraint = 0;
left_foot_constraint = 1;
for i_time = 1 : number_of_time_steps
    kinematic_tree.externalTorques = torque_unconstrained_trajectories(i_time, :)';

    % define the constraint matrix
    [A, ADot] = createConstraintMatrix_pointConstraints(kinematic_tree, right_foot_constraint, left_foot_constraint);
    A = kinematic_tree.bodyJacobians{3};
    ADot = kinematic_tree.bodyJacobianTemporalDerivatives{3};

    vector_to_center = kinematic_tree.endEffectorPositions{3} - constraint_sphere_center;
    q = vector_to_center;
    d_q_norm_by_d_p = q' * 1 / norm(q);
    d_p_by_d_theta = kinematic_tree.endEffectorJacobians{3};
    J_p = d_p_by_d_theta;
    A = d_q_norm_by_d_p * d_p_by_d_theta;
    JDot = kinematic_tree.endEffectorJacobianTemporalDerivatives{3};
    A = J;
    ADot = d_q_norm_by_d_p * JDot;

    % remove the workless part from the applied torques
    external_torques = kinematic_tree.externalTorques;
    k_constraints_full = rank(A);
    [~, ~, W_constraints] = svd(A);                 % use singular value decomposition to get the basis
    B_constraints = W_constraints(:, 1:k_constraints_full);                         % E_rng contains the base vectors of the range space - the torque combinations that do no work
    C_constraints = W_constraints(:, k_constraints_full+1:end);                     % E_nul contains the base vectors of the null space - the torque combinations that do work
    P_constraint_full = B_constraints*B_constraints';                           % K_rng projects onto the range space
    P_constraint_free = C_constraints*C_constraints';                           % K_nul projects onto the null space
    external_torques_workpure = P_constraint_free * external_torques;
    external_torques_workless = P_constraint_full * external_torques;
    external_torques_check = external_torques_workpure + external_torques_workless;

    B_virtualspace = [eye(6); zeros(kinematic_tree.numberOfJoints - 6, 6)];             % E_joint contains the base vectors of the joint torque space
    span_virtualspace_or_constraint_free = [B_virtualspace C_constraints];
    k_virtualspace_or_constraint_free = rank(span_virtualspace_or_constraint_free);
    span_virtualspace_or_constraint_full = [B_virtualspace B_constraints];
    k_virtualspace_or_constraint_full = rank(span_virtualspace_or_constraint_full);

    B_jointspace = [zeros(6, 6); eye(18, 6)];             % E_joint contains the base vectors of the joint torque space
    span_jointspace_or_constraint_free = [B_jointspace C_constraints];
    k_jointspace_or_constraint_free = rank(span_virtualspace_or_constraint_free);
    span_jointspace_or_constraint_full = [B_jointspace B_constraints];
    k_jointspace_or_constraint_full = rank(span_virtualspace_or_constraint_full);

    [~, d, W_virtualspace_or_constraint_free] = svd(span_virtualspace_or_constraint_free');
    B_virtualspace_or_constraint_free = W_virtualspace_or_constraint_free(:, 1:k_virtualspace_or_constraint_free);
    C_virtualspace_or_constraint_free = W_virtualspace_or_constraint_free(:, k_virtualspace_or_constraint_free+1:end);

    constraints_check = P_constraint_full * C_virtualspace_or_constraint_free; % should not be 0 if C_... pushes against the constraints
    work_check = P_constraint_free * C_virtualspace_or_constraint_free; % should not be 0 if C_... does some work
    B_jointspace_and_constraint_full = C_virtualspace_or_constraint_free; % this is the direction of joint space that does no work

    % remove the part that does no work from the joint torques
    projector = eye(24, 24) - (B_jointspace_and_constraint_full*B_jointspace_and_constraint_full');
    external_torques_new = external_torques - (B_jointspace_and_constraint_full*B_jointspace_and_constraint_full')*external_torques;



    kinematic_tree.externalTorques = external_torques_new;










    M = kinematic_tree.inertiaMatrix;
    lambda = (A*M^(-1)*A')^(-1) ...
        * (A*M^(-1)*(kinematic_tree.externalTorques - kinematic_tree.coriolisMatrix*kinematic_tree.jointVelocities - kinematic_tree.gravitationalTorqueMatrix) + ADot*kinematic_tree.jointVelocities);
    kinematic_tree.constraintTorques = A'*lambda;
    constraint_torques = A'*lambda;

    kinematic_tree.calculateAccelerationsFromExternalTorques;


    % check
%         lambda_reconstructed = inverseDynamicsDataPoint(kinematic_tree, A, ADot, virtual_joints);




    accelerations_check = kinematic_tree.jointAccelerations;
    kinematic_tree.jointVelocities = kinematic_tree.jointVelocities + time_step*kinematic_tree.jointAccelerations;
    kinematic_tree.jointAngles = kinematic_tree.jointAngles + time_step*kinematic_tree.jointVelocities + time_step^2*kinematic_tree.jointAccelerations;
    kinematic_tree.updateInternals;



    timeseries_inverse_dynamics_joint_angle(i_time, :) = kinematic_tree.jointAngles';
    timeseries_inverse_dynamics_joint_velocity(i_time, :) = kinematic_tree.jointVelocities';
    timeseries_inverse_dynamics_joint_acceleration(i_time, :) = kinematic_tree.jointAccelerations';
    timeseries_inverse_dynamics_joint_torques_applied(i_time, :) = kinematic_tree.externalTorques';
    timeseries_inverse_dynamics_constraint_torques_applied(i_time, :) = kinematic_tree.constraintTorques';
    timeseries_inverse_dynamics_end_effector_positions(i_time, 1:3) = kinematic_tree.endEffectorPositions{1}';
    timeseries_inverse_dynamics_end_effector_positions(i_time, 4:6) = kinematic_tree.endEffectorPositions{2}';
    timeseries_inverse_dynamics_contact_forces_applied(i_time, 1 : length(lambda)) = lambda;

    display_step = 5;
    if (i_time / display_step) == floor(i_time / display_step)
        disp([num2str(i_time) '(' num2str(number_of_time_steps) ')']);
    end
end











% create marker trajectories
marker_trajectories = zeros(length(time_mocap), size(kinematic_tree.exportMarkerPositions, 2));
for i_time = 1 : length(time_mocap)
    theta = torque_unconstrained_trajectories(i_time, :)';
    kinematic_tree.jointAngles = theta;
    kinematic_tree.updateConfiguration;
    
    marker_positions = kinematic_tree.exportMarkerPositions;
    marker_trajectories(i_time, :) = marker_positions;
end

% add noise
noise_strength = 0.002;
marker_trajectories = marker_trajectories + randn(size(marker_trajectories)) * noise_strength;

joint_angle_trajectories_simulated = torque_unconstrained_trajectories;

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







