%     This file is part of the KinematicChain library
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

% check out constraints

real_joint_left_pos = [-2; 0; 2];
real_joint_right_pos = [2; 0; 2];
real_joint_left_axis_tree = [0; 1; 0];
real_joint_right_axis = [0; 1; 0];
virtual_joint_tree_pos = [0; 0; 2];
virtual_joint_one_axis = [1; 0; 0];
virtual_joint_two_axis = [0; 0; 1];
virtual_joint_three_axis = [0; 1; 0];
center_link_position = [0; 0; 0;];
left_link_position = [-2; 0; 0;];
right_link_position = [2; 0; 0;];

joint_positions_tree = {virtual_joint_tree_pos, virtual_joint_tree_pos, virtual_joint_tree_pos, real_joint_left_pos, real_joint_right_pos};
joint_axes_tree = {virtual_joint_one_axis, virtual_joint_two_axis, virtual_joint_three_axis, real_joint_left_axis_tree, real_joint_right_axis};
end_effectors_tree_pos = {[-1; 0; 0]; [3; 0; 0];};
link_positions_tree = {virtual_joint_tree_pos, virtual_joint_tree_pos, center_link_position, left_link_position, right_link_position};
link_masses_tree = [0 0 2 1 1];
link_moments_of_inertia_tree = [0 0 0; 0 0 0; 5 5 5; 1 1 1; 1 1 1] * 0.1;
branch_matrix_tree = ...
  [ ...
    1 1 1 1 0; ...
    1 1 1 0 1; ...
  ];


link_orientations = {eye(3), eye(3), eye(3), eye(3), eye(3)};
joint_types = [2 2 1 1 1];

kinematic_tree = GeneralKinematicTree ...
( ...
  joint_positions_tree, ...
  joint_axes_tree, ...
  joint_types, ...
  branch_matrix_tree, ...
  end_effectors_tree_pos, ...
  link_positions_tree, ...
  link_orientations, ...
  link_masses_tree, ...
  link_moments_of_inertia_tree ...
);



%% check out instantaneous torques 
M = kinematic_tree.inertiaMatrix;
C = kinematic_tree.coriolisMatrix;
N = kinematic_tree.gravitationalTorqueMatrix;
A = [kinematic_tree.endEffectorJacobians{1}([1 3], :); kinematic_tree.endEffectorJacobians{2}([1 3], :)];
A_dot = [kinematic_tree.endEffectorJacobianTemporalDerivatives{1}([1 3], :); kinematic_tree.endEffectorJacobianTemporalDerivatives{2}([1 3], :)];
P = eye(kinematic_tree.numberOfJoints) - A' * (A * M^(-1) * A')^(-1) * A * M^(-1);
k_v = 3;
number_of_joints = kinematic_tree.numberOfJoints;


theta_dot = kinematic_tree.jointVelocities;
desired_joint_accelerations = zeros(size(kinematic_tree.jointAngles));
b_acc = (M*desired_joint_accelerations + P*C*theta_dot + P*N + A'*(A*M^(-1)*A')^(-1)*(A_dot * theta_dot));

%% solution one - only contact constraints
H = ...
  [ ...
    P; ...
  ];
b = ...
  [ ...
    b_acc; ...
  ];
T = pinv(H) * b;
lambda = (A*M^(-1)*A')^(-1) ...
    * (A*M^(-1)*(T - C*theta_dot - N) + A_dot*theta_dot);
kinematic_tree.externalTorques = T;
kinematic_tree.constraintTorques = A'*lambda;
kinematic_tree.calculateAccelerationsFromExternalTorques;
accelerations_check = kinematic_tree.jointAccelerations;

% check for workless torques
k_constraints_full = rank(A);
[~, ~, W_constraints] = svd(A);                                             % use singular value decomposition to get the basis
B_constraints = W_constraints(:, 1:k_constraints_full);                     % E_rng contains the base vectors of the range space - the torque combinations that do no work
C_constraints = W_constraints(:, k_constraints_full+1:end);                 % E_nul contains the base vectors of the null space - the torque combinations that do work
P_constraint_full = B_constraints*B_constraints';                           % K_rng projects onto the range space
P_constraint_free = C_constraints*C_constraints';                           % K_nul projects onto the null space
T_workpure = P_constraint_free * T;
T_workless = P_constraint_full * T;
T_check = T_workpure + T_workless;

% check for workless torques - but only those that also have no virtual torques
B = [eye(k_v) zeros(k_v, number_of_joints - k_v)];
k_c = rank(A);
[~, ~, V_c] = svd(A);
C_c = V_c(:, k_c+1:end); % C_c spans the null space of A
k_w = rank([B' C_c]);
[~, ~, V_w] = svd([B' C_c]');
B_w = V_w(:, 1:k_w);                                                        % B_w columns span the space of torque vectors that do work, or have a virtual component, or both
C_w = V_w(:, k_w+1:end);                                                    % C_w columns span the space of torque vectors that do no work and have no virtual component                        
P_workpure = B_w*B_w';                                                      % P_workpure projects onto the space spanned by B_w columns
P_workless = C_w*C_w';                                                      % P_workless projects onto the space spanned by C_w columns
T_workpure = P_workpure * T;
T_workless = P_workless * T;
T_check = T_workpure + T_workless;


ground_reaction_wrench = calculateInstantaneousGroundReactionWrench(kinematic_tree, kinematic_tree.constraintTorques, eye(4));

%% solution two - contact constraints and no torques in the free DoFs
H = ...
  [ ...
    P; ...
    eye(k_v, number_of_joints); ...
  ];
b = ...
  [ ...
    b_acc; ...
    zeros(k_v, 1); ...
  ];

T_2 = pinv(H) * b;
lambda = (A*M^(-1)*A')^(-1) ...
    * (A*M^(-1)*(T_2 - C*theta_dot - N) + A_dot*theta_dot);
kinematic_tree.externalTorques = T_2;
kinematic_tree.constraintTorques = A'*lambda;
kinematic_tree.calculateAccelerationsFromExternalTorques;
accelerations_check = kinematic_tree.jointAccelerations;

% check for workless torques
k_constraints_full = rank(A);
[~, ~, W_constraints] = svd(A);                                             % use singular value decomposition to get the basis
B_constraints = W_constraints(:, 1:k_constraints_full);                     % E_rng contains the base vectors of the range space - the torque combinations that do no work
C_constraints = W_constraints(:, k_constraints_full+1:end);                 % E_nul contains the base vectors of the null space - the torque combinations that do work
P_constraint_full = B_constraints*B_constraints';                           % K_rng projects onto the range space
P_constraint_free = C_constraints*C_constraints';                           % K_nul projects onto the null space
T_workpure_2 = P_constraint_free * T_2;
T_workless_2 = P_constraint_full * T_2;
T_check_2 = T_workpure_2 + T_workless_2;

% check for workless torques - but only those that also have no virtual torques
B = [eye(k_v) zeros(k_v, number_of_joints - k_v)];
k_c = rank(A);
[~, ~, V_c] = svd(A);
C_c = V_c(:, k_c+1:end); % C_c spans the null space of A
k_w = rank([B' C_c]);
[~, ~, V_w] = svd([B' C_c]');
B_w = V_w(:, 1:k_w);                                                        % B_w columns span the space of torque vectors that do work, or have a virtual component, or both
C_w = V_w(:, k_w+1:end);                                                    % C_w columns span the space of torque vectors that do no work and have no virtual component                        
P_workpure = B_w*B_w';                                                      % P_workpure projects onto the space spanned by B_w columns
P_workless = C_w*C_w';                                                      % P_workless projects onto the space spanned by C_w columns
T_workpure_2 = P_workpure * T_2;
T_workless_2 = P_workless * T_2;
T_check_2 = T_workpure_2 + T_workless_2;

% check what happens if I remove workless torques
lambda = (A*M^(-1)*A')^(-1) ...
    * (A*M^(-1)*(T_workpure_2 - C*theta_dot - N) + A_dot*theta_dot);
kinematic_tree.externalTorques = T_workpure_2;
kinematic_tree.constraintTorques = A'*lambda;
kinematic_tree.calculateAccelerationsFromExternalTorques;
accelerations_check = kinematic_tree.jointAccelerations;


%% solution three - contact constraints, no workless torques

B = [eye(k_v) zeros(k_v, number_of_joints - k_v)];
k_c = rank(A);
[~, ~, V_c] = svd(A);
C_c = V_c(:, k_c+1:end); % C_c spans the null space of A
k_w = rank([B' C_c]);
[~, ~, V_w] = svd([B' C_c]');
C_w = V_w(:, k_w+1:end);                        
% combine
H = ...
  [ ...
    P; ...
    C_w'; ...
  ];

b = ...
  [ ...
    b_acc; ...
    zeros(number_of_joints - k_w, 1); ...
  ];
T_3 = pinv(H) * b;
lambda = (A*M^(-1)*A')^(-1) ...
    * (A*M^(-1)*(T_3 - C*theta_dot - N) + A_dot*theta_dot);
kinematic_tree.externalTorques = T_3;
kinematic_tree.constraintTorques = A'*lambda;
kinematic_tree.calculateAccelerationsFromExternalTorques;
accelerations_check = kinematic_tree.jointAccelerations;

% check for workless torques
k_constraints_full = rank(A);
[~, ~, W_constraints] = svd(A);                                             % use singular value decomposition to get the basis
B_constraints = W_constraints(:, 1:k_constraints_full);                     % E_rng contains the base vectors of the range space - the torque combinations that do no work
C_constraints = W_constraints(:, k_constraints_full+1:end);                 % E_nul contains the base vectors of the null space - the torque combinations that do work
P_constraint_full = B_constraints*B_constraints';                           % P_constraint_full projects onto the range space
P_constraint_free = C_constraints*C_constraints';                           % P_constraint_free projects onto the null space
T_workpure_3 = P_constraint_free * T_3;
T_workless_3 = P_constraint_full * T_3;
T_check_3 = T_workpure_3 + T_workless_3;

% check for workless torques - but only those that also have no virtual torques
B = [eye(k_v) zeros(k_v, number_of_joints - k_v)];
k_c = rank(A);
[~, ~, V_c] = svd(A);
C_c = V_c(:, k_c+1:end); % C_c spans the null space of A
k_w = rank([B' C_c]);
[~, ~, V_w] = svd([B' C_c]');
B_w = V_w(:, 1:k_w);                                                        % B_w columns span the space of torque vectors that do work, or have a virtual component, or both
C_w = V_w(:, k_w+1:end);                                                    % C_w columns span the space of torque vectors that do no work and have no virtual component                        
P_workpure = B_w*B_w';                                                      % P_workpure projects onto the space spanned by B_w columns
P_workless = C_w*C_w';                                                      % P_workless projects onto the space spanned by C_w columns
T_workpure_3 = P_workpure * T_3;
T_workless_3 = P_workless * T_3;
T_check_3 = T_workpure_3 + T_workless_3;

% check what happens if I remove workless torques
lambda = (A*M^(-1)*A')^(-1) ...
    * (A*M^(-1)*(T_workpure_3 - C*theta_dot - N) + A_dot*theta_dot);
kinematic_tree.externalTorques = T_workpure_3;
kinematic_tree.constraintTorques = A'*lambda;
kinematic_tree.calculateAccelerationsFromExternalTorques;
accelerations_check = kinematic_tree.jointAccelerations;


%% solution four - contact constraints, no torques in the free DoFs and no workless torques
B = eye(k_v, number_of_joints);
k_c = rank(A);
[~, ~, V_c] = svd(A);
C_c = V_c(:, k_c+1:end); % C_c spans the null space of A
k_w = rank([B' C_c]);
[~, ~, V_w] = svd([B' C_c]');
C_w = V_w(:, k_w+1:end);                        
% combine
H = ...
  [ ...
    P; ...
    B; ...
    C_w'; ...
  ];


% this does not work, apparently I have either virtual torques or workless torques

b = ...
  [ ...
    b_acc; ...
    zeros(k_v, 1); ...
    zeros(number_of_joints - k_w, 1); ...
  ];


T_4 = pinv(H) * b;
lambda = (A*M^(-1)*A')^(-1) ...
    * (A*M^(-1)*(T_4 - C*theta_dot - N) + A_dot*theta_dot);
kinematic_tree.externalTorques = T_4;
kinematic_tree.constraintTorques = A'*lambda;
kinematic_tree.calculateAccelerationsFromExternalTorques;
accelerations_check = kinematic_tree.jointAccelerations;

% check for workless torques
k_constraints_full = rank(A);
[~, ~, W_constraints] = svd(A);                                             % use singular value decomposition to get the basis
B_constraints = W_constraints(:, 1:k_constraints_full);                     % E_rng contains the base vectors of the range space - the torque combinations that do no work
C_constraints = W_constraints(:, k_constraints_full+1:end);                 % E_nul contains the base vectors of the null space - the torque combinations that do work
P_constraint_full = B_constraints*B_constraints';                           % P_constraint_full projects onto the range space
P_constraint_free = C_constraints*C_constraints';                           % P_constraint_free projects onto the null space
T_workpure_4 = P_constraint_free * T_4;
T_workless_4 = P_constraint_full * T_4;
T_check_4 = T_workpure_4 + T_workless_4;

% check for workless torques - but only those that also have no virtual torques
B = [eye(k_v) zeros(k_v, number_of_joints - k_v)];
k_c = rank(A);
[~, ~, V_c] = svd(A);
C_c = V_c(:, k_c+1:end); % C_c spans the null space of A
k_w = rank([B' C_c]);
[~, ~, V_w] = svd([B' C_c]');
B_w = V_w(:, 1:k_w);                                                        % B_w columns span the space of torque vectors that do work, or have a virtual component, or both
C_w = V_w(:, k_w+1:end);                                                    % C_w columns span the space of torque vectors that do no work and have no virtual component                        
P_workpure = B_w*B_w';                                                      % P_workpure projects onto the space spanned by B_w columns
P_workless = C_w*C_w';                                                      % P_workless projects onto the space spanned by C_w columns
T_workpure_4 = P_workpure * T_4;
T_workless_4 = P_workless * T_4;
T_check_4 = T_workpure_4 + T_workless_4;


% check what happens if I remove workless torques
lambda = (A*M^(-1)*A')^(-1) ...
    * (A*M^(-1)*(T_workpure_4 - C*theta_dot - N) + A_dot*theta_dot);
kinematic_tree.externalTorques = T_workpure_4;
kinematic_tree.constraintTorques = A'*lambda;
kinematic_tree.calculateAccelerationsFromExternalTorques;
accelerations_check = kinematic_tree.jointAccelerations;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I was here before the weekend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return







% do a longer simulation

scene_bound = 3.5*[-1 1; -1 1; -1 1];

stick_figure_tree = KinematicTreeStickFigure(kinematic_tree, scene_bound);
% stick_figure_tree = KinematicTreeStickFigure(kinematic_tree, scene_bound);
% stick_figure_tree.showLinkMassEllipsoids = true;
position = get(stick_figure_tree.sceneFigure, 'Position');
position(1) = position(1) + position(3);
set(stick_figure_tree.sceneFigure, 'Position', position);

stick_figure_tree.setLinkPlotsLinewidth(5);
stick_figure_tree.setLinkPlotsColor([1 0.7 0]);
stick_figure_tree.update();

kinematic_tree.externalTorques = [0; 0; 0; 2; 0] * 1;
time_step = 0.005;
counter = 0;
total_time = 10;
time = time_step : time_step : total_time;
number_of_time_steps = length(time);

joint_angles_tree = zeros(number_of_time_steps, 5);
joint_velocities_tree = zeros(number_of_time_steps, 5);
joint_accelerations_tree = zeros(number_of_time_steps, 5);
joint_torques_tree = zeros(number_of_time_steps, 5);
contact_forces_tree = zeros(number_of_time_steps, 4);
for i_time = 1 : length(time);
    % apply constraints
    A_tree = [kinematic_tree.endEffectorJacobians{1}([1 3], :); kinematic_tree.endEffectorJacobians{2}([1 3], :)];
    A_treeDot = [kinematic_tree.endEffectorJacobianTemporalDerivatives{1}([1 3], :); kinematic_tree.endEffectorJacobianTemporalDerivatives{2}([1 3], :)];
    M_tree = kinematic_tree.inertiaMatrix;
    lambda_tree = (A_tree*M_tree^(-1)*A_tree')^(-1) ...
        * (A_tree*M_tree^(-1)*(kinematic_tree.externalTorques - kinematic_tree.coriolisMatrix*kinematic_tree.jointVelocities - kinematic_tree.gravitationalTorqueMatrix) + A_treeDot*kinematic_tree.jointVelocities);
    kinematic_tree.constraintTorques = A_tree'*lambda_tree;
    constraint_torques_tree = A_tree'*lambda_tree;
    
    % make Euler step
    kinematic_tree.calculateAccelerationsFromExternalTorques;
    theta_two_dot = kinematic_tree.jointAccelerations;
    kinematic_tree.jointVelocities = kinematic_tree.jointVelocities + time_step*kinematic_tree.jointAccelerations;
    kinematic_tree.jointAngles = kinematic_tree.jointAngles + time_step*kinematic_tree.jointVelocities + time_step^2*kinematic_tree.jointAccelerations;
    kinematic_tree.updateInternals;
    
    % store data
    joint_angles_tree(i_time, :) = kinematic_tree.jointAngles;
    joint_velocities_tree(i_time, :) = kinematic_tree.jointVelocities;
    joint_accelerations_tree(i_time, :) = kinematic_tree.jointAccelerations;
    joint_torques_tree(i_time, :) = kinematic_tree.externalTorques;
    contact_forces_tree(i_time, :) = lambda_tree;

    counter = counter+1;
    if counter == 5
        counter = 0;
        stick_figure_tree.update;
        drawnow;
    end
end

figure; axes; hold on; title('joint angles');
plot(time, joint_angles_tree(:, 1), 'r:', 'linewidth', 2, 'displayname', 'tree 1');
plot(time, joint_angles_tree(:, 2), 'g:', 'linewidth', 2, 'displayname', 'tree 2');
plot(time, joint_angles_tree(:, 3), 'b:', 'linewidth', 2, 'displayname', 'tree 3');
plot(time, joint_angles_tree(:, 4), 'c:', 'linewidth', 2, 'displayname', 'tree 4');
plot(time, joint_angles_tree(:, 5), 'm:', 'linewidth', 2, 'displayname', 'tree 5');
legend('toggle')


figure; axes; hold on; title('contact forces');
plot(time, contact_forces_tree(:, 1), 'r:', 'linewidth', 1, 'displayname', 'tree 1');
plot(time, contact_forces_tree(:, 2), 'g:', 'linewidth', 1, 'displayname', 'tree 1');
plot(time, contact_forces_tree(:, 3), 'b:', 'linewidth', 1, 'displayname', 'tree 1');


