
% load('subjectModel.mat')
% kinematic_tree.updateInternals;

number_of_joints = kinematic_tree.numberOfJoints;
k_v = 6;
M = kinematic_tree.inertiaMatrix;
C = kinematic_tree.coriolisMatrix;
N = kinematic_tree.gravitationalTorqueMatrix;
A = [kinematic_tree.bodyJacobians{3}; kinematic_tree.bodyJacobians{6}];
A_dot = [kinematic_tree.bodyJacobianTemporalDerivatives{3}; kinematic_tree.bodyJacobianTemporalDerivatives{6}];
P = eye(kinematic_tree.numberOfJoints) - A' * (A * M^(-1) * A')^(-1) * A * M^(-1);

% find applied torques that result in zero acceleration
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
[~, ~, W_constraints] = svd(A);                 % use singular value decomposition to get the basis
B_constraints = W_constraints(:, 1:k_constraints_full);                         % E_rng contains the base vectors of the range space - the torque combinations that do no work
C_constraints = W_constraints(:, k_constraints_full+1:end);                     % E_nul contains the base vectors of the null space - the torque combinations that do work
P_constraint_full = B_constraints*B_constraints';                           % K_rng projects onto the range space
P_constraint_free = C_constraints*C_constraints';                           % K_nul projects onto the null space
T_workpure = P_constraint_free * T;
T_workless = P_constraint_full * T;
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

T = pinv(H) * b;
lambda = (A*M^(-1)*A')^(-1) ...
    * (A*M^(-1)*(T - C*theta_dot - N) + A_dot*theta_dot);
kinematic_tree.externalTorques = T;
kinematic_tree.constraintTorques = A'*lambda;
kinematic_tree.calculateAccelerationsFromExternalTorques;
accelerations_check = kinematic_tree.jointAccelerations;

% check for workless torques
k_constraints_full = rank(A);
[~, ~, W_constraints] = svd(A);                 % use singular value decomposition to get the basis
B_constraints = W_constraints(:, 1:k_constraints_full);                         % E_rng contains the base vectors of the range space - the torque combinations that do no work
C_constraints = W_constraints(:, k_constraints_full+1:end);                     % E_nul contains the base vectors of the null space - the torque combinations that do work
P_constraint_full = B_constraints*B_constraints';                           % K_rng projects onto the range space
P_constraint_free = C_constraints*C_constraints';                           % K_nul projects onto the null space
T_workpure = P_constraint_free * T;
T_workless = P_constraint_full * T;
T_check = T_workpure + T_workless;

%% solution three - contact constraints, no torques in the free DoFs and no workless torques
% this doesn't seem to work at the moment

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
    eye(k_v, number_of_joints); ...
    C_w'; ...
  ];

b = ...
  [ ...
    b_acc; ...
    zeros(k_v, 1); ...
    zeros(number_of_joints - k_w, 1); ...
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
[~, ~, W_constraints] = svd(A);                 % use singular value decomposition to get the basis
B_constraints = W_constraints(:, 1:k_constraints_full);                         % E_rng contains the base vectors of the range space - the torque combinations that do no work
C_constraints = W_constraints(:, k_constraints_full+1:end);                     % E_nul contains the base vectors of the null space - the torque combinations that do work
P_constraint_full = B_constraints*B_constraints';                           % K_rng projects onto the range space
P_constraint_free = C_constraints*C_constraints';                           % K_nul projects onto the null space
T_workpure = P_constraint_free * T;
T_workless = P_constraint_full * T;
T_check = T_workpure + T_workless;







return







%% old stuff, copied from somewhere else but not needed anymore I believe


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












% check
%         lambda_reconstructed = inverseDynamicsDataPoint(kinematic_tree, A, ADot, virtual_joints);




accelerations_check = kinematic_tree.jointAccelerations;
kinematic_tree.jointVelocities = kinematic_tree.jointVelocities + time_step*kinematic_tree.jointAccelerations;
kinematic_tree.jointAngles = kinematic_tree.jointAngles + time_step*kinematic_tree.jointVelocities + time_step^2*kinematic_tree.jointAccelerations;
kinematic_tree.updateInternals;






