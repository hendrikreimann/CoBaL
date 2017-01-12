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


% test calculation of the ground reaction forces and torques


% load and initialize
cd /Users/reimajbi/TempleDrive/20160517_HR_balanceBeam/processed

load subjectInfo.mat;
model_file_name = makeFileName(date, subject_id, 'model');
load(model_file_name);

plant.jointVelocities = rand(plant.numberOfJoints, 1);
plant.updateInternals;

% calculate constraints
A_left = plant.bodyJacobians{3};
A_left_dot = plant.bodyJacobianTemporalDerivatives{3};
A_right = plant.bodyJacobians{6};
A_right_dot = plant.bodyJacobianTemporalDerivatives{6};

% A_left = [];
% A_left_dot = [];
% A_right = [];
% A_right_dot = [];

number_of_left_foot_constraints = size(A_left, 1);
number_of_right_foot_constraints = size(A_right, 1);

A = [A_left; A_right];
A_dot = [A_left_dot; A_right_dot];

% try hinge constraints
[A, A_dot] = createConstraintMatrix_hingeConstraints(plant, 1, 2);
number_of_left_foot_constraints = 5;
number_of_right_foot_constraints = 5;
A_left = A(1 : number_of_left_foot_constraints, :);
A_left_dot = A_dot(1 : number_of_left_foot_constraints, :);
A_right = A(number_of_left_foot_constraints + 1 : end, :);
A_right_dot = A_dot(number_of_left_foot_constraints + 1 : end, :);

M = plant.inertiaMatrix;
C = plant.coriolisMatrix;
N = plant.gravitationalTorqueMatrix;

% find joint torques for zero accelerations
P = eye(plant.numberOfJoints) - A' * (A * M^(-1) * A')^(-1) * A * M^(-1);

% this next part is for zero workless torques
virtual_joints = 1 : 6;
k_c = rank(A);
k_v = length(virtual_joints);
[~, ~, V_c] = svd(A);
C_c = V_c(:, k_c+1:end); % C_c spans the null space of A
B_v = [eye(k_v); zeros(plant.numberOfJoints - k_v, k_v)];
k_w = rank([B_v C_c]);
[~, ~, V_w] = svd([B_v C_c]');
C_w = V_w(:, k_w+1:end);
% D = [B_v C_w];


theta_two_dot_des = zeros(plant.numberOfJoints, 1);
H = ...
    [ ...
      P; ...                                                    % desired joint acceleration
      B_v'; ...                                                 % torques in the free body dofs
      C_w'; ...                                                 % workless torques
    ]; ...       
    
T = pinv ...
  ( ...
    [ ...
      P; ...                                                    % desired joint acceleration
      B_v'; ...                                                 % torques in the free body dofs
      C_w'; ...                                                 % workless torques
    ] ...       
  ) * ...
  [ ...
    P*C*plant.jointVelocities + P*N + M*theta_two_dot_des; ...  % desired joint acceleration
    zeros(size(B_v, 2), 1); ...                                            % torques in free dofs
    zeros(size(C_w, 2), 1); ...                                            % workless torques
  ];

% try applying the torque that keeps the body still for only the right foot constrained, but with both feet constrained
% T = T_right_only;
% this works, the accelerations are still zero, as they should be
% but why do I not get a torque that actually gives me zero accelerations for both feet constrained?

% try this without constraining the free DoFs
% T = pinv ...
%   ( ...
%     P ...                                                    % desired joint acceleration
%   ) * ...
%     P*C*plant.jointVelocities + P*N + M*theta_two_dot_des; ...  % desired joint acceleration


% apply constraints and check resulting joint accelerations
lambda = (A*M^(-1)*A')^(-1) ...
    * (A*M^(-1)*(T - C*plant.jointVelocities - N) + A_dot*plant.jointVelocities);
constraint_torques = A'*lambda;

theta_two_dot_res = M^(-1) * (T - C*plant.jointVelocities - N - constraint_torques);

unconstrained_equation_of_motion_check_left = M*theta_two_dot_res + C*plant.jointVelocities + N + constraint_torques;
unconstrained_equation_of_motion_check_right = T;
constrained_equation_of_motion_check_left = M*theta_two_dot_res + P*(C*plant.jointVelocities + N);
constrained_equation_of_motion_check_right = P*T;


% there's some problem here, the resulting accelerations are not 0 when constraining both feet

% maybe the matrix I am pseudo-inverting doesn't have full rank?
H = [ ...
      P; ...                                                    % desired joint acceleration
      [eye(6, 6), zeros(6, plant.numberOfJoints - 6)] ...       % torques in the free body dofs
    ];
% yes, that is the right. With both feet constrained, rank(P) = 27, with one feet constrained rank(P) = 32

% I suspect that this is because I have too many degrees of freedom. I can probably constrain this by not permitting
% workless joints. But how do I do that?    

% return

% calculate ground reaction wrenches
lambda_left = [lambda(1:number_of_left_foot_constraints); zeros(number_of_right_foot_constraints, 1)];
lambda_right = [zeros(number_of_left_foot_constraints, 1); lambda(number_of_left_foot_constraints+1 : number_of_left_foot_constraints+number_of_right_foot_constraints);];

constraint_torque_left = A' * lambda_left;
constraint_torque_right = A' * lambda_right;
constraint_torque_left_alt = A_left' * lambda(1:number_of_left_foot_constraints);
constraint_torque_right_alt = A_right' * lambda(number_of_left_foot_constraints+1 : number_of_left_foot_constraints+number_of_right_foot_constraints);

left_ground_reaction_wrench = calculateInstantaneousGroundReactionWrench(plant, constraint_torque_left, eye(4, 4));
right_ground_reaction_wrench = calculateInstantaneousGroundReactionWrench(plant, constraint_torque_right, eye(4, 4));
combined_ground_reaction_wrench = left_ground_reaction_wrench + right_ground_reaction_wrench;
combined_ground_reaction_wrench = [left_ground_reaction_wrench; right_ground_reaction_wrench];

left_ground_reaction_wrench_alt = calculateInstantaneousGroundReactionWrench_alt(plant, constraint_torque_left_alt, eye(4, 4), 12);
right_ground_reaction_wrench_alt = calculateInstantaneousGroundReactionWrench_alt(plant, constraint_torque_right_alt, eye(4, 4), 18);
added_ground_reaction_wrench_alt = left_ground_reaction_wrench_alt + right_ground_reaction_wrench_alt;
% combined_ground_reaction_wrench_alt = right_ground_reaction_wrench_alt;

left_ground_reaction_wrench = pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 12)') * constraint_torque_left;
right_ground_reaction_wrench = pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 18)') * constraint_torque_right;



% inverse dynamics
P = eye(plant.numberOfJoints) - A' * (A * M^(-1) * A')^(-1) * A * M^(-1);
virtual_joints = 1 : 6;
k_c = rank(A);
k_v = length(virtual_joints);
[~, ~, V_c] = svd(A);
C_c = V_c(:, k_c+1:end); % C_c spans the null space of A
B_v = [eye(k_v); zeros(plant.numberOfJoints - k_v, k_v)];
k_w = rank([B_v C_c]);
[~, ~, V_w] = svd([B_v C_c]');
C_w = V_w(:, k_w+1:end);
D = [B_v C_w];


% do not use the "no workless torques"-assumption
% D = B_v;

% set up equation matrices
R_left = pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 12)');
H_1_left = R_left * A' * (A*M^(-1)*A')^(-1)*A*M^(-1);
b_1_left = left_ground_reaction_wrench_alt - R_left * A' * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*plant.jointVelocities - N) + A_dot*plant.jointVelocities);

R_right = pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 18)');
H_1_right = R_right * A' * (A*M^(-1)*A')^(-1)*A*M^(-1);
b_1_right = right_ground_reaction_wrench_alt - R_right * A' * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*plant.jointVelocities - N) + A_dot*plant.jointVelocities);

% R_both = [pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 12)'); pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 18)')];
% H_1_both = R_both * A' * (A*M^(-1)*A')^(-1)*A*M^(-1);
% b_1_both = combined_ground_reaction_wrench - R_both * A' * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*plant.jointVelocities - N) + A_dot*plant.jointVelocities);
H_1_both = [H_1_left; H_1_right];
b_1_both = [b_1_left; b_1_right];

% same equation, just using A_left and A_right instead of full A
% R_left = pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 12)');
% H_1_left = R_left * A_left' * (A_left*M^(-1)*A_left')^(-1)*A_left*M^(-1);
% b_1_left = left_ground_reaction_wrench_alt - R_left * A_left' * (A_left*M^(-1)*A_left')^(-1) * (A_left*M^(-1)*(- C*plant.jointVelocities - N) + A_left_dot*plant.jointVelocities);
% R_right = pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 18)');
% H_1_right = R_right * A_right' * (A_right*M^(-1)*A_right')^(-1)*A_right*M^(-1);
% b_1_right = right_ground_reaction_wrench_alt - R_right * A_right' * (A_right*M^(-1)*A_right')^(-1) * (A_right*M^(-1)*(- C*plant.jointVelocities - N) + A_right_dot*plant.jointVelocities);

% GRF equations with appropriate projections
Q_left = [eye(number_of_left_foot_constraints) zeros(number_of_right_foot_constraints, number_of_left_foot_constraints); zeros(number_of_left_foot_constraints, number_of_right_foot_constraints) zeros(number_of_right_foot_constraints)];
R_left = pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 12)');
H_1_left = R_left * A' * Q_left * (A*M^(-1)*A')^(-1)*A*M^(-1);
b_1_left = left_ground_reaction_wrench_alt - R_left * A' * Q_left * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*plant.jointVelocities - N) + A_dot*plant.jointVelocities);

Q_right = [zeros(number_of_left_foot_constraints) zeros(number_of_right_foot_constraints, number_of_left_foot_constraints); zeros(number_of_left_foot_constraints, number_of_right_foot_constraints) eye(number_of_right_foot_constraints)];
R_right = pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 18)');
H_1_right = R_right * A' * Q_right * (A*M^(-1)*A')^(-1)*A*M^(-1);
b_1_right = right_ground_reaction_wrench_alt - R_right * A' * Q_right * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*plant.jointVelocities - N) + A_dot*plant.jointVelocities);

% GRF equation with only one forceplate for two feet
H_1_both = H_1_left + H_1_right;
b_1_both = b_1_left + b_1_right;


% measured accelerations are generated
H_2 = M^(-1)*P;
b_2 = theta_two_dot_res + M^(-1)*P*(C*plant.jointVelocities + N) + M^(-1)*A'*(A*M^(-1)*A')^(-1)*A_dot*plant.jointVelocities;
% no virtual and no workless torques
H_3 = D';
b_3 = zeros(size(D, 2), 1);

H = [H_2; H_3];
b = [b_2; b_3];

% H = [H_1_left; H_2; H_3];
% b = [b_1_left; b_2; b_3];

H = [H_1_left; H_1_right; H_2; H_3];
b = [b_1_left; b_1_right; b_2; b_3];


F = pinv(H) * b;

TF = [T F];
TmF = T - F;
% works for one contact foot, i.e. F = T, the joint torques calculated by inverse dynamics are the same as I applied

% test the GRF equations for the solution without GRF data
check_left = H_1_left * T;
b_1_left;
check_right = H_1_right * T;



% these equations do not hold. 
% Why not? check that from the beginning

left_ground_reaction_wrench = R_left * constraint_torque_left;
combined_ground_reaction_wrench = [R_left zeros(size(R_left)); zeros(size(R_left)) R_right] * [constraint_torque_left; constraint_torque_right];



left_ground_reaction_wrench_check = R_left * A' * Q_left * lambda; 

check_left_side = R_left * A' * Q_left * (A*M^(-1)*A')^(-1) * A * M^(-1) * T;
check_right_side = left_ground_reaction_wrench_alt - R_left * A' * Q_left * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*plant.jointVelocities - N) + A_dot*plant.jointVelocities);

return






