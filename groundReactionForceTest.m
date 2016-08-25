
% test calculation of the ground reaction forces and torques


% load and initialize
cd /Users/reimajbi/TempleDrive/20160517_HR_balanceBeam/processed

load subjectInfo.mat;
model_file_name = makeFileName(date, subject_id, 'model');
load(model_file_name);

% plant.updateInternals;


% calculate constraints
A_left = plant.bodyJacobians{3};
A_left_dot = plant.bodyJacobianTemporalDerivatives{3};
A_right = [];%plant.bodyJacobians{6};
A_right_dot = [];%plant.bodyJacobianTemporalDerivatives{6};
number_of_left_foot_constraints = size(A_left, 1);
number_of_right_foot_constraints = size(A_right, 1);

A = [A_left; A_right];
A_dot = [A_left_dot; A_right_dot];
% A = A_left;
% A_dot = A_left_dot;


M = plant.inertiaMatrix;
C = plant.coriolisMatrix;
N = plant.gravitationalTorqueMatrix;

% find joint torques for zero accelerations
P = eye(plant.numberOfJoints) - A' * (A * M^(-1) * A')^(-1) * A * M^(-1);
theta_two_dot_des = zeros(plant.numberOfJoints, 1);
T = pinv ...
  ( ...
    [ ...
      P; ...                                                    % desired joint acceleration
      [eye(6, 6), zeros(6, plant.numberOfJoints - 6)] ...       % torques in the free body dofs
    ] ...       
  ) * ...
  [ ...
    P*C*plant.jointVelocities + P*N + M*theta_two_dot_des; ...  % desired joint acceleration
    zeros(6, 1) ...                                             % torques in free dofs
  ];

% apply constraints and check resulting joint accelerations
lambda = (A*M^(-1)*A')^(-1) ...
    * (A*M^(-1)*(T - C*plant.jointVelocities - N) + A_dot*plant.jointVelocities);
constraint_torques = A'*lambda;

theta_two_dot_res = M^(-1) * (T - C*plant.jointVelocities - N - constraint_torques);


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
combined_ground_reaction_wrench_alt = left_ground_reaction_wrench_alt;% + right_ground_reaction_wrench_alt;




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

% set up equation matrices
R_left = pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 12)');
H_1_left = R_left * A' * (A*M^(-1)*A')^(-1)*A*M^(-1);
b_1_left = left_ground_reaction_wrench_alt - R_left * A' * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*plant.jointVelocities - N) + A_dot*plant.jointVelocities);

% R_right = pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 18)');
% H_1_right = R_right * A' * (A*M^(-1)*A')^(-1)*A*M^(-1);
% b_1_right = right_ground_reaction_wrench_alt - R_right * A' * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*plant.jointVelocities - N) + A_dot*plant.jointVelocities);

% R_both = [pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 12)'); pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 18)')];
% H_1_both = R_both * A' * (A*M^(-1)*A')^(-1)*A*M^(-1);
% b_1_both = combined_ground_reaction_wrench - R_both * A' * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*plant.jointVelocities - N) + A_dot*plant.jointVelocities);


H_2 = M^(-1)*P;
b_2 = theta_two_dot_res + M^(-1)*P*(C*plant.jointVelocities + N) + M^(-1)*A'*(A*M^(-1)*A')^(-1)*A_dot*plant.jointVelocities;
H_3 = D';
b_3 = zeros(size(D, 2), 1);

H = [H_1_left; H_2; H_3];
b = [b_1_left; b_2; b_3];

% H = [H_1_both; H_2; H_3];
% b = [b_1_both; b_2; b_3];

H = [H_2; H_3];
b = [b_2; b_3];

F = pinv(H) * b;

% works for one contact foot, i.e. F = T, the joint torques calculated by inverse dynamics are the same as I applied



return






