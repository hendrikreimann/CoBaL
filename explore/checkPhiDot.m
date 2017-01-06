% test some differentiation

%% create plant
jointPositions = ...
{ ...
  [0; 0; 0], ...
  [0; 0; 0], ...
  [0; 0; 1], ...
  [0; 0; 1] ...
};
jointAxes = ...
{ ...
  [0; 0; 1], ...
  [1; 0; 0], ...
  [0; 0; 1], ...
  [1; 0; 0] ...
};

plant = GeneralKinematicTree ...
( ...
  jointPositions, ...
  jointAxes, ...
  [1 1 1 1], ...
  [1 1 1 1], ...
  {[eye(3) [0; 0; 2]; [0 0 0 1]]}, ...
  {[0; 0; 0.5]; [0; 0; 0.5]; [0; 0; 1.5]; [0; 0; 1.5]}, ...
  {eye(3), eye(3), eye(3), eye(3)}, ...
  {zeros(6, 6), zeros(6, 6), zeros(6, 6), eye(6)} ...
);
plant.updateInternals();
number_of_joints = plant.numberOfJoints;


%% set configuration and update

theta_0 = [pi/2; pi/2; 0; -pi/2];
thetaDot = [.4; .6; .2; .4];
thetaTwoDot = [3; 2; 1; 3];
theta_0 = randn(number_of_joints, 1);
thetaDot = randn(number_of_joints, 1);
% thetaTwoDot = randn(numberOfJoints, 1);
% thetaTwoDot = zeros(numberOfJoints, 1);

% update plant
plant.jointAngles = theta_0;
plant.jointVelocities = thetaDot;
plant.updateInternals;

p_local = [0; 0; 1; 1];
p = (plant.endEffectorTransformations{1}*p_local);
p = p(1:3);
k = [0; 0; 1];
kNorm = norm(k);
w_local = [0; 1; 0; 0];
w = (plant.endEffectorTransformations{1}*w_local);
w = w(1:3);

%% calculate derivatives analytically
J_p = plant.calculateArbitraryPointJacobian(p_local, number_of_joints, 'local');
J_w = plant.calculateArbitraryPointJacobian(w_local, number_of_joints, 'local');

cosPhi = w'*k / norm(k);
phi = acos(cosPhi);


% calculate analytical derivatives
wDot = J_w * thetaDot;
wDot_ana = wDot;
phiDot_ana = (-wDot'*k/kNorm^2) / (kNorm^2 - (w'*k)^2)^(0.5);

d_phi_by_d_w = -1 / (1 - cosPhi^2)^(0.5) * k' * 1/kNorm;
d_w_by_d_theta = J_w;
d_phi_by_d_theta_ana = d_phi_by_d_w * d_w_by_d_theta;
phiDot_ana = d_phi_by_d_theta_ana * thetaDot;

%% calculate numerical derivative of basic variables
delta_t = 1e-8;

theta_new = theta_0 + delta_t*thetaDot;

plant.jointAngles = theta_new;
plant.updateInternals;

% update plant with new values
p_new = eye(3, 4)*(plant.endEffectorTransformations{1}*p_local);
w_new = eye(3, 4)*(plant.endEffectorTransformations{1}*w_local);
J_p_new = plant.calculateArbitraryPointJacobian(p_local, number_of_joints, 'local');
J_w_new = plant.calculateArbitraryPointJacobian(w_local, number_of_joints, 'local');

wDot_new = J_w_new * thetaDot;
cosPhi_new = w_new'*k / norm(k);
phi_new = acos(cosPhi_new);

% calculate numerical derivatives
cosPhiDot_num = (cosPhi_new - cosPhi) / delta_t;
phiDot_num = (phi_new - phi) / delta_t;
wDot_num = (w_new - w) / delta_t;



% disp(['phiDot_ana = ' num2str(phiDot_ana) ', phiDot_num = ' num2str(phiDot_num)]);

%% calculate joint derivatives
delta_theta = 1e-8;
d_phi_by_d_theta_num = zeros(1, number_of_joints);
for i_joint = 1 : number_of_joints
    plant.jointAngles = theta_0;
    plant.jointAngles(i_joint) = plant.jointAngles(i_joint) + delta_theta;
    plant.updateKinematics;
    
    w_new = eye(3, 4) * (plant.endEffectorTransformations{1}*w_local);
    cosPhi_new = w_new'*k / norm(k);
    phi_new = acos(cosPhi_new);
    
    delta_phi = phi_new - phi;
    d_phi_by_d_theta_num(i_joint) = delta_phi / delta_theta;
end

disp(['d_phi_by_d_theta_ana = ' num2str(d_phi_by_d_theta_ana)]);
disp(['d_phi_by_d_theta_num = ' num2str(d_phi_by_d_theta_num)]);
disp(['phiDot_ana = ' num2str(phiDot_ana)]);
disp(['phiDot_num = ' num2str(phiDot_num)]);




























