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



shoulder = [0; 0; 0];
elbow_reference = [1; 0; 0];
wrist_reference = [2; 0; 0];
marker_reference = [1; 0; 1];

joint_positions = {shoulder, shoulder, shoulder, elbow_reference};
joint_axes = {[1; 0; 0], [0; 1; 0], [0; 0; 1], [0; 0; 1]};
joint_types = [1 1 1 1];
end_effectors = {marker_reference, wrist_reference};

link_positions = {[0; 0; 0]; [0; 0; 0]; [0; 0; 0.5]; [0; 0; 1.5];};
branch_matrix = [1 1 1 0; 1 1 1 1]; % first eef = marker, second eef = wrist
link_orientations = {eye(3), eye(3), eye(3), eye(3)};

degrees_of_freedom = length(joint_axes);
link_masses = ones(degrees_of_freedom, 1);
link_moments_of_inertia = [[1 5 5] ; repmat([1 1 .5], degrees_of_freedom-1, 1)] * 0.01;

link_orientations = {eye(3), eye(3), eye(3), eye(3)};
generalized_link_inertia_matrices = {eye(6), eye(6), eye(6), eye(6)};



test_arm = GeneralKinematicTree ...
( ...
  joint_positions, ...
  joint_axes, ...
  joint_types, ...
  branch_matrix, ...
  end_effectors, ...
  link_positions, ...
  link_orientations, ...
  generalized_link_inertia_matrices ...
);


% elbow angle should be between 0 and pi, then this is a right arm
test_arm.jointAngles = rand(4, 1);
test_arm.jointAngles = [0; 0; 0; -pi/2];

test_arm.updateInternals();

marker_current = test_arm.endEffectorTransformations{1}(1:3, 4);
wrist_current = test_arm.endEffectorTransformations{2}(1:3, 4);

if test_arm.jointAngles(4) >= 0
    direction = 'positive';
else
    direction = 'negative';
end
elbow_current_actual = test_arm.jointTransformations{4}(1:3, 4)
elbow_current_found = findHingeJointCenter(shoulder, wrist_current, marker_current, norm(elbow_reference - marker_reference), norm(elbow_reference - shoulder), norm(elbow_reference - wrist_reference), direction)

return