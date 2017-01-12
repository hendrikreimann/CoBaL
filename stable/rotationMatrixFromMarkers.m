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

function R = estimateRotation(X, Y)
%
% ---------------------------------------------------------------------
%
% estimateRotation.m
%
% estimates the rotation matrix between two configurations of a rigid body from a set of points in both configurations,
% using Soderkvist & Wedin "Determining the movements of the skeleton using well-configured markers", J. Biomech
% 26:1473-1477.
% 
% Hendrik Reimann, 2011
% Institut fuer Neuroinformatik
% Ruhr-University Bochum
% hendrikreimann@gmail.com
%
% ---------------------------------------------------------------------
%
% Input
% X = matrix of n points on the rigid body before the rotation
%     3 x n    
% Y = matrix of n points on the rigid body after the rotation
%     3 x n    
%
% Output
% R = estimated rotation matrix
%     3 x 3
%
% ---------------------------------------------------------------------

% prune the entries containing NaNs

nanEntries = any(isnan(X) + isnan(Y), 1);
X(:, nanEntries) = [];
Y(:, nanEntries) = [];

if size(X, 2) < 3
    R = NaN * ones(3, 3);
else
    % calculate rigid body centers
    X_center = mean(X, 2);
    Y_center = mean(Y, 2);

    % calculate the vectors to each point 
    X_local = X - X_center * ones(1, size(X, 2));
    Y_local = Y - Y_center * ones(1, size(Y, 2));

    % calculate the rotation
    C = Y_local * X_local';
    [P, ~, Q] = svd(C);
    R = (P * diag([1 1 det(P*Q')], 0) * Q');
end
