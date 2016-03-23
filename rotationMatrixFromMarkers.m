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
