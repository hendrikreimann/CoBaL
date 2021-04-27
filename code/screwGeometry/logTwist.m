%     This file is part of the ScrewGeometry library
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


function [twist, angle] = logTwist(transformationMatrix, useNegativeTheta)
% Logarithm map SE(3) -> se(3), using Rodriguez' formula
%
% This is the inverse of the exp: se(3) -> SE(3)
%
% [xi, theta] = LOGTWIST(R) transforms a 4x4 rigid transformation matrix
%               into the twist xi and the angle of rotation theta. The 
%               angle of rotation is chosen to be positive.
% [xi, theta] = LOGTWIST(R, useNegativeTheta) chooses the angle of
%               rotation theta to be negative.
%
% See also LOGAXIS, EXPTWIST

    if (size(transformationMatrix, 1) ~= 4) || (size(transformationMatrix, 2) ~= 4)
        error('transformationMatrix must be a 3x3 matrix.')
    end
    
    if nargin < 2
        useNegativeTheta = false;
    end

    % extract rotation and translation
    rotation = transformationMatrix(1:3, 1:3);
    translation = transformationMatrix(1:3, 4);

    % calculate exponential coordinates for R = exp(\wedge \omega \theta)
    [omega, angle] = logAxis(rotation, useNegativeTheta);

    % calculate v
    omega_wedge = wedgeAxis(omega);
    if isZero(norm(translation))
        v = zeros(3, 1);
    else
        v = ((eye(3, 3) - rotation) * omega_wedge + omega * omega' * angle)^(-1) * translation;
    end

    % concatenate omega and v to get twist coordinates
    twist = [v; omega];
end
