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


function [axis, angle] = logAxis(rotationMatrix, useNegativeTheta)
% Logarithm map SO(3) -> so(3), using Rodriguez' formula
%
% This is the inverse of the exp: so(3) -> SO(3)
%
% [omega, theta] = LOGAXIS(R) transforms a 3x3 rotation matrix into the
%                  axis of rotation omega and the angle of rotation theta. 
%                  The angle of rotation is chosen to be positive.
% [omega, theta] = LOGAXIS(R, useNegativeTheta) chooses the angle of
%                  rotation theta to be negative.
%
% See also EXPAXIS

    if (size(rotationMatrix, 1) ~= 3) || (size(rotationMatrix, 2) ~= 3)
        error('rotationMatrix must be a 3x3 matrix.')
    end

    axis = zeros(3, 1);
    if nargin < 2
        useNegativeTheta = false;
    end
    try
        rotation_matrix_trace = trace(rotationMatrix);
        % calculate rotation angle
        angle = acos((rotation_matrix_trace - 1.0) / 2.0);
        if (rotation_matrix_trace <= -1.0)
            % capture numeric failures for rotation_matrix_trace very close to one
            angle = pi;
        end
        if useNegativeTheta
            angle = 2.0 * pi - angle;
        end    
        % calculate axis of rotation
        sin_theta = sin(angle);
        cos_theta = cos(angle);
        axis(1, 1) = 1 / (2 * sin_theta) * (rotationMatrix(3, 2) - rotationMatrix(2, 3));
        axis(2, 1) = 1 / (2 * sin_theta) * (rotationMatrix(1, 3) - rotationMatrix(3, 1));
        axis(3, 1) = 1 / (2 * sin_theta) * (rotationMatrix(2, 1) - rotationMatrix(1, 2));
        axis = axis * (1 / norm(axis));

        if (isZero(angle - pi))
            % easy way of calculating axis fails for numerical reasons, choose different formula, depending on signs
            if ((rotationMatrix(2, 1) > 0) && (rotationMatrix(3, 1) > 0) && (rotationMatrix(3, 2) > 0))

              axis(1, 1) = + sqrt((rotationMatrix(1, 1) - cos_theta) / (1-cos_theta));
              axis(2, 1) = + sqrt((rotationMatrix(2, 2) - cos_theta) / (1-cos_theta));
              axis(3, 1) = + sqrt((rotationMatrix(3, 3) - cos_theta) / (1-cos_theta));
            end
            if ((rotationMatrix(2, 1) > 0) && (rotationMatrix(3, 1) > 0) && (rotationMatrix(3, 2) < 0))

              axis(1, 1) = + sqrt((rotationMatrix(1, 1) - cos_theta) / (1-cos_theta));
              axis(2, 1) = + sqrt((rotationMatrix(2, 2) - cos_theta) / (1-cos_theta));
              axis(3, 1) = - sqrt((rotationMatrix(3, 3) - cos_theta) / (1-cos_theta));
            end
            if ((rotationMatrix(2, 1) < 0) && (rotationMatrix(3, 1) > 0) && (rotationMatrix(3, 2) < 0))

              axis(1, 1) = + sqrt((rotationMatrix(1, 1) - cos_theta) / (1-cos_theta));
              axis(2, 1) = - sqrt((rotationMatrix(2, 2) - cos_theta) / (1-cos_theta));
              axis(3, 1) = + sqrt((rotationMatrix(3, 3) - cos_theta) / (1-cos_theta));
            end
            if ((rotationMatrix(2, 1) < 0) && (rotationMatrix(3, 1) < 0) && (rotationMatrix(3, 2) > 0))

              axis(1, 1) = + sqrt((rotationMatrix(1, 1) - cos_theta) / (1-cos_theta));
              axis(2, 1) = - sqrt((rotationMatrix(2, 2) - cos_theta) / (1-cos_theta));
              axis(3, 1) = - sqrt((rotationMatrix(3, 3) - cos_theta) / (1-cos_theta));
            end
        end
        angle = normalizeAngle(angle);
    catch exception
        rethrow(exception);
    end


