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


function twist = generateTwistCoordinates(supportPoint, axis, type)
% Generate a twist representation from support point and axis
%
% xi = GENERATETWISTCOORDINATES(p, omega, type) generates a 6x1 coordinate
%      representation of a twist corresponding to rotation around the axis
%      with direction omega through the support point p for type one, or
%      corresponding to translation in the direction given by axis
%
% See also EXPTWIST

    if nargin < 3
        type = 1;
    end
    if type == 1
        if (size(supportPoint, 1) ~= 3) || (size(supportPoint, 2) ~= 1)
            error('supportPoint must be a 3x1 vector.')
        end
        if (size(axis, 1) ~= 3) || (size(axis, 2) ~= 1)
            error('axis must be a 3x1 vector.')
        end
        twist = [-cross(axis, supportPoint); axis * 1 / norm(axis)];
    elseif type == 2
        if (size(axis, 1) ~= 3) || (size(axis, 2) ~= 1)
            error('axis must be a 3x1 vector.')
        end
        twist = [axis; zeros(3, 1)];
    else
        error('incorrect joint type specified - please use "1" for revolute joints, "2" for prismatic joints');
        
    end
end