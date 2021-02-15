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


function rotationMatrix = expAxis(axis, angle)
% Exponential map so(3) -> SO(3), using Rodriguez' formula
%
% R = EXPAXIS(omega, theta) transforms a normed axis omega and angle theta
%     into the 3x3 matrix representing rotation by theta around omega.
%
% See also LOGAXIS, EXPTWIST

    if (size(axis, 1) ~= 3) || (size(axis, 2) ~= 1)
        error('omega must be a 3x1 vector.')
    end

    omega_wedge = wedgeAxis(axis);
    rotationMatrix = eye(3, 3) + (omega_wedge * sin(angle)) + (omega_wedge * omega_wedge * (1 - cos(angle)));
end