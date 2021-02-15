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


function skewMatrix = wedgeAxis(axis)
% Wedge-operator R^3 --> so(3)
%
% transforms a vector to the skew-symmetric matrix representing the cross
% product with that vector
%
% omegaHat = WEDGEAXIS(omega) transforms a 3d axis vector omega into a 
%                        3x3 skew-symmetric matrix
%
% See also VEEAXIS, EXPAXIS


    skewMatrix = zeros(3, 3);
    try
        skewMatrix(1, 1) = 0;
        skewMatrix(1, 2) = -axis(3);
        skewMatrix(1, 3) = axis(2);
        skewMatrix(2, 1) = axis(3);
        skewMatrix(2, 2) = 0;
        skewMatrix(2, 3) = -axis(1);
        skewMatrix(3, 1) = -axis(2);
        skewMatrix(3, 2) = axis(1);
        skewMatrix(3, 3) = 0;
    catch exception
        if length(axis) ~= 3
            error('axis must be a vector of length 3')
        else
            rethrow(exception)
        end
    end
end
