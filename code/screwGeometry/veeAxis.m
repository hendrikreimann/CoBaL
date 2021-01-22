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


function axis = veeAxis(matrix, ignoreSkewRequirement)
% Vee-operator so(3) --> R^3
%
% transforms a skew-symmetric matrix representing the cross
% product with a vector to that vector
%
% omega = VEEAXIS(omegaHat) transforms a 3x3 skew-symmetric matrix into a 
%                           3d vector omega
%
% See also WEDGEAXIS
    axis = zeros(3, 1);
    if isZero(matrix + matrix') || ignoreSkewRequirement
        try
            axis(1, 1) = matrix(3, 2);
            axis(2, 1) = matrix(1, 3);
            axis(3, 1) = matrix(2, 1);
        catch exception
            if (size(matrix, 1) ~= 3) || (size(matrix, 2) ~= 3)
                error('matrix must be 3x3')
            else
                rethrow(exception);
            end
        end
    else        
        error('matrix must be skew-symmetric')
    end
end
