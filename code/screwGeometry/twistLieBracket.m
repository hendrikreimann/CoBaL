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


function xi_lie = twistLieBracket(xi_1, xi_2)
% calculates the Lie bracket of two twists
%
% xi_lie = TWISTLIEBRACKET(xi_1, xi_2) calculates the Lie bracket [xi_1, xi_2] 
%          of two twists. The result is given in 6x1 twist coordinates, the
%          inputs can be given in 6x1 twist coordinates or 4x4 twist matrices.
%
% See also WEDGETWIST, VEETWIST

    if (size(xi_1, 1) == 6) || (size(xi_1, 2) == 1)
        xi_1 = wedgeTwist(xi_1);
    end
    if (size(xi_1, 1) ~= 4) || (size(xi_1, 2) ~= 4)
        error('xi_1 must be a 6x1 twist coordinate vector or a 4x4 twist.')
    end
    if (size(xi_2, 1) == 6) || (size(xi_2, 2) == 1)
        xi_2 = wedgeTwist(xi_2);
    end
    if (size(xi_2, 1) ~= 4) || (size(xi_2, 2) ~= 4)
        error('xi_2 must be a 6x1 twist coordinate vector or a 4x4 twist.')
    end

    xi_lie = veeTwist(xi_1 * xi_2 - xi_2 * xi_1);

end