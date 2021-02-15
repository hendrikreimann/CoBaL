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


function inverse = invertAdjointTransformation(adjoint)
% Calculates the inverse of the adjoint of a rigid transformation
%
% B = INVERTADJOINTTRANSFORMATION(A) inverts the 6x6 adjoint
%           transformation A.
%
% This function is faster than directly inverting the 6x6 matrix A.
%
% See also RIGIDTOADJOINTTRANSFORMATION, ADJOINTTORIGIDTRANSFORMATION

    inverse = [adjoint(1:3, 1:3)', -adjoint(1:3, 1:3)' * adjoint(1:3, 4:6) * adjoint(1:3, 1:3)'; zeros(3, 3) adjoint(1:3, 1:3)'];
end
