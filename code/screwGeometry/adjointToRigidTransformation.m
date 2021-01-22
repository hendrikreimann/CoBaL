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


function rigidTransformation = adjointToRigidTransformation(adjoint)
% Inverse of the adjoint operator for rigid transformations
%
% g = ADJOINTTORIGIDTRANSFORMATION(g) transforms a 6x6 adjoint
%           transformation Ad_g into the corresponding 4x4 rigid  
%           transformation g
%
% See also RIGIDTOADJOINTTRANSFORMATION, INVERTADJOINTTRANSFORMATION
    rotation = adjoint(1:3, 1:3);
    position = veeAxis(adjoint(1:3, 4:6) * rotation^(-1));

    rigidTransformation = [rotation position; 0 0 0 1];
end