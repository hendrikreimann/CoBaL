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


function adjoint = rigidToAdjointTransformation(rigidTransformation)
% Adjoint operator for rigid transformations
%
% adjoint = RIGIDTOADJOINTTRANSFORMATION(g) transforms a 4x4 rigid 
%           transformation g into the corresponding 6x6 adjoint 
%           transformation Ad_g
%
% See also ADJOINTTORIGIDTRANSFORMATION, INVERTADJOINTTRANSFORMATION


  % extract componenets
  rotation = rigidTransformation(1:3, 1:3);
  translation = rigidTransformation(1:3, 4);
  translation_wedge_times_rotation = wedgeAxis(translation) * rotation;

  % concatenate to form adjoint
  adjoint = [rotation translation_wedge_times_rotation; zeros(3, 3) rotation];
end

