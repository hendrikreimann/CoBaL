%     This file is part of the CoBaL code base
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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.% compare the kinematic tree against the kinematic chain

function colorStruct = createColorStruct
    colorStruct = struct;
    
    colorStruct.left_heel_pos = [0.2081    0.1663    0.5292];
    colorStruct.left_heel_vel = [0.0582    0.4677    0.8589];
    colorStruct.left_heel_acc = [0.0282    0.6663    0.7574];
    colorStruct.left_toes_pos = [0.4783    0.7489    0.4877];
    colorStruct.left_toes_vel = [0.9264    0.7256    0.2996];
    colorStruct.left_toes_acc = [0.9763    0.9831    0.0538];
    colorStruct.right_heel_pos = [0.2081    0.1663    0.5292];
    colorStruct.right_heel_vel = [0.0582    0.4677    0.8589];
    colorStruct.right_heel_acc = [0.0282    0.6663    0.7574];
    colorStruct.right_toes_pos = [0.4783    0.7489    0.4877];
    colorStruct.right_toes_vel = [0.9264    0.7256    0.2996];
    colorStruct.right_toes_acc = [0.9763    0.9831    0.0538];

    colorStruct.left_touchdown = [0 0.4470 0.7410];
    colorStruct.left_pushoff = [0.8500 0.3250 0.0980];
    colorStruct.right_touchdown = [0.9290 0.6940 0.1250];
    colorStruct.right_pushoff = [0.4940 0.1840 0.5560];
    
end