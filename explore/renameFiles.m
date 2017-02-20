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


new_code = 'ZRUI';

% find files
clear file_name_list;
data_dir = dir;
[file_name_list{1:length(data_dir)}] = deal(data_dir.name);

for i_file = 4 : length(file_name_list)
    data_file_name = file_name_list{i_file};
    
    file_name_split = strsplit(data_file_name, '_');
    file_name_split{3} = new_code;
    new_file_name = [file_name_split{1} '_' file_name_split{2} '_' new_code];
    for i_step = 4 : length(file_name_split)
        new_file_name = [new_file_name '_' file_name_split{i_step}];
    end
    
    movefile(data_file_name, new_file_name)
end


