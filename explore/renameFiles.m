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

% this script renames files according to the CoBaL file naming scheme
% it will go through all files in the current folder and replace part of
% the name with string specified in "new_code"
target_folder = 'analysis';

% 1 = date, 2 = subjectID, 3 = trial type, 4 = trial number, 5 = data type
spot_to_replace = 5;
code_to_replace = 'Pevents';
new_code = 'events';

% find files
clear file_name_list;
data_dir = dir([pwd filesep target_folder]);
[file_name_list{1:length(data_dir)}] = deal(data_dir.name);

for i_file = 1 : length(file_name_list)
    data_file_name = file_name_list{i_file};
    file_name_split_dot = strsplit(data_file_name, '.');
    data_file_name_body = file_name_split_dot{1};
    if length(file_name_split_dot) > 1
        data_file_name_type = file_name_split_dot{2};
    end
    
    file_name_split = strsplit(data_file_name_body, '_');
    if length(file_name_split) >= spot_to_replace && strcmp(file_name_split{spot_to_replace}, code_to_replace)
        file_name_split{spot_to_replace} = new_code;
        new_file_name_body = file_name_split{1};
        for i_step = 2 : length(file_name_split)
            new_file_name_body = [new_file_name_body '_' file_name_split{i_step}];
        end
        new_file_name = [new_file_name_body '.' data_file_name_type];
        
        % rename the file
        old_file_name_with_path = [pwd filesep target_folder filesep data_file_name];
        new_file_name_with_path = [pwd filesep target_folder filesep new_file_name];
        movefile(old_file_name_with_path, new_file_name_with_path);
        disp(['Moved file ' old_file_name_with_path ' --> ' new_file_name_with_path])
    end
    
    
end


