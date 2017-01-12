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


% rename all the qtm files in the current directory to the naming scheme used by CoBaL

% hard-code information here
date = '20160314';
subject_id = 'SS';
trial_type = 'walking';

% prepare folder
if ~exist('misnamed', 'dir')
    mkdir('misnamed')
end

% find files
clear file_name_list;
data_dir = dir('*.qtm');
[file_name_list{1:length(data_dir)}] = deal(data_dir.name);

number_of_header_lines = 26;
trial_number_line = 25;
for i_file = 1 : length(file_name_list)
    data_file_name = file_name_list{i_file};
    
    % extract trial number
    file_name_split_dot = strsplit(data_file_name, '.');
%     file_name_split = strsplit(file_name_split_dot{1}, '_');
%     trial_number = str2num(file_name_split{end});
    
    trial_number = str2num(file_name_split_dot{1}(end-3:end));

    % save file
    new_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'qualisysData.qtm');
    while exist(['renamed' filesep new_file_name], 'file')
        disp([new_file_name ' already exists, adding ".new"']);
        new_file_name = [new_file_name '.new'];
    end
    movefile(data_file_name, ['misnamed' filesep data_file_name])
    copyfile(['misnamed' filesep data_file_name], new_file_name)
end


