% this script loads the settings file and returns a struct
% input: 
% file studySettings.txt in study root folder (one up from subject folder)
%
% output: 
% struct containing the loaded settings

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


function settings = loadStudySettings()
    settings = struct();

    % open file
    study_settings_file = ['..' filesep 'studySettings.txt'];
    if ~exist(study_settings_file, 'file')
        fileID = fopen(study_settings_file, 'w');
        fprintf(fileID,'example settings file, edit this later\n');
        fclose(fileID);
        
        disp('Failed to load "studySettings.txt", a sample file has been created. Please edit it with your study information or copy the correct file.')
        return
    end

    fileID = fopen(study_settings_file, 'r');
    
    % for now, just read line by line and assume each one is a variable. Later I should add something like
    % readSettingsBlock()
    
    % read text to cell
    text_line = fgetl(fileID);
    text_cell = {};
    while ischar(text_line)
        text_cell = [text_cell, text_line];
        text_line = fgetl(fileID);
    end
    fclose(fileID);
    
    % store in struct
    for i_line = 1 : length(text_cell)
        text_line = text_cell{i_line};
        % parse this line
        line_split = strsplit(text_line, ':');
        variable_name = strrep(line_split{1}, ' ', '_');;
        variable_value_string = strrep(line_split{2}, ' ', '');
        variable_value_cell = strsplit(variable_value_string, ',');
        if length(variable_value_cell) == 1
            variable_value = variable_value_cell{1};
        else
            variable_value = variable_value_cell;
        end
        
        evalstring = ['settings.' variable_name ' = variable_value;'];
        eval(evalstring);
    end 
    
end






















