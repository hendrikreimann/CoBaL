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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% this script loads the settings file and returns a struct
% input: 
% file studySettings.txt in study root folder (one up from subject folder)
%
% output: 
% struct containing the loaded settings


function settings = loadSettingsFile(filename)
    settings = struct();

    % open file
    if ~exist(filename, 'file')
        fileID = fopen(filename, 'w');
        fprintf(fileID,'example settings file, edit this later\n');
        fclose(fileID);
        
        error(['Failed to load "' filename '", a sample file has been created. Please edit it with your study information or copy the correct file.'])
    end

    fileID = fopen(filename, 'r');
    
    % for now, just read line by line and assume each one is a variable. Later I should add something like
    % readSettingsBlock()
    
    % read text to cell
    text_line = fgetl(fileID);
    text_cell = {};
    while ischar(text_line)
        text_cell = [text_cell; text_line]; %#ok<AGROW>
        text_line = fgetl(fileID);
    end
    fclose(fileID);
    
    % extract data and store in settings struct
    while ~isempty(text_cell)
        [text_cell, settings] = parseNextBlock(text_cell, settings);
    end
end


function [text_cell, settings] = parseNextBlock(text_cell, settings)
    % get first line of remaining text
    text_line = text_cell{1};
    
    if isempty(text_line)
        % empty line, remove line
        text_cell = text_cell(2:end);
        return
    end
    
    if ~any(text_line ~= ' ')
        % only spaces, remove line
        text_cell = text_cell(2:end);
        return
    end
    
    if length(text_line) >= 2 && strcmp(text_line(1:2), '//')
        % comment, remove line
        text_cell = text_cell(2:end);
        return
    end

    if (length(text_cell) > 1) && (text_cell{2}(1) == '{')
        % this is the beginning of a block
        line_split = strsplit(text_line, ':');
        variable_name = strrep(line_split{1}, ' ', '_');
        
        % get data
        block_end_line_index = find(strcmp(text_cell, '}'), 1, 'first');
        variable_data_lines = text_cell(3 : block_end_line_index-1);
        variable_value = {};
        for i_line = 1 : length(variable_data_lines)
            this_line_text = variable_data_lines{i_line};
            while ~isempty(this_line_text) && this_line_text(1) == ' '
                this_line_text(1) = [];
            end
            this_line_text = strrep(this_line_text, ', ', ',');
            this_line_cell = strsplit(this_line_text, ',');
            variable_value(i_line, :) = this_line_cell; %#ok<AGROW>
        end
                
        % try to transform this into a double array
        variable_array = zeros(size(variable_value)) * NaN;
        for i_row = 1 : size(variable_value, 1)
            for i_col = 1 : size(variable_value, 2)
                if all(ismember(variable_value{i_row, i_col}, '0123456789+-.eEdD'))
                    variable_array(i_row, i_col) = str2num(variable_value{i_row, i_col}); %#ok<ST2NM>
                end
            end
        end
        if ~any(isnan(variable_array))
            variable_value = variable_array; %#ok<NASGU>
        end
        
        % store in struct
        evalstring = ['settings.' variable_name ' = variable_value;'];
        eval(evalstring);
        
        % remove parsed line
        text_cell = text_cell(block_end_line_index+1:end);
        
        return
    end
    
    % parse first line as a single entry=
    line_split = strsplit(text_line, ':');
    variable_name = strrep(line_split{1}, ' ', '_');
    variable_value_string = line_split{2};
    while ~isempty(variable_value_string) && variable_value_string(1) == ' '
        variable_value_string(1) = [];
    end
    variable_value_string = strrep(variable_value_string, ', ', ',');
    variable_value_cell = strsplit(variable_value_string, ',');
    if length(variable_value_cell) == 1
        variable_value = variable_value_cell{1};

        % try to transform to a single double
        if ~isempty(str2num(variable_value)) %#ok<ST2NM>
            variable_value = str2num(variable_value); %#ok<ST2NM,NASGU>
        end
    else
        variable_value = variable_value_cell;

        % try to transform to a double array
        variable_value_array = zeros(size(variable_value)) * NaN;
        for i_entry = 1 : length(variable_value)
            if ~isempty(str2num(variable_value{i_entry})) %#ok<ST2NM>
                variable_value_array(i_entry) = str2num(variable_value{i_entry}); %#ok<ST2NM>
            end
        end
        if ~any(isnan(variable_value_array))
            variable_value = variable_value_array; %#ok<NASGU>
        end
    end
    
    % add to settings
    evalstring = ['settings.' variable_name ' = variable_value;'];
    eval(evalstring);
    
    % remove parsed line
    text_cell = text_cell(2:end);
end



















