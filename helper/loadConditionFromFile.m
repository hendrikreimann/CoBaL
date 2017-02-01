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

% This loads the condition for a single type and trial from a .csv file.

function condition = loadConditionFromFile(filename, condition_label, trial_number)
    % read file
    fileID = fopen(filename, 'r');
    header_line = fgetl(fileID);
    text_cell = {};
    text_line = fgetl(fileID);
    while ischar(text_line)
        text_cell = [text_cell; text_line]; %#ok<AGROW>
        text_line = fgetl(fileID);
    end
    fclose(fileID);
    
    % transform to arrays
    header = strsplit(strrep(header_line, ' ', ''), ',');
    condition_header = header(2 : end);
    number_of_trials = size(text_cell, 1);
    number_of_conditions = size(header, 2) - 1;
    trial_list = zeros(number_of_trials, 1) * NaN;
    condition_cell = cell(number_of_trials, number_of_conditions);
    for i_trial = 1 : number_of_trials
        text_line = text_cell{i_trial};
        line_split = strsplit(text_line, ',');
        trial_list(i_trial) = str2num(line_split{1});
        condition_cell(i_trial, :) = line_split(2:end);
    end
    
    % get condition
    trial_row = find(trial_list == trial_number, 1, 'first');
    if isempty(trial_row)
        condition = [];
    else
        condition = condition_cell{trial_row, strcmp(condition_header, condition_label)};
    end
end