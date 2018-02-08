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

% input
% name: name of the variable containing the data
% time: name of the variable containing the corresponding time vector
% labels: name of the variable containing the corresponding labels
% file: identifier of the file in which this information is saved
% folder: folder containing the file

% function addAvailableData(name, time, sampling_rate, labels, folder, file_name)
function addAvailableData_new(name, time, sampling_rate, labels, directions, folder, file_name)
    [date, subject_id, trial_type, trial_number] = getFileParameters(file_name);
    
    % load list of already available variables if existing
    available_variables = {};
    variable_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'availableVariables.mat');
    
    if exist(['analysis' filesep variable_file_name], 'file')
        load(['analysis' filesep variable_file_name]);
    end

    % add new entry
    new_entry = {name, time, sampling_rate, labels, directions, folder, file_name};
    if any(any(strcmp(name, available_variables)))
        available_variables(strcmp(name, available_variables(:, 1)), :) = new_entry;
    else
        available_variables(end+1, :) = new_entry;
    end
    
    % save
    save(['analysis' filesep variable_file_name], 'available_variables');
end