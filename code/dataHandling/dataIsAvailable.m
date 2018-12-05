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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.% input

% input
% name: name of the variable containing the data

function is_available = dataIsAvailable(date, subject_id, trial_type, trial_number, data_name)
    % load list of available variables
    info_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'availableVariables.mat');
    load(['analysis' filesep info_file_name], 'available_variables');
    
    available_variable_names = available_variables(:, 1);
    
    % check availability
    is_available = false;
    if any(strcmp(data_name, available_variable_names))
        is_available = true;
    end
end
    
    