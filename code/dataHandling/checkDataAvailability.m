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

function available = checkDataAvailability(date, subject_id, trial_type, trial_number, data_name)
    available = false;
    
    % load list of available variables
    info_file_name = ['analysis' filesep makeFileName(date, subject_id, trial_type, trial_number, 'availableVariables.mat')];
    if ~exist(info_file_name, 'file')
        available = false;
        return
    end
    load(info_file_name, 'available_variables');
    
    % check availability
    if size(available_variables, 2) > 0 && any(strcmp(data_name, available_variables(:, 1)))
        available = true;
    end
    
end