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

% TODO: this should be merged into the WalkingData object

function [data, time, sampling_rate, labels, success] = loadData(date, subject_id, trial_type, trial_number, data_name, optional)
    if nargin < 6
        optional = 'required';
    end
    
    % load list of available variables
    info_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'availableVariables.mat');
    load(['analysis' filesep info_file_name], 'available_variables');
    
    % load requested data
    if ~any(any(strcmp(data_name, available_variables)))
        if ~strcmp(optional, 'optional')
            % data is not optional, throw error
            error(['Required data "' data_name '" not available.'])
        end
        data = [];
        time = [];
        sampling_rate = [];
        labels = [];
        success = false;
        return
    end
    
    % load data
    row = find(strcmp(available_variables(:, 1), data_name));
    data_time = available_variables{row, 2};
    data_sampling_rate = available_variables{row, 3};
    data_labels = available_variables{row, 4};
    data_folder = available_variables{row, 5};
    data_file_name = available_variables{row, 6};
    loaded_data = load([data_folder filesep data_file_name], data_name, data_time, data_sampling_rate, data_labels);

    % hand over
    eval(['data = loaded_data.' data_name ';'])
    eval(['time = loaded_data.' data_time ';'])
    eval(['sampling_rate = loaded_data.' data_sampling_rate ';'])
    eval(['labels = loaded_data.' data_labels ';'])
    success = true;
end