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

function [data, time, sampling_rate, labels, directions, success, data_file_name] = loadData(date, subject_id, trial_type, trial_number, data_name, optional)
    if nargin < 6
        optional = 'required';
    end
    
    % load list of available variables
    info_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'availableVariables.mat');
    load(['analysis' filesep info_file_name], 'available_variables');
    
    % check if this is a compound
    if any(data_name==':')
        data_label = data_name;
        this_variable_split = strsplit(data_name, ':');
        data_name = [this_variable_split{1} '_trajectories'];
        data_sublabel = this_variable_split{2};
        is_compound = true;
    else
        is_compound = false;
    end
    
    
    % load requested data
    if ~any(any(strcmp(data_name, available_variables)))
        if ~strcmp(optional, 'optional')
            % data is not optional, throw error
            error(['Required data "' data_name '" not available.'])
        end
        % this was optional, so return empty arrays
        data = [];
        time = [];
        sampling_rate = [];
        labels = [];
        directions = [];
        success = false;
        return
    end
    
    % load data
    row = find(strcmp(available_variables(:, 1), data_name));
    data_time = available_variables{row, 2};
    data_sampling_rate = available_variables{row, 3};
    data_labels = available_variables{row, 4};
    if strcmp(data_labels(1), '_')
        load_data_labels = true;
        data_labels_variable = data_labels(2:end);
    else
        load_data_labels = false;
    end
    data_directions = available_variables{row, 5};
    if strcmp(data_directions(1), '_')
        load_data_directions = true;
        data_directions_variable = data_directions(2:end);
    else
        load_data_directions = false;
    end
    data_folder = available_variables{row, 6};
    data_file_name = available_variables{row, 7};
    
    if load_data_labels && load_data_directions
        loaded_data = load([data_folder filesep data_file_name], data_name, data_time, data_sampling_rate, data_labels_variable, data_directions_variable);
    end
    if load_data_labels && ~load_data_directions
        loaded_data = load([data_folder filesep data_file_name], data_name, data_time, data_sampling_rate, data_labels_variable);
    end
    if ~load_data_labels && load_data_directions
        loaded_data = load([data_folder filesep data_file_name], data_name, data_time, data_sampling_rate, data_directions_variable);
    end
    if ~load_data_labels && ~load_data_directions
        loaded_data = load([data_folder filesep data_file_name], data_name, data_time, data_sampling_rate);
    end

    % hand over to output arguments
    eval(['data = loaded_data.' data_name ';'])
    eval(['time = loaded_data.' data_time ';'])
    eval(['sampling_rate = loaded_data.' data_sampling_rate ';'])
    
    if load_data_labels
        eval(['labels = loaded_data.' data_labels_variable ';'])
    else
        labels = data_labels;
    end
    if load_data_directions
        eval(['directions = loaded_data.' data_directions_variable ';'])
    else
        directions = data_directions;
    end

    % 
    if is_compound
        data = data(:, strcmp(labels, data_sublabel)); %#ok<NODEF>
        directions = directions(:, strcmp(labels, data_sublabel)); %#ok<NODEF>
        labels = data_label;
    end
    
%     if isempty(data_labels)
%         labels = data_name;
%     else
%         eval(['labels = loaded_data.' data_labels ';'])
%     end
    success = true;
end