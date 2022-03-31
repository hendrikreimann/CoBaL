%     This file is part of the CoBaL code base
%     Copyright (C) 2021 Hendrik Reimann <hendrikreimann@gmail.com>
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

function exportLinearModelGroupResults(varargin)
    % load settings
    linear_model_settings = loadSettingsFromFile('linearModel');
    
    % load data
    model_data_filename = ['groupResults' filesep 'linearModelResults.mat'];
    model_data = load(model_data_filename);
    model_list = linear_model_settings.getTable('export_table');
    
    % process
    for i_model = 1 : size(model_list, 1)
        this_model_label = model_list.label{i_model};
        this_model_index = findModelIndex(model_data.model_results, this_model_label);
        this_model_data = model_data.model_results(this_model_index);
        this_model_filename = ['groupResults' filesep model_list.filename{i_model}];
        
        exportThisModelData(this_model_data, this_model_label, this_model_filename);
    end

end

function exportThisModelData(model_data, label, filename)
    if strcmp(model_data.type, 'discrete')
        writetable(model_data.data, [filename '.csv']);
    elseif strcmp(model_data.type, 'continuous')
        warning(['Exporting model "' label '" failed, export not yet implemented for continuous data models'])
    end
end

function index = findModelIndex(model_data, requested_label)
    % loop through models to find the requested one
    index = 0;
    for i_model = 1 : size(model_data, 1)
        if strcmp(model_data(i_model).label, requested_label)
            index = i_model;
        end
    end
    
    if index == 0
        error(['Model "' requested_label '" not available'])
    end
end
