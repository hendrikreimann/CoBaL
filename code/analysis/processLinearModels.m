
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

% This function uses the previously calculated stretch variables to process analysis variables. For all stretch variables,
% the response is calculated, i.e. the difference from the control mean.

function processLinearModels(varargin)
    study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');

    % load settings and existing results
    band_labels = study_settings.get('band_labels');
    results_file_name = ['results' filesep makeFileName(collection_date, subject_id, 'results')];
    data = load(results_file_name);
    
    variable_table = study_settings.getTable('linear_model_variables', 1);
    step_time_index_in_saved_data = find(strcmp(data.stretch_names_session, 'step_time'), 1, 'first');
    this_step_time_data = data.stretch_data_session{step_time_index_in_saved_data};
    
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    number_of_time_steps_per_continuous_stretch = (number_of_time_steps_normalized-1) * data.bands_per_stretch + 1;
    number_of_stretches = numel(data.time_list_session);
    
    for i_variable = 1 : size(variable_table, 1)
        % get data
        variable_A_name = variable_table.variable_A_name{i_variable};
        variable_A_type = variable_table.variable_A_type{i_variable};
        variable_A_time = variable_table.variable_A_time{i_variable};
        variable_A_time_type = variable_table.variable_A_time_type{i_variable};
        variable_B_name = variable_table.variable_B_name{i_variable};
        variable_B_type = variable_table.variable_B_type{i_variable};
        variable_B_time = variable_table.variable_B_time{i_variable};
        variable_B_time_type = variable_table.variable_B_time_type{i_variable};
        
        % pick data depending on source specification
        variable_A_data_source = data.([variable_A_type '_data_session']);
        variable_A_names_source = data.([variable_A_type '_names_session']);
        variable_A_directions_source = data.([variable_A_type '_directions_session']);
        variable_B_data_source = data.([variable_B_type '_data_session']);
        variable_B_names_source = data.([variable_B_type '_names_session']);
        variable_B_directions_source = data.([variable_B_type '_directions_session']);
        
        % extract
        variable_A_data_all = variable_A_data_source{strcmp(variable_A_names_source, variable_A_name)};
        variable_B_data_all = variable_B_data_source{strcmp(variable_B_names_source, variable_B_name)};
        variable_A_directions = variable_A_directions_source(strcmp(variable_A_names_source, variable_A_name), :);
        variable_B_directions = variable_B_directions_source(strcmp(variable_B_names_source, variable_B_name), :);
        
        % figure out time
        if strcmp(variable_A_time_type, 'percentage')
            % haven't checked whether this case works yet
            variable_A_indices = ones(size(this_step_time_data)) * str2double(start_info);
        elseif strcmp(variable_A_time_type, 'stretch')
            start_data_time_within_band = data.stretch_data_session{strcmp(data.stretch_names_session, variable_A_time)};
            start_data_ratio = start_data_time_within_band ./ this_step_time_data;
            variable_A_indices_within_band = round(start_data_ratio * 100);
            variable_A_indices = variable_A_indices_within_band + (0:data.bands_per_stretch-1)' * (number_of_time_steps_normalized-1);
        elseif strcmp(variable_A_time_type, 'discrete')
            
        else
            error(['Unrecognized time type "' variable_A_time_type '" specified for variable "' variable_A_name]);
        end
        if strcmp(variable_B_time_type, 'percentage')
            % haven't checked whether this case works yet
            variable_B_indices = ones(size(this_step_time_data)) * str2double(start_info);
        elseif strcmp(variable_B_time_type, 'stretch')
            start_data_time_within_band = data.stretch_data_session{strcmp(data.stretch_names_session, variable_B_time)};
            start_data_ratio = start_data_time_within_band ./ this_step_time_data;
            variable_B_indices = round(start_data_ratio * 100);
        elseif strcmp(variable_B_time_type, 'discrete')
            variable_B_indices = ones(data.bands_per_stretch, number_of_time_steps_per_continuous_stretch) + (0:data.bands_per_stretch-1)';
        else
            error(['Unrecognized time type "' variable_B_time_type '" specified for variable "' variable_B_name]);
        end
        
        % extract data
        variable_A_data = zeros(data.bands_per_stretch, number_of_stretches);
        variable_B_data = zeros(data.bands_per_stretch, number_of_stretches);
        for i_band = 1 : data.bands_per_stretch
            for i_stretch = 1 : number_of_stretches
                variable_A_index = variable_A_indices(i_band, i_stretch);
                variable_A_data(i_band, i_stretch) = variable_A_data_all(variable_A_index, i_stretch);
                variable_B_index = variable_B_indices(i_band, i_stretch);
                variable_B_data(i_band, i_stretch) = variable_B_data_all(variable_B_index, i_stretch);
            end
        end
        
        % linear fits
        for i_band = 1 : data.bands_per_stretch
            predictor = variable_A_data(i_band, :)';
            outcome = variable_B_data(i_band, :)';
            [fit_object, fit_stats] = fit(predictor, outcome, 'poly1'); %#ok<ASGLU>
            
            [R,P] = corrcoef(predictor, outcome);
            
            prediction = fit_object(predictor);
            error = outcome - prediction;
            R_square = 1 - sum((error).^2)/sum((outcome - mean(outcome)).^2);
            
            figure; 
            plot(fit_object, variable_A_data(i_band, :)', variable_B_data(i_band, :)')
            title([band_labels{i_band} ' - $R^2$ = ' num2str(R_square) ', slope = ' num2str(fit_object.p1)]);
            
            xlabel(strrep(variable_A_name, '_', ' '));
            ylabel(strrep(variable_B_name, '_', ' '));
        end
    end
    
      
end









