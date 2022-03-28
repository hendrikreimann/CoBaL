%     This file is part of the CoBaL code base
%     Copyright (C) 2022 Hendrik Reimann <hendrikreimann@gmail.com>
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

function collectLinearModelGroupResults(subjects, varargin)
    % parse input parameters
%     parser = inputParser;
%     parser.KeepUnmatched = true;
%     addParameter(parser, 'subjects', {})
%     parse(parser, varargin{:})
%     arguments.subjects = parser.Results.subjects;

    % load settings
    study_settings = loadSettingsFromFile('study');
    linear_model_settings = loadSettingsFromFile('linearModel');
%     condition_to_compare = linear_model_settings.get('condition_to_compare');
    
    % define subjects
    if nargin < 1
        subjects = {'S1', 'S2'};
    end
    
    number_of_subjects = length(subjects);

    % load data
    model_data = cell(number_of_subjects, 1);
    for i_subject = 1 : number_of_subjects
        this_data_folder_path = subjects{i_subject};
        subject_settings = loadSettingsFromFile('subject', this_data_folder_path);
        collection_date = subject_settings.get('collection_date');
        subject_id = subject_settings.get('subject_id');
        
        model_file_name = [this_data_folder_path filesep 'results' filesep collection_date '_' subject_id '_linearModels.mat'];
        model_data{i_subject} = load(model_file_name);
    end
    model_list = linear_model_settings.getTable('models');
    
    % process
    number_of_models = size(model_list, 1);
%     model_results_header = {'label', 'type', 'results'};
%     model_results = cell(number_of_models, 3);
    model_results = [];
for i_model = 1 : number_of_models
        this_model = model_list(i_model, :);
        this_model_results = struct;
        this_model_results.type = this_model.type{1};
        this_model_results.label = this_model.label{1};
        this_model_results.headers = [];
        
        if strcmp(this_model.type{1}, 'discrete')
            this_model_results.data = table;
        elseif strcmp(this_model.type{1}, 'continuous')
            this_model_results.data = {};
        end
        
        % find the requested model in the data
        for i_subject = 1 : number_of_subjects
            % extract model data for this subject
            this_subject_model_data = model_data{i_subject};
            this_model_index = findModelIndex(this_subject_model_data, this_model);
            this_subject_this_model_data  = this_subject_model_data.linear_model_results{this_model_index, strcmp(this_subject_model_data.linear_model_results_header, 'results')};

            % build data table
            if strcmp(this_model.type{1}, 'discrete')
                this_subject_data_table = buildDataTableForThisSubject(this_subject_this_model_data);
            elseif strcmp(this_model.type{1}, 'continuous')
                [this_subject_data_table, this_subject_data_table_headers] = buildDataCellForThisSubject(this_subject_this_model_data);
                this_model_results.headers = this_subject_data_table_headers;
            end
            
            this_model_results.data = [this_model_results.data; this_subject_data_table];
        end
        
        % store
        model_results = [model_results; this_model_results]; %#ok<AGROW>
    end
    
    % save
    save('linearModelResults.mat', 'model_results');
    
end

function index = findModelIndex(data, model)
    this_model_label = model.label{1};
    
    label_column = strcmp(data.linear_model_results_header, 'label');
    index = strcmp(data.linear_model_results(:, label_column), this_model_label);

    if all(index==0)
        error(['Model not available'])
    end
    
% below is old code from before using a label for each model
%     % extract info
%     requested_outcome_variable = model.outcome_variable_name{1};
%     requested_predictor_variable_list_name = model.predictor_variable_list{1};
%     
%     % get list of predictor variables
%     requested_predictor_variable_list = settings.get(requested_predictor_variable_list_name, 1);
%     if isempty(requested_predictor_variable_list)
%         % failed to load list from settings, assume this is already a variable name
%         requested_predictor_variable_list = {requested_predictor_variable_list_name};
%     end
% 
%     predictor_column = strcmp(data.linear_model_results_header, 'predictor_variables');
%     outcome_column = strcmp(data.linear_model_results_header, 'outcome_variable');
%     
%     % loop through models to find the requested one
%     index = 0;
%     for i_model = 1 : size(data.linear_model_results, 1)
%         this_model_outcome_variable = data.linear_model_results{i_model, outcome_column};
%         this_model_predictor_variables = data.linear_model_results{i_model, predictor_column};
%         
%         if isequal(this_model_outcome_variable, requested_outcome_variable) ...
%             && isequal(this_model_predictor_variables, requested_predictor_variable_list)
%             index = i_model;
%         end
%     end
%     
%     if index == 0
%         error(['Model not available'])
%     end
end

function data_table = buildDataTableForThisSubject(model_data)
    % transform conditions from cell array to table
    condition_table = model_data.row_info;
    condition_table_headers = model_data.row_info_headers;
    data_table = cell2table(condition_table, 'VariableNames', condition_table_headers);
    
    % add R^2 column
    r_square_table = model_data.R_square;
    data_table = addvars(data_table, r_square_table, 'NewVariableNames', 'R_square');
    
    % add slope columns
    slope_data = model_data.slope;
    outcome_name = model_data.names.outcome;
    predictor_names = model_data.names.predictors;
    number_of_data_points = length(slope_data);
    for i_predictor = 1 : length(predictor_names)
        % extract data
        this_predictor_slope_data = zeros(number_of_data_points, 1);
        for i_data_point = 1 : number_of_data_points
            this_predictor_slope_data(i_data_point) = slope_data{i_data_point}(:, i_predictor);
        end
        
        this_predictor_header = ['d_' outcome_name '_by_d_' predictor_names{i_predictor}];
        
        data_table = addvars(data_table, this_predictor_slope_data, 'NewVariableNames', this_predictor_header);
    end
end

function [data_table, data_table_headers] = buildDataCellForThisSubject(model_data)
    % take condition table as starting point for data table
    data_table_headers = model_data.row_info_headers;
    data_table = model_data.row_info;
    
    % add R^2 column
    r_square_table = model_data.R_square;
    data_table_headers = [data_table_headers, 'R_square'];
    data_table = [data_table r_square_table];
    
    % add slope columns
    slope_data = model_data.slope;
    outcome_name = model_data.names.outcome;
    predictor_names = model_data.names.predictors;
    number_of_data_points = length(slope_data);
    for i_predictor = 1 : length(predictor_names)
        % extract data
        this_predictor_slope_data = cell(number_of_data_points, 1);
        for i_data_point = 1 : number_of_data_points
            this_predictor_slope_data{i_data_point} = slope_data{i_data_point}(:, i_predictor);
        end
        data_table = [data_table this_predictor_slope_data]; %#ok<AGROW>
        
        this_predictor_header = ['d_' outcome_name '_by_d_' predictor_names{i_predictor}];
        data_table_headers = [data_table_headers, this_predictor_header]; %#ok<AGROW>
    end
end

