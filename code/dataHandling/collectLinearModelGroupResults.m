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

function collectLinearModelGroupResults(varargin)
    % parse input parameters
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'subjects', {})
    parse(parser, varargin{:})
    arguments.subjects = parser.Results.subjects;

    % load settings
    linear_model_settings = loadSettingsFromFile('linearModel');
    
    % define subjects
    [data_folder_list, subjects] = determineDataStructure(arguments.subjects);
    number_of_subjects = length(subjects);

    % load data
    model_data = cell(number_of_subjects, 1);
    for i_subject = 1 : number_of_subjects
        this_data_folder_path = data_folder_list{i_subject};
        subject_settings = loadSettingsFromFile('subject', this_data_folder_path);
        collection_date = subject_settings.get('collection_date');
        subject_id = subject_settings.get('subject_id');
        
        model_file_name = [this_data_folder_path filesep 'results' filesep collection_date '_' subject_id '_linearModels.mat'];
        model_data{i_subject} = load(model_file_name);
    end
    model_list_first_order = linear_model_settings.getTable('models');
    model_list_second_order = linear_model_settings.getTable('models_second_order');
    model_list = [model_list_first_order; model_list_second_order];
    
    % process
    number_of_models = size(model_list, 1);
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
                this_subject_data_table = buildDataTableForThisSubject(this_subject_this_model_data, this_model_results.label);
            elseif strcmp(this_model.type{1}, 'continuous')
                [this_subject_data_table, this_subject_data_table_headers] = buildDataCellForThisSubject(this_subject_this_model_data, this_model_results.label);
                this_model_results.headers = this_subject_data_table_headers;
            end
            
            this_model_results.data = [this_model_results.data; this_subject_data_table];
        end
        
        % store
        model_results = [model_results; this_model_results]; %#ok<AGROW>
    end
    
    % save
    if ~directoryExists('groupResults')
        mkdir('groupResults')
    end
    filename = ['groupResults' filesep 'linearModelResults.mat'];
    save(filename, 'model_results');
    
end

function index = findModelIndex(data, model)
    this_model_label = model.label{1};
    
    label_column = strcmp(data.linear_model_results_header, 'label');
    index = strcmp(data.linear_model_results(:, label_column), this_model_label);

    if all(index==0)
        error('Model not available')
    end
end

function data_table = buildDataTableForThisSubject(model_data, model_label)
    % transform conditions from cell array to table
    condition_table = model_data.row_info;
    condition_table_headers = model_data.row_info_headers;
    data_table = cell2table(condition_table, 'VariableNames', condition_table_headers);
    
    % add variance measures
    r_square_table = model_data.R_square;
    ss_residual_table = model_data.SS_residual;
    ss_total_table = model_data.SS_total;
    ss_difference_table = model_data.SS_difference;
    data_table = addvars(data_table, r_square_table, 'NewVariableNames', 'R_square');
    data_table = addvars(data_table, ss_residual_table, 'NewVariableNames', 'SS_residual');
    data_table = addvars(data_table, ss_total_table, 'NewVariableNames', 'SS_total');
    data_table = addvars(data_table, ss_difference_table, 'NewVariableNames', 'SS_difference');
    
    % add slope columns
    slope_data = model_data.slope;
    predictor_names = model_data.names.predictors;
    number_of_data_points = length(slope_data);
    for i_predictor = 1 : length(predictor_names)
        % extract data
        this_predictor_slope_data = zeros(number_of_data_points, 1);
        for i_data_point = 1 : number_of_data_points
            this_predictor_slope_data(i_data_point) = slope_data{i_data_point}(:, i_predictor);
        end
        this_predictor_header = [model_label '_slope_' num2str(i_predictor)];
        
        data_table = addvars(data_table, this_predictor_slope_data, 'NewVariableNames', this_predictor_header);
    end
    
    % add covariate columns
    covariate_data = model_data.covariate_means;
    covariate_names = model_data.names.covariates;
    for i_covariate = 1 : length(covariate_names)
        % extract data
        this_covariate_slope_data = zeros(number_of_data_points, 1);
        for i_data_point = 1 : number_of_data_points
            this_covariate_slope_data(i_data_point) = covariate_data{i_data_point}(:, i_covariate);
        end
        this_covariate_header = covariate_names{i_covariate};
        
        data_table = addvars(data_table, this_covariate_slope_data, 'NewVariableNames', this_covariate_header);
    end
end

function [data_table, data_table_headers] = buildDataCellForThisSubject(model_data, model_label)
    % take condition table as starting point for data table
    data_table_headers = model_data.row_info_headers;
    data_table = model_data.row_info;
    
    % add R^2 column
    r_square_table = model_data.R_square;
    ss_residual_table = model_data.SS_residual;
    ss_total_table = model_data.SS_total;
    ss_difference_table = model_data.SS_difference;
    data_table_headers = [data_table_headers, 'R_square', 'SS_residual', 'SS_total', 'SS_difference'];
    data_table = [data_table r_square_table ss_residual_table ss_total_table ss_difference_table];
    
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
        
%         this_predictor_header = ['d_' outcome_name '_by_d_' predictor_names{i_predictor}];
        this_predictor_header = [model_label '_slope_' num2str(i_predictor)];
        data_table_headers = [data_table_headers, this_predictor_header]; %#ok<AGROW>
    end
end

