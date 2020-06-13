%     This file is part of the CoBaL code base
%     Copyright (C) 2020 Hendrik Reimann <hendrikreimann@gmail.com>
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

classdef StretchDataCustodian < handle
    properties
        data;
        
        root;
        source;
        session_folder_list;
        number_of_source_sessions;
        study_settings;
    end
    methods
        function this = StretchDataCustodian(root, source, subjects)
            this.root = root;
            this.source = source;
            this.loadData(subjects);
        end
        function loadData(this, subjects)
            % load all data under the root directory
            this.study_settings = loadSettingsFromFile('study');
            
            % declare variables
            file_label = ['results' this.source];

            % load data
            this.session_folder_list = determineDataStructure(subjects);
            this.number_of_source_sessions = length(this.session_folder_list);
            this.data = cell(this.number_of_source_sessions, 1);

            for i_folder = 1 : this.number_of_source_sessions
                % get information
                this_data_folder_path = this.session_folder_list{i_folder};
                subject_settings = loadSettingsFromFile('subject', this_data_folder_path);
                collection_date = subject_settings.get('collection_date');
                subject_id = subject_settings.get('subject_id');

                % find results file
                results_file_candidate_analysis = [this_data_folder_path filesep 'analysis' filesep makeFileName(collection_date, subject_id, file_label) '.mat'];
                results_file_candidate_subject = [this_data_folder_path filesep makeFileName(collection_date, subject_id, file_label) '.mat'];
                results_file_candidate_results = [this_data_folder_path filesep 'results' filesep  makeFileName(collection_date, subject_id, file_label) '.mat'];
                if exist(results_file_candidate_analysis, 'file')
                    results_file_name = results_file_candidate_analysis;
                end    
                if exist(results_file_candidate_subject, 'file')
                    results_file_name = results_file_candidate_subject;
                end    
                if exist(results_file_candidate_results, 'file')
                    results_file_name = results_file_candidate_results;
                end    

                % load data
                disp(['loading data from ' results_file_name])
                loaded_data = load(results_file_name);
                this.data{i_folder} = loaded_data;
                
            end
        end
        function data = getData(this, variable_names, variable_types)
            data = struct;

            number_of_variables_to_get = length(variable_names);
            conditions_settings = this.study_settings.get('conditions');
            condition_labels = conditions_settings(:, 1)';
            condition_source_variables = conditions_settings(:, 2)';
            number_of_condition_labels = length(condition_labels);
            
            data.condition_data = {};
            data.variable_data = cell(number_of_variables_to_get, 1);
            data.directions = cell(number_of_variables_to_get, 2);

            data.step_time_data = [];
            pushoff_time_data = [];
            data.bands_per_stretch = [];
            
            % extract data from individual
            for i_folder = 1 : this.number_of_source_sessions
                this_session_data = this.data{i_folder};
                number_of_stretches_this_session = length(this_session_data.time_list_session);
                bands_per_stretch_this_session = this_session_data.bands_per_stretch;

                % transform conditions into cell array
                conditions_session = this_session_data.conditions_session;
                condition_array_session = cell(number_of_stretches_this_session, number_of_condition_labels);
                for i_condition = 1 : number_of_condition_labels
                    condition_array_session(:, i_condition) = conditions_session.(condition_source_variables{i_condition});
                end

                data.condition_data = [data.condition_data; condition_array_session];
%                 origin_trial_list_all = [origin_trial_list_all; this_session_data.origin_trial_list_session]; %#ok<AGROW>
%                 origin_start_time_list_all = [origin_start_time_list_all; this_session_data.origin_start_time_list_session]; %#ok<AGROW>
%                 origin_end_time_list_all = [origin_end_time_list_all; this_session_data.origin_end_time_list_session]; %#ok<AGROW>

                % extract data
                for i_variable = 1 : number_of_variables_to_get

                    this_variable_name = variable_names{i_variable};
                    this_variable_type = variable_types{i_variable};
                    this_variable_source_index = find(strcmp(this_session_data.([this_variable_type '_names_session']), this_variable_name), 1, 'first');
                    if isempty(this_variable_source_index)
                        error(['Variable not found: ' this_variable_name])
                    end
                    this_variable_data = this_session_data.([this_variable_type '_data_session']){this_variable_source_index};
                    this_variable_directions = this_session_data.([this_variable_type '_directions_session'])(this_variable_source_index, :);

                    % store
                    data.variable_data{i_variable} = [data.variable_data{i_variable} this_variable_data];
                    data.directions(i_variable, :) = this_variable_directions;
                end

                % get time variables
                step_time_available = 0;
                if isfield(this_session_data, 'stretch_names_session') && any(find(strcmp(this_session_data.stretch_names_session, 'step_time')))
                    index_in_saved_data = find(strcmp(this_session_data.stretch_names_session, 'step_time'), 1, 'first');
                    this_step_time_data = this_session_data.stretch_data_session{index_in_saved_data};
                    data.step_time_data = [data.step_time_data this_step_time_data];
                    step_time_available = 1;
                end
                pushoff_time_available = 0;
                if isfield(this_session_data, 'stretch_names_session') && any(find(strcmp(this_session_data.stretch_names_session, 'pushoff_time')))
                    index_in_saved_data = find(strcmp(this_session_data.stretch_names_session, 'pushoff_time'), 1, 'first');
                    this_pushoff_time_data = this_session_data.stretch_data_session{index_in_saved_data};
                    pushoff_time_data = [pushoff_time_data this_pushoff_time_data]; %#ok<AGROW>
                    pushoff_time_available = 1;
                end
                if isempty(data.bands_per_stretch)
                    data.bands_per_stretch = bands_per_stretch_this_session;
                else
                    if data.bands_per_stretch ~= bands_per_stretch_this_session
                       warning('Different sessions have different numbers of bands per stretch') 
                    end
                end
                
            end
            % calculate mean pushoff index
            if step_time_available && pushoff_time_available
                pushoff_time_ratio = pushoff_time_data ./ data.step_time_data;
                mean_pushoff_ratio = mean(pushoff_time_ratio, 2);
                data.pushoff_index = round(mean_pushoff_ratio * 100);
            end

        end
    end
end








