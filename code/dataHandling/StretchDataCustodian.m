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
        bands_per_stretch;
        normalized_time;
        
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
            this.bands_per_stretch = [];

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
                
                bands_per_stretch_this_session = loaded_data.bands_per_stretch;
                if isempty(this.bands_per_stretch)
                    this.bands_per_stretch = bands_per_stretch_this_session;
                else
                    if this.bands_per_stretch ~= bands_per_stretch_this_session
                       warning('Different sessions have different numbers of bands per stretch') 
                    end
                end
            end
            
            number_of_time_steps = this.study_settings.get('number_of_time_steps_normalized');
            this.normalized_time = 1 : (number_of_time_steps-1)*this.bands_per_stretch+1;
            
        end
        function condition_data = getConditionData(this)
            conditions_settings = this.study_settings.get('conditions');
            condition_labels = conditions_settings(:, 1)';
            condition_source_variables = conditions_settings(:, 2)';
            number_of_condition_labels = length(condition_labels);
            condition_data = {};
            for i_folder = 1 : this.number_of_source_sessions
                this_session_data = this.data{i_folder};
                
                % transform conditions into cell array
                number_of_stretches_this_session = length(this_session_data.time_list_session);
                conditions_session = this_session_data.conditions_session;
                condition_array_session = cell(number_of_stretches_this_session, number_of_condition_labels);
                for i_condition = 1 : number_of_condition_labels
                    condition_array_session(:, i_condition) = conditions_session.(condition_source_variables{i_condition});
                end

                condition_data = [condition_data; condition_array_session]; %#ok<AGROW>
            end            
        end
        function [subjects, trials, start_times, end_times] = getOriginData(this)
            % create containers
            subjects = {};
            trials = [];
            start_times = [];
            end_times = [];
            
            % go through source folders, grab data and collect in containers
            for i_folder = 1 : this.number_of_source_sessions
                this_session_data = this.data{i_folder};
                
                this_session_subject_data = this_session_data.conditions_session.subject_list;
                subjects = [subjects; this_session_subject_data]; %#ok<AGROW>
                trials = [trials; this_session_data.origin_trial_list_session]; %#ok<AGROW>
                start_times = [start_times; this_session_data.origin_start_time_list_session]; %#ok<AGROW>
                end_times = [end_times; this_session_data.origin_end_time_list_session]; %#ok<AGROW>
            end           
        end
        
        function data = getData(this, variable_names, variable_types)
            % normalize data format
            if ~iscell(variable_names)
                variable_names = {variable_names};
            end
            if ~iscell(variable_types)
                variable_types = {variable_types};
            end
            
            data = struct;
            data.variable_names = variable_names;
            data.variable_types = variable_types;

            number_of_variables_to_get = length(variable_names);
            data.variable_data = cell(number_of_variables_to_get, 1);
            data.directions = cell(number_of_variables_to_get, 2);
            
            % extract data from individual
            for i_folder = 1 : this.number_of_source_sessions
                this_session_data = this.data{i_folder};

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

                
            end
            if number_of_variables_to_get == 1
                data.variable_names = data.variable_names{1};
                data.variable_types = data.variable_types{1};
                data.variable_data = data.variable_data{1};
            end
        end
        function time = getNormalizedTime(this)
            time = this.normalized_time;
        end
    end
end








