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

% note: the inversion could be automated more elegantly, but that's for a later date

function processAnalysisVariables(varargin)
    load('subjectInfo.mat', 'date', 'subject_id');
    % load settings and existing results
    study_settings_file = '';
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    results_file_name = ['analysis' filesep makeFileName(date, subject_id, 'results')];
    loaded_data = load(results_file_name);
    
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    bands_per_stretch = loaded_data.bands_per_stretch;
    stretch_names_session = loaded_data.stretch_names_session;
    stretch_data_session = loaded_data.stretch_data_session;
    stretch_directions_session = loaded_data.stretch_directions_session;
    
    number_of_stretch_variables = length(loaded_data.stretch_names_session);
    number_of_stretches = size(loaded_data.stretch_data_session{1}, 2); %#ok<*USENS>
    
    % make condition data tables
    conditions_settings = study_settings.get('conditions');
    condition_labels = conditions_settings(:, 1)';
    condition_source_variables = conditions_settings(:, 2)';
    number_of_condition_labels = length(condition_labels);
    conditions_session = loaded_data.conditions_session;
    condition_data_all = cell(number_of_stretches, number_of_condition_labels);
    for i_condition = 1 : number_of_condition_labels
        condition_data_all(:, i_condition) = conditions_session.(condition_source_variables{i_condition});
    end
    labels_to_ignore = study_settings.get('conditions_to_ignore');
    levels_to_remove = study_settings.get('levels_to_remove');
    [condition_combination_labels, ~, condition_combinations_control] = determineConditionCombinations(condition_data_all, conditions_settings, labels_to_ignore, levels_to_remove);
    condition_combinations_control_unique = table2cell(unique(cell2table(condition_combinations_control), 'rows'));
    
    if isfield(loaded_data, 'analysis_data_session')
        analysis_data_session = loaded_data.analysis_data_session;
        analysis_directions_session = loaded_data.analysis_directions_session;
        analysis_names_session = loaded_data.analysis_names_session;
    else
        analysis_data_session = {};
        analysis_directions_session = {};
        analysis_names_session = {};
    end
    
    %% calculate response (i.e. difference from control mean)
    response_data_session = {};
    response_directions_session = loaded_data.stretch_directions_session;
    response_names_session = loaded_data.stretch_names_session;
    if ~isempty(condition_combinations_control)
        % prepare container
        response_data_session = cell(size(loaded_data.stretch_data_session));
        for i_variable = 1 : number_of_stretch_variables
            response_data_session{i_variable} = zeros(size(loaded_data.stretch_data_session{i_variable}));
        end        
        
        % go stretch by stretch
        for i_stretch = 1 : number_of_stretches
            % extract this stretches relevant conditions
            this_stretch_condition_string = cell(1, length(condition_combination_labels));
            for i_label = 1 : length(condition_combination_labels)
                this_label = condition_combination_labels{i_label};
                this_label_source_variable = condition_source_variables{strcmp(condition_labels, this_label)};
                this_label_condition_list = conditions_session.(this_label_source_variable);
                this_stretch_condition_string{i_label} = this_label_condition_list{i_stretch};
            end
            
            % determine applicable control condition index
            if strcmp(study_settings.get('experimental_paradigm'), 'Vision') || strcmp(study_settings.get('experimental_paradigm'), 'GVS') || strcmp(study_settings.get('experimental_paradigm'), 'GVS_old')
                if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
                    applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_LEFT'));
                end
                if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
                    applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT'));
                end
            end
            if strcmp(study_settings.get('experimental_paradigm'), 'CadenceGVS')
                if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'cadence')}, '80BPM') && strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
                    applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'cadence')), '80BPM') & strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_LEFT'));
                elseif strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'cadence')}, '80BPM') && strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
                    applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'cadence')), '80BPM') & strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT'));
                end
                if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'cadence')}, '110BPM') && strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
                    applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'cadence')), '110BPM') & strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_LEFT'));
                elseif strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'cadence')}, '110BPM') && strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
                     applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'cadence')), '110BPM') & strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT'));
                end
            end
            if strcmp(study_settings.get('stimulus_condition'), 'VISUAL') || strcmp(study_settings.get('stimulus_condition'), 'GVS')
                if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'stance_foot')}, 'STANCE_LEFT')
                    applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'stance_foot')), 'STANCE_LEFT'));
                end
                if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'stance_foot')}, 'STANCE_RIGHT')
                    applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'stance_foot')), 'STANCE_RIGHT'));
                end
            end
            
            % determine indicator for control
            control_condition_indicator = true(number_of_stretches, 1);
            for i_label = 1 : length(condition_combination_labels)
                this_label = condition_combination_labels{i_label};
                this_label_list = condition_data_all(:, strcmp(conditions_settings(:, 1), this_label));
                this_label_indicator = strcmp(this_label_list, condition_combinations_control_unique(applicable_control_condition_index, i_label));
                control_condition_indicator = control_condition_indicator .* this_label_indicator;
            end        
            control_condition_indicator = logical(control_condition_indicator);            
            
            % calculate responses
            for i_variable = 1 : number_of_stretch_variables
                % calculate control mean
                data_this_variable = loaded_data.stretch_data_session{i_variable};
                this_condition_control_data = data_this_variable(:, control_condition_indicator);
                % why are there nans in the control data??
                if any(any(isnan(this_condition_control_data)))
                    [rows, col_to_remove] = find(isnan(this_condition_control_data));
                    col_to_remove = unique(col_to_remove);
                    this_condition_control_data(:, col_to_remove) = []; 
                end
                this_condition_control_mean = mean(this_condition_control_data, 2);
                
                % calculate response
                response_data_session{i_variable}(:, i_stretch) = loaded_data.stretch_data_session{i_variable}(:, i_stretch) - this_condition_control_mean;
            end
            
        end

    end

    %% load necessary info
%     variables_to_integrate = study_settings.get('analysis_variables_from_integration');
    step_time_index_in_saved_data = find(strcmp(loaded_data.stretch_names_session, 'step_time'), 1, 'first');
    this_step_time_data = loaded_data.stretch_data_session{step_time_index_in_saved_data};
    pushoff_time_index_in_saved_data = find(strcmp(loaded_data.stretch_names_session, 'pushoff_time'), 1, 'first');
    if ~isempty(pushoff_time_index_in_saved_data)
        this_pushoff_time_data = loaded_data.stretch_data_session{pushoff_time_index_in_saved_data};
    else
        this_pushoff_time_data = [];
    end
        
    %% calculate band end variables
    variables_step_end = study_settings.get('analysis_variables_from_band_end');
    for i_variable = 1 : size(variables_step_end, 1)
        this_variable_name = variables_step_end{i_variable, 1};
        this_variable_source_name = variables_step_end{i_variable, 2};
        this_variable_source_type = variables_step_end{i_variable, 3};
        % pick data depending on source specification
        eval(['data_source = ' this_variable_source_type '_data_session;']);
        eval(['names_source = ' this_variable_source_type '_names_session;']);
        eval(['directions_source = ' this_variable_source_type '_directions_session;']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        new_variable_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        step_end_data = zeros(bands_per_stretch, number_of_stretches);
        for i_band = 1 : bands_per_stretch
            [~, end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
            step_end_data(i_band, :) = this_variable_source_data(end_index, :);
        end
        % store
        [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
            addOrReplaceResultsData ...
              ( ...
                analysis_data_session, analysis_names_session, analysis_directions_session, ...
                step_end_data, this_variable_name, new_variable_directions ...
              );
    end
    
    %% calculate variables referenced by stretch
    variables_referenced_by_stretch = study_settings.get('analysis_variables_referenced_by_stretch');
    for i_variable = 1 : size(variables_referenced_by_stretch, 1)
        this_variable_name = variables_referenced_by_stretch{i_variable, 1};
        this_variable_source_name = variables_referenced_by_stretch{i_variable, 2};
        this_variable_source_type = variables_referenced_by_stretch{i_variable, 3};
        this_variable_reference_band_index = str2num(variables_referenced_by_stretch{i_variable, 4});
        this_variable_reference_point_percentage_within_band = str2num(variables_referenced_by_stretch{i_variable, 5});
        
        % pick data depending on source specification
        eval(['data_source = ' this_variable_source_type '_data_session;']);
        eval(['names_source = ' this_variable_source_type '_names_session;']);
        eval(['directions_source = ' this_variable_source_type '_directions_session;']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        new_variable_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        new_variable_data = zeros(size(this_variable_source_data));
        reference_index_within_band = round((number_of_time_steps_normalized-1)*this_variable_reference_point_percentage_within_band * 1/100) + 1;
        relevant_band_start_index = getBandIndices(this_variable_reference_band_index, number_of_time_steps_normalized);
        reference_index = relevant_band_start_index + reference_index_within_band - 1;
        
        for i_stretch = 1 : size(new_variable_data, 2)
            new_variable_data(:, i_stretch) = this_variable_source_data(:, i_stretch) - this_variable_source_data(reference_index, i_stretch);
        end
        % store
        [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
            addOrReplaceResultsData ...
              ( ...
                analysis_data_session, analysis_names_session, analysis_directions_session, ...
                new_variable_data, this_variable_name, new_variable_directions ...
              );
    end
    
    %% calculate integrated variables
    variables_to_integrate_header = study_settings.get('analysis_variables_from_integration_header');
    variables_to_integrate = study_settings.get('analysis_variables_from_integration');
    names_source = response_names_session;
    directions_source = response_directions_session;
    for i_variable = 1 : size(variables_to_integrate, 1)
        this_variable_name = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'new_variable_name')};
        this_variable_source_name = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'source_variable_name')};
        this_variable_source_type = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'source_variable_type')};
        start_info = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'start')};
        start_variable_source_type = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'start_variable_type')};
        end_info = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'end')};
        end_variable_source_type = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'end_variable_type')};
        % pick data depending on source specification
        eval(['data_source = ' this_variable_source_type '_data_session;']);
        eval(['names_source = ' this_variable_source_type '_names_session;']);
        eval(['directions_source = ' this_variable_source_type '_directions_session;']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        this_variable_source_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        
        % determine start and end of integration in percent of band time
        if strcmp(start_variable_source_type, 'percentage')
            start_data_percent = ones(size(this_step_time_data)) * str2num(start_info);
        else
            eval(['start_data_source = ' start_variable_source_type '_data_session;']);
            eval(['start_names_source = ' start_variable_source_type '_names_session;']);
            start_data_time_within_band = start_data_source{strcmp(start_names_source, start_info)};
            start_data_ratio = start_data_time_within_band ./ this_step_time_data;
            start_data_percent = round(start_data_ratio * 100);
        end
        if strcmp(end_variable_source_type, 'percentage')
            end_data_percent = ones(size(this_step_time_data)) * str2num(end_info);
        else
            eval(['end_data_source = ' end_variable_source_type '_data_session;']);
            eval(['end_names_source = ' end_variable_source_type '_names_session;']);
            end_data_time_within_band = end_data_source{strcmp(end_names_source, end_info)};
            end_data_ratio = end_data_time_within_band ./ this_step_time_data;
            end_data_percent = round(end_data_ratio * 100);
        end
        
        % integrate
        integrated_data = zeros(bands_per_stretch, number_of_stretches);
        for i_stretch = 1 : number_of_stretches
            for i_band = 1 : bands_per_stretch
                [band_start_index, band_end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
                this_band_time_full = linspace(0, this_step_time_data(i_band), 100);
                this_band_data_full = this_variable_source_data(band_start_index : band_end_index, i_stretch);

                range = start_data_percent(i_band, i_stretch) : end_data_percent(i_band, i_stretch);
                this_band_time_range = this_band_time_full(range);
                this_band_data_range = this_band_data_full(range);
                
                % integrate
                this_band_data_integrated = cumtrapz(this_band_time_range, this_band_data_range);
                integrated_data(i_band, i_stretch) = this_band_data_integrated(end);
            end
        end        
        
        % store
        [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
            addOrReplaceResultsData ...
              ( ...
                analysis_data_session, analysis_names_session, analysis_directions_session, ...
                integrated_data, this_variable_name, this_variable_source_directions ...
              );

    end
    
    %% calculate inversion variables
    inversion_variables = study_settings.get('inversion_variables');
    for i_variable = 1 : size(inversion_variables, 1)
        % get data
        this_variable_name = inversion_variables{i_variable, 1};
        this_variable_source_name = inversion_variables{i_variable, 2};
        this_variable_source_type = inversion_variables{i_variable, 3};
        % pick data depending on source specification
        eval(['data_source = ' this_variable_source_type '_data_session;']);
        eval(['names_source = ' this_variable_source_type '_names_session;']);
        eval(['directions_source = ' this_variable_source_type '_directions_session;']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        this_variable_source_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        new_variable_directions = inversion_variables(i_variable, 6:7);
        
        relevant_condition = inversion_variables{i_variable, 4};
        inversion_table = study_settings.get(inversion_variables{i_variable, 5});
        
        % go through levels and invert
        this_variable_data = this_variable_source_data;
        level_list = conditions_session.(condition_source_variables{strcmp(condition_labels, relevant_condition)});
        for i_level = 1 : size(inversion_table, 1)
            % determine whether this has to be inverted
            this_level_direction_map = inversion_table(i_level, 2:3);
            if strcmp(this_level_direction_map{1}, this_variable_source_directions{1}) && strcmp(this_level_direction_map{2}, this_variable_source_directions{2})
                % directions of the new variable and the source variable are the same, no need to invert here
                sign_this_level = 1;
            elseif strcmp(this_level_direction_map{1}, this_variable_source_directions{2}) && strcmp(this_level_direction_map{2}, this_variable_source_directions{1})
                % positive direction for new variable is negative for source variable, and vice versa, so we need to invert data for this level
                sign_this_level = -1;
            else
                error(['Trying to invert variable ' this_variable_name ', but direction labels do not match.'])
            end
            
            % get matches
            label_this_level = inversion_table{i_level, 1};
            match_this_level = strcmp(level_list, label_this_level);
            
            % invert
            this_variable_data(:, match_this_level) = sign_this_level * this_variable_data(:, match_this_level);
        end

        % store
        [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
            addOrReplaceResultsData ...
              ( ...
                analysis_data_session, analysis_names_session, analysis_directions_session, ...
                this_variable_data, this_variable_name, new_variable_directions ...
              );
    end
    
    %% gather variables that are selected from different sources depending on condition
    % TODO: deal with bands
    % TODO: deal with directions
    variables_to_select = study_settings.get('analysis_variables_from_selection');
    for i_variable = 1 : size(variables_to_select, 1)
        % get signs
        if strcmp(variables_to_select{i_variable, 5}, '+')
            sign_trigger_left = 1;
        elseif strcmp(variables_to_select{i_variable, 5}, '-')
            sign_trigger_left = -1;
        else
            error('Sign must be either "+" or "-"')
        end
        if strcmp(variables_to_select{i_variable, 6}, '+')
            sign_trigger_right = 1;
        elseif strcmp(variables_to_select{i_variable, 6}, '-')
            sign_trigger_right = -1;
        else
            error('Sign must be either "+" or "-"')
        end
        
        % get data
        this_variable_name = variables_to_select{i_variable, 1};
        this_variable_source_name_triggerLeft = variables_to_select{i_variable, 3};
        this_variable_source_name_triggerRight = variables_to_select{i_variable, 4};
        this_variable_source_type = variables_to_select{i_variable, 2};
        if strcmp(this_variable_source_type, 'response')
            this_variable_source_index_triggerLeft = find(strcmp(response_names_session, this_variable_source_name_triggerLeft), 1, 'first');
            this_variable_source_index_triggerRight = find(strcmp(response_names_session, this_variable_source_name_triggerRight), 1, 'first');
            if isempty(this_variable_source_index_triggerLeft)
                error(['Data not found: ' this_variable_source_name_triggerLeft])
            end
            if isempty(this_variable_source_index_triggerRight)
                error(['Data not found: ' this_variable_source_name_triggerRight])
            end
            this_variable_source_data_triggerLeft = response_data_session{this_variable_source_index_triggerLeft};
            this_variable_source_data_triggerRight = response_data_session{this_variable_source_index_triggerRight};
            new_variable_directions = response_directions_session(strcmp(response_names_session, this_variable_source_name_triggerLeft), :);
        end
        if strcmp(this_variable_source_type, 'analysis')
            this_variable_source_index_triggerLeft = find(strcmp(analysis_names_session, this_variable_source_name_triggerLeft), 1, 'first');
            this_variable_source_index_triggerRight = find(strcmp(analysis_names_session, this_variable_source_name_triggerRight), 1, 'first');
            if isempty(this_variable_source_index_triggerLeft)
                error(['Data not found: ' this_variable_source_name_triggerLeft])
            end
            if isempty(this_variable_source_index_triggerRight)
                error(['Data not found: ' this_variable_source_name_triggerRight])
            end
            this_variable_source_data_triggerLeft = analysis_data_session{this_variable_source_index_triggerLeft};
            this_variable_source_data_triggerRight = analysis_data_session{this_variable_source_index_triggerRight};
            new_variable_directions = analysis_directions_session(strcmp(analysis_names_session, this_variable_source_name_triggerLeft), :);
        end
        
        % select
        this_variable_data = zeros(size(this_variable_source_data_triggerLeft));
        stance_foot_list = conditions_session.(condition_source_variables{strcmp(condition_labels, 'stance_foot')});
        index_list = conditions_session.(condition_source_variables{strcmp(condition_labels, 'index')});
        for i_stretch = 1 : number_of_stretches
            if ...
              (strcmp(stance_foot_list{i_stretch}, 'STANCE_LEFT') && strcmp(index_list{i_stretch}, 'ONE')) || ...
              (strcmp(stance_foot_list{i_stretch}, 'STANCE_RIGHT') && strcmp(index_list{i_stretch}, 'TWO')) || ...
              (strcmp(stance_foot_list{i_stretch}, 'STANCE_LEFT') && strcmp(index_list{i_stretch}, 'THREE')) || ...
              (strcmp(stance_foot_list{i_stretch}, 'STANCE_RIGHT') && strcmp(index_list{i_stretch}, 'FOUR'))
                this_variable_data(:, i_stretch) = sign_trigger_left * this_variable_source_data_triggerLeft(:, i_stretch);
            end
            if ...
              (strcmp(stance_foot_list{i_stretch}, 'STANCE_RIGHT') && strcmp(index_list{i_stretch}, 'ONE')) || ...
              (strcmp(stance_foot_list{i_stretch}, 'STANCE_LEFT') && strcmp(index_list{i_stretch}, 'TWO')) || ...
              (strcmp(stance_foot_list{i_stretch}, 'STANCE_RIGHT') && strcmp(index_list{i_stretch}, 'THREE')) || ...
              (strcmp(stance_foot_list{i_stretch}, 'STANCE_LEFT') && strcmp(index_list{i_stretch}, 'FOUR'))
                this_variable_data(:, i_stretch) = sign_trigger_right * this_variable_source_data_triggerRight(:, i_stretch);
            end
        end
        
        % store
        [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
            addOrReplaceResultsData ...
              ( ...
                analysis_data_session, analysis_names_session, analysis_directions_session, ...
                this_variable_data, this_variable_name, new_variable_directions ...
              );
        % TODO: check whether the directions are actually still correct here
    end
    
    %% process variables where something specific happens for each variable
    % TO DO: automate the source type and data extraction
    special_variables_to_calculate = study_settings.get('analysis_variables_special');
    for i_variable = 1:size(special_variables_to_calculate, 1)
        this_variable_name = special_variables_to_calculate{i_variable, 1};
        this_variable_source_name = special_variables_to_calculate{i_variable, 2};
        this_variable_source_type = special_variables_to_calculate{i_variable, 3};
        
          % pick data depending on source specification
        eval(['data_source = ' this_variable_source_type '_data_session;']);
        eval(['names_source = ' this_variable_source_type '_names_session;']);
        eval(['directions_source = ' this_variable_source_type '_directions_session;']);
%         this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        this_variable_source_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        new_variable_directions = this_variable_source_directions;
          
        if strcmp(this_variable_name, 'step_symmetry_index')
            this_variable_source_index = find(strcmp(stretch_names_session, this_variable_source_name), 1, 'first');
            this_variable_source_data = stretch_data_session{this_variable_source_index};

            average_step_time = mean(reshape(this_variable_source_data, 1, length(this_variable_source_data)*2));
            for i_stretch = 1:  length(this_variable_source_data)
                % find the left and right stance data
                if strcmp(condition_data_all(i_stretch,3), 'STANCE_LEFT')
                    left_step_index = 2;
                    right_step_index = 1;
                else
                    left_step_index = 1;
                    right_step_index = 2;
                end
                this_left_step_time = this_variable_source_data(left_step_index,i_stretch);
                this_right_step_time = this_variable_source_data(right_step_index,i_stretch);
                this_variable_data(1:2,i_stretch) = (this_left_step_time - this_right_step_time) / average_step_time;
            end
        end
        if strcmp(this_variable_name, 'com_from_com_init_x')
            % should we be looking at response or stretch variable here??
            this_variable_source_index = find(strcmp(response_names_session, this_variable_source_name), 1, 'first');
            this_variable_source_data = response_data_session{this_variable_source_index};  
            this_variable_data = [];
            for i_stretch = 1:length(this_variable_source_data)
                this_variable_data(:,i_stretch) = this_variable_source_data(:,i_stretch) - this_variable_source_data(1,i_stretch);
            end
        end
       if strcmp(this_variable_name, 'trigger_leg_ankle_dorsiflexion_inverted_max')
            this_variable_source_index = find(strcmp(analysis_names_session, this_variable_source_name), 1, 'first');
            this_variable_source_data = analysis_data_session{this_variable_source_index};
            
            
            for i_stretch = 1:length(this_variable_source_data)
                % create time
                this_stretch_time_full = linspace(0, this_step_time_data(i_stretch), 100);
                this_stretch_data_full = this_variable_source_data(:, i_stretch);

                % interpolate double stance to 100 data points
                this_stretch_time_double = linspace(0, this_pushoff_time_data(i_stretch), 100);       
                this_stretch_data_double = interp1(this_stretch_time_full, this_stretch_data_full, this_stretch_time_double);


                this_dorsi_angle_value = max(findpeaks(this_stretch_data_double));
                if ~isempty(this_dorsi_angle_value)
                    this_variable_data(i_stretch) = this_dorsi_angle_value;
                else
                    this_variable_data(i_stretch) = NaN;
                end
            end
       end
       if strcmp(this_variable_name, 'cop_from_com_x_integrated_twice')
            this_variable_source_index = find(strcmp(response_names_session, this_variable_source_name), 1, 'first');
            this_variable_source_data = response_data_session{this_variable_source_index};
            number_of_stretches = size(this_variable_source_data, 2);

            for i_stretch = 1 : number_of_stretches
                % get data for full step
                this_stretch_time_full = linspace(0, this_step_time_data(i_stretch), 100);
                this_stretch_data_full = this_variable_source_data(:, i_stretch);

                % interpolate single stance to 100 data points
                this_stretch_time_single = linspace(this_pushoff_time_data(i_stretch), this_step_time_data(i_stretch), 100);
                this_stretch_data_single = interp1(this_stretch_time_full, this_stretch_data_full, this_stretch_time_single);

                % integrate data in single stance
                this_stretch_data_single_integrated = cumtrapz(this_stretch_time_single, this_stretch_data_single);
                this_stretch_data_single_integrated_twice = cumtrapz(this_stretch_time_single, this_stretch_data_single_integrated);
                this_variable_data(i_stretch) = this_stretch_data_single_integrated_twice(end);
            end 
       end
              % store
        [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
            addOrReplaceResultsData ...
              ( ...
                analysis_data_session, analysis_names_session, analysis_directions_session, ...
                this_variable_data, this_variable_name, new_variable_directions ...
              );
    end
    
    %% calculate variables from extrema
    variables_from_extrema = study_settings.get('analysis_variables_from_extrema');
    for i_variable = 1 : size(variables_from_extrema, 1)
        % get data
        this_variable_name = variables_from_extrema{i_variable, 1};
        this_variable_source_name = variables_from_extrema{i_variable, 2};
        this_variable_source_type = variables_from_extrema{i_variable, 3};
        this_variable_extremum_type = variables_from_extrema{i_variable, 4};
        % pick data depending on source specification
        eval(['data_source = ' this_variable_source_type '_data_session;']);
        eval(['names_source = ' this_variable_source_type '_names_session;']);
        eval(['directions_source = ' this_variable_source_type '_directions_session;']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        new_variable_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        
        extrema_data = zeros(bands_per_stretch, number_of_stretches);
        for i_band = 1 : bands_per_stretch
            [band_start_index, band_end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
            data_this_band = this_variable_source_data(band_start_index : band_end_index, :);
            if strcmp(this_variable_extremum_type, 'min')
                extrema_data(i_band, :) = min(data_this_band);
            end
            if strcmp(this_variable_extremum_type, 'max')
                extrema_data(i_band, :) = max(data_this_band);
            end
            if ~strcmp(this_variable_extremum_type, 'min') && ~strcmp(this_variable_extremum_type, 'max')
                error(['"' this_variable_extremum_type '" is not a valid type for variables_from_extrema. Acceptable types are "min" or "max".']);
            end
        end
        % store
        [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
            addOrReplaceResultsData ...
              ( ...
                analysis_data_session, analysis_names_session, analysis_directions_session, ...
                extrema_data, this_variable_name, new_variable_directions ...
              );
    end

    %% calculate variables from inversion
    % TODO: deal with bands
    variables_to_invert = study_settings.get('analysis_variables_from_inversion');
    for i_variable = 1 : size(variables_to_invert, 1)
        % get data
        this_variable_name = variables_to_invert{i_variable, 1};
        this_variable_source_name = variables_to_invert{i_variable, 2};
        this_variable_source_type = variables_to_invert{i_variable, 3};
        % pick data depending on source specification
        eval(['data_source = ' this_variable_source_type '_data_session;']);
        eval(['names_source = ' this_variable_source_type '_names_session;']);
        eval(['directions_source = ' this_variable_source_type '_directions_session;']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        new_variable_directions = variables_to_invert(i_variable, 4:5);
        
        relevant_condition = variables_to_invert{i_variable, 6};
        condition_sign_map = reshape(variables_to_invert(i_variable, 7:end), 2, (size(variables_to_invert, 2)-6)/2)';
        
        % go through levels and invert
        this_variable_data = this_variable_source_data;
        level_list = conditions_session.(condition_source_variables{strcmp(condition_labels, relevant_condition)});
        for i_level = 1 : size(condition_sign_map, 1)
            % get sign
            if strcmp(condition_sign_map{i_level, 2}, '+')
                sign_this_level = 1;
            elseif strcmp(condition_sign_map{i_level, 2}, '-')
                sign_this_level = -1;
            else
                error('Sign must be either "+" or "-"')
            end
            
            % get matches
            label_this_level = condition_sign_map{i_level, 1};
            match_this_level = strcmp(level_list, label_this_level);
            
            % invert
            this_variable_data(:, match_this_level) = sign_this_level * this_variable_data(:, match_this_level);
        end

        
        
        % store
        [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
            addOrReplaceResultsData ...
              ( ...
                analysis_data_session, analysis_names_session, analysis_directions_session, ...
                this_variable_data, this_variable_name, new_variable_directions ...
              );
    end
    
    variables_to_affected_side = study_settings.get('variables_to_affected_side');
    for i_variable = 1 : size(variables_to_affected_side, 1)
        % get data
        this_affected_variable_source_name = variables_to_affected_side{i_variable, 1};
        this_unaffected_variable_source_name = variables_to_affected_side{i_variable, 2};
        source_affected_side_info =  variables_to_affected_side{i_variable, 3};
        left_sided_variable_name = variables_to_affected_side{i_variable, 4};
        right_sided_variable_name = variables_to_affected_side{i_variable, 5};
        this_variable_source_type = variables_to_affected_side{i_variable, 6};
        
        % may get an error here with AS data.. probably need to add
        % arm_angle and phase variables to new code structure (add
        % directions)
        
        eval(['data_source = ' this_variable_source_type '_data_session;']);
        eval(['names_source = ' this_variable_source_type '_names_session;']);
        eval(['directions_source = ' this_variable_source_type '_directions_session;']);
              
        this_affected_side_info = conditions_session.condition_affectedSide_list; % check this
        
        % only need to check one index for this type of variable
        if strcmp(this_affected_side_info{1}, 'L')
            this_affected_variable_source_name = left_sided_variable_name;
            this_unaffected_variable_source_name = right_sided_variable_name;
            this_affected_variable_data = data_source{strcmp(names_source, this_affected_variable_source_name)};
            this_unaffected_variable_data = data_source{strcmp(names_source, this_unaffected_variable_source_name)};
                        
        elseif strcmp(this_affected_side_info{1}, 'R')
            this_affected_variable_source_name = right_sided_variable_name;
            this_unaffected_variable_source_name = left_sided_variable_name;
            this_affected_variable_data = data_source{strcmp(names_source, this_affected_variable_source_name)};
            this_unaffected_variable_data = data_source{strcmp(names_source, this_unaffected_variable_source_name)};
        else
            Warning('Either the variable specificed in studySettings.txt cannot be processed here or the affectedSide info is innappropriate')
        end
        
        % only need to take 1 side, assuming the variables are the same
        % type
        new_variable_directions = directions_source(strcmp(names_source, left_sided_variable_name), :);
        
         % store affected
        [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
            addOrReplaceResultsData ...
              ( ...
                analysis_data_session, analysis_names_session, analysis_directions_session, ...
                this_affected_variable_data, 'affected_arm_angle', new_variable_directions ...
              );
         % store unaffected
        [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
            addOrReplaceResultsData ...
              ( ...
                analysis_data_session, analysis_names_session, analysis_directions_session, ...
                this_unaffected_variable_data, 'unaffected_arm_angle', new_variable_directions ...
              ); 
    end
    
    
    
    
    
%% LEGACY CODE
    
    %% calculate step end variables
    variables_step_end = study_settings.get('analysis_variables_from_step_end');
    names_source = response_names_session;
    directions_source = response_directions_session;
    for i_variable = 1 : size(variables_step_end, 1)
        this_variable_name = variables_step_end{i_variable, 1};
        this_variable_source_name = variables_step_end{i_variable, 2};
        this_variable_response_data = response_data_session{strcmp(loaded_data.stretch_names_session, this_variable_source_name)};
        new_variable_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        step_end_data = this_variable_response_data(end, :);
        
        % store
        [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
            addOrReplaceResultsData ...
              ( ...
                analysis_data_session, analysis_names_session, analysis_directions_session, ...
                step_end_data, this_variable_name, new_variable_directions ...
              );
    end
        
    %% gather variables with inversion by perturbation
    % THIS IS LEGACY CODE
    % used this for the Vision experiment, it doesn't deal with bands or directions
    % use the general solution for variables_to_invert instead
    variables_to_invert = study_settings.get('analysis_variables_from_inversion_by_perturbation');
    for i_variable = 1 : size(variables_to_invert, 1)
        warning('analysis_variables_from_inversion_by_perturbation is in the process of being phased out, look for another solution.')
        % get data
        this_variable_name = variables_to_invert{i_variable, 1};
        this_variable_source_name = variables_to_invert{i_variable, 2};
        this_variable_source_type = variables_to_invert{i_variable, 3};
        if strcmp(this_variable_source_type, 'response')
            this_variable_source_index = find(strcmp(response_names_session, this_variable_source_name), 1, 'first');
            if isempty(this_variable_source_index)
                error(['Data not found: ' this_variable_source_name])
            end
            this_variable_source_data = response_data_session{this_variable_source_index};
            new_variable_directions = response_directions_session(strcmp(response_names_session, this_variable_source_name), :);
        end
        if strcmp(this_variable_source_type, 'analysis')
            this_variable_source_index = find(strcmp(analysis_names_session, this_variable_source_name), 1, 'first');
            if isempty(this_variable_source_index)
                error(['Data not found: ' this_variable_source_name])
            end
            this_variable_source_data = analysis_data_session{this_variable_source_index};
            new_variable_directions = analysis_directions_session(strcmp(analysis_names_session, this_variable_source_name), :);
        end
 
        
        % get signs
        if strcmp(variables_to_invert{i_variable, 4}, '+')
            sign_illusion_left = 1;
        elseif strcmp(variables_to_invert{i_variable, 4}, '-')
            sign_illusion_left = -1;
        else
            error('Sign must be either "+" or "-"')
        end
        if strcmp(variables_to_invert{i_variable, 5}, '+')
            sign_illusion_right = 1;
        elseif strcmp(variables_to_invert{i_variable, 5}, '-')
            sign_illusion_right = -1;
        else
            error('Sign must be either "+" or "-"')
        end
        
        % invert
        this_variable_data = this_variable_source_data;
        for i_stretch = 1 : number_of_stretches
            perturbation_list = conditions_session.(condition_source_variables{strcmp(condition_labels, 'perturbation')});
            if strcmp(perturbation_list{i_stretch}, 'ILLUSION_LEFT')
                this_variable_data(:, i_stretch) = sign_illusion_left * this_variable_source_data(:, i_stretch);
            elseif strcmp(perturbation_list{i_stretch}, 'ILLUSION_RIGHT')
                this_variable_data(:, i_stretch) = sign_illusion_right * this_variable_source_data(:, i_stretch);
            end
%             if strcmp(loaded_data.condition_perturbation_list_session{i_stretch}, 'ILLUSION_LEFT')
%                 this_variable_data(:, i_stretch) = sign_illusion_left * this_variable_source_data(:, i_stretch);
%             elseif strcmp(loaded_data.condition_perturbation_list_session{i_stretch}, 'ILLUSION_RIGHT')
%                 this_variable_data(:, i_stretch) = sign_illusion_right * this_variable_source_data(:, i_stretch);
%             end
        end
        
        % store
        [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
            addOrReplaceResultsData ...
              ( ...
                analysis_data_session, analysis_names_session, analysis_directions_session, ...
                this_variable_data, this_variable_name, new_variable_directions ...
              );
        % TODO: check whether the directions are actually still correct here
    end
    
    %% gather variables with inversion by direction
    % THIS IS LEGACY CODE
    % used this for the Vision experiment, it doesn't deal with bands or directions
    variables_to_invert = study_settings.get('analysis_variables_from_inversion_by_direction');
    for i_variable = 1 : size(variables_to_invert, 1)
        warning(['analysis_variables_from_inversion_by_direction is in the process of being phased out, look for another solution '])
        % get data
        this_variable_name = variables_to_invert{i_variable, 1};
        this_variable_source_name = variables_to_invert{i_variable, 2};
        this_variable_source_type = variables_to_invert{i_variable, 3};
        if strcmp(this_variable_source_type, 'response')
            this_variable_source_index = find(strcmp(response_names_session, this_variable_source_name), 1, 'first');
            this_variable_source_data = response_data_session{this_variable_source_index};
            new_variable_directions = response_directions_session(strcmp(response_names_session, this_variable_source_name), :);
        end
        if strcmp(this_variable_source_type, 'analysis')
            this_variable_source_index = find(strcmp(analysis_names_session, this_variable_source_name), 1, 'first');
            this_variable_source_data = analysis_data_session{this_variable_source_index};
            new_variable_directions = analysis_directions_session(strcmp(analysis_names_session, this_variable_source_name), :);
        end
        
        % get signs
        if strcmp(variables_to_invert{i_variable, 4}, '+')
            sign_illusion_left = 1;
        elseif strcmp(variables_to_invert{i_variable, 4}, '-')
            sign_illusion_left = -1;
        else
            error('Sign must be either "+" or "-"')
        end
        if strcmp(variables_to_invert{i_variable, 5}, '+')
            sign_illusion_right = 1;
        elseif strcmp(variables_to_invert{i_variable, 5}, '-')
            sign_illusion_right = -1;
        else
            error('Sign must be either "+" or "-"')
        end
        
        % invert
        this_variable_data = this_variable_source_data;
        for i_stretch = 1 : number_of_stretches
            direction_list = conditions_session.(condition_source_variables{strcmp(condition_labels, 'direction')});
            if strcmp(direction_list{i_stretch}, 'TOWARDS')
                this_variable_data(:, i_stretch) = sign_illusion_left * this_variable_source_data(:, i_stretch);
            elseif strcmp(direction_list{i_stretch}, 'AWAY')
                this_variable_data(:, i_stretch) = sign_illusion_right * this_variable_source_data(:, i_stretch);
            end
            
%             if strcmp(condition_direction_list_session{i_stretch}, 'TOWARDS')
%                 this_variable_data(:, i_stretch) = sign_illusion_left * this_variable_source_data(:, i_stretch);
%             elseif strcmp(condition_direction_list_session{i_stretch}, 'AWAY')
%                 this_variable_data(:, i_stretch) = sign_illusion_right * this_variable_source_data(:, i_stretch);
%             end
        end
        
        % store
        [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
            addOrReplaceResultsData ...
              ( ...
                analysis_data_session, analysis_names_session, analysis_directions_session, ...
                this_variable_data, this_variable_name, new_variable_directions ...
              );
        % TODO: check whether the directions are actually still correct here
    end

    %% calculate integrated variables - old
    % TODO: deal with bands
    % TODO: deal with directions .. check this
    variables_to_integrate = study_settings.get('analysis_variables_from_integration_old');
    names_source = response_names_session;
    directions_source = response_directions_session;
    for i_variable = 1 : size(variables_to_integrate, 1)
        this_variable_name = variables_to_integrate{i_variable, 1};
        this_variable_source_name = variables_to_integrate{i_variable, 2};
        this_variable_response_data = response_data_session{strcmp(loaded_data.stretch_names_session, this_variable_source_name)};
        number_of_stretches = size(this_variable_response_data, 2);
        new_variable_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        integrated_data = zeros(1, number_of_stretches);
        for i_stretch = 1 : number_of_stretches
            % get data for full step
            this_stretch_time_full = linspace(0, this_step_time_data(i_stretch), 100);
            this_stretch_data_full = this_variable_response_data(:, i_stretch);
            
            % interpolate single stance to 100 data points
            this_stretch_time_single = linspace(this_pushoff_time_data(i_stretch), this_step_time_data(i_stretch), 100);
            this_stretch_data_single = interp1(this_stretch_time_full, this_stretch_data_full, this_stretch_time_single);
            
            % integrate data in single stance
            this_stretch_data_single_integrated = cumtrapz(this_stretch_time_single, this_stretch_data_single);
            integrated_data(i_stretch) = this_stretch_data_single_integrated(end);
            
        end
        
        % store
        [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
            addOrReplaceResultsData ...
              ( ...
                analysis_data_session, analysis_names_session, analysis_directions_session, ...
                integrated_data, this_variable_name, new_variable_directions ...
              );

    end

    
    %% save data
    variables_to_save = loaded_data;
    variables_to_save.response_data_session = response_data_session;
    variables_to_save.response_directions_session = response_directions_session;
    variables_to_save.response_names_session = response_names_session;
    variables_to_save.analysis_data_session = analysis_data_session;
    variables_to_save.analysis_directions_session = analysis_directions_session;
    variables_to_save.analysis_names_session = analysis_names_session;
    save(results_file_name, '-struct', 'variables_to_save');    

    

end













