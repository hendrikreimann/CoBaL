
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

function processAnalysisVariables(varargin)
    study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');

    % load settings and existing results
    results_file_name = ['results' filesep makeFileName(collection_date, subject_id, 'results')];
    data = load(results_file_name);
    
    % make condition data tables
    conditions.settings_table = study_settings.get('conditions');
    conditions.factor_labels = conditions.settings_table(:, 1)';
    conditions.source_variables = conditions.settings_table(:, 2)';
    conditions.number_of_factor_labels = length(conditions.factor_labels);
    conditions.conditions_session = data.conditions_session;
    
    % get table with analysis information
    analysis_table = getAnalysisTable(study_settings);
    
    % analyze
    for i_row = 1 : height(analysis_table)
        this_action = analysis_table{i_row, 'action'};
        this_settings_table_name = analysis_table{i_row, 'settings_table'};
        this_settings_table = study_settings.get(this_settings_table_name{1}, 1);
        settings_table_header_name = analysis_table{i_row, 'settings_table_header'};
        this_settings_table_header = study_settings.get(settings_table_header_name{1}, 1);
        
        if strcmp(this_action, 'calculate stimulus response')
            data = calculateStimulusResponse(this_settings_table, this_settings_table_header, study_settings, data, conditions);
        end
        if strcmp(this_action, 'integrate over time')
            data = calculateIntegratedVariables(this_settings_table, this_settings_table_header, study_settings, data);
        end
        if strcmp(this_action, 'integrate over range')
            data = integrateOverRange(this_settings_table, this_settings_table_header, study_settings, data);
        end
        if strcmp(this_action, 'select from multiple variables by condition')
            data = calculateSelectionVariables(this_settings_table, this_settings_table_header, study_settings, data, conditions);
        end
        if strcmp(this_action, 'invert by condition')
            data = calculateInversionVariables(this_settings_table, this_settings_table_header, study_settings, data, conditions);
        end
        if strcmp(this_action, 'calculate mean over time')
            data = calculateMeanVariables(this_settings_table, this_settings_table_header, study_settings, data);
        end
        if strcmp(this_action, 'calculate rms over time')
            data = calculateRmsVariables(this_settings_table, this_settings_table_header, study_settings, data);
        end
        if strcmp(this_action, 'take value at band end')
            data = calculateBandEndVariables(this_settings_table, this_settings_table_header, study_settings, data);
        end
        if strcmp(this_action, 'take value at point given by percentage within band')
            data = calculateBandPercentVariables(this_settings_table, this_settings_table_header, study_settings, data);
        end
        if strcmp(this_action, 'take value at point given by absolute time within band')
            data = calculateTimePointVariables(this_settings_table, this_settings_table_header, study_settings, data);
        end
        if strcmp(this_action, 'take value at point given by absolute time across range')
            data = calculateTimePointRangeVariables(this_settings_table, this_settings_table_header, study_settings, data);
        end
        if strcmp(this_action, 'take extremum within whole band')
            data = calculateExtremaVariables(this_settings_table, this_settings_table_header, study_settings, data);
        end
        if strcmp(this_action, 'take extremum over range') 
            data = calculateExtremaOverRangeVariables(this_settings_table, this_settings_table_header, data);
        end
        if strcmp(this_action, 'average across range') 
            data = calculateAverageAcrossRangeVariables(this_settings_table, this_settings_table_header, data);
        end
        if strcmp(this_action, 'combine two variables')
            data = combineTwoVariables(this_settings_table, this_settings_table_header, study_settings, data);
        end
        if strcmp(this_action, 'multiply two variables')
            data = multiplyTwoVariables(this_settings_table, this_settings_table_header, study_settings, data);
        end
        
        
        
        
    end
    
    % save results
    save(results_file_name, '-struct', 'data');    
end

function data = calculateStimulusResponse(response_variables, response_variables_header, study_settings, data, conditions)
    % make condition data tables
    conditions_settings = study_settings.get('conditions');
    condition_labels = conditions_settings(:, 1)';
    condition_source_variables = conditions_settings(:, 2)';
    number_of_condition_labels = length(condition_labels);
    conditions_session = conditions.conditions_session;
    number_of_stretches = size(data.time_list_session, 1); %#ok<*USENS>
    condition_data_all = cell(number_of_stretches, number_of_condition_labels);
    for i_condition = 1 : number_of_condition_labels
        condition_data_all(:, i_condition) = conditions_session.(condition_source_variables{i_condition});
    end
    labels_to_ignore = study_settings.get('conditions_to_ignore');
    levels_to_remove = study_settings.get('levels_to_remove', 1);
    comparisons = struct;
    [comparisons.combination_labels, comparisons.condition_combinations_control] ...
        = determineConditionCombinations ...
            (condition_data_all, conditions_settings, labels_to_ignore, levels_to_remove, 'control');
    comparisons.condition_combinations_control_unique ...
        = table2cell(unique(cell2table(comparisons.condition_combinations_control), 'rows'));

    %% calculate response (i.e. difference from control mean)
    if ~isempty(comparisons.condition_combinations_control)
        % prepare containers
        number_of_response_variables = size(response_variables, 1);
        response_data = cell(number_of_response_variables, 1);
        response_directions = cell(number_of_response_variables, 2);
        response_names = cell(number_of_response_variables, 1);
        for i_variable = 1 : number_of_response_variables
            this_variable_source_name = response_variables{i_variable, strcmp(response_variables_header, 'source_variable_name')};
            this_variable_source_type = response_variables{i_variable, strcmp(response_variables_header, 'source_type')};
            
            % pick data depending on source specification
            this_variable_source_index = find(strcmp(data.([this_variable_source_type '_names_session']), this_variable_source_name), 1, 'first');
            this_variable_source_data = data.([this_variable_source_type '_data_session']){this_variable_source_index};
            
            response_data{i_variable} = zeros(size(this_variable_source_data));
        end        
        
        % go stretch by stretch
        for i_stretch = 1 : number_of_stretches
            % extract this stretches relevant conditions
            this_stretch_condition_string = cell(1, length(comparisons.combination_labels));
            for i_label = 1 : length(comparisons.combination_labels)
                this_label = comparisons.combination_labels{i_label};
                this_label_source_variable = condition_source_variables{strcmp(condition_labels, this_label)};
                this_label_condition_list = conditions_session.(this_label_source_variable);
                this_stretch_condition_string{i_label} = this_label_condition_list{i_stretch};
            end
            
            % determine applicable control condition index
            applicable_control_condition_index ...
                = determineControlConditionIndex(study_settings, comparisons, this_stretch_condition_string);
            
            % determine indicator for control
            control_condition_indicator = true(number_of_stretches, 1);
            for i_label = 1 : length(comparisons.combination_labels)
                this_label = comparisons.combination_labels{i_label};
                this_label_list = condition_data_all(:, strcmp(conditions_settings(:, 1), this_label));
                this_label_control ...
                    = comparisons.condition_combinations_control_unique(applicable_control_condition_index, i_label);
                this_label_indicator = strcmp(this_label_list, this_label_control);
                control_condition_indicator = control_condition_indicator .* this_label_indicator;
            end        
            control_condition_indicator = logical(control_condition_indicator);            
            
            % calculate responses
            for i_variable = 1 : number_of_response_variables
                this_variable_name = response_variables{i_variable, strcmp(response_variables_header, 'new_variable_name')};
                this_variable_source_name = response_variables{i_variable, strcmp(response_variables_header, 'source_variable_name')};
                this_variable_source_type = response_variables{i_variable, strcmp(response_variables_header, 'source_type')};

                % pick data depending on source specification
                this_variable_source_index = find(strcmp(data.([this_variable_source_type '_names_session']), this_variable_source_name), 1, 'first');
                this_variable_source_data = data.([this_variable_source_type '_data_session']){this_variable_source_index};
                this_variable_source_directions = data.([this_variable_source_type '_directions_session'])(this_variable_source_index, :);
                
                % calculate control mean
                this_condition_control_data = this_variable_source_data(:, control_condition_indicator);
                this_condition_control_mean = mean(this_condition_control_data, 2);
                
                % calculate response
                response_data{i_variable}(:, i_stretch) ...
                    = this_variable_source_data(:, i_stretch) - this_condition_control_mean;
                response_directions(i_variable, :) = this_variable_source_directions;
                response_names{i_variable} = this_variable_name;
            end
            
        end

    end

    % store
    for i_variable = 1 : number_of_response_variables
        new_data = struct;
        new_data.data = response_data{i_variable};
        new_data.directions = response_directions(i_variable, :);
        new_data.name = response_names{i_variable};
        data = addOrReplaceResultsData(data, new_data, 'response');
    end
end



function data = calculateInversionVariables(inversion_variables, inversion_variables_header, study_settings, data, conditions)
    for i_variable = 1 : size(inversion_variables, 1)
        % get data
        this_variable_name = inversion_variables{i_variable, strcmp(inversion_variables_header, 'new_variable_name')};
        this_variable_source_name = inversion_variables{i_variable, strcmp(inversion_variables_header, 'source_variable_name')};
        this_variable_source_type = inversion_variables{i_variable, strcmp(inversion_variables_header, 'source_type')};
        relevant_condition = inversion_variables{i_variable, strcmp(inversion_variables_header, 'relevant_condition')};
        inversion_table = study_settings.get(inversion_variables{i_variable, strcmp(inversion_variables_header, 'information_table')});
        new_variable_direction_pos = inversion_variables(i_variable, strcmp(inversion_variables_header, 'direction_label_positive'));
        new_variable_direction_neg = inversion_variables(i_variable, strcmp(inversion_variables_header, 'direction_label_negative'));
        new_variable_directions = [new_variable_direction_pos, new_variable_direction_neg];
        
        % pick data depending on source specification
        this_variable_source_index = find(strcmp(data.([this_variable_source_type '_names_session']), this_variable_source_name), 1, 'first');
        this_variable_source_data = data.([this_variable_source_type '_data_session']){this_variable_source_index};
        this_variable_source_directions = data.([this_variable_source_type '_directions_session'])(this_variable_source_index, :);
        
        % go through levels and invert
        this_variable_data = this_variable_source_data;
        level_list = conditions.conditions_session.(conditions.source_variables{strcmp(conditions.factor_labels, relevant_condition)});
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
        new_data = struct;
        new_data.data = this_variable_data;
        new_data.directions = new_variable_directions;
        new_data.name = this_variable_name;
        
        if strcmp(this_variable_source_type, 'range')
            data = addOrReplaceResultsData(data, new_data, 'range');
        else
            data = addOrReplaceResultsData(data, new_data, 'analysis');
        end
    end
end

function data = calculateSelectionVariables(selection_variables, selection_variables_header, study_settings, data, conditions)
    for i_variable = 1 : size(selection_variables, 1)
        % get data
        this_variable_name = selection_variables{i_variable, strcmp(selection_variables_header, 'new_variable_name')};
        this_variable_source_type = selection_variables{i_variable, strcmp(selection_variables_header, 'source_type')};
        this_variable_relevant_condition = selection_variables{i_variable, strcmp(selection_variables_header, 'relevant_condition')};
        this_variable_information_table = selection_variables{i_variable, strcmp(selection_variables_header, 'information_table')};
        % pick data depending on source specification
        data_source = data.([this_variable_source_type '_data_session']);
        names_source = data.([this_variable_source_type '_names_session']);
        directions_source = data.([this_variable_source_type '_directions_session']);
        selection_table = study_settings.get(this_variable_information_table);
        
        % go through levels and select
        this_variable_data = [];
        new_variable_directions = [];
        level_list = conditions.conditions_session.(conditions.source_variables{strcmp(conditions.factor_labels, this_variable_relevant_condition)});
        for i_level = 1 : size(selection_table, 1)
            % extract level info
            label_this_level = selection_table(i_level, 1);
            data_source_this_level = selection_table(i_level, 2);
            % get data
            source_data_this_level = data_source{strcmp(names_source, data_source_this_level)};
            source_directions_this_level = directions_source(strcmp(names_source, data_source_this_level), :);
            % extract data
            match_this_level = strcmp(level_list, label_this_level);
            if isempty(this_variable_data)
                this_variable_data = source_data_this_level;
                this_variable_data(:, ~match_this_level) = NaN;
            else
                this_variable_data(:, match_this_level) = source_data_this_level(:, match_this_level); %#ok<AGROW>
            end
            if isempty(new_variable_directions)
                new_variable_directions = source_directions_this_level;
            else
                if any(~(strcmp(new_variable_directions, source_directions_this_level)))
                    error('Trying to combine variables with different directions')
                end
            end
        end

        % store
        new_data = struct;
        new_data.data = this_variable_data;
        new_data.directions = new_variable_directions;
        new_data.name = this_variable_name;
        data = addOrReplaceResultsData(data, new_data, 'analysis');
    end    

end

function data = calculateIntegratedVariables(variables_to_integrate, variables_to_integrate_header, study_settings, data)
    step_time_index_in_saved_data = find(strcmp(data.stretch_names_session, 'step_time'), 1, 'first');
    this_step_time_data = data.stretch_data_session{step_time_index_in_saved_data};
    number_of_stretches = size(data.stretch_data_session{1}, 2);
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');

    for i_variable = 1 : size(variables_to_integrate, 1)
        this_variable_name = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'new_variable_name')};
        this_variable_source_name = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'source_variable_name')};
        this_variable_source_type = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'source_variable_type')};
        start_info = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'start')};
        start_variable_source_type = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'start_variable_type')};
        end_info = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'end')};
        end_variable_source_type = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'end_variable_type')};
        % pick data depending on source specification
        data_source = data.([this_variable_source_type '_data_session']);
        names_source = data.([this_variable_source_type '_names_session']);
        directions_source = data.([this_variable_source_type '_directions_session']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        this_variable_source_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        
        % determine start and end of integration in percent of band time
        if strcmp(start_variable_source_type, 'percentage')
            start_data_percent = ones(size(this_step_time_data)) * str2double(start_info);
        else
            start_data_source = data.([start_variable_source_type '_data_session']);
            start_names_source = data.([start_variable_source_type '_names_session']);
            
            start_data_time_within_band = start_data_source{strcmp(start_names_source, start_info)};
            start_data_ratio = start_data_time_within_band ./ this_step_time_data;
            start_data_percent = round(start_data_ratio * 100);
			% HR (02-OCT-2019) -- this is a dirty hack. Leaving it in because it's for TF's data. Do not do things in this way in the future!
%             if strcmp(subject_id,'OLR02')
%                 start_data_percent(4,128) = 30;
%             end
%             if strcmp(subject_id,'OLR03') 
%                 start_data_percent(4,10) = 30;
%             end
%             if strcmp(subject_id,'OLR05') 
%                 start_data_percent(3,162) = 30;
%             end
%             if strcmp(subject_id,'OLR06') 
%                 start_data_percent(4,79) = 30;
%             end
%             if strcmp(subject_id,'OLR06') 
%                 end_data_percent(4,108) = 30;
%             end
%             if strcmp(subject_id,'OLR06') 
%                 start_data_percent(4,133) = 30;
%             end
%             if strcmp(subject_id,'OLR08') 
%                 start_data_percent(4,19) = 30;
%             end
%             
%             
%             if strcmp(subject_id, 'CAD09')
%                 start_data_percent(1,202) = 30;
%             end
%             if strcmp(subject_id, 'CAD15')
%                 start_data_percent(1,188) = 30;
%             end
%             if strcmp(subject_id, 'CAD16')
%                 start_data_percent(2,132) = 30;
%             end
%             if strcmp(subject_id, 'CAD19')
%                 start_data_percent(2,48) = 30;
%             end
%             
%             if strcmp(subject_id, 'CAD21')
%                 start_data_percent(1,97) = 30;
%                 start_data_percent(1,27) = 30;
%             end
        end
        if strcmp(end_variable_source_type, 'percentage')
            end_data_percent = ones(size(this_step_time_data)) * str2double(end_info);
        else
            end_data_source = data.([end_variable_source_type '_data_session']);
            end_names_source = data.([end_variable_source_type '_names_session']);
            end_data_time_within_band = end_data_source{strcmp(end_names_source, end_info)};
            end_data_ratio = end_data_time_within_band ./ this_step_time_data;
            end_data_percent = round(end_data_ratio * 100);
			% HR (02-OCT-2019) -- this is a dirty hack. Leaving it in because it's for TF's data. Do not do things in this way in the future!
%             if strcmp(subject_id,'OLR02') 
%                 end_data_percent(4,128) = 30;
%             end
%             if strcmp(subject_id,'OLR03') 
%                 end_data_percent(4,10) = 30;
%             end
%             if strcmp(subject_id,'OLR05') 
%                 end_data_percent(3,162) = 30;
%             end
%             if strcmp(subject_id,'OLR06') 
%                 end_data_percent(4,79) = 30;
%             end
%             if strcmp(subject_id,'OLR06') 
%                 end_data_percent(4,27) = 30;
%             end
%             if strcmp(subject_id,'OLR06') 
%                 end_data_percent(4,133) = 30;
%             end
%             if strcmp(subject_id,'OLR08') 
%                 end_data_percent(4,19) = 30;
%             end
%             
%             
%             if strcmp(subject_id, 'CAD09')
%                 end_data_percent(1,202) = 30;
%             end
%             if strcmp(subject_id, 'CAD15')
%                 end_data_percent(1,188) = 30;
%             end
%             if strcmp(subject_id, 'CAD16')
%                 end_data_percent(2,132) = 30;
%             end
%             if strcmp(subject_id, 'CAD19')
%                 end_data_percent(2,48) = 30;
%             end
%             if strcmp(subject_id, 'CAD21')
%                 end_data_percent(1,97) = 30;
%                 end_data_percent(1,27) = 30;
%             end
        end
        
        % integrate
        integrated_data = zeros(data.bands_per_stretch, number_of_stretches);
        for i_stretch = 1 : number_of_stretches
            for i_band = 1 : data.bands_per_stretch
                [band_start_index, band_end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
                this_band_time_full = linspace(0, this_step_time_data(i_band), number_of_time_steps_normalized);
                this_band_data_full = this_variable_source_data(band_start_index : band_end_index, i_stretch);
                range = (start_data_percent(i_band, i_stretch) : end_data_percent(i_band, i_stretch)) + 1;
                if isempty(range)
                    integrated_data(i_band, i_stretch) = 0;
                else
                    this_band_time_range = this_band_time_full(range);
                    this_band_data_range = this_band_data_full(range);
                    % integrate
                    this_band_data_integrated = cumtrapz(this_band_time_range, this_band_data_range);           
                    integrated_data(i_band, i_stretch) = this_band_data_integrated(end);
                end

            end
        end        
        
        % store
        new_data = struct;
        new_data.data = integrated_data;
        new_data.directions = this_variable_source_directions;
        new_data.name = this_variable_name;
        data = addOrReplaceResultsData(data, new_data, 'analysis');

    end

end

function data = integrateOverRange(variables_to_integrate, variables_to_integrate_header, study_settings, data)
    step_time_index_in_saved_data = find(strcmp(data.stretch_names_session, 'step_time'), 1, 'first');
    this_step_time_data = data.stretch_data_session{step_time_index_in_saved_data};
    number_of_stretches = size(data.stretch_data_session{1}, 2);
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');

    for i_variable = 1 : size(variables_to_integrate, 1)
        this_variable_name = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'new_variable_name')};
        this_variable_source_name = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'source_variable_name')};
        this_variable_source_type = variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'source_variable_type')};
        
        % pick data depending on source specification
        data_source = data.([this_variable_source_type '_data_session']);
        names_source = data.([this_variable_source_type '_names_session']);
        directions_source = data.([this_variable_source_type '_directions_session']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        this_variable_source_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        
        % determine integration range
        start_data_percent = str2num(variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'start_percent')});
        end_data_percent = str2num(variables_to_integrate{i_variable, strcmp(variables_to_integrate_header, 'end_percent')});
        start_index = round(start_data_percent/100 * (size(this_variable_source_data, 1) - 1)) + 1;
        end_index = round(end_data_percent/100 * (size(this_variable_source_data, 1) - 1)) + 1;
        
        % integrate
        integrated_data = zeros(1, number_of_stretches);
        for i_stretch = 1 : number_of_stretches
            % create time vector
            time_stretch = zeros(size(this_variable_source_data, 1), 1);
            for i_band = 1 : data.bands_per_stretch
                [band_start_index, band_end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
                this_band_time = linspace(sum(this_step_time_data(1:i_band-1)), sum(this_step_time_data(1:i_band)), number_of_time_steps_normalized);
                time_stretch(band_start_index : band_end_index) = this_band_time;
            end
            
            % integrate
            this_stretch_data = this_variable_source_data(start_index:end_index, i_stretch);
            integrated_data(i_stretch) = trapz(time_stretch(start_index:end_index), this_stretch_data);
        end        
        
        % store
        new_data = struct;
        new_data.data = integrated_data;
        new_data.directions = this_variable_source_directions;
        new_data.name = this_variable_name;
        data = addOrReplaceResultsData(data, new_data, 'range');

    end

end


function data = calculateMeanVariables(variables_mean, variables_mean_header, study_settings, data)
    step_time_index_in_saved_data = find(strcmp(data.stretch_names_session, 'step_time'), 1, 'first');
    this_step_time_data = data.stretch_data_session{step_time_index_in_saved_data};
    number_of_stretches = size(data.stretch_data_session{1}, 2);
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    for i_variable = 1 : size(variables_mean, 1)
        this_variable_name = variables_mean{i_variable, strcmp(variables_mean_header, 'new_variable_name')};
        this_variable_source_name = variables_mean{i_variable, strcmp(variables_mean_header, 'source_variable_name')};
        this_variable_source_type = variables_mean{i_variable, strcmp(variables_mean_header, 'source_variable_type')};
        start_info = variables_mean{i_variable, strcmp(variables_mean_header, 'start')};
        start_variable_source_type = variables_mean{i_variable, strcmp(variables_mean_header, 'start_variable_type')};
        end_info = variables_mean{i_variable, strcmp(variables_mean_header, 'end')};
        end_variable_source_type = variables_mean{i_variable, strcmp(variables_mean_header, 'end_variable_type')};

        % pick data depending on source specification
        data_source = data.([this_variable_source_type '_data_session']);
        names_source = data.([this_variable_source_type '_names_session']);
        directions_source = data.([this_variable_source_type '_directions_session']);
        
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        this_variable_source_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        
        % determine start and end of integration in percent of band time
        if strcmp(start_variable_source_type, 'percentage')
            start_data_percent = ones(size(this_step_time_data)) * str2double(start_info);
        else
            start_data_source = data.([start_variable_source_type '_data_session']);
            start_names_source = data.([start_variable_source_type '_names_session']);
            start_data_time_within_band = start_data_source{strcmp(start_names_source, start_info)};
            start_data_ratio = start_data_time_within_band ./ this_step_time_data;
            start_data_percent = round(start_data_ratio * 100);
        end
        if strcmp(end_variable_source_type, 'percentage')
            end_data_percent = ones(size(this_step_time_data)) * str2double(end_info);
        else
            end_data_source = data.([end_variable_source_type '_data_session']);
            end_names_source = data.([end_variable_source_type '_names_session']);
            end_data_time_within_band = end_data_source{strcmp(end_names_source, end_info)};
            end_data_ratio = end_data_time_within_band ./ this_step_time_data;
            end_data_percent = round(end_data_ratio * 100);
        end
        
        % mean
        mean_data = zeros(data.bands_per_stretch, number_of_stretches);
        for i_stretch = 1 : number_of_stretches
            for i_band = 1 : data.bands_per_stretch
                [band_start_index, band_end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
                this_band_data_full = this_variable_source_data(band_start_index : band_end_index, i_stretch);
                range = (start_data_percent(i_band, i_stretch) : end_data_percent(i_band, i_stretch)) + 1;
                if isempty(range)
                    mean_data(i_band, i_stretch) = 0;
                else
                    this_band_data_range = this_band_data_full(range);
                    % mean
                     this_band_data_mean = mean(this_band_data_range);           
                     mean_data(i_band, i_stretch) = this_band_data_mean(end);
                end

            end
        end        
        
        % store
        new_data = struct;
        new_data.data = mean_data;
        new_data.directions = this_variable_source_directions;
        new_data.name = this_variable_name;
        data = addOrReplaceResultsData(data, new_data, 'analysis');
          
    end
end

function data = calculateRmsVariables(variables_rms, variables_rms_header, study_settings, data)
    step_time_index_in_saved_data = find(strcmp(data.stretch_names_session, 'step_time'), 1, 'first');
    this_step_time_data = data.stretch_data_session{step_time_index_in_saved_data};
    number_of_stretches = size(data.stretch_data_session{1}, 2);
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    for i_variable = 1 : size(variables_rms, 1)
        this_variable_name = variables_rms{i_variable, strcmp(variables_rms_header, 'new_variable_name')};
        this_variable_source_name = variables_rms{i_variable, strcmp(variables_rms_header, 'source_variable_name')};
        this_variable_source_type = variables_rms{i_variable, strcmp(variables_rms_header, 'source_variable_type')};
        start_info = variables_rms{i_variable, strcmp(variables_rms_header, 'start')};
        start_variable_source_type = variables_rms{i_variable, strcmp(variables_rms_header, 'start_variable_type')};
        end_info = variables_rms{i_variable, strcmp(variables_rms_header, 'end')};
        end_variable_source_type = variables_rms{i_variable, strcmp(variables_rms_header, 'end_variable_type')};

        % pick data depending on source specification
        data_source = data.([this_variable_source_type '_data_session']);
        names_source = data.([this_variable_source_type '_names_session']);
        directions_source = data.([this_variable_source_type '_directions_session']);
        
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        this_variable_source_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        
        % determine start and end of integration in percent of band time
        if strcmp(start_variable_source_type, 'percentage')
            start_data_percent = ones(size(this_step_time_data)) * str2double(start_info);
        else
            start_data_source = data.([start_variable_source_type '_data_session']);
            start_names_source = data.([start_variable_source_type '_names_session']);
            start_data_time_within_band = start_data_source{strcmp(start_names_source, start_info)};
            start_data_ratio = start_data_time_within_band ./ this_step_time_data;
            start_data_percent = round(start_data_ratio * 100);
        end
        if strcmp(end_variable_source_type, 'percentage')
            end_data_percent = ones(size(this_step_time_data)) * str2double(end_info);
        else
            end_data_source = data.([end_variable_source_type '_data_session']);
            end_names_source = data.([end_variable_source_type '_names_session']);
            end_data_time_within_band = end_data_source{strcmp(end_names_source, end_info)};
            end_data_ratio = end_data_time_within_band ./ this_step_time_data;
            end_data_percent = round(end_data_ratio * 100);
        end
        
        % rms
        rms_data = zeros(data.bands_per_stretch, number_of_stretches);
        for i_stretch = 1 : number_of_stretches
            for i_band = 1 : data.bands_per_stretch
                [band_start_index, band_end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
                this_band_data_full = this_variable_source_data(band_start_index : band_end_index, i_stretch);
                range = (start_data_percent(i_band, i_stretch) : end_data_percent(i_band, i_stretch)) + 1;
                if isempty(range)
                    rms_data(i_band, i_stretch) = 0;
                else
                    this_band_data_range = this_band_data_full(range);
                    % rms
                     this_band_data_rms = rms(this_band_data_range);           
                     rms_data(i_band, i_stretch) = this_band_data_rms(end);
                end

            end
        end        
        
        % store
        new_data = struct;
        new_data.data = rms_data;
        new_data.directions = this_variable_source_directions;
        new_data.name = this_variable_name;
        data = addOrReplaceResultsData(data, new_data, 'analysis');
          
    end
end

function data = calculateBandEndVariables(variables_step_end, variables_step_end_header, study_settings, data)
    number_of_stretches = size(data.stretch_data_session{1}, 2);
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    for i_variable = 1 : size(variables_step_end, 1)
        this_variable_name = variables_step_end{i_variable, strcmp(variables_step_end_header, 'new_variable_name')};
        this_variable_source_name = variables_step_end{i_variable, strcmp(variables_step_end_header, 'source_variable_name')};
        this_variable_source_type = variables_step_end{i_variable, strcmp(variables_step_end_header, 'source_type')};
        
        % pick data depending on source specification
        data_source = data.([this_variable_source_type '_data_session']);
        names_source = data.([this_variable_source_type '_names_session']);
        directions_source = data.([this_variable_source_type '_directions_session']);
        
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        new_variable_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        step_end_data = zeros(data.bands_per_stretch, number_of_stretches);
        for i_band = 1 : data.bands_per_stretch
            [~, end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
            step_end_data(i_band, :) = this_variable_source_data(end_index, :);
        end
        
        % store
        new_data = struct;
        new_data.data = step_end_data;
        new_data.directions = new_variable_directions;
        new_data.name = this_variable_name;
        data = addOrReplaceResultsData(data, new_data, 'analysis');
    end

end

function data = calculateBandPercentVariables(variables_band_percent, variables_band_percent_header, study_settings, data)
    number_of_stretches = size(data.stretch_data_session{1}, 2);
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    for i_variable = 1 : size(variables_band_percent, 1)
        this_variable_name = variables_band_percent{i_variable, strcmp(variables_band_percent_header, 'new_variable_name')};
        this_variable_source_name = variables_band_percent{i_variable, strcmp(variables_band_percent_header, 'source_variable_name')};
        this_variable_source_type = variables_band_percent{i_variable, strcmp(variables_band_percent_header, 'source_variable_type')};
        this_variable_source_percent = str2double(variables_band_percent{i_variable, strcmp(variables_band_percent_header, 'percentage')});

        % pick data depending on source specification
        data_source = data.([this_variable_source_type '_data_session']);
        names_source = data.([this_variable_source_type '_names_session']);
        directions_source = data.([this_variable_source_type '_directions_session']);
        
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        new_variable_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        band_percent_data = zeros(data.bands_per_stretch, number_of_stretches);
        for i_band = 1 : data.bands_per_stretch
            [start_index, end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
            this_index = start_index + this_variable_source_percent/100 * (end_index - start_index);
            band_percent_data(i_band, :) = this_variable_source_data(this_index, :);
        end
        
        % store
        new_data = struct;
        new_data.data = band_percent_data;
        new_data.directions = new_variable_directions;
        new_data.name = this_variable_name;
        data = addOrReplaceResultsData(data, new_data, 'analysis');
    end

end

function data = calculateTimePointVariables(variables_time_point, variables_time_point_header, study_settings, data)
    step_time_index_in_saved_data = find(strcmp(data.stretch_names_session, 'step_time'), 1, 'first');
    this_step_time_data = data.stretch_data_session{step_time_index_in_saved_data};
    number_of_stretches = size(data.stretch_data_session{1}, 2);
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    for i_variable = 1 : size(variables_time_point, 1)
        this_variable_name = variables_time_point{i_variable, strcmp(variables_time_point_header, 'new_variable_name')};
        this_variable_source_name = variables_time_point{i_variable, strcmp(variables_time_point_header, 'source_variable_name')};
        this_variable_source_type = variables_time_point{i_variable, strcmp(variables_time_point_header, 'source_variable_type')};
        time_point = variables_time_point{i_variable, strcmp(variables_time_point_header, 'time_point')};
        time_point_source_type = variables_time_point{i_variable, strcmp(variables_time_point_header, 'time_point_specifier')};
        % pick data depending on source specification
        data_source = data.([this_variable_source_type '_data_session']);
        names_source = data.([this_variable_source_type '_names_session']);
        directions_source = data.([this_variable_source_type '_directions_session']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        new_variable_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        
        % determine time point in percent of band time
        if strcmp(time_point_source_type, 'percentage')
            this_time_point_percent = ones(size(this_step_time_data)) * str2double(time_point);
        else
            time_point_data_source = data.([time_point_source_type '_data_session']);
            time_point_names_source = data.([time_point_source_type '_names_session']);
            
            time_point_within_band = time_point_data_source{strcmp(time_point_names_source, time_point)};
            time_point_ratio = time_point_within_band ./ this_step_time_data;
            this_time_point_percent = round(time_point_ratio * 100);
        end
        
        % get the data at this time point
        time_point_data = zeros(data.bands_per_stretch, number_of_stretches);
        for i_stretch = 1 : number_of_stretches
            for i_band = 1 : data.bands_per_stretch
                [start_index, end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
                this_index = round(start_index + this_time_point_percent(i_band, i_stretch)/100 * (end_index - start_index));
                time_point_data(i_band, i_stretch) = this_variable_source_data(this_index, i_stretch);
            end
        end
        
        % store
        new_data = struct;
        new_data.data = time_point_data;
        new_data.directions = new_variable_directions;
        new_data.name = this_variable_name;
        data = addOrReplaceResultsData(data, new_data, 'analysis');
    end

end

function data = calculateTimePointRangeVariables(variables_time_point, variables_time_point_header, study_settings, data)
    number_of_stretches = size(data.stretch_data_session{1}, 2);
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    for i_variable = 1 : size(variables_time_point, 1)
        this_variable_name = variables_time_point{i_variable, strcmp(variables_time_point_header, 'new_variable_name')};
        this_variable_source_name = variables_time_point{i_variable, strcmp(variables_time_point_header, 'source_variable_name')};
        this_variable_source_type = variables_time_point{i_variable, strcmp(variables_time_point_header, 'source_variable_type')};
        time_point = variables_time_point{i_variable, strcmp(variables_time_point_header, 'time_point')};
        time_point_source_type = variables_time_point{i_variable, strcmp(variables_time_point_header, 'time_point_specifier')};
        % pick data depending on source specification
        data_source = data.([this_variable_source_type '_data_session']);
        names_source = data.([this_variable_source_type '_names_session']);
        directions_source = data.([this_variable_source_type '_directions_session']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        new_variable_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        
        % determine time point in percent of band time
        if strcmp(time_point_source_type, 'range')
            time_point_data_source = data.([time_point_source_type '_data_session']);
            time_point_names_source = data.([time_point_source_type '_names_session']);
            time_point_data_index = strcmp(time_point_names_source, time_point);
            time_point_data = time_point_data_source{time_point_data_index};
        else
            error(['Unknown time point source type "' time_point_source_type '"'])
        end
        
        % get the data at this time point
        new_variable_data = zeros(1, number_of_stretches);
        for i_stretch = 1 : number_of_stretches
            % recreate time vector for normalized stretch data from stretch_times
            this_stretch_times = data.stretch_times(:, i_stretch);
            this_stretch_time = zeros((number_of_time_steps_normalized-1) * data.bands_per_stretch + 1, 1);
            for i_band = 1 : data.bands_per_stretch
                % make time vector for this band
                this_band_start_time = this_stretch_times(i_band);
                this_band_end_time = this_stretch_times(i_band+1);
                this_band_time = linspace(this_band_start_time, this_band_end_time, number_of_time_steps_normalized);
                
                % store time vector for this band at proper location
                [start_index, end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
                this_stretch_time(start_index : end_index) = this_band_time;
            end
            
            % find index for specified time point
            this_time_point = time_point_data(i_stretch);
            specified_time_point_index = findClosestIndex(this_time_point, this_stretch_time);
            
            % extract data
            new_variable_data(i_stretch) = this_variable_source_data(specified_time_point_index, i_stretch);
        end
        
        % store
        new_data = struct;
        new_data.data = new_variable_data;
        new_data.directions = new_variable_directions;
        new_data.name = this_variable_name;
        data = addOrReplaceResultsData(data, new_data, 'range');
    end

end

function data = calculateExtremaVariables(variables_from_extrema, variables_from_extrema_header, study_settings, data)
    number_of_stretches = size(data.stretch_data_session{1}, 2);
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    for i_variable = 1 : size(variables_from_extrema, 1)
        % get data
        this_variable_name = variables_from_extrema{i_variable, strcmp(variables_from_extrema_header, 'new_variable_name')};
        this_variable_source_name = variables_from_extrema{i_variable, strcmp(variables_from_extrema_header, 'source_variable_name')};
        this_variable_source_type = variables_from_extrema{i_variable, strcmp(variables_from_extrema_header, 'source_type')};
        this_variable_extremum_type = variables_from_extrema{i_variable, strcmp(variables_from_extrema_header, 'extremum_type')};

        % pick data depending on source specification
        data_source = data.([this_variable_source_type '_data_session']);
        names_source = data.([this_variable_source_type '_names_session']);
        directions_source = data.([this_variable_source_type '_directions_session']);
        
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        new_variable_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        
        extrema_data = zeros(data.bands_per_stretch, number_of_stretches);
        for i_band = 1 : data.bands_per_stretch
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
        new_data = struct;
        new_data.data = extrema_data;
        new_data.directions = new_variable_directions;
        new_data.name = this_variable_name;
        data = addOrReplaceResultsData(data, new_data, 'analysis');
    end
end

function data = calculateExtremaOverRangeVariables(variables_from_extrema_range, variables_from_extrema_range_header, data)
    for i_variable = 1 : size(variables_from_extrema_range, 1)
        % get data
        this_variable_name = variables_from_extrema_range{i_variable, strcmp(variables_from_extrema_range_header, 'new_variable_name')};
        this_variable_source_name = variables_from_extrema_range{i_variable, strcmp(variables_from_extrema_range_header, 'source_variable_name')};
        this_variable_source_type = variables_from_extrema_range{i_variable, strcmp(variables_from_extrema_range_header, 'source_type')};
        this_variable_extremum_type = variables_from_extrema_range{i_variable, strcmp(variables_from_extrema_range_header, 'extremum_type')};

        % pick data depending on source specification
        data_source = data.([this_variable_source_type '_data_session']);
        names_source = data.([this_variable_source_type '_names_session']);
        directions_source = data.([this_variable_source_type '_directions_session']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        new_variable_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        
        if strcmp(this_variable_extremum_type, 'min')
            extrema_data = min(this_variable_source_data);
        end
        if strcmp(this_variable_extremum_type, 'max')
            extrema_data = max(this_variable_source_data);
        end
        if ~strcmp(this_variable_extremum_type, 'min') && ~strcmp(this_variable_extremum_type, 'max')
            error(['"' this_variable_extremum_type '" is not a valid type for variables_from_extrema. Acceptable types are "min" or "max".']);
        end

        % store
        new_data = struct;
        new_data.data = extrema_data;
        new_data.directions = new_variable_directions;
        new_data.name = this_variable_name;
        data = addOrReplaceResultsData(data, new_data, 'range');
    end
end


function data = calculateAverageAcrossRangeVariables(variables_from_average_across_range, variables_from_average_across_range_header, data)
    for i_variable = 1 : size(variables_from_average_across_range, 1)
        % get data
        this_variable_name = variables_from_average_across_range{i_variable, strcmp(variables_from_average_across_range_header, 'new_variable_name')};
        this_variable_source_name = variables_from_average_across_range{i_variable, strcmp(variables_from_average_across_range_header, 'source_variable_name')};
        this_variable_source_type = variables_from_average_across_range{i_variable, strcmp(variables_from_average_across_range_header, 'source_type')};

        % pick data depending on source specification
        data_source = data.([this_variable_source_type '_data_session']);
        names_source = data.([this_variable_source_type '_names_session']);
        directions_source = data.([this_variable_source_type '_directions_session']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        new_variable_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        
        mean_data = mean(this_variable_source_data);

        % store
        new_data = struct;
        new_data.data = mean_data;
        new_data.directions = new_variable_directions;
        new_data.name = this_variable_name;
        data = addOrReplaceResultsData(data, new_data, 'range');
    end
end


function data = combineTwoVariables(variable_table, table_header, study_settings, data)
    for i_variable = 1 : size(variable_table, 1)
        % get data
        this_variable_name = variable_table{i_variable, strcmp(table_header, 'new_variable_name')};
        variable_A_name = variable_table{i_variable, strcmp(table_header, 'variable_A_name')};
        variable_A_type = variable_table{i_variable, strcmp(table_header, 'variable_A_type')};
        variable_A_gain = str2double(variable_table{i_variable, strcmp(table_header, 'variable_A_gain')});
        variable_B_name = variable_table{i_variable, strcmp(table_header, 'variable_B_name')};
        variable_B_type = variable_table{i_variable, strcmp(table_header, 'variable_B_type')};
        variable_B_gain = str2double(variable_table{i_variable, strcmp(table_header, 'variable_B_gain')});
        offset = str2double(variable_table{i_variable, strcmp(table_header, 'offset')});

        % pick data depending on source specification
        variable_A_data_source = data.([variable_A_type '_data_session']);
        variable_A_names_source = data.([variable_A_type '_names_session']);
        variable_A_directions_source = data.([variable_A_type '_directions_session']);
        variable_B_data_source = data.([variable_B_type '_data_session']);
        variable_B_names_source = data.([variable_B_type '_names_session']);
        variable_B_directions_source = data.([variable_B_type '_directions_session']);
        
        % extract
        variable_A_data = variable_A_data_source{strcmp(variable_A_names_source, variable_A_name)};
        variable_B_data = variable_B_data_source{strcmp(variable_B_names_source, variable_B_name)};
        variable_A_directions = variable_A_directions_source(strcmp(variable_A_names_source, variable_A_name), :);
        variable_B_directions = variable_B_directions_source(strcmp(variable_B_names_source, variable_B_name), :);

        % compare directions
        if ~strcmp(variable_A_directions{1}, variable_B_directions{1})
            error ...
              ( ...
                [ ...
                  'Positive direction labels "' variable_A_directions{1} '" ' ...
                  'in variable "' variable_A_name '" ' ...
                  'and "' variable_B_directions{1} '" ' ...
                  'in variable "' variable_B_name '" ' ...
                  'do not match' ...
                ] ...
              );
        end
        if ~strcmp(variable_A_directions{2}, variable_B_directions{2})
            error ...
              ( ...
                [ ...
                  'Positive direction labels "' variable_A_directions{2} '" ' ...
                  'in variable "' variable_A_name '" ' ...
                  'and "' variable_B_directions{2} '" ' ...
                  'in variable "' variable_B_name '" ' ...
                  'do not match' ...
                ] ...
              );
        end
        combined_variable_directions = variable_A_directions;
        if ~isequal(size(variable_A_data), size(variable_B_data))
            % one variable is discrete and the other one continuous, deal with this
            if size(variable_A_data, 1) < size(variable_B_data, 1)
                % variable A is the discrete one
                discrete_variable_data = variable_A_data;
                continuous_variable_data = variable_B_data;
            else
                % variable B is the discrete one
                discrete_variable_data = variable_B_data;
                continuous_variable_data = variable_A_data;
            end
            
            number_of_bands = data.bands_per_stretch;
            number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
            discrete_variable_data_extended = zeros(size(continuous_variable_data)) * NaN;
            for i_band = 1 : number_of_bands
                [band_start_index, band_end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
                discrete_variable_data_extended(band_start_index+1 : band_end_index-1, :) = repmat(discrete_variable_data(i_band, :), number_of_time_steps_normalized-2, 1);
            end
            
            if size(variable_A_data, 1) < size(variable_B_data, 1)
                % variable A is the discrete one
                variable_A_data = discrete_variable_data_extended;
            else
                % variable B is the discrete one
                variable_B_data = discrete_variable_data_extended;
            end
            
        end
        combined_variable_data = variable_A_data * variable_A_gain + variable_B_data * variable_B_gain + offset;
        
        % store
        new_data = struct;
        new_data.data = combined_variable_data;
        new_data.directions = combined_variable_directions;
        new_data.name = this_variable_name;
        
        if strcmp(variable_A_type, 'range') && strcmp(variable_A_type, 'range')
            data = addOrReplaceResultsData(data, new_data, 'range');
        else
            data = addOrReplaceResultsData(data, new_data, 'analysis');
        end
        
    end
end

function data = multiplyTwoVariables(variable_table, table_header, study_settings, data)
    for i_variable = 1 : size(variable_table, 1)
        % get data
        this_variable_name = variable_table{i_variable, strcmp(table_header, 'new_variable_name')};
        variable_A_name = variable_table{i_variable, strcmp(table_header, 'variable_A_name')};
        variable_A_type = variable_table{i_variable, strcmp(table_header, 'variable_A_type')};
        variable_A_exponent = str2double(variable_table{i_variable, strcmp(table_header, 'variable_A_exponent')});
        variable_B_name = variable_table{i_variable, strcmp(table_header, 'variable_B_name')};
        variable_B_type = variable_table{i_variable, strcmp(table_header, 'variable_B_type')};
        variable_B_exponent = str2double(variable_table{i_variable, strcmp(table_header, 'variable_B_exponent')});
        factor = str2double(variable_table{i_variable, strcmp(table_header, 'factor')});
        directions_source = variable_table{i_variable, strcmp(table_header, 'directions_source')};

        % pick data depending on source specification
        variable_A_data_source = data.([variable_A_type '_data_session']);
        variable_A_names_source = data.([variable_A_type '_names_session']);
        variable_A_directions_source = data.([variable_A_type '_directions_session']);
        variable_B_data_source = data.([variable_B_type '_data_session']);
        variable_B_names_source = data.([variable_B_type '_names_session']);
        variable_B_directions_source = data.([variable_B_type '_directions_session']);
        
        % extract
        variable_A_data = variable_A_data_source{strcmp(variable_A_names_source, variable_A_name)};
        variable_B_data = variable_B_data_source{strcmp(variable_B_names_source, variable_B_name)};
        variable_A_directions = variable_A_directions_source(strcmp(variable_A_names_source, variable_A_name), :);
        variable_B_directions = variable_B_directions_source(strcmp(variable_B_names_source, variable_B_name), :);

        % get directions
        if strcmp(directions_source, variable_A_name)
            combined_variable_directions = variable_A_directions;            
        elseif strcmp(directions_source, variable_B_name)
            combined_variable_directions = variable_B_directions;            
        else
            error(['Directions for variable "' this_variable_name '" must come from one of the two source variables.']);
        end
                
        if ~isequal(size(variable_A_data), size(variable_B_data))
            % one variable is discrete and the other one continuous, deal with this
            if size(variable_A_data, 1) < size(variable_B_data, 1)
                % variable A is the discrete one
                discrete_variable_data = variable_A_data;
                continuous_variable_data = variable_B_data;
            else
                % variable B is the discrete one
                discrete_variable_data = variable_B_data;
                continuous_variable_data = variable_A_data;
            end
            
            number_of_bands = data.bands_per_stretch;
            number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
            discrete_variable_data_extended = zeros(size(continuous_variable_data)) * NaN;
            for i_band = 1 : number_of_bands
                [band_start_index, band_end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
                discrete_variable_data_extended(band_start_index+1 : band_end_index-1, :) = repmat(discrete_variable_data(i_band, :), number_of_time_steps_normalized-2, 1);
            end
            
            if size(variable_A_data, 1) < size(variable_B_data, 1)
                % variable A is the discrete one
                variable_A_data = discrete_variable_data_extended;
            else
                % variable B is the discrete one
                variable_B_data = discrete_variable_data_extended;
            end
            
        end
        combined_variable_data = variable_A_data.^variable_A_exponent .* variable_B_data.^variable_B_exponent * factor;
        
        % store
        new_data = struct;
        new_data.data = combined_variable_data;
        new_data.directions = combined_variable_directions;
        new_data.name = this_variable_name;
        
        if strcmp(variable_A_type, 'range') && strcmp(variable_A_type, 'range')
            data = addOrReplaceResultsData(data, new_data, 'range');
        else
            data = addOrReplaceResultsData(data, new_data, 'analysis');
        end
        
    end
end

function applicable_control_condition_index = determineControlConditionIndex(study_settings, comparisons, this_stretch_condition_string)
    experimental_paradigm = study_settings.get('experimental_paradigm');
    paradigms_with_intermittent_perturbation = ...
      { ...
        'Vision', 'GVS', 'GVS_old', 'Vision_old', 'CadenceGVS', 'FatigueGVS', 'OculusLaneRestriction', ...
        'CognitiveLoadVision', 'CognitiveLoadGvs', 'GvsOverground' ...
      };
    paradigms_with_stochastic_resonance = {'SR_VisualStim', 'nGVS_Vision'};
  
    % define the factors for which the levels have to match to determine the control condition for this condition
    if any(strcmp(experimental_paradigm, paradigms_with_stochastic_resonance))
        if strcmp(experimental_paradigm, 'SR_VisualStim')
            relevant_factors_for_control = {'stim_amplitude', 'trigger_foot'};
        end
        if strcmp(experimental_paradigm, 'nGVS_Vision')
            relevant_factors_for_control = {'ngvs_settings', 'trigger_foot'};
        end
    end

    if any(strcmp(experimental_paradigm, paradigms_with_intermittent_perturbation))
        relevant_factors_for_control = {'trigger_foot'}; % control only differs by trigger foot in this paradigm
        if strcmp(experimental_paradigm, 'CadenceGVS')
            relevant_factors_for_control = [relevant_factors_for_control, 'cadence'];
        end
        if strcmp(experimental_paradigm, 'FatigueGVS')
            relevant_factors_for_control = [relevant_factors_for_control, 'fatigue'];
        end
        if strcmp(experimental_paradigm, 'OculusLaneRestriction')
            relevant_factors_for_control = [relevant_factors_for_control, 'zone_side'];
        end
        if strcmp(experimental_paradigm, 'CognitiveLoadVision') || strcmp(experimental_paradigm, 'CognitiveLoadGvs')
            relevant_factors_for_control = [relevant_factors_for_control, 'cognitive_load'];
        end
    end
    
    % go through each factor and find the matches in the control
    control_row_indicator = true(size(comparisons.condition_combinations_control_unique, 1), 1);
    for i_factor = 1 : length(relevant_factors_for_control)
        this_factor_label = relevant_factors_for_control{i_factor};
        this_factor_column = strcmp(comparisons.combination_labels, this_factor_label);
        this_factor_this_level = this_stretch_condition_string(this_factor_column);
        this_factor_control_levels = comparisons.condition_combinations_control_unique(:, this_factor_column);
        this_factor_candidate_rows ...
            = strcmp ...
              ( ...
                this_factor_control_levels, ...
                this_factor_this_level ...
              );
        % keep only rows that match previous ones and this one
        control_row_indicator = control_row_indicator & this_factor_candidate_rows;
    end

    % transform the resulting row into an index
    applicable_control_condition_index = find(control_row_indicator);

end




