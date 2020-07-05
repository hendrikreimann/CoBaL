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
        
        if strcmp(this_action, 'integrate over time')
            data = calculateIntegratedVariables(this_settings_table, this_settings_table_header, study_settings, data);
        end
        if strcmp(this_action, 'select from multiple variables by condition')
            data = calculateSelectionVariables(this_settings_table, this_settings_table_header, study_settings, data, conditions);
        end
        if strcmp(this_action, 'invert by condition')
            data = calculateInversionVariables(this_settings_table, this_settings_table_header, study_settings, data, conditions);
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
        if strcmp(this_action, 'take extremum within whole band')
    data = calculateExtremaVariables(study_settings, data);
        end
        if strcmp(this_action, 'take extremum over time interval within band')
    data = calculateExtremaOverRangeVariables(study_settings, data);
        end
    end
    
    
%     % calculate inversion variables
%     data = calculateInversionVariables(study_settings, data, conditions);
%         
%     % calculate selection variables
%     data = calculateSelectionVariables(study_settings, data, conditions);
%     
%     % calculate integrated variables
%     data = calculateIntegratedVariables(study_settings, data);
%     
%     % calculate rms variables
%     data = calculateRmsVariables(study_settings, data);
%     
%     % calculate band variables
%     data = calculateBandEndVariables(study_settings, data);
%     data = calculateBandPercentVariables(study_settings, data);
%     data = calculateTimePointVariables(study_settings, data);
%     
%     % calculate extrema variables
%     data = calculateExtremaVariables(study_settings, data);
%     data = calculateExtremaOverRangeVariables(study_settings, data);
    
    % save results
    save(results_file_name, '-struct', 'data');    
end

function analysis_table = getAnalysisTable(study_settings)
    analysis_table_body = study_settings.get('analysis_table', 1);
    analysis_table_header = study_settings.get('analysis_table_header', 1);
    
    analysis_table = cell2table(analysis_table_body);
    analysis_table.Properties.VariableNames = analysis_table_header;
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
        data = addOrReplaceResultsData(data, new_data, 'analysis');
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

function data = calculateRmsVariables(variables_to_rms, variables_to_rms_header, study_settings, data)
    step_time_index_in_saved_data = find(strcmp(data.stretch_names_session, 'step_time'), 1, 'first');
    this_step_time_data = data.stretch_data_session{step_time_index_in_saved_data};
    number_of_stretches = size(data.stretch_data_session{1}, 2);
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    for i_variable = 1 : size(variables_to_rms, 1)
        this_variable_name = variables_to_rms{i_variable, strcmp(variables_to_rms_header, 'new_variable_name')};
        this_variable_source_name = variables_to_rms{i_variable, strcmp(variables_to_rms_header, 'source_variable_name')};
        this_variable_source_type = variables_to_rms{i_variable, strcmp(variables_to_rms_header, 'source_variable_type')};
        start_info = variables_to_rms{i_variable, strcmp(variables_to_rms_header, 'start')};
        start_variable_source_type = variables_to_rms{i_variable, strcmp(variables_to_rms_header, 'start_variable_type')};
        end_info = variables_to_rms{i_variable, strcmp(variables_to_rms_header, 'end')};
        end_variable_source_type = variables_to_rms{i_variable, strcmp(variables_to_rms_header, 'end_variable_type')};

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
%         this_variable_name = variables_step_end{i_variable, 1};
%         this_variable_source_name = variables_step_end{i_variable, 2};
%         this_variable_source_type = variables_step_end{i_variable, 3};
        
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

function data = calculateExtremaVariables(study_settings, data)
    variables_from_extrema = study_settings.get('analysis_variables_from_extrema', 1);
    number_of_stretches = size(data.stretch_data_session{1}, 2);
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    for i_variable = 1 : size(variables_from_extrema, 1)
        % get data
        this_variable_name = variables_from_extrema{i_variable, 1};
        this_variable_source_name = variables_from_extrema{i_variable, 2};
        this_variable_source_type = variables_from_extrema{i_variable, 3};
        this_variable_extremum_type = variables_from_extrema{i_variable, 4};

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

function data = calculateExtremaOverRangeVariables(study_settings, data)
    variables_from_extrema_range = study_settings.get('analysis_variables_from_extrema_range', 1);
    for i_variable = 1 : size(variables_from_extrema_range, 1)
        % get data
        this_variable_name = variables_from_extrema_range{i_variable, 1};
        this_variable_source_name = variables_from_extrema_range{i_variable, 2};
        this_variable_source_type = variables_from_extrema_range{i_variable, 3};
        this_variable_extremum_type = variables_from_extrema_range{i_variable, 4};

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






