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
    loaded_data = load(results_file_name);
    
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    bands_per_stretch = loaded_data.bands_per_stretch;
%     stretch_names_session = loaded_data.stretch_names_session;
%     stretch_data_session = loaded_data.stretch_data_session;
%     stretch_directions_session = loaded_data.stretch_directions_session;
%     
%     number_of_stretch_variables = length(loaded_data.stretch_names_session);
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
    
    if isfield(loaded_data, 'response_data_session')
        response_data_session = loaded_data.response_data_session;
        response_directions_session = loaded_data.response_directions_session;
        response_names_session = loaded_data.response_names_session;
    end
    if isfield(loaded_data, 'analysis_data_session')
        analysis_data.data = loaded_data.analysis_data_session;
        analysis_data.directions = loaded_data.analysis_directions_session;
        analysis_data.names = loaded_data.analysis_names_session;
    else
        analysis_data.data = {};
        analysis_data.directions = {};
        analysis_data.names = {};
    end
    if isfield(loaded_data, 'range_data_session')
        range_data.data = loaded_data.range_data_session;
        range_data.directions = loaded_data.range_directions_session;
        range_data.names = loaded_data.range_names_session;
    else
        range_data.data = {};
        range_data.directions = {};
        range_data.names = {};
    end
    
    %% load necessary info
    step_time_index_in_saved_data = find(strcmp(loaded_data.stretch_names_session, 'step_time'), 1, 'first');
    this_step_time_data = loaded_data.stretch_data_session{step_time_index_in_saved_data};
    pushoff_time_index_in_saved_data = find(strcmp(loaded_data.stretch_names_session, 'pushoff_time'), 1, 'first');
    if ~isempty(pushoff_time_index_in_saved_data)
        this_pushoff_time_data = loaded_data.stretch_data_session{pushoff_time_index_in_saved_data};
    else
        this_pushoff_time_data = [];
    end
        
    %% calculate inversion variables
    inversion_variables = study_settings.get('inversion_variables', 1);
    for i_variable = 1 : size(inversion_variables, 1)
        % get data
        this_variable_name = inversion_variables{i_variable, 1};
        this_variable_source_name = inversion_variables{i_variable, 2};
        this_variable_source_type = inversion_variables{i_variable, 3};
        
        % pick data depending on source specification
        this_variable_source_index = find(strcmp(loaded_data.([this_variable_source_type '_names_session']), this_variable_source_name), 1, 'first');
        this_variable_source_data = loaded_data.([this_variable_source_type '_data_session']){this_variable_source_index};
        this_variable_source_directions = loaded_data.([this_variable_source_type '_directions_session'])(this_variable_source_index, :);
%         eval(['data_source = loaded_data.' this_variable_source_type '_data_session;']);
%         eval(['names_source = loaded_data.' this_variable_source_type '_names_session;']);
%         eval(['directions_source = loaded_data.' this_variable_source_type '_directions_session;']);
%         this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
%         this_variable_source_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
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
        analysis_data = addOrReplaceResultsData(analysis_data, this_variable_data, this_variable_name, new_variable_directions);
    end
    
    %% calculate selection variables
    selection_variables_header = study_settings.get('selection_variables_header', 1);
    selection_variables = study_settings.get('selection_variables', 1);
    for i_variable = 1 : size(selection_variables, 1)
        % get data
        this_variable_name = selection_variables{i_variable, strcmp(selection_variables_header, 'new_variable_name')};
        this_variable_source_type = selection_variables{i_variable, strcmp(selection_variables_header, 'source_type')};
        this_variable_relevant_condition = selection_variables{i_variable, strcmp(selection_variables_header, 'relevant_condition')};
        this_variable_information_table = selection_variables{i_variable, strcmp(selection_variables_header, 'information_table')};
        % pick data depending on source specification
        eval(['data_source = loaded_data.' this_variable_source_type '_data_session;']);
        eval(['names_source = loaded_data.' this_variable_source_type '_names_session;']);
        eval(['directions_source = loaded_data.' this_variable_source_type '_directions_session;']);
        
        selection_table = study_settings.get(this_variable_information_table);
        
        
        
        
        
        
        
%         this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
%         this_variable_source_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
%         new_variable_directions = selection_variables(i_variable, 6:7);
        
        
        % go through levels and select
        this_variable_data = [];
        new_variable_directions = [];
        level_list = conditions_session.(condition_source_variables{strcmp(condition_labels, this_variable_relevant_condition)});
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
                this_variable_data(:, match_this_level) = source_data_this_level(:, match_this_level);
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
        analysis_data = addOrReplaceResultsData(analysis_data, this_variable_data, this_variable_name, new_variable_directions);
    end    
    
    %% calculate integrated variables
    variables_to_integrate_header = study_settings.get('analysis_variables_from_integration_header', 1);
    variables_to_integrate = study_settings.get('analysis_variables_from_integration', 1);
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
        eval(['data_source = loaded_data.' this_variable_source_type '_data_session;']);
        eval(['names_source = loaded_data.' this_variable_source_type '_names_session;']);
        eval(['directions_source = loaded_data.' this_variable_source_type '_directions_session;']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        this_variable_source_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        
        % determine start and end of integration in percent of band time
        if strcmp(start_variable_source_type, 'percentage')
            start_data_percent = ones(size(this_step_time_data)) * str2num(start_info);
        else
            eval(['start_data_source = loaded_data.' start_variable_source_type '_data_session;']);
            eval(['start_names_source = loaded_data.' start_variable_source_type '_names_session;']);
            start_data_time_within_band = start_data_source{strcmp(start_names_source, start_info)};
            % find a way to take the step time per condition (maybe only
            % for CadenceGVS) and then perform integration on only
            % stretches within each condition
            % HR (09-OCT-2019) -- the above comment was from TF, but I think this is not required
            start_data_ratio = start_data_time_within_band ./ this_step_time_data;
            start_data_percent = round(start_data_ratio * 100);
			% HR (02-OCT-2019) -- this is a dirty hack. Leaving it in because it's for TF's data. Do not do things in this way in the future!
            if strcmp(subject_id,'OLR02')
                start_data_percent(4,128) = 30;
            end
            if strcmp(subject_id,'OLR03') 
                start_data_percent(4,10) = 30;
            end
            if strcmp(subject_id,'OLR05') 
                start_data_percent(3,162) = 30;
            end
            if strcmp(subject_id,'OLR06') 
                start_data_percent(4,79) = 30;
            end
            if strcmp(subject_id,'OLR06') 
                end_data_percent(4,108) = 30;
            end
            if strcmp(subject_id,'OLR06') 
                start_data_percent(4,133) = 30;
            end
            if strcmp(subject_id,'OLR08') 
                start_data_percent(4,19) = 30;
            end
            
            
            if strcmp(subject_id, 'CAD09')
                start_data_percent(1,202) = 30;
            end
            if strcmp(subject_id, 'CAD15')
                start_data_percent(1,188) = 30;
            end
            if strcmp(subject_id, 'CAD16')
                start_data_percent(2,132) = 30;
            end
            if strcmp(subject_id, 'CAD19')
                start_data_percent(2,48) = 30;
            end
            
            if strcmp(subject_id, 'CAD21')
                start_data_percent(1,97) = 30;
                start_data_percent(1,27) = 30;
            end
        end
        if strcmp(end_variable_source_type, 'percentage')
            end_data_percent = ones(size(this_step_time_data)) * str2num(end_info);
        else
            eval(['end_data_source = loaded_data.' end_variable_source_type '_data_session;']);
            eval(['end_names_source = loaded_data.' end_variable_source_type '_names_session;']);
            end_data_time_within_band = end_data_source{strcmp(end_names_source, end_info)};
            end_data_ratio = end_data_time_within_band ./ this_step_time_data;
            end_data_percent = round(end_data_ratio * 100);
			% HR (02-OCT-2019) -- this is a dirty hack. Leaving it in because it's for TF's data. Do not do things in this way in the future!
            if strcmp(subject_id,'OLR02') 
                end_data_percent(4,128) = 30;
            end
            if strcmp(subject_id,'OLR03') 
                end_data_percent(4,10) = 30;
            end
            if strcmp(subject_id,'OLR05') 
                end_data_percent(3,162) = 30;
            end
            if strcmp(subject_id,'OLR06') 
                end_data_percent(4,79) = 30;
            end
            if strcmp(subject_id,'OLR06') 
                end_data_percent(4,27) = 30;
            end
            if strcmp(subject_id,'OLR06') 
                end_data_percent(4,133) = 30;
            end
            if strcmp(subject_id,'OLR08') 
                end_data_percent(4,19) = 30;
            end
            
            
            if strcmp(subject_id, 'CAD09')
                end_data_percent(1,202) = 30;
            end
            if strcmp(subject_id, 'CAD15')
                end_data_percent(1,188) = 30;
            end
            if strcmp(subject_id, 'CAD16')
                end_data_percent(2,132) = 30;
            end
            if strcmp(subject_id, 'CAD19')
                end_data_percent(2,48) = 30;
            end
            if strcmp(subject_id, 'CAD21')
                end_data_percent(1,97) = 30;
                end_data_percent(1,27) = 30;
            end
        end
        
        % integrate
        integrated_data = zeros(bands_per_stretch, number_of_stretches);
        for i_stretch = 1 : number_of_stretches
            for i_band = 1 : bands_per_stretch
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
        analysis_data = addOrReplaceResultsData(analysis_data, integrated_data, this_variable_name, this_variable_source_directions);

    end
        
    %% calculate rms variables
    variables_to_rms_header = study_settings.get('analysis_variables_from_rms_header', 1);
    variables_to_rms = study_settings.get('analysis_variables_from_rms', 1);
    names_source = response_names_session;
    directions_source = response_directions_session;
    for i_variable = 1 : size(variables_to_rms, 1)
        this_variable_name = variables_to_rms{i_variable, strcmp(variables_to_rms_header, 'new_variable_name')};
        this_variable_source_name = variables_to_rms{i_variable, strcmp(variables_to_rms_header, 'source_variable_name')};
        this_variable_source_type = variables_to_rms{i_variable, strcmp(variables_to_rms_header, 'source_variable_type')};
        start_info = variables_to_rms{i_variable, strcmp(variables_to_rms_header, 'start')};
        start_variable_source_type = variables_to_rms{i_variable, strcmp(variables_to_rms_header, 'start_variable_type')};
        end_info = variables_to_rms{i_variable, strcmp(variables_to_rms_header, 'end')};
        end_variable_source_type = variables_to_rms{i_variable, strcmp(variables_to_rms_header, 'end_variable_type')};
        % pick data depending on source specification
        eval(['data_source = loaded_data.' this_variable_source_type '_data_session;']);
        eval(['names_source = loaded_data.' this_variable_source_type '_names_session;']);
        eval(['directions_source = loaded_data.' this_variable_source_type '_directions_session;']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        this_variable_source_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        
        % determine start and end of integration in percent of band time
        if strcmp(start_variable_source_type, 'percentage')
            start_data_percent = ones(size(this_step_time_data)) * str2num(start_info);
        else
            eval(['start_data_source = loaded_data.' start_variable_source_type '_data_session;']);
            eval(['start_names_source = loaded_data.' start_variable_source_type '_names_session;']);
            start_data_time_within_band = start_data_source{strcmp(start_names_source, start_info)};
            % find a way to take the step time per condition (maybe only
            % for CadenceGVS) and then perform integration on only
            % stretches within each condition
            % HR (09-OCT-2019) -- the above comment was from TF, but I think this is not required
            start_data_ratio = start_data_time_within_band ./ this_step_time_data;
            start_data_percent = round(start_data_ratio * 100);
			% HR (02-OCT-2019) -- this is a dirty hack. Leaving it in because it's for TF's data. Do not do things in this way in the future!
            if strcmp(subject_id,'OLR02')
                start_data_percent(4,128) = 30;
            end
            if strcmp(subject_id,'OLR03') 
                start_data_percent(4,10) = 30;
            end
            if strcmp(subject_id,'OLR05') 
                start_data_percent(3,162) = 30;
            end
            if strcmp(subject_id,'OLR06') 
                start_data_percent(4,79) = 30;
            end
            if strcmp(subject_id,'OLR06') 
                end_data_percent(4,108) = 30;
            end
            if strcmp(subject_id,'OLR06') 
                start_data_percent(4,133) = 30;
            end
            if strcmp(subject_id,'OLR08') 
                start_data_percent(4,19) = 30;
            end
            
            
            if strcmp(subject_id, 'CAD09')
                start_data_percent(1,202) = 30;
            end
            if strcmp(subject_id, 'CAD15')
                start_data_percent(1,188) = 30;
            end
            if strcmp(subject_id, 'CAD16')
                start_data_percent(2,132) = 30;
            end
            if strcmp(subject_id, 'CAD19')
                start_data_percent(2,48) = 30;
            end
            
            if strcmp(subject_id, 'CAD21')
                start_data_percent(1,97) = 30;
                start_data_percent(1,27) = 30;
            end
        end
        if strcmp(end_variable_source_type, 'percentage')
            end_data_percent = ones(size(this_step_time_data)) * str2num(end_info);
        else
            eval(['end_data_source = loaded_data.' end_variable_source_type '_data_session;']);
            eval(['end_names_source = loaded_data.' end_variable_source_type '_names_session;']);
            end_data_time_within_band = end_data_source{strcmp(end_names_source, end_info)};
            end_data_ratio = end_data_time_within_band ./ this_step_time_data;
            end_data_percent = round(end_data_ratio * 100);
			% HR (02-OCT-2019) -- this is a dirty hack. Leaving it in because it's for TF's data. Do not do things in this way in the future!
            if strcmp(subject_id,'OLR02') 
                end_data_percent(4,128) = 30;
            end
            if strcmp(subject_id,'OLR03') 
                end_data_percent(4,10) = 30;
            end
            if strcmp(subject_id,'OLR05') 
                end_data_percent(3,162) = 30;
            end
            if strcmp(subject_id,'OLR06') 
                end_data_percent(4,79) = 30;
            end
            if strcmp(subject_id,'OLR06') 
                end_data_percent(4,27) = 30;
            end
            if strcmp(subject_id,'OLR06') 
                end_data_percent(4,133) = 30;
            end
            if strcmp(subject_id,'OLR08') 
                end_data_percent(4,19) = 30;
            end
            
            
            if strcmp(subject_id, 'CAD09')
                end_data_percent(1,202) = 30;
            end
            if strcmp(subject_id, 'CAD15')
                end_data_percent(1,188) = 30;
            end
            if strcmp(subject_id, 'CAD16')
                end_data_percent(2,132) = 30;
            end
            if strcmp(subject_id, 'CAD19')
                end_data_percent(2,48) = 30;
            end
            if strcmp(subject_id, 'CAD21')
                end_data_percent(1,97) = 30;
                end_data_percent(1,27) = 30;
            end
        end
        
        % rms
        rms_data = zeros(bands_per_stretch, number_of_stretches);
        for i_stretch = 1 : number_of_stretches
            for i_band = 1 : bands_per_stretch
                [band_start_index, band_end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
                this_band_time_full = linspace(0, this_step_time_data(i_band), number_of_time_steps_normalized);
                this_band_data_full = this_variable_source_data(band_start_index : band_end_index, i_stretch);
                range = (start_data_percent(i_band, i_stretch) : end_data_percent(i_band, i_stretch)) + 1;
                if isempty(range)
                    rms_data(i_band, i_stretch) = 0;
                else
                    this_band_time_range = this_band_time_full(range);
                    this_band_data_range = this_band_data_full(range);
                    % rms
%                    this_band_data_rms = cumtrapz(this_band_time_range, this_band_data_range);           
 %                   rms_data(i_band, i_stretch) = this_band_data_rms(end);
                     this_band_data_rms = rms(this_band_data_range);           
                     rms_data(i_band, i_stretch) = this_band_data_rms(end);
                end

            end
        end        
        
        % store
        analysis_data = addOrReplaceResultsData(analysis_data, rms_data, this_variable_name, this_variable_source_directions);
          
    end
    
    %% calculate band end variables
    variables_step_end = study_settings.get('analysis_variables_from_band_end', 1);
    for i_variable = 1 : size(variables_step_end, 1)
        this_variable_name = variables_step_end{i_variable, 1};
        this_variable_source_name = variables_step_end{i_variable, 2};
        this_variable_source_type = variables_step_end{i_variable, 3};
        % pick data depending on source specification
        eval(['data_source = loaded_data.' this_variable_source_type '_data_session;']);
        eval(['names_source = loaded_data.' this_variable_source_type '_names_session;']);
        eval(['directions_source = loaded_data.' this_variable_source_type '_directions_session;']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        new_variable_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        step_end_data = zeros(bands_per_stretch, number_of_stretches);
        for i_band = 1 : bands_per_stretch
            [~, end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
            step_end_data(i_band, :) = this_variable_source_data(end_index, :);
        end
        
        % store
        analysis_data = addOrReplaceResultsData(analysis_data, step_end_data, this_variable_name, new_variable_directions);
    end
    
    %% calculate band percent variables
    variables_band_percent_header = study_settings.get('analysis_variables_from_band_percent_header', 1);
    variables_band_percent = study_settings.get('analysis_variables_from_band_percent', 1);
    for i_variable = 1 : size(variables_band_percent, 1)
        this_variable_name = variables_band_percent{i_variable, strcmp(variables_band_percent_header, 'new_variable_name')};
        this_variable_source_name = variables_band_percent{i_variable, strcmp(variables_band_percent_header, 'source_variable_name')};
        this_variable_source_type = variables_band_percent{i_variable, strcmp(variables_band_percent_header, 'source_variable_type')};
        this_variable_source_percent = str2num(variables_band_percent{i_variable, strcmp(variables_band_percent_header, 'percentage')});
        % pick data depending on source specification
        eval(['data_source = loaded_data.' this_variable_source_type '_data_session;']);
        eval(['names_source = loaded_data.' this_variable_source_type '_names_session;']);
        eval(['directions_source = loaded_data.' this_variable_source_type '_directions_session;']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        new_variable_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        band_percent_data = zeros(bands_per_stretch, number_of_stretches);
        for i_band = 1 : bands_per_stretch
            [start_index, end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
            this_index = start_index + this_variable_source_percent/100 * (end_index - start_index);
            band_percent_data(i_band, :) = this_variable_source_data(this_index, :);
        end
        
        % store
        analysis_data = addOrReplaceResultsData(analysis_data, band_percent_data, this_variable_name, new_variable_directions);
    end
    
    %% calculate variables referenced by stretch
    variables_referenced_by_stretch = study_settings.get('analysis_variables_referenced_by_stretch', 1);
    for i_variable = 1 : size(variables_referenced_by_stretch, 1)
        this_variable_name = variables_referenced_by_stretch{i_variable, 1};
        this_variable_source_name = variables_referenced_by_stretch{i_variable, 2};
        this_variable_source_type = variables_referenced_by_stretch{i_variable, 3};
        this_variable_reference_band_index = str2num(variables_referenced_by_stretch{i_variable, 4});
        this_variable_reference_point_percentage_within_band = str2num(variables_referenced_by_stretch{i_variable, 5});
        
        % pick data depending on source specification
        eval(['data_source = loaded_data.' this_variable_source_type '_data_session;']);
        eval(['names_source = loaded_data.' this_variable_source_type '_names_session;']);
        eval(['directions_source = loaded_data.' this_variable_source_type '_directions_session;']);
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
        analysis_data = addOrReplaceResultsData(analysis_data, new_variable_data, this_variable_name, new_variable_directions);
    end
        
    %% process variables where something specific happens for each variable
    special_variables_to_calculate = study_settings.get('analysis_variables_special', 1);
    for i_variable = 1:size(special_variables_to_calculate, 1)
        this_variable_name = special_variables_to_calculate{i_variable, 1};
        this_variable_source_name = special_variables_to_calculate{i_variable, 2};
        this_variable_source_type = special_variables_to_calculate{i_variable, 3};
        
          % pick data depending on source specification
        eval(['data_source = loaded_data.' this_variable_source_type '_data_session;']);
        eval(['names_source = loaded_data.' this_variable_source_type '_names_session;']);
        eval(['directions_source = loaded_data.' this_variable_source_type '_directions_session;']);
        this_variable_source_directions = directions_source(strcmp(names_source, this_variable_source_name), :);
        new_variable_directions = this_variable_source_directions;
        % HR: something went wrong here during merge on 5. December, 2018, fix later

        % TO DO: need to automate the data type loading and storage
        if false  
%         if strcmp(this_variable_name, 'step_symmetry_index')
%             % TO DO: adjust for bands
%             this_variable_source_index = find(strcmp(stretch_names_session, this_variable_source_name), 1, 'first');
%             this_variable_source_data = stretch_data_session{this_variable_source_index};
% 
%             average_step_time = mean(reshape(this_variable_source_data, 1, length(this_variable_source_data)*2));
%             for i_stretch = 1:  length(this_variable_source_data)
%                 % find the left and right stance data
%                 if strcmp(condition_data_all(i_stretch,3), 'STANCE_LEFT')
%                     left_step_index = 2;
%                     right_step_index = 1;
%                 else
%                     left_step_index = 1;
%                     right_step_index = 2;
%                 end
%                 this_left_step_time = this_variable_source_data(left_step_index,i_stretch);
%                 this_right_step_time = this_variable_source_data(right_step_index,i_stretch);
%                 this_variable_data(1:2,i_stretch) = (this_left_step_time - this_right_step_time) / average_step_time;
%             end
%         end
%         if strcmp(this_variable_name, 'com_from_com_init_x')
%             % TO DO: adjust for bands
%             % should we be looking at response or stretch variable here??
%             this_variable_source_index = find(strcmp(response_names_session, this_variable_source_name), 1, 'first');
%             this_variable_source_data = response_data_session{this_variable_source_index};  
%             this_variable_data = [];
%             for i_stretch = 1:length(this_variable_source_data)
%                 this_variable_data(:,i_stretch) = this_variable_source_data(:,i_stretch) - this_variable_source_data(1,i_stretch);
%             end
%         end
%        if strcmp(this_variable_name, 'trigger_leg_ankle_dorsiflexion_inverted_max')
%            % TO DO: adjust for bands
%             this_variable_source_index = find(strcmp(analysis_names_session, this_variable_source_name), 1, 'first');
%             this_variable_source_data = analysis_data_session{this_variable_source_index};
%             
%             
%             for i_stretch = 1:length(this_variable_source_data)
%                 % create time
%                 this_stretch_time_full = linspace(0, this_step_time_data(i_stretch), 100);
%                 this_stretch_data_full = this_variable_source_data(:, i_stretch);
% 
%                 % interpolate double stance to 100 data points
%                 this_stretch_time_double = linspace(0, this_pushoff_time_data(i_stretch), 100);       
%                 this_stretch_data_double = interp1(this_stretch_time_full, this_stretch_data_full, this_stretch_time_double);
% 
% 
%                 this_dorsi_angle_value = max(findpeaks(this_stretch_data_double));
%                 if ~isempty(this_dorsi_angle_value)
%                     this_variable_data(i_stretch) = this_dorsi_angle_value;
%                 else
%                     this_variable_data(i_stretch) = NaN;
%                 end
%             end
%        end
%        if strcmp(this_variable_name, 'cop_from_com_x_integrated_twice')
%            % TO DO: adjust for bands
%             this_variable_source_index = find(strcmp(response_names_session, this_variable_source_name), 1, 'first');
%             this_variable_source_data = response_data_session{this_variable_source_index};
%             number_of_stretches = size(this_variable_source_data, 2);
% 
%             for i_stretch = 1 : number_of_stretches
%                 % get data for full step
%                 this_stretch_time_full = linspace(0, this_step_time_data(i_stretch), 100);
%                 this_stretch_data_full = this_variable_source_data(:, i_stretch);
% 
%                 % interpolate single stance to 100 data points
%                 this_stretch_time_single = linspace(this_pushoff_time_data(i_stretch), this_step_time_data(i_stretch), 100);
%                 this_stretch_data_single = interp1(this_stretch_time_full, this_stretch_data_full, this_stretch_time_single);
% 
%                 % integrate data in single stance
%                 this_stretch_data_single_integrated = cumtrapz(this_stretch_time_single, this_stretch_data_single);
%                 this_stretch_data_single_integrated_twice = cumtrapz(this_stretch_time_single, this_stretch_data_single_integrated);
%                 this_variable_data(i_stretch) = this_stretch_data_single_integrated_twice(end);
%             end 
%             
%        end
        end
        if strcmp(this_variable_name, 'trigger_leg_ankle_dorsiflexion_dsmid') || strcmp(this_variable_name, 'com_x_vel_sym_doublestance_mid')
            this_variable_source_index = find(strcmp(analysis_names_session, this_variable_source_name), 1, 'first');
            this_variable_source_data = analysis_data_session{this_variable_source_index};
            number_of_stretches = size(this_variable_source_data, 2);
            
            start_data_percent = 1;
            
            end_data_time_within_band = this_pushoff_time_data;
            end_data_ratio = end_data_time_within_band ./ this_step_time_data;
            end_data_percent = round(end_data_ratio * 100);

            this_variable_data = zeros(bands_per_stretch, number_of_stretches);
            for i_stretch = 1 : number_of_stretches
                for i_band = 1 : bands_per_stretch
                    [band_start_index, band_end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
                    this_band_time_full = linspace(0, this_step_time_data(i_band), 100);
                    this_band_data_full = this_variable_source_data(band_start_index : band_end_index, i_stretch);

                    double_stance_mid_index = round((end_data_percent(i_band, i_stretch) - start_data_percent)/2);

                    this_variable_data(i_band, i_stretch) = this_band_data_full(double_stance_mid_index);

                end
            end
        end
       
       if strcmp(this_variable_name,'com_x_inverted_pushoff_end') || strcmp(this_variable_name,'com_x_vel_inverted_pushoff_end') || ... 
               strcmp(this_variable_name,'com_x_sym_pushoff_end') || strcmp(this_variable_name,'com_x_vel_sym_pushoff_end') 
           
           this_variable_source_index = find(strcmp(analysis_names_session, this_variable_source_name), 1, 'first');
           this_variable_source_data = analysis_data_session{this_variable_source_index};
           number_of_stretches = size(this_variable_source_data, 2);

           end_data_time_within_band = this_pushoff_time_data;
           end_data_ratio = end_data_time_within_band ./ this_step_time_data;
           end_data_percent = round(end_data_ratio * 100);
           
           if strcmp(subject_id,'OLR02') 
                end_data_percent(4,128) = 30;
            end
            if strcmp(subject_id,'OLR03') 
                end_data_percent(4,10) = 30;
            end
            if strcmp(subject_id,'OLR05') 
                end_data_percent(3,162) = 30;
            end
            if strcmp(subject_id,'OLR06') 
                end_data_percent(4,79) = 30;
            end
            if strcmp(subject_id,'OLR06') 
                end_data_percent(4,27) = 30;
            end
            if strcmp(subject_id,'OLR06') 
                end_data_percent(4,133) = 30;
            end
            if strcmp(subject_id,'OLR08') 
                end_data_percent(4,19) = 30;
            end
            if strcmp(subject_id, 'CAD09')
                end_data_percent(1,202) = 30;
            end
            if strcmp(subject_id, 'CAD15')
                end_data_percent(1,188) = 30;
            end
            if strcmp(subject_id, 'CAD16')
                end_data_percent(2,132) = 30;
            end
            if strcmp(subject_id, 'CAD19')
                end_data_percent(2,48) = 30;
            end
            if strcmp(subject_id, 'CAD21')
                end_data_percent(1,97) = 30;
                end_data_percent(1,27) = 30;
            end
           
           
            for i_stretch = 1 : number_of_stretches
                for i_band = 1 : bands_per_stretch
                    [band_start_index, band_end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
            %                    this_band_time_full = linspace(0, this_step_time_data(i_band), 100);
                    this_band_data_full = this_variable_source_data(band_start_index : band_end_index, i_stretch);

                    this_variable_data(i_band, i_stretch) = this_band_data_full(end_data_percent(i_band, i_stretch));

                end
            end
        end
       if strcmp(this_variable_name, 'cop_from_com_x_integrated_twice')
            this_variable_source_index = find(strcmp(response_names_session, this_variable_source_name), 1, 'first');
            this_variable_source_data = response_data_session{this_variable_source_index};
            number_of_stretches = size(this_variable_source_data, 2);

            this_variable_source_index = find(strcmp(response_names_session, this_variable_source_name), 1, 'first');
            this_variable_source_data = response_data_session{this_variable_source_index};
            number_of_stretches = size(this_variable_source_data, 2);
       end
       
       if strcmp(this_variable_name,'com_x_inverted_band_end') || strcmp(this_variable_name,'com_x_vel_inverted_band_end') || strcmp(this_variable_name,'fy_band_end') || ...
               strcmp(this_variable_name,'com_x_sym_band_end') || strcmp(this_variable_name,'com_x_vel_sym_band_end')
           
           if strcmp(this_variable_name,'com_x_inverted_band_end') || strcmp(this_variable_name,'com_x_vel_inverted_band_end') || ...
                   strcmp(this_variable_name,'com_x_sym_band_end') || strcmp(this_variable_name,'com_x_vel_sym_band_end') 
               this_variable_source_index = find(strcmp(analysis_names_session, this_variable_source_name), 1, 'first');
               this_variable_source_data = analysis_data_session{this_variable_source_index};
           else
               this_variable_source_index = find(strcmp(response_names_session, this_variable_source_name), 1, 'first');
               this_variable_source_data = response_data_session{this_variable_source_index};
           end
%            end_data_time_within_band = this_pushoff_time_data;
%            end_data_ratio = end_data_time_within_band ./ this_step_time_data;
%            end_data_percent = round(end_data_ratio * 100);
            for i_stretch = 1 : number_of_stretches
               for i_band = 1 : bands_per_stretch
                   [band_start_index, band_end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
%                    this_band_time_full = linspace(0, this_step_time_data(i_band), 100);
                   this_band_data_full = this_variable_source_data(band_start_index : band_end_index, i_stretch);
                   
                   this_variable_data(i_band, i_stretch) = this_band_data_full(end);
                    
               end
           end
       end
       
       if strcmp(this_variable_name,'com_x_sym_1sec') || strcmp(this_variable_name,'com_x_vel_sym_1sec')
           this_variable_source_index = find(strcmp(analysis_names_session, this_variable_source_name), 1, 'first');
           this_variable_source_data = analysis_data_session{this_variable_source_index};
           
           number_of_stretches = size(this_variable_source_data, 2);
           
           % need to assign end data ratio
           % 1) take first band and determine ratio step time band one - 1
           % (total time)
           % what to do if time greater than 1 falls outside stride?
           end_data_time = 1;  
           
           
%            end_data_ratio = end_data_time_within_band ./ this_step_time_data;
%           
%            end_data_percent = round(end_data_ratio * 100);
           
            for i_stretch = 1 : number_of_stretches
               i_band_one = 1;
               i_band_two = 2;
                   [band_one_start_index, band_one_end_index] = getBandIndices(i_band_one, number_of_time_steps_normalized);
                   [band_two_start_index, band_two_end_index] = getBandIndices(i_band_two, number_of_time_steps_normalized);
                   band_one_data_time = this_step_time_data(i_band_one, i_stretch);
                   band_two_data_time = this_step_time_data(i_band_two, i_stretch);
                   
                   band_two_data_full = this_variable_source_data(band_two_start_index : band_two_end_index, i_stretch);
                   
                   end_data_ratio_two = (end_data_time - band_one_data_time) / this_step_time_data(i_band_two, i_stretch);
                   end_data_percent_two = round(end_data_ratio_two * 100);
                   if end_data_percent_two > 100
                       end_data_percent_two = 100;
                   end
                   this_variable_data(:, i_stretch) = band_two_data_full(end_data_percent_two);
           end
       end
       
        % store
        analysis_data = addOrReplaceResultsData(analysis_data, this_variable_data, this_variable_name, new_variable_directions);
        
    end
    
    %% calculate variables from extrema
    variables_from_extrema = study_settings.get('analysis_variables_from_extrema', 1);
    for i_variable = 1 : size(variables_from_extrema, 1)
        % get data
        this_variable_name = variables_from_extrema{i_variable, 1};
        this_variable_source_name = variables_from_extrema{i_variable, 2};
        this_variable_source_type = variables_from_extrema{i_variable, 3};
        this_variable_extremum_type = variables_from_extrema{i_variable, 4};
        % pick data depending on source specification
        eval(['data_source = loaded_data.' this_variable_source_type '_data_session;']);
        eval(['names_source = loaded_data.' this_variable_source_type '_names_session;']);
        eval(['directions_source = loaded_data.' this_variable_source_type '_directions_session;']);
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
        analysis_data = addOrReplaceResultsData(analysis_data, extrema_data, this_variable_name, new_variable_directions);
    end
    
    %% calculate variables from extrema for range
    variables_from_extrema_range = study_settings.get('analysis_variables_from_extrema_range', 1);
    for i_variable = 1 : size(variables_from_extrema_range, 1)
        % get data
        this_variable_name = variables_from_extrema_range{i_variable, 1};
        this_variable_source_name = variables_from_extrema_range{i_variable, 2};
        this_variable_source_type = variables_from_extrema_range{i_variable, 3};
        this_variable_extremum_type = variables_from_extrema_range{i_variable, 4};
        % pick data depending on source specification
        eval(['data_source = loaded_data.' this_variable_source_type '_data_session;']);
        eval(['names_source = loaded_data.' this_variable_source_type '_names_session;']);
        eval(['directions_source = loaded_data.' this_variable_source_type '_directions_session;']);
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
        range_data = addOrReplaceResultsData(range_data, extrema_data, this_variable_name, new_variable_directions);
    end
    
return
    
    %% save data
    variables_to_save = loaded_data;

    variables_to_save.analysis_data_session = analysis_data.data;
    variables_to_save.analysis_directions_session = analysis_data.directions;
    variables_to_save.analysis_names_session = analysis_data.names;
    variables_to_save.range_data_session = range_data.data;
    variables_to_save.range_directions_session = range_data.directions;
    variables_to_save.range_names_session = range_data.names;
    
    save(results_file_name, '-struct', 'variables_to_save');    

    

end













