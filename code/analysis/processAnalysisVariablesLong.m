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

function processAnalysisVariablesLong(varargin)
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
    results_file_name = ['analysis' filesep makeFileName(date, subject_id, 'longStretchResults')];
    loaded_data = load(results_file_name);
    
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    bands_per_stretch = loaded_data.bands_per_long_stretch;
    stretch_names_session = loaded_data.long_stretch_data_labels_session;
    stretch_data_session = loaded_data.long_stretch_data_session;
    response_names_session = loaded_data.long_stretch_data_labels_session;
    response_data_session = loaded_data.long_stretch_response_data_session;
    
    number_of_stretch_variables = length(stretch_names_session);
    number_of_stretches = size(stretch_data_session{1}, 2); %#ok<*USENS>
    
    % make condition data tables
    conditions_settings = study_settings.get('conditions');
    condition_labels = conditions_settings(:, 1)';
    condition_source_variables = conditions_settings(:, 2)';
    number_of_condition_labels = length(condition_labels);
    conditions_session = loaded_data.long_stretch_conditions_session;
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
    
    %% calculate variables from inversion
    variables_to_invert = study_settings.get('inversion_variables_long');
    for i_variable = 1 : size(variables_to_invert, 1)
        % get data
        this_variable_name = variables_to_invert{i_variable, 1};
        this_variable_source_name = variables_to_invert{i_variable, 2};
        this_variable_source_type = variables_to_invert{i_variable, 3};
        % pick data depending on source specification
        eval(['data_source = ' this_variable_source_type '_data_session;']);
        eval(['names_source = ' this_variable_source_type '_names_session;']);
        this_variable_source_data = data_source{strcmp(names_source, this_variable_source_name)};
        
        relevant_condition = variables_to_invert{i_variable, 4};
        inversion_table = study_settings.get(variables_to_invert{i_variable, 5});
        
        % go through levels and invert
        this_variable_data = this_variable_source_data;
        level_list = conditions_session.(condition_source_variables{strcmp(condition_labels, relevant_condition)});
        for i_level = 1 : size(inversion_table, 1)
            sign_this_level = str2double(inversion_table(i_level, 2));
            
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
                this_variable_data, this_variable_name, {'~', '~'} ...
              );
        
    end
    
    
    %% calculate variables from extrema
    variables_from_extrema = study_settings.get('analysis_variables_from_long_extrema');
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
        
        if strcmp(this_variable_extremum_type, 'min')
            extrema_data = min(this_variable_source_data);
        end
        if strcmp(this_variable_extremum_type, 'max')
            extrema_data = max(this_variable_source_data);
        end
        if ~strcmp(this_variable_extremum_type, 'min') && ~strcmp(this_variable_extremum_type, 'max')
            error(['"' this_variable_extremum_type '" is not a valid type for variables_from_extrema. Acceptable types are "min" or "max".']);
        end
     
          [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
            addOrReplaceResultsData ...
              ( ...
                analysis_data_session, analysis_names_session, analysis_directions_session, ...
                extrema_data, this_variable_name, new_variable_directions ...
              );
    end


    
    %% save data
    variables_to_save = loaded_data;
    variables_to_save.analysis_data_session = analysis_data_session;
    variables_to_save.analysis_directions_session = analysis_directions_session;
    variables_to_save.analysis_names_session = analysis_names_session;
    variables_to_save.analysis_directions_session = analysis_directions_session;
    save(results_file_name, '-struct', 'variables_to_save');    

    

end













