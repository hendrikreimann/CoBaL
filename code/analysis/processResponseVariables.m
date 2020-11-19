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

function processResponseVariables(varargin)
    % load settings and existing results
    study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
    results_file_name = ['results' filesep makeFileName(collection_date, subject_id, 'results')];
    loaded_data = load(results_file_name);
    
    % get some numbers
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
    levels_to_remove = study_settings.get('levels_to_remove', 1);
    comparisons = struct;
    [comparisons.combination_labels, comparisons.condition_combinations_control] ...
        = determineConditionCombinations ...
            (condition_data_all, conditions_settings, labels_to_ignore, levels_to_remove, 'control');
    comparisons.condition_combinations_control_unique ...
        = table2cell(unique(cell2table(comparisons.condition_combinations_control), 'rows'));
    
    %% calculate response (i.e. difference from control mean)
    response_data_session = {};
    response_directions_session = loaded_data.stretch_directions_session;
    response_names_session = loaded_data.stretch_names_session;
    if ~isempty(comparisons.condition_combinations_control)
        % prepare container
        response_data_session = cell(size(loaded_data.stretch_data_session));
        for i_variable = 1 : number_of_stretch_variables
            response_data_session{i_variable} = zeros(size(loaded_data.stretch_data_session{i_variable}));
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
            for i_variable = 1 : number_of_stretch_variables
                % calculate control mean
                data_this_variable = loaded_data.stretch_data_session{i_variable};
                this_condition_control_data = data_this_variable(:, control_condition_indicator);
                this_condition_control_mean = mean(this_condition_control_data, 2);
                
                % calculate response
                response_data_session{i_variable}(:, i_stretch) ...
                    = loaded_data.stretch_data_session{i_variable}(:, i_stretch) - this_condition_control_mean;
            end
            
        end

    end

    %% save data
    variables_to_save = loaded_data;
    variables_to_save.response_data_session = response_data_session;
    variables_to_save.response_directions_session = response_directions_session;
    variables_to_save.response_names_session = response_names_session;
    save(results_file_name, '-struct', 'variables_to_save');    

end



function applicable_control_condition_index = determineControlConditionIndex(study_settings, comparisons, this_stretch_condition_string)
    experimental_paradigm = study_settings.get('experimental_paradigm');
    paradigms_with_intermittent_perturbation = {'Vision', 'GVS', 'GVS_old', 'Vision_old'};
    paradigms_with_stochastic_resonance = {'SR_VisualStim', 'nGVS_Vision'};
  
    if any(strcmp(experimental_paradigm, paradigms_with_stochastic_resonance))
        % HR 2020-06-04: this is supposed to be the general approach, but
        % only testing this for SR_VisualStim for now
        if strcmp(study_settings.get('experimental_paradigm', 1), 'SR_VisualStim')
            relevant_factors_for_control = {'stim_amplitude', 'trigger_foot'};
        end
        if strcmp(study_settings.get('experimental_paradigm', 1), 'nGVS_Vision')
            relevant_factors_for_control = {'ngvs_settings', 'trigger_foot'};
        end
        control_row_indicator = true(size(comparisons.condition_combinations_control_unique, 1), 1);
        for i_factor = 1 : length(relevant_factors_for_control)
            this_factor_label = relevant_factors_for_control{i_factor};
            this_factor_column = strcmp(comparisons.combination_labels, this_factor_label);
            this_factor_this_level = this_stretch_condition_string(this_factor_column);
            this_factor_control_levels = comparisons.condition_combinations_control_unique(:, this_factor_column);
            this_factor_candidate_rows = strcmp(this_factor_control_levels, this_factor_this_level);
            control_row_indicator = control_row_indicator & this_factor_candidate_rows;
        end
        applicable_control_condition_index = find(control_row_indicator);
    end

    if any(strcmp(experimental_paradigm, paradigms_with_intermittent_perturbation))
        % define the factors for which the levels have to match to determine the control condition for this condition
        relevant_factors_for_control = {'trigger_foot'}; % control only differs by trigger foot in this paradigm
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
            control_row_indicator = control_row_indicator & this_factor_candidate_rows;
        end
        applicable_control_condition_index = find(control_row_indicator);
    end
    if strcmp(study_settings.get('experimental_paradigm', 1), 'CadenceGVS')
        if strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'cadence')}, '80BPM') && strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
            applicable_control_condition_index = find(strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'cadence')), '80BPM') & strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_LEFT'));
        elseif strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'cadence')}, '80BPM') && strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
            applicable_control_condition_index = find(strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'cadence')), '80BPM') & strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT'));
        end
        if strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'cadence')}, '110BPM') && strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
            applicable_control_condition_index = find(strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'cadence')), '110BPM') & strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_LEFT'));
        elseif strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'cadence')}, '110BPM') && strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
             applicable_control_condition_index = find(strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'cadence')), '110BPM') & strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT'));
        end
    end
    if strcmp(study_settings.get('experimental_paradigm', 1), 'FatigueGVS')
        if strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'fatigue')}, 'UNFATIGUED') ...
                && strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
            applicable_control_condition_index = ...
                find ...
                  ( ...
                    strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'fatigue')), 'UNFATIGUED') ...
                    & strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_LEFT') ...
                  );
        elseif strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'fatigue')}, 'UNFATIGUED') ...
                && strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
            applicable_control_condition_index = ...
                find ...
                  ( ...
                    strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'fatigue')), 'UNFATIGUED') ...
                    & strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT') ...
                  );
        end
        if strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'fatigue')}, 'FATIGUED') ...
            && strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
            applicable_control_condition_index = ...
                find ...
                  ( ...
                    strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'fatigue')), 'FATIGUED') ...
                    & strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_LEFT') ...
                  );
        elseif strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'fatigue')}, 'FATIGUED') ...
                && strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
             applicable_control_condition_index = ...
                 find ...
                   ( ...
                     strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'fatigue')), 'FATIGUED') ...
                     & strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT') ...
                   );
        end
    end
    if strcmp(study_settings.get('experimental_paradigm', 1), 'CognitiveLoadVision') || strcmp(study_settings.get('experimental_paradigm', 1), 'CognitiveLoadGvs')
        if strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'cognitive_load')}, 'NO_LOAD') ...
                && strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
            applicable_control_condition_index = ...
                find ...
                  ( ...
                    strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'cognitive_load')), 'NO_LOAD') ...
                    & strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_LEFT') ...
                  );
        elseif strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'cognitive_load')}, 'NO_LOAD') ...
                && strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
            applicable_control_condition_index = ...
                find ...
                  ( ...
                    strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'cognitive_load')), 'NO_LOAD') ...
                    & strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT') ...
                  );
        end
        if strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'cognitive_load')}, 'BACK_7') ...
            && strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
            applicable_control_condition_index = ...
                find ...
                  ( ...
                    strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'cognitive_load')), 'BACK_7') ...
                    & strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_LEFT') ...
                  );
        elseif strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'cognitive_load')}, 'BACK_7') ...
                && strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
             applicable_control_condition_index = ...
                 find ...
                   ( ...
                     strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'cognitive_load')), 'BACK_7') ...
                     & strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT') ...
                   );
        end
    end
    if strcmp(study_settings.get('experimental_paradigm', 1), 'GvsOverground')
         if strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
             applicable_control_condition_index = find(strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_LEFT'));
         end
         if strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
             applicable_control_condition_index = find(strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT'));
         end
    end       
    if strcmp(study_settings.get('experimental_paradigm'), 'OculusLaneRestriction')
        if strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'zone_side')}, 'STIM_ZONE_LEFT') && strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
            applicable_control_condition_index = find(strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'zone_side')), 'STIM_ZONE_LEFT') & strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_LEFT'));
        elseif strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'zone_side')}, 'STIM_ZONE_LEFT') && strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
            applicable_control_condition_index = find(strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'zone_side')), 'STIM_ZONE_LEFT') & strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT'));
        elseif strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'zone_side')}, 'STIM_ZONE_RIGHT') && strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
            applicable_control_condition_index = find(strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'zone_side')), 'STIM_ZONE_RIGHT') & strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_LEFT'));
        elseif strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'zone_side')}, 'STIM_ZONE_RIGHT') && strcmp(this_stretch_condition_string{strcmp(comparisons.combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
            applicable_control_condition_index = find(strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'zone_side')), 'STIM_ZONE_RIGHT') & strcmp(comparisons.condition_combinations_control_unique(:, strcmp(comparisons.combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT'));
        end
    end

end







