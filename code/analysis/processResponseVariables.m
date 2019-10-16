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
    levels_to_remove = study_settings.get('levels_to_remove', 1);
    [condition_combination_labels, ~, condition_combinations_control] = determineConditionCombinations(condition_data_all, conditions_settings, labels_to_ignore, levels_to_remove);
    condition_combinations_control_unique = table2cell(unique(cell2table(condition_combinations_control), 'rows'));
    
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
            if strcmp(study_settings.get('experimental_paradigm', 1), 'Vision') || strcmp(study_settings.get('experimental_paradigm', 1), 'GVS') || strcmp(study_settings.get('experimental_paradigm', 1), 'GVS_old')
                 if exist('affected_side')
                     if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'affected_side')}, 'TRIGGER_AFFECTED')
                         applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'affected_side')), 'TRIGGER_AFFECTED'));
                     end
                     if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'affected_side')}, 'TRIGGER_UNAFFECTED')
                         applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'affected_side')), 'TRIGGER_UNAFFECTED'));
                     end
                 else
                     if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
                         applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_LEFT'));
                     end
                     if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
                         applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT'));
                     end
                 end
            end       
            if strcmp(study_settings.get('experimental_paradigm', 1), 'CadenceGVS')
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
            if strcmp(study_settings.get('stimulus_condition', 1), 'VISUAL') || strcmp(study_settings.get('stimulus_condition', 1), 'GVS')
                if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
                    applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_LEFT'));
                end
                if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
                    applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT'));
                end
            end
            if strcmp(study_settings.get('experimental_paradigm', 1), 'FatigueGVS')
                if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'fatigue')}, 'UNFATIGUED') ...
                        && strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
                    applicable_control_condition_index = ...
                        find ...
                          ( ...
                            strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'fatigue')), 'UNFATIGUED') ...
                            & strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_LEFT') ...
                          );
                elseif strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'fatigue')}, 'UNFATIGUED') ...
                        && strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
                    applicable_control_condition_index = ...
                        find ...
                          ( ...
                            strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'fatigue')), 'UNFATIGUED') ...
                            & strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT') ...
                          );
                end
                if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'fatigue')}, 'FATIGUED') ...
                    && strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
                    applicable_control_condition_index = ...
                        find ...
                          ( ...
                            strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'fatigue')), 'FATIGUED') ...
                            & strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_LEFT') ...
                          );
                elseif strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'fatigue')}, 'FATIGUED') ...
                        && strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
                     applicable_control_condition_index = ...
                         find ...
                           ( ...
                             strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'fatigue')), 'FATIGUED') ...
                             & strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT') ...
                           );
                end
            end
            if strcmp(study_settings.get('experimental_paradigm', 1), 'CognitiveLoadVision') || strcmp(study_settings.get('experimental_paradigm', 1), 'CognitiveLoadGvs')
                if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'cognitive_load')}, 'NO_LOAD') ...
                        && strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
                    applicable_control_condition_index = ...
                        find ...
                          ( ...
                            strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'cognitive_load')), 'NO_LOAD') ...
                            & strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_LEFT') ...
                          );
                elseif strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'cognitive_load')}, 'NO_LOAD') ...
                        && strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
                    applicable_control_condition_index = ...
                        find ...
                          ( ...
                            strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'cognitive_load')), 'NO_LOAD') ...
                            & strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT') ...
                          );
                end
                if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'cognitive_load')}, 'BACK_7') ...
                    && strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
                    applicable_control_condition_index = ...
                        find ...
                          ( ...
                            strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'cognitive_load')), 'BACK_7') ...
                            & strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_LEFT') ...
                          );
                elseif strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'cognitive_load')}, 'BACK_7') ...
                        && strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
                     applicable_control_condition_index = ...
                         find ...
                           ( ...
                             strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'cognitive_load')), 'BACK_7') ...
                             & strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT') ...
                           );
                end
            end
            if strcmp(study_settings.get('experimental_paradigm', 1), 'GvsOverground')
                 if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
                     applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_LEFT'));
                 end
                 if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
                     applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT'));
                 end
            end       
            if strcmp(study_settings.get('experimental_paradigm'), 'OculusLaneRestriction')
                if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'zone_side')}, 'STIM_ZONE_LEFT') && strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
                    applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'zone_side')), 'STIM_ZONE_LEFT') & strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_LEFT'));
                elseif strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'zone_side')}, 'STIM_ZONE_LEFT') && strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
                    applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'zone_side')), 'STIM_ZONE_LEFT') & strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT'));
                elseif strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'zone_side')}, 'STIM_ZONE_RIGHT') && strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
                    applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'zone_side')), 'STIM_ZONE_RIGHT') & strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_LEFT'));
                elseif strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'zone_side')}, 'STIM_ZONE_RIGHT') && strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
                    applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'zone_side')), 'STIM_ZONE_RIGHT') & strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT'));
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
                
                this_condition_control_mean = mean(this_condition_control_data, 2);
                
                % calculate response
                response_data_session{i_variable}(:, i_stretch) = loaded_data.stretch_data_session{i_variable}(:, i_stretch) - this_condition_control_mean;
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











