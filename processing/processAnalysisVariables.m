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

% process stretches


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
    load(results_file_name);
    
    conditions_control = study_settings.get('conditions_control');
    number_of_stretch_variables = length(stretch_names_session);
    
    analysis_data_session = {};
    analysis_names_session = {};
    
    %% calculate response (i.e. difference from control mean)
    response_data_session = {};
    response_names_session = stretch_names_session;
    if ~isempty(conditions_control)
        % prepare container
        response_data_session = cell(size(stretch_data_session)); %#ok<USENS>
        for i_variable = 1 : number_of_stretch_variables
            response_data_session{i_variable} = zeros(size(stretch_data_session{i_variable}));
        end        
        
        % new code - go stretch by stretch
        number_of_stretches = size(stretch_data_session{1}, 2);
        for i_stretch = 1 : number_of_stretches
            % determine conditions and applicable control
            this_stretch_condition_string = ...
              { ...
                condition_stance_foot_list_session{i_stretch}, ...
                condition_perturbation_list_session{i_stretch}, ...
                condition_delay_list_session{i_stretch}, ...
                condition_index_list_session{i_stretch}, ...
                condition_experimental_list_session{i_stretch}, ...
                condition_stimulus_list_session{i_stretch}, ...
                condition_day_list_session{i_stretch} ...
              };  %#ok<USENS>
            applicable_control_condition_index = findApplicableControlConditionIndex(this_stretch_condition_string, conditions_control);
            applicable_control_condition_labels = conditions_control(applicable_control_condition_index, :);
            
            % determine indicator for control
            stance_foot_indicator = strcmp(condition_stance_foot_list_session, applicable_control_condition_labels{1});
            perturbation_indicator = strcmp(condition_perturbation_list_session, applicable_control_condition_labels{2});
            delay_indicator = strcmp(condition_delay_list_session, applicable_control_condition_labels{3});
            index_indicator = strcmp(condition_index_list_session, applicable_control_condition_labels{4});
            experimental_indicator = strcmp(condition_experimental_list_session, applicable_control_condition_labels{5});
            stimulus_indicator = strcmp(condition_stimulus_list_session, applicable_control_condition_labels{6});
            day_indicator = strcmp(condition_day_list_session, applicable_control_condition_labels{7});
            this_condition_control_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
            
            % calculate responses
            for i_variable = 1 : number_of_stretch_variables
                % calculate control mean
                data_this_variable = stretch_data_session{i_variable};
                this_condition_control_data = data_this_variable(:, this_condition_control_indicator);
                this_condition_control_mean = mean(this_condition_control_data, 2);
                
                % calculate response
                response_data_session{i_variable}(:, i_stretch) = stretch_data_session{i_variable}(:, i_stretch) - this_condition_control_mean;
            end
            
        end

    end

    %% calculate integrated variables
    variables_to_integrate = study_settings.get('analysis_variables_from_integration');
    step_time_index_in_saved_data = find(strcmp(stretch_names_session, 'step_time'), 1, 'first');
    this_step_time_data = stretch_data_session{step_time_index_in_saved_data};
    pushoff_time_index_in_saved_data = find(strcmp(stretch_names_session, 'pushoff_time'), 1, 'first');
    this_pushoff_time_data = stretch_data_session{pushoff_time_index_in_saved_data};
    for i_variable = 1 : size(variables_to_integrate, 1)
        this_variable_name = variables_to_integrate{i_variable, 1};
        this_variable_source_name = variables_to_integrate{i_variable, 2};
        this_variable_response_data = response_data_session{strcmp(stretch_names_session, this_variable_source_name)};
        number_of_stretches = size(this_variable_response_data, 2);
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
        
        % add to analyzed variables
        analysis_data_session = [analysis_data_session; integrated_data]; %#ok<AGROW>
        analysis_names_session = [analysis_names_session; this_variable_name]; %#ok<AGROW>
    end
    
    %% calculate step end variables
    variables_step_end = study_settings.get('analysis_variables_from_step_end');
    for i_variable = 1 : size(variables_step_end, 1)
        this_variable_name = variables_step_end{i_variable, 1};
        this_variable_source_name = variables_step_end{i_variable, 2};
        this_variable_response_data = response_data_session{strcmp(stretch_names_session, this_variable_source_name)};
        step_end_data = this_variable_response_data(end, :);
        
        % add to analyzed variables
        analysis_data_session = [analysis_data_session; step_end_data]; %#ok<AGROW>
        analysis_names_session = [analysis_names_session; this_variable_name]; %#ok<AGROW>
    end




    %% save data
    save ...
      ( ...
        results_file_name, ...
        'stretch_data_session', ...
        'stretch_names_session', ...
        'response_data_session', ...
        'response_names_session', ...
        'analysis_data_session', ...
        'analysis_names_session', ...
        'condition_stance_foot_list_session', ...
        'condition_perturbation_list_session', ...
        'condition_delay_list_session', ...
        'condition_index_list_session', ...
        'condition_experimental_list_session', ...
        'condition_stimulus_list_session', ...
        'condition_day_list_session', ...
        'origin_trial_list_session', ...
        'origin_start_time_list_session', ...
        'origin_end_time_list_session', ...
        'time_list_session' ...
      )
end













