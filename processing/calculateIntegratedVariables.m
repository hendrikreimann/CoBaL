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

% input
% ... results.mat for each subject
% 

 %#ok<*AGROW>
 
function calculateIntegratedVariables(varargin)
    
    %% load and extract data
    study_settings_file = '';
    if exist('studySettings.txt', 'file')
        study_settings_file = 'studySettings.txt';
    end    
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    
    load('subjectInfo.mat', 'date', 'subject_id');
    load(['analysis' filesep date '_' subject_id '_results.mat']);

    variables_to_integrate = study_settings.get('variables_to_integrate');
    step_time_index_in_saved_data = find(strcmp(variable_names_session, 'step_time'), 1, 'first'); %#ok<NODEF>
    this_step_time_data = variable_data_session{step_time_index_in_saved_data}; %#ok<NODEF>
    pushoff_time_index_in_saved_data = find(strcmp(variable_names_session, 'pushoff_time'), 1, 'first');
    this_pushoff_time_data = variable_data_session{pushoff_time_index_in_saved_data};
    for i_variable = 1 : length(variables_to_integrate)
        this_variable_name = variables_to_integrate{i_variable};
        this_variable_response_data = response_data_session{strcmp(variable_names_session, this_variable_name)};
        number_of_stretches = size(this_variable_response_data, 2);
        integrated_data = zeros(1, number_of_stretches);
        for i_stretch = 1 : number_of_stretches
            % get data for full step
            this_stretch_time_full = linspace(0, this_step_time_data(i_stretch), 100);
            this_stretch_data_full = this_variable_response_data(:, i_stretch);
            
            % interpolate single stance to 100 data points
            this_stretch_time_single = linspace(this_pushoff_time_data(i_stretch), this_step_time_data(i_stretch), 100);
%             this_stretch_time_single = linspace(0.4, 0.5, 50);
            this_stretch_data_single = interp1(this_stretch_time_full, this_stretch_data_full, this_stretch_time_single);
            
            % integrate data in single stance
            this_stretch_data_single_integrated = cumtrapz(this_stretch_time_single, this_stretch_data_single);
            integrated_data(i_stretch) = this_stretch_data_single_integrated(end);
            
        end
        
        % add to saved data
        integrated_variable_name = [this_variable_name '_single_stance_integrated'];
        integrated_variable_index = find(strcmp(variable_names_session, integrated_variable_name));
        if isempty(integrated_variable_index)
            % stimulus_response_x doesn't exist yet, so add it to end
            variable_data_session = [variable_data_session; NaN];
            response_data_session = [response_data_session; integrated_data];
            variable_names_session = [variable_names_session; integrated_variable_name];
        end
        if ~isempty(integrated_variable_index)
            % stimulus_response_x exists, so overwrite this entry
            variable_data_session{integrated_variable_index} = NaN;
            response_data_session{integrated_variable_index} = integrated_data;
        end        
    end

    results_file_name = ['analysis' filesep makeFileName(date, subject_id, 'results')];
    save ...
      ( ...
        results_file_name, ...
        'variable_data_session', ...
        'response_data_session', ...
        'integrated_data_session', ...
        'variable_names_session', ...
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












