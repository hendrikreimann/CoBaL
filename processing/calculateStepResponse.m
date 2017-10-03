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

function calculateStepResponse(varargin)
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;
    
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

    conditions_control = study_settings.get('conditions_control');
    conditions_to_analyze = study_settings.get('conditions_to_analyze');
    step_placement_x_data = variable_data_session{strcmp(variable_names_session, 'step_placement_x')};
    mpsis_x_data = variable_data_session{strcmp(variable_names_session, 'mpsis_x')};
    com_x_data = variable_data_session{strcmp(variable_names_session, 'com_x')};
    mpsis_x_vel_data = deriveByTime(mpsis_x_data, 0.01); % XXX ACHTUNG: this is a hack because I don't have the data yet
    com_x_vel_data = variable_data_session{strcmp(variable_names_session, 'com_x_vel')};
%     lheelx_data = variable_data_session{strcmp(variable_names_session, 'lheel_x')};
%     rheelx_data = variable_data_session{strcmp(variable_names_session, 'rheel_x')};
    lanklex_data = variable_data_session{strcmp(variable_names_session, 'lankle_x')};
    ranklex_data = variable_data_session{strcmp(variable_names_session, 'rankle_x')};
    cop_x_data = variable_data_session{strcmp(variable_names_session, 'cop_x')};
    midstance_index_data = ones(1, size(mpsis_x_data, 2)) * 65; % XXX ACHTUNG: this is a hack because I don't have the midstance index data yet
    midstance_index_data = variable_data_session{strcmp(variable_names_session, 'midstance_index')};
    
    stimulus_response_x_data = zeros(1, size(mpsis_x_data, 2));

    %% estimate step response model parameters
    Jacobians = cell(1, size(conditions_control, 1));
    correlations_c = cell(1, size(conditions_control, 1));
    correlations_p = cell(1, size(conditions_control, 1));
    step_placement_x_means = zeros(1, size(conditions_control, 1));
    com_from_ankle_x_midstance_means = zeros(1, size(conditions_control, 1));
    com_x_vel_midstance_means = zeros(1, size(conditions_control, 1));
    
    for i_condition = 1 : size(conditions_control, 1)
        % determine stance foot
        stance_ankle_x_data = [];
        if strcmp(conditions_control{i_condition, 1}, 'STANCE_LEFT')
            stance_ankle_x_data = lanklex_data;
        end
        if strcmp(conditions_control{i_condition, 1}, 'STANCE_RIGHT')
            stance_ankle_x_data = ranklex_data;
        end
        
        % get control indicators
        stance_foot_indicator = strcmp(condition_stance_foot_list_session, conditions_control{i_condition, 1});
        perturbation_indicator = strcmp(condition_perturbation_list_session, conditions_control{i_condition, 2});
        delay_indicator = strcmp(condition_delay_list_session, conditions_control{i_condition, 3});
        index_indicator = strcmp(condition_index_list_session, conditions_control{i_condition, 4});
        experimental_indicator = strcmp(condition_experimental_list_session, conditions_control{i_condition, 5});
        stimulus_indicator = strcmp(condition_stimulus_list_session, conditions_control{i_condition, 6});
        day_indicator = strcmp(condition_day_list_session, conditions_control{i_condition, 7});
        this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
        
        % extract condition data
        step_placement_x_this_condition = step_placement_x_data(:, this_condition_indicator);
        mpsis_x_this_condition = mpsis_x_data(:, this_condition_indicator);
        com_x_this_condition = com_x_data(:, this_condition_indicator);
        mpsis_x_vel_this_condition = mpsis_x_vel_data(:, this_condition_indicator);
        com_x_vel_this_condition = com_x_vel_data(:, this_condition_indicator);
        stance_ankle_x_this_condition = stance_ankle_x_data(:, this_condition_indicator);
        cop_x_this_condition = cop_x_data(:, this_condition_indicator);
        
        % extract midstance data
        midstance_index_data_this_condition = midstance_index_data(this_condition_indicator);
        midstance_sub_indices = sub2ind(size(mpsis_x_this_condition), midstance_index_data_this_condition, 1:length(midstance_index_data_this_condition));
        
        mpsis_x_midstance = mpsis_x_this_condition(midstance_sub_indices);
        com_x_midstance = com_x_this_condition(midstance_sub_indices);
        stance_ankle_x_midstance = stance_ankle_x_this_condition(midstance_sub_indices);
        cop_x_midstance = cop_x_this_condition(midstance_sub_indices);
        mpsis_x_vel_midstance = mpsis_x_vel_this_condition(midstance_sub_indices);
        com_x_vel_midstance = com_x_vel_this_condition(midstance_sub_indices);
        
        mpsis_from_ankle_x_midstance = mpsis_x_midstance - stance_ankle_x_midstance;
        com_from_ankle_x_midstance = com_x_midstance - stance_ankle_x_midstance;
        mpsis_from_cop_x_midstance = mpsis_x_midstance - cop_x_midstance;
        com_from_cop_x_midstance = com_x_midstance - cop_x_midstance;
        
        % calculate means
        step_placement_x_means(i_condition) = mean(step_placement_x_this_condition);
        com_from_ankle_x_midstance_means(i_condition) = mean(com_from_ankle_x_midstance);
        com_x_vel_midstance_means(i_condition) = mean(com_x_vel_midstance);
        
        % calculate delta from mean
        step_placement_x_this_condition_delta = step_placement_x_this_condition - mean(step_placement_x_this_condition);
        mpsis_from_ankle_x_midstance_delta = mpsis_from_ankle_x_midstance - mean(mpsis_from_ankle_x_midstance);
        com_from_ankle_x_midstance_delta = com_from_ankle_x_midstance - mean(com_from_ankle_x_midstance);
        mpsis_from_cop_x_midstance_delta = mpsis_from_cop_x_midstance - mean(mpsis_from_cop_x_midstance);
        com_from_cop_x_midstance_delta = com_from_cop_x_midstance - mean(com_from_cop_x_midstance);
        mpsis_x_vel_midstance_delta = mpsis_x_vel_midstance - mean(mpsis_x_vel_midstance);
        com_x_vel_midstance_delta = com_x_vel_midstance - mean(com_x_vel_midstance);
        
        % calculate regression coefficients for com from ankle
        input_matrix = [com_from_ankle_x_midstance_delta; com_x_vel_midstance_delta];
        output_matrix = step_placement_x_this_condition_delta;
        Jacobian = output_matrix * pinv(input_matrix);
        [correlation_c, correlation_p] = corr(input_matrix', output_matrix');
        Jacobians{i_condition} = Jacobian;
        correlations_c{i_condition} = correlation_c;
        correlations_p{i_condition} = correlation_p;

        if visualize
            % calculate regression coefficients for mpsis from ankle
            input_matrix = [mpsis_from_ankle_x_midstance_delta; mpsis_x_vel_midstance_delta];
            output_matrix = step_placement_x_this_condition_delta;
            Jacobian = output_matrix * pinv(input_matrix);
            [correlation_c, correlation_p] = corr(input_matrix', output_matrix');

            figure; plot(mpsis_from_ankle_x_midstance_delta, step_placement_x_this_condition_delta, 'x')
            title(['mpsis from ankle - pos, J = ' num2str(Jacobian(1)) ', c = ' num2str(correlation_c(1)) ', p = ' num2str(correlation_p(1))]); axis equal
            figure; plot(mpsis_x_vel_midstance_delta, step_placement_x_this_condition_delta, 'x')
            title(['mpsis - vel, J = ' num2str(Jacobian(2)) ', c = ' num2str(correlation_c(2)) ', p = ' num2str(correlation_p(2))]); axis equal

            % calculate regression coefficients for com from ankle
            input_matrix = [com_from_ankle_x_midstance_delta; com_x_vel_midstance_delta];
            output_matrix = step_placement_x_this_condition_delta;
            Jacobian = output_matrix * pinv(input_matrix);
            [correlation_c, correlation_p] = corr(input_matrix', output_matrix');

            figure; plot(com_from_ankle_x_midstance_delta, step_placement_x_this_condition_delta, 'x')
            title(['com from ankle - pos, J = ' num2str(Jacobian(1)) ', c = ' num2str(correlation_c(1)) ', p = ' num2str(correlation_p(1))]); axis equal
            figure; plot(com_x_vel_midstance_delta, step_placement_x_this_condition_delta, 'x')
            title(['com - vel, J = ' num2str(Jacobian(2)) ', c = ' num2str(correlation_c(2)) ', p = ' num2str(correlation_p(2))]); axis equal

            % calculate regression coefficients for mpsis from cop
            input_matrix = [mpsis_from_cop_x_midstance_delta; mpsis_x_vel_midstance_delta];
            output_matrix = step_placement_x_this_condition_delta;
            Jacobian = output_matrix * pinv(input_matrix);
            [correlation_c, correlation_p] = corr(input_matrix', output_matrix');

            figure; plot(mpsis_from_cop_x_midstance_delta, step_placement_x_this_condition_delta, 'x')
            title(['mpsis from cop - pos, J = ' num2str(Jacobian(1)) ', c = ' num2str(correlation_c(1)) ', p = ' num2str(correlation_p(1))]); axis equal
            figure; plot(mpsis_x_vel_midstance_delta, step_placement_x_this_condition_delta, 'x')
            title(['mpsis from cop - vel, J = ' num2str(Jacobian(2)) ', c = ' num2str(correlation_c(2)) ', p = ' num2str(correlation_p(2))]); axis equal

            % calculate regression coefficients for com from cop
            input_matrix = [com_from_cop_x_midstance_delta; com_x_vel_midstance_delta];
            output_matrix = step_placement_x_this_condition_delta;
            Jacobian = output_matrix * pinv(input_matrix);
            [correlation_c, correlation_p] = corr(input_matrix', output_matrix');

            figure; plot(com_from_cop_x_midstance_delta, step_placement_x_this_condition_delta, 'x')
            title(['com from cop - pos, J = ' num2str(Jacobian(1)) ', c = ' num2str(correlation_c(1)) ', p = ' num2str(correlation_p(1))]); axis equal
            figure; plot(com_x_vel_midstance_delta, step_placement_x_this_condition_delta, 'x')
            title(['com - vel, J = ' num2str(Jacobian(2)) ', c = ' num2str(correlation_c(2)) ', p = ' num2str(correlation_p(2))]); axis equal
        end
        
        
    end
    
    %% estimate stimulus response
    all_conditions = [conditions_to_analyze; conditions_control];
    for i_condition = 1 : size(all_conditions, 1)
        
        % determine stance foot
        stance_ankle_x_data = [];
        if strcmp(all_conditions{i_condition, 1}, 'STANCE_LEFT')
            stance_ankle_x_data = lanklex_data;
            applicable_control_condition = find(strcmp(conditions_control(:, 1), 'STANCE_LEFT'));
        end
        if strcmp(all_conditions{i_condition, 1}, 'STANCE_RIGHT')
            stance_ankle_x_data = ranklex_data;
            applicable_control_condition = find(strcmp(conditions_control(:, 1), 'STANCE_RIGHT'));
        end
        
        % get control indicators
        stance_foot_indicator = strcmp(condition_stance_foot_list_session, all_conditions{i_condition, 1});
        perturbation_indicator = strcmp(condition_perturbation_list_session, all_conditions{i_condition, 2});
        delay_indicator = strcmp(condition_delay_list_session, all_conditions{i_condition, 3});
        index_indicator = strcmp(condition_index_list_session, all_conditions{i_condition, 4});
        experimental_indicator = strcmp(condition_experimental_list_session, all_conditions{i_condition, 5});
        stimulus_indicator = strcmp(condition_stimulus_list_session, all_conditions{i_condition, 6});
        day_indicator = strcmp(condition_day_list_session, all_conditions{i_condition, 7});
        this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
        
        % extract condition data
        step_placement_x_this_condition = step_placement_x_data(:, this_condition_indicator);
        mpsis_x_this_condition = mpsis_x_data(:, this_condition_indicator);
        com_x_this_condition = com_x_data(:, this_condition_indicator);
        com_x_vel_this_condition = com_x_vel_data(:, this_condition_indicator);
        stance_ankle_x_this_condition = stance_ankle_x_data(:, this_condition_indicator);
        
        % extract midstance data
        midstance_index_data_this_condition = midstance_index_data(this_condition_indicator);
        midstance_sub_indices = sub2ind(size(mpsis_x_this_condition), midstance_index_data_this_condition, 1:length(midstance_index_data_this_condition));
        
        com_x_midstance = com_x_this_condition(midstance_sub_indices);
        stance_ankle_x_midstance = stance_ankle_x_this_condition(midstance_sub_indices);
        com_x_vel_midstance = com_x_vel_this_condition(midstance_sub_indices);
        
        com_from_ankle_x_midstance = com_x_midstance - stance_ankle_x_midstance;
        
        % calculate delta from control mean
        step_placement_x_this_condition_delta = step_placement_x_this_condition - step_placement_x_means(applicable_control_condition);
        com_from_ankle_x_midstance_delta = com_from_ankle_x_midstance - com_from_ankle_x_midstance_means(applicable_control_condition);
        com_x_vel_midstance_delta = com_x_vel_midstance - com_x_vel_midstance_means(applicable_control_condition); 
        
        % calculate and remove expected part
        step_response_x = step_placement_x_this_condition_delta;
        Jacobian = Jacobians{applicable_control_condition};
        expected_response_x = Jacobian(1) * com_from_ankle_x_midstance_delta + Jacobian(2) * com_x_vel_midstance_delta;
        stimulus_response_x = step_response_x - expected_response_x;
        
        stimulus_response_x_data(this_condition_indicator) = stimulus_response_x;
        
%         hold on; plot(stimulus_response, 'x')
%         mean(stimulus_response)
%         mean(step_response)
    end
    
    
    % save data
    stimulus_response_index = find(strcmp(variable_names_session, 'stimulus_response_x'));
    if isempty(stimulus_response_index)
        % stimulus_response_x doesn't exist yet, so add it to end
        variable_data_session = [variable_data_session; stimulus_response_x_data];
        response_data_session = [response_data_session; stimulus_response_x_data];
        variable_names_session = [variable_names_session; 'stimulus_response_x'];
    end
    if ~isempty(stimulus_response_index)
        % stimulus_response_x exists, so overwrite this entry
        variable_data_session{stimulus_response_index} = stimulus_response_x_data;
        response_data_session{stimulus_response_index} = stimulus_response_x_data;
        variable_names_session{stimulus_response_index} = 'stimulus_response_x';
    end
    
    results_file_name = ['analysis' filesep makeFileName(date, subject_id, 'results')];
    save ...
      ( ...
        results_file_name, ...
        'variable_data_session', ...
        'response_data_session', ...
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












