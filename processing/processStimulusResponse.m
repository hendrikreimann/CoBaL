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

function processStimulusResponse(varargin)
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
    
    results_file_name = ['analysis' filesep makeFileName(date, subject_id, 'results')];
    loaded_data = load(results_file_name);
    
    % find relevant conditions
    conditions_session = loaded_data.conditions_session;
    conditions_settings = study_settings.get('conditions');
    condition_labels = conditions_settings(:, 1)';
    condition_source_variables = conditions_settings(:, 2)';
    number_of_condition_labels = length(condition_labels);
    
    % extract data
    step_placement_x_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'step_placement_x')};
    step_placement_x_directions = loaded_data.stretch_directions_session(strcmp(loaded_data.stretch_names_session, 'step_placement_x'), :);
    mpsis_x_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'mpsis_x')};
    com_x_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'com_x')};
    mpsis_x_vel_data = deriveByTime(mpsis_x_data, 0.01); % XXX ACHTUNG: this is a hack because I don't have the data yet
    com_x_vel_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'com_x_vel')};
    lanklex_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'lankle_x')};
    ranklex_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'rankle_x')};
    cop_x_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'cop_x')};
    midstance_index_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'midstance_index')};
    
    stimulus_response_x_data = zeros(1, size(mpsis_x_data, 2));
    number_of_stretches_session = length(step_placement_x_data);
    
    % make condition data tables
    condition_data_all = cell(number_of_stretches_session, number_of_condition_labels);
    for i_condition = 1 : number_of_condition_labels
        condition_data_all(:, i_condition) = conditions_session.(condition_source_variables{i_condition});
    end
    labels_to_ignore = study_settings.get('conditions_to_ignore');
    levels_to_remove = study_settings.get('levels_to_remove');
    [condition_combination_labels, condition_combinations_stimulus, condition_combinations_control] = determineConditionCombinations(condition_data_all, conditions_settings, labels_to_ignore, levels_to_remove);
    condition_combinations_control_unique = table2cell(unique(cell2table(condition_combinations_control), 'rows'));
    
    

    %% estimate step response model parameters
    number_of_conditions_control = size(condition_combinations_control_unique, 1);
    Jacobians = cell(1, number_of_conditions_control);
    correlations_c = cell(1, number_of_conditions_control);
    correlations_p = cell(1, number_of_conditions_control);
    step_placement_x_means = zeros(1, number_of_conditions_control);
    com_from_ankle_x_midstance_means = zeros(1, number_of_conditions_control);
    com_x_vel_midstance_means = zeros(1, number_of_conditions_control);
    
    for i_condition = 1 : number_of_conditions_control
        % determine stance foot
        stance_ankle_x_data = [];
        if strcmp(condition_combinations_control_unique{i_condition, strcmp(condition_combination_labels, 'stance_foot')}, 'STANCE_LEFT')
            stance_ankle_x_data = lanklex_data;
        end
        if strcmp(condition_combinations_control_unique{i_condition, strcmp(condition_combination_labels, 'stance_foot')}, 'STANCE_RIGHT')
            stance_ankle_x_data = ranklex_data;
        end
        
        % get control indicators
        this_condition_indicator = true(number_of_stretches_session, 1);
        for i_label = 1 : length(condition_combination_labels)
            this_label = condition_combination_labels{i_label};
            this_label_list = condition_data_all(:, strcmp(conditions_settings(:, 1), this_label));
            this_label_indicator = strcmp(this_label_list, condition_combinations_control_unique(i_condition, i_label));
            this_condition_indicator = this_condition_indicator .* this_label_indicator;
        end        
        this_condition_indicator = logical(this_condition_indicator);
        
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

%             figure; plot(mpsis_from_ankle_x_midstance, step_placement_x_this_condition, 'x'); hold on, plot(0, 0);
%             title('mpsis from ankle - pos'); axis equal
%             figure; plot(mpsis_x_vel_midstance_delta, step_placement_x_this_condition, 'x'); hold on, plot(0, 0);
%             title('mpsis from ankle - vel'); axis equal
%             
%             figure; plot(mpsis_from_ankle_x_midstance_delta, step_placement_x_this_condition_delta, 'x')
%             title(['Delta mpsis from ankle - pos, J = ' num2str(Jacobian(1)) ', c = ' num2str(correlation_c(1)) ', p = ' num2str(correlation_p(1))]); axis equal
%             figure; plot(mpsis_x_vel_midstance_delta, step_placement_x_this_condition_delta, 'x')
%             title(['Delta mpsis - vel, J = ' num2str(Jacobian(2)) ', c = ' num2str(correlation_c(2)) ', p = ' num2str(correlation_p(2))]); axis equal

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
    all_conditions = [condition_combinations_stimulus; condition_combinations_control_unique];
    for i_condition = 1 : size(all_conditions, 1)
        
        % determine stance foot
        stance_ankle_x_data = [];
        if strcmp(all_conditions{i_condition, strcmp(condition_combination_labels, 'stance_foot')}, 'STANCE_LEFT')
            stance_ankle_x_data = lanklex_data;
            applicable_control_condition = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'stance_foot')), 'STANCE_LEFT'));
        end
        if strcmp(all_conditions{i_condition, strcmp(condition_combination_labels, 'stance_foot')}, 'STANCE_RIGHT')
            stance_ankle_x_data = ranklex_data;
            applicable_control_condition = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'stance_foot')), 'STANCE_RIGHT'));
        end
        
        % get condition indicator
        this_condition_indicator = true(number_of_stretches_session, 1);
        for i_label = 1 : length(condition_combination_labels)
            this_label = condition_combination_labels{i_label};
            this_label_list = condition_data_all(:, strcmp(conditions_settings(:, 1), this_label));
            this_label_indicator = strcmp(this_label_list, all_conditions(i_condition, i_label));
            this_condition_indicator = this_condition_indicator .* this_label_indicator;
        end        
        this_condition_indicator = logical(this_condition_indicator);
        
%         stance_foot_indicator = strcmp(loaded_data.condition_stance_foot_list_session, all_conditions{i_condition, 1});
%         perturbation_indicator = strcmp(loaded_data.condition_perturbation_list_session, all_conditions{i_condition, 2});
%         delay_indicator = strcmp(loaded_data.condition_delay_list_session, all_conditions{i_condition, 3});
%         index_indicator = strcmp(loaded_data.condition_index_list_session, all_conditions{i_condition, 4});
%         experimental_indicator = strcmp(loaded_data.condition_experimental_list_session, all_conditions{i_condition, 5});
%         stimulus_indicator = strcmp(loaded_data.condition_stimulus_list_session, all_conditions{i_condition, 6});
%         day_indicator = strcmp(loaded_data.condition_day_list_session, all_conditions{i_condition, 7});
%         this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
        
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
    
    % save
    if isfield(loaded_data, 'analysis_data_session')
        analysis_data_session = loaded_data.analysis_data_session;
        analysis_names_session = loaded_data.analysis_names_session;
        analysis_directions_session = loaded_data.analysis_directions_session;
    else
        analysis_data_session = {};
        analysis_names_session = {};
        analysis_directions_session = {};
    end
%     [analysis_data_session, analysis_names_session] = addOrOverwriteData(analysis_data_session, analysis_names_session, stimulus_response_x_data, 'stimulus_response_x');
    [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
        addOrOverwriteResultsData ...
          ( ...
            analysis_data_session, ...
            analysis_names_session, ...
            analysis_directions_session, ...
            stimulus_response_x_data, ...
            'stimulus_response_x',...
            step_placement_x_directions ...
          );
    
    variables_to_save = loaded_data;
    variables_to_save.analysis_data_session = analysis_data_session;
    variables_to_save.analysis_names_session = analysis_names_session;
    variables_to_save.analysis_directions_session = analysis_directions_session;
    save(results_file_name, '-struct', 'variables_to_save');
        
    
end












