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
    
    number_of_stretches = length(loaded_data.time_list_session);
    stance_feet = {'STANCE_LEFT', 'STANCE_RIGHT'};
    
    % find relevant conditions
    conditions_session = loaded_data.conditions_session;
    conditions_settings = study_settings.get('conditions');
    condition_labels = conditions_settings(:, 1)';
    condition_source_variables = conditions_settings(:, 2)';
    number_of_condition_labels = length(condition_labels);
    
    % extract data
    step_placement_x_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'step_placement_x')};
    step_placement_x_directions = loaded_data.stretch_directions_session(strcmp(loaded_data.stretch_names_session, 'step_placement_x'), :);
    com_x_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'com_x')};
    com_x_vel_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'com_x_vel')};
    lankle_x_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'lankle_x')};
    rankle_x_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'rankle_x')};
    midstance_index_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'midstance_index')};
    number_of_bands_per_stretch = size(step_placement_x_data, 1);
    
    % make condition data tables
    condition_data_all = cell(number_of_stretches, number_of_condition_labels);
    for i_condition = 1 : number_of_condition_labels
        condition_data_all(:, i_condition) = conditions_session.(condition_source_variables{i_condition});
    end
    labels_to_ignore = study_settings.get('conditions_to_ignore');
    levels_to_remove = study_settings.get('levels_to_remove');
    [condition_combination_labels, condition_combinations_stimulus, condition_combinations_control] = determineConditionCombinations(condition_data_all, conditions_settings, labels_to_ignore, levels_to_remove);
    condition_combinations_control_unique = table2cell(unique(cell2table(condition_combinations_control), 'rows'));
    
    %% gather data for unperturbed steps from all control stretches
    number_of_conditions_control = size(condition_combinations_control_unique, 1);
    com_x_pos_midstance_control_data = [];
    com_x_vel_midstance_control_data = [];
    lankle_x_midstance_control_data = [];
    rankle_x_midstance_control_data = [];
    step_placement_x_control_data = [];
    stance_foot_control_data = {};
    
    for i_condition = 1 : number_of_conditions_control
        % get indicator for this control condition
        this_condition_combination = condition_combinations_control_unique(i_condition, :);
        this_condition_indicator = true(number_of_stretches, 1);
        for i_label = 1 : length(condition_combination_labels)
            this_label = condition_combination_labels{i_label};
            this_label_list = condition_data_all(:, strcmp(conditions_settings(:, 1), this_label));
            this_label_indicator = strcmp(this_label_list, this_condition_combination{i_label});
            this_condition_indicator = this_condition_indicator .* this_label_indicator;
        end
        this_condition_indices = find(this_condition_indicator);
        
        % get data
        for i_stretch = 1 : length(this_condition_indices)
            % extract data for this stretch
            this_stretch_stance_foot_data = conditions_session.stance_foot_data(this_condition_indices(i_stretch), :);
            this_condition_midstance_index_data = midstance_index_data(:, this_condition_indices(i_stretch));
            this_stretch_com_x_pos_data = com_x_data(:, this_condition_indices(i_stretch));
            this_stretch_com_x_vel_data = com_x_vel_data(:, this_condition_indices(i_stretch));
            this_stretch_lankle_x_data = lankle_x_data(:, this_condition_indices(i_stretch));
            this_stretch_rankle_x_data = rankle_x_data(:, this_condition_indices(i_stretch));
            this_stretch_step_placement_x_data = step_placement_x_data(:, this_condition_indices(i_stretch));
            
            % go through bands and store data if the band is a swing
            for i_band = 1 : number_of_bands_per_stretch
                if ~isnan(this_stretch_step_placement_x_data(i_band))
                    stance_foot_control_data = [stance_foot_control_data; this_stretch_stance_foot_data{i_band}]; %#ok<AGROW>
                    com_x_pos_midstance_control_data = [com_x_pos_midstance_control_data; this_stretch_com_x_pos_data(this_condition_midstance_index_data(i_band))]; %#ok<AGROW>
                    com_x_vel_midstance_control_data = [com_x_vel_midstance_control_data; this_stretch_com_x_vel_data(this_condition_midstance_index_data(i_band))]; %#ok<AGROW>
                    lankle_x_midstance_control_data = [lankle_x_midstance_control_data; this_stretch_lankle_x_data(this_condition_midstance_index_data(i_band))]; %#ok<AGROW>
                    rankle_x_midstance_control_data = [rankle_x_midstance_control_data; this_stretch_rankle_x_data(this_condition_midstance_index_data(i_band))]; %#ok<AGROW>
                    step_placement_x_control_data = [step_placement_x_control_data; this_stretch_step_placement_x_data(i_band)]; %#ok<AGROW>
                end
            end
        end
    end    
%     figure; axes; hold on; set(gca, 'ylim', [-0.1, 0.1])
%     plot(com_x_pos_midstance_control_data(strcmp(stance_foot_control_data, 'STANCE_RIGHT')), 'rx')
%     plot(com_x_pos_midstance_control_data(strcmp(stance_foot_control_data, 'STANCE_LEFT')), 'gx')
% 
%     figure; axes; hold on; set(gca, 'ylim', [-0.1, 0.1])
%     plot(com_x_vel_midstance_control_data(strcmp(stance_foot_control_data, 'STANCE_RIGHT')), 'rx')
%     plot(com_x_vel_midstance_control_data(strcmp(stance_foot_control_data, 'STANCE_LEFT')), 'gx')
% 
%     figure; axes; hold on; set(gca, 'ylim', [-0.25, 0.25])
%     plot(step_placement_x_control_data(strcmp(stance_foot_control_data, 'STANCE_RIGHT')), 'rx')
%     plot(step_placement_x_control_data(strcmp(stance_foot_control_data, 'STANCE_LEFT')), 'gx')
% 
%     figure; axes; hold on; set(gca, 'ylim', [-0.25, 0.25])
%     plot(lankle_x_midstance_control_data(strcmp(stance_foot_control_data, 'STANCE_RIGHT')), 'rx')
%     plot(lankle_x_midstance_control_data(strcmp(stance_foot_control_data, 'STANCE_LEFT')), 'gx')
% 
%     figure; axes; hold on; set(gca, 'ylim', [-0.25, 0.25])
%     plot(rankle_x_midstance_control_data(strcmp(stance_foot_control_data, 'STANCE_RIGHT')), 'rx')
%     plot(rankle_x_midstance_control_data(strcmp(stance_foot_control_data, 'STANCE_LEFT')), 'gx')

    %% fit linear model
    stance_styles = {'gx', 'rx'};
    linear_models = cell(size(stance_feet));
    com_from_ankle_means = zeros(size(stance_feet));
    com_vel_means = zeros(size(stance_feet));
    foot_placement_means = zeros(size(stance_feet));
    for i_stance = 1 : length(stance_feet)
        % get data for this stance foot
        this_stance_foot_indicator = strcmp(stance_foot_control_data, stance_feet{i_stance});
        this_stance_com_x_pos_midstance_data = com_x_pos_midstance_control_data(this_stance_foot_indicator);
        this_stance_com_x_vel_midstance_data = com_x_vel_midstance_control_data(this_stance_foot_indicator);
        this_stance_step_placement_x_data = step_placement_x_control_data(this_stance_foot_indicator);
        if strcmp(stance_feet{i_stance}, 'STANCE_LEFT')
            this_stance_foot_ankle_x_data = lankle_x_midstance_control_data(this_stance_foot_indicator);
        end
        if strcmp(stance_feet{i_stance}, 'STANCE_RIGHT')
            this_stance_foot_ankle_x_data = rankle_x_midstance_control_data(this_stance_foot_indicator);
        end
        this_stance_com_from_stance_ankle_data = this_stance_com_x_pos_midstance_data - this_stance_foot_ankle_x_data;
        
        % calculate and remove means
        com_from_ankle_means(i_stance) = mean(this_stance_com_from_stance_ankle_data);
        com_vel_means(i_stance) = mean(this_stance_com_x_vel_midstance_data);
        foot_placement_means(i_stance) = mean(this_stance_step_placement_x_data);
        this_stance_com_from_stance_ankle_data_mean_free = this_stance_com_from_stance_ankle_data - mean(this_stance_com_from_stance_ankle_data);
        this_stance_com_x_vel_midstance_data_mean_free = this_stance_com_x_vel_midstance_data - mean(this_stance_com_x_vel_midstance_data);
        this_stance_step_placement_x_data_mean_free = this_stance_step_placement_x_data - mean(this_stance_step_placement_x_data);
        
        % remove nans
        
        
        % fit regression model
%         input_matrix = [this_stance_com_from_stance_ankle_data_mean_free'; this_stance_com_x_vel_midstance_data_mean_free'];
%         output_matrix = this_stance_step_placement_x_data_mean_free';
%         Jacobian = output_matrix * pinv(input_matrix);
%         [correlation_c, correlation_p] = corr(input_matrix', output_matrix');
%         Jacobians{i_stance} = Jacobian;
%         correlations_c{i_stance} = correlation_c;
%         correlations_p{i_stance} = correlation_p;
        [fit_object, fit_stats] = fit([this_stance_com_from_stance_ankle_data_mean_free, this_stance_com_x_vel_midstance_data_mean_free], this_stance_step_placement_x_data_mean_free, 'poly11');
        linear_models{i_stance} = [fit_object.p10, fit_object.p01];

%         figure; axes; hold on; set(gca, 'ylim', [-0.1, 0.1])
%         plot(this_stance_com_from_stance_ankle_data, 'bx')
%         plot(this_stance_com_from_stance_ankle_data_mean_free, stance_styles{i_stance})
%         
%         figure; axes; hold on; set(gca, 'ylim', [-0.1, 0.1])
%         plot(this_stance_com_x_vel_midstance_data_mean_free, 'bx')
%         plot(this_stance_com_x_vel_midstance_data, stance_styles{i_stance})
%         
%         figure; axes; hold on; set(gca, 'ylim', [-0.25, 0.25])
%         plot(this_stance_step_placement_x_data_mean_free, 'bx')
%         plot(this_stance_step_placement_x_data, stance_styles{i_stance})
        
%         figure; axes; hold on; 
%         plot(this_stance_com_from_stance_ankle_data, 'bx')
%         plot(this_stance_com_x_pos_midstance_data, 'mo')
%         plot(this_stance_foot_ankle_x_data, 'co')
        

        if visualize
            figure; 
            plot(fit_object, [this_stance_com_from_stance_ankle_data_mean_free, this_stance_com_x_vel_midstance_data_mean_free], this_stance_step_placement_x_data_mean_free)
            xlabel('\Delta CoM from stance ankle'); ylabel('\Delta CoM vel'); zlabel('\Delta foot placement')
            title(['slopes = [' num2str(linear_model_slopes(1)) ', ' num2str(linear_model_slopes(2)) '], r^2 = ' num2str(fit_stats.rsquare)]);
            
%             figure; hold on; 
%             plot(this_stance_com_from_stance_ankle_data_mean_free, this_stance_step_placement_x_data_mean_free, 'o')
%             plot([-0.04, 0.04], [0, 0], 'color', [1 1 1]*0.5);
%             plot([0, 0], [-0.15, 0.15], 'color', [1 1 1]*0.5);
%             xlabel('\Delta CoM at midstance')
%             ylabel('\Delta foot placement')
%             title(['J = ' num2str(Jacobian(1)) ', c = ' num2str(correlation_c(1)) ', p = ' num2str(correlation_p(1))]);
%             title(['com from ankle - pos, J = ' num2str(Jacobian(1)) ', c = ' num2str(correlation_c(1)) ', p = ' num2str(correlation_p(1))]);
%             set(gca, 'xlim', [-0.04, 0.04], 'ylim', [-0.15, 0.15])
%             
%             figure; plot(this_stance_com_x_vel_midstance_data_mean_free, this_stance_step_placement_x_data_mean_free, 'x')
%             title(['com - vel, J = ' num2str(Jacobian(2)) ', c = ' num2str(correlation_c(2)) ', p = ' num2str(correlation_p(2))]); axis equal
%             
%             figure; plot3(this_stance_com_from_stance_ankle_data_mean_free, this_stance_com_x_vel_midstance_data_mean_free, this_stance_step_placement_x_data_mean_free, 'x')
        end    

    end

    %% for each data point, calculate difference from prediction
    stimulus_response_x_data = zeros(size(step_placement_x_data)) * NaN;
    stimulus_response_x_directions = step_placement_x_directions;
    for i_stretch = 1 : number_of_stretches
        this_stretch_stance_foot_data = conditions_session.stance_foot_data(i_stretch, :);
        this_condition_midstance_index_data = midstance_index_data(:, i_stretch);
        this_stretch_com_x_pos_data = com_x_data(:, i_stretch);
        this_stretch_com_x_vel_data = com_x_vel_data(:, i_stretch);
        this_stretch_lankle_x_data = lankle_x_data(:, i_stretch);
        this_stretch_rankle_x_data = rankle_x_data(:, i_stretch);
        this_stretch_step_placement_x_data = step_placement_x_data(:, i_stretch);

        % go through bands
        for i_band = 1 : number_of_bands_per_stretch
            if ~isnan(this_stretch_step_placement_x_data(i_band))
                % get data
                foot_index = find(strcmp(stance_feet, this_stretch_stance_foot_data{i_band}));
                com_x_pos_midstance_data = this_stretch_com_x_pos_data(this_condition_midstance_index_data(i_band));
                com_x_vel_midstance_data = this_stretch_com_x_vel_data(this_condition_midstance_index_data(i_band));
                lankle_x_midstance_data = this_stretch_lankle_x_data(this_condition_midstance_index_data(i_band));
                rankle_x_midstance_data = this_stretch_rankle_x_data(this_condition_midstance_index_data(i_band));
                step_placement_x_here = this_stretch_step_placement_x_data(i_band);
                if strcmp(this_stretch_stance_foot_data{i_band}, 'STANCE_LEFT')
                    this_stance_foot_ankle_x_data = lankle_x_midstance_data;
                end
                if strcmp(this_stretch_stance_foot_data{i_band}, 'STANCE_RIGHT')
                    this_stance_foot_ankle_x_data = rankle_x_midstance_data;
                end
                com_from_stance_ankle_data = com_x_pos_midstance_data - this_stance_foot_ankle_x_data;
                
                % calculated deltas
                com_from_stance_ankle_delta = com_from_stance_ankle_data - com_from_ankle_means(foot_index);
                com_vel_delta = com_x_vel_midstance_data - com_vel_means(foot_index);
                foot_placement_delta = step_placement_x_here - foot_placement_means(foot_index);
                model_slopes = linear_models{i_stance};
                
                predicted_foot_placement_change = com_from_stance_ankle_delta * model_slopes(1) + com_vel_delta * model_slopes(2);
                actual_foot_placement_change = foot_placement_delta;
                stimulus_response = actual_foot_placement_change - predicted_foot_placement_change;

                stimulus_response_x_data(i_band, i_stretch) = stimulus_response;
                
                
            end
        end
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
    [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
        addOrReplaceResultsData ...
          ( ...
            analysis_data_session, ...
            analysis_names_session, ...
            analysis_directions_session, ...
            stimulus_response_x_data, ...
            'stimulus_response_x',...
            stimulus_response_x_directions ...
          );
    
    variables_to_save = loaded_data;
    variables_to_save.analysis_data_session = analysis_data_session;
    variables_to_save.analysis_names_session = analysis_names_session;
    variables_to_save.analysis_directions_session = analysis_directions_session;
    save(results_file_name, '-struct', 'variables_to_save');

    return % stuff below this point is old
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    %% estimate step response model parameters
    number_of_conditions_control = size(condition_combinations_control_unique, 1);
    Jacobians = cell(1, number_of_conditions_control);
    correlations_c = cell(1, number_of_conditions_control);
    correlations_p = cell(1, number_of_conditions_control);
    step_placement_x_means = zeros(1, number_of_conditions_control);
    com_from_ankle_x_midstance_means = zeros(1, number_of_conditions_control);
    com_x_vel_midstance_means = zeros(1, number_of_conditions_control);
    
    % TODO: make this work for bands
    for i_condition = 1 : number_of_conditions_control
        % get relevant control condition
        this_condition_combination = condition_combinations_control_unique(i_condition, :);
        this_condition_indicator = true(number_of_stretches, 1);
        for i_label = 1 : length(condition_combination_labels)
            this_label = condition_combination_labels{i_label};
            this_label_list = condition_data_all(:, strcmp(conditions_settings(:, 1), this_label));
            this_label_indicator = strcmp(this_label_list, this_condition_combination{i_label});
            this_condition_indicator = this_condition_indicator .* this_label_indicator;
        end
        condition_representant_index = find(this_condition_indicator, 1);
        this_condition_indicator = logical(this_condition_indicator);            
%         this_condition_stance_foot_data = conditions_session.stance_foot_data(condition_representant_index, :);
        this_condition_stance_foot_data = conditions_session.condition_stance_foot_list(condition_representant_index, :);
        
        % determine stance ankle data
        stance_ankle_x_data = zeros(size(lankle_x_data));
        if strcmp(this_condition_stance_foot_data, 'STANCE_LEFT')
            stance_ankle_x_data = lankle_x_data;
        end
        if strcmp(this_condition_stance_foot_data, 'STANCE_RIGHT')
            stance_ankle_x_data = rankle_x_data;
        end
        
        % extract condition data
        step_placement_x_this_condition = step_placement_x_data(:, this_condition_indicator);
        mpsis_x_this_condition = mpsis_x_data(:, this_condition_indicator);
        com_x_this_condition = com_x_data(:, this_condition_indicator);
        mpsis_x_vel_this_condition = mpsis_x_vel_data(:, this_condition_indicator);
        com_x_vel_this_condition = com_x_vel_data(:, this_condition_indicator);
        stance_ankle_x_this_condition = stance_ankle_x_data(:, this_condition_indicator);
        cop_x_this_condition = cop_x_data(:, this_condition_indicator);
        
        
        
        
        
        
        
        
        % determine stance foot
        stance_ankle_x_data = [];
        if strcmp(condition_combinations_control_unique{i_condition, strcmp(condition_combination_labels, 'stance_foot')}, 'STANCE_LEFT')
            stance_ankle_x_data = lankle_x_data;
        end
        if strcmp(condition_combinations_control_unique{i_condition, strcmp(condition_combination_labels, 'stance_foot')}, 'STANCE_RIGHT')
            stance_ankle_x_data = rankle_x_data;
        end
        
        % get condition indicators
        this_condition_indicator = true(number_of_stretches, 1);
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

            figure; hold on;
            plot(com_from_ankle_x_midstance_delta, step_placement_x_this_condition_delta, 'o')
            plot([-0.04, 0.04], [0, 0], 'color', [1 1 1]*0.5);
            plot([0, 0], [-0.15, 0.15], 'color', [1 1 1]*0.5);
            xlabel('\Delta CoM at midstance')
            ylabel('\Delta foot placement')
            title(['J = ' num2str(Jacobian(1)) ', c = ' num2str(correlation_c(1)) ', p = ' num2str(correlation_p(1))]);
%             title(['com from ankle - pos, J = ' num2str(Jacobian(1)) ', c = ' num2str(correlation_c(1)) ', p = ' num2str(correlation_p(1))]);
            set(gca, 'xlim', [-0.04, 0.04], 'ylim', [-0.15, 0.15])
            
            
            axis equal
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
            stance_ankle_x_data = lankle_x_data;
            applicable_control_condition = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'stance_foot')), 'STANCE_LEFT'));
        end
        if strcmp(all_conditions{i_condition, strcmp(condition_combination_labels, 'stance_foot')}, 'STANCE_RIGHT')
            stance_ankle_x_data = rankle_x_data;
            applicable_control_condition = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'stance_foot')), 'STANCE_RIGHT'));
        end
        
        % get condition indicator
        this_condition_indicator = true(number_of_stretches, 1);
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
    
    return % temporary so I won't accidentally save over data
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
    [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
        addOrReplaceResultsData ...
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












