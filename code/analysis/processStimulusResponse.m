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
    fx_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'fx')};
    fx_directions = loaded_data.stretch_directions_session(strcmp(loaded_data.stretch_names_session, 'fx'), :);
    lankle_angle_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'joint_angle:left_ankle_eversion')};
    rankle_angle_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'joint_angle:right_ankle_eversion')};
    
    cop_com_int_data = loaded_data.stretch_data_session{strcmp(loaded_data.analysis_names_session, 'cop_from_com_x_integrated')};
    cop_com_int_directions = loaded_data.stretch_directions_session(strcmp(loaded_data.analysis_names_session, 'cop_from_com_x_integrated'), :);
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
    com_x_pos_end_control_data = [];
    com_x_vel_end_control_data = [];
    com_x_pos_initial_control_data = [];
    com_x_vel_initial_control_data = [];
    lankle_x_midstance_control_data = [];
    rankle_x_midstance_control_data = [];
    lankle_x_end_control_data = [];
    rankle_x_end_control_data = [];
    lankle_x_initial_control_data = [];
    rankle_x_initial_control_data = [];
    lankle_angle_end_control_data = [];
    rankle_angle_end_control_data = [];
    step_placement_x_control_data = [];
    fx_control_data = [];
    cop_com_int_control_data = [];
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
            this_stretch_fx_data = fx_data(:, this_condition_indices(i_stretch)); 
            this_stretch_lankle_angle_data = lankle_angle_data(:, this_condition_indices(i_stretch));
            this_stretch_rankle_angle_data = -rankle_angle_data(:, this_condition_indices(i_stretch));
%             this_stretch_cop_com_data = cop_com_data(:, this_condition_indices(i_stretch));
            this_stretch_cop_com_int_data = cop_com_int_data(:, this_condition_indices(i_stretch));
            this_stretch_step_placement_x_data = step_placement_x_data(:, this_condition_indices(i_stretch));
            
            % go through bands and store data if the band is a swing
            for i_band = 1 : number_of_bands_per_stretch
                if ~isnan(this_stretch_step_placement_x_data(i_band))
                    stance_foot_control_data = [stance_foot_control_data; this_stretch_stance_foot_data{i_band}]; %#ok<AGROW>
                    com_x_pos_midstance_control_data = [com_x_pos_midstance_control_data; this_stretch_com_x_pos_data(this_condition_midstance_index_data(i_band))]; %#ok<AGROW>
                    com_x_vel_midstance_control_data = [com_x_vel_midstance_control_data; this_stretch_com_x_vel_data(this_condition_midstance_index_data(i_band))]; %#ok<AGROW>
                    com_x_pos_end_control_data = [com_x_pos_end_control_data; this_stretch_com_x_pos_data(100*i_band - 1)]; %#ok<AGROW>
                    com_x_vel_end_control_data = [com_x_vel_end_control_data; this_stretch_com_x_vel_data(100*i_band - 1)]; %#ok<AGROW>
                    com_x_pos_initial_control_data = [com_x_pos_initial_control_data; this_stretch_com_x_pos_data(100*i_band - 99)]; %#ok<AGROW>
                    com_x_vel_initial_control_data = [com_x_vel_initial_control_data; this_stretch_com_x_vel_data(100*i_band - 99)]; %#ok<AGROW>
                    lankle_x_midstance_control_data = [lankle_x_midstance_control_data; this_stretch_lankle_x_data(this_condition_midstance_index_data(i_band))]; %#ok<AGROW>
                    rankle_x_midstance_control_data = [rankle_x_midstance_control_data; this_stretch_rankle_x_data(this_condition_midstance_index_data(i_band))]; %#ok<AGROW>
                    fx_control_data = [fx_control_data; this_stretch_fx_data(end)]; %#ok<AGROW>
                    lankle_x_end_control_data = [lankle_x_end_control_data; this_stretch_lankle_x_data(100*i_band - 1)]; %#ok<AGROW>
                    rankle_x_end_control_data = [rankle_x_end_control_data; this_stretch_rankle_x_data(100*i_band - 1)]; %#ok<AGROW>
                    lankle_x_initial_control_data = [lankle_angle_end_control_data; this_stretch_lankle_x_data(100*i_band - 99)]; %#ok<AGROW>
                    rankle_x_initial_control_data = [rankle_angle_end_control_data; this_stretch_rankle_x_data(100*i_band - 99)]; %#ok<AGROW>
                    lankle_angle_end_control_data = [lankle_angle_end_control_data; this_stretch_lankle_angle_data(100*i_band - 5)]; %#ok<AGROW>
                    rankle_angle_end_control_data = [rankle_angle_end_control_data; this_stretch_rankle_angle_data(100*i_band - 5)]; %#ok<AGROW>
                    cop_com_int_control_data = [cop_com_int_control_data; this_stretch_cop_com_int_data(i_band)]; %#ok<AGROW>
                    step_placement_x_control_data = [step_placement_x_control_data; this_stretch_step_placement_x_data(i_band)]; %#ok<AGROW>
                end
            end
        end
    end
    
    %% Linear Model
    stance_styles = {'gx', 'rx'};
    step_linear_models = cell(size(stance_feet));
    com_from_ankle_midstance_means = zeros(size(stance_feet));
    com_vel_midstance_means = zeros(size(stance_feet));
    com_from_ankle_end_means = zeros(size(stance_feet));
    com_x_vel_initial_means = zeros(size(stance_feet));
    com_from_ankle_initial_means = zeros(size(stance_feet));
    cop_com_int_data_means = zeros(size(stance_feet));
    foot_placement_means = zeros(size(stance_feet));
    
    for i_stance = 1 : length(stance_feet)
        % get data for this stance foot
        this_stance_foot_indicator = strcmp(stance_foot_control_data, stance_feet{i_stance});
        this_stance_com_x_pos_midstance_data = com_x_pos_midstance_control_data(this_stance_foot_indicator);
        this_stance_com_x_vel_midstance_data = com_x_vel_midstance_control_data(this_stance_foot_indicator);
        this_stance_com_x_pos_end_data = com_x_pos_end_control_data(this_stance_foot_indicator);
        this_stance_com_x_pos_initial_data = com_x_pos_initial_control_data(this_stance_foot_indicator);
        this_stance_cop_com_int_data = cop_com_int_control_data(this_stance_foot_indicator);
        this_stance_com_x_vel_initial_data  = com_x_vel_initial_control_data(this_stance_foot_indicator);    
        this_stance_step_placement_x_data = step_placement_x_control_data(this_stance_foot_indicator);
        this_stance_fx_data = fx_control_data(this_stance_foot_indicator);

        if strcmp(stance_feet{i_stance}, 'STANCE_LEFT')
            this_stance_foot_ankle_x_data_midstance = lankle_x_midstance_control_data(this_stance_foot_indicator);
        end
        if strcmp(stance_feet{i_stance}, 'STANCE_RIGHT')
            this_stance_foot_ankle_x_data_midstance = rankle_x_midstance_control_data(this_stance_foot_indicator);
        end
        this_stance_com_from_stance_ankle_data_midstance = this_stance_com_x_pos_midstance_data - this_stance_foot_ankle_x_data_midstance;
        
        if strcmp(stance_feet{i_stance}, 'STANCE_LEFT')
            this_stance_foot_ankle_x_data_end = lankle_x_end_control_data(this_stance_foot_indicator);
        end
        if strcmp(stance_feet{i_stance}, 'STANCE_RIGHT')
            this_stance_foot_ankle_x_data_end = rankle_x_end_control_data(this_stance_foot_indicator);
        end
        this_stance_com_from_stance_ankle_data_end = this_stance_com_x_pos_end_data - this_stance_foot_ankle_x_data_end;
        
        if strcmp(stance_feet{i_stance}, 'STANCE_LEFT')
            this_stance_foot_ankle_x_data_initial = lankle_x_initial_control_data(this_stance_foot_indicator);
        end
        if strcmp(stance_feet{i_stance}, 'STANCE_RIGHT')
            this_stance_foot_ankle_x_data_initial = rankle_x_initial_control_data(this_stance_foot_indicator);
        end
        this_stance_com_from_stance_ankle_data_initial = this_stance_com_x_pos_initial_data - this_stance_foot_ankle_x_data_initial;
        
        % calculate and remove means
        com_from_ankle_midstance_means(i_stance) = mean(this_stance_com_from_stance_ankle_data_midstance);
        com_from_ankle_end_means(i_stance) = mean(this_stance_com_from_stance_ankle_data_end);
        com_vel_midstance_means(i_stance) = mean(this_stance_com_x_vel_midstance_data);
        foot_placement_means(i_stance) = mean(this_stance_step_placement_x_data);
        
        
%         fx_means(i_stance) = mean(this_stance_fx_data);
        
        com_from_ankle_initial_means(i_stance) = mean(this_stance_com_from_stance_ankle_data_initial);
        
        com_x_vel_initial_means(i_stance) = mean(this_stance_com_x_vel_initial_data);
        
        cop_com_int_data_means(i_stance) = mean(this_stance_cop_com_int_data);
        
        this_stance_com_from_stance_ankle_data_mean_free_midstance = this_stance_com_from_stance_ankle_data_midstance - com_from_ankle_midstance_means(i_stance);
        this_stance_com_x_vel_midstance_data_mean_free_midstance = this_stance_com_x_vel_midstance_data - com_vel_midstance_means(i_stance);
        this_stance_com_x_from_stance_ankle_data_mean_free_end = this_stance_com_from_stance_ankle_data_end - com_from_ankle_end_means(i_stance);
        
        this_stance_com_from_stance_ankle_data_mean_free_initial = this_stance_com_from_stance_ankle_data_initial - com_from_ankle_initial_means(i_stance);
        
        this_stance_cop_com_x_int_data_mean_free = this_stance_cop_com_int_data - cop_com_int_data_means(i_stance);
        
        this_stance_com_x_vel_data_mean_free_initial = this_stance_com_x_vel_initial_data - com_x_vel_initial_means(i_stance);
        this_stance_step_placement_x_data_mean_free = this_stance_step_placement_x_data - foot_placement_means(i_stance);
        
        if strcmp(stance_feet{i_stance}, 'STANCE_LEFT')
            this_stance_foot_ankle_angle_data_end = lankle_angle_end_control_data(this_stance_foot_indicator);
        end
        if strcmp(stance_feet{i_stance}, 'STANCE_RIGHT')
            this_stance_foot_ankle_angle_data_end = rankle_angle_end_control_data(this_stance_foot_indicator);
        end
        
        ankle_angle_data_end_means(i_stance) = mean(this_stance_foot_ankle_angle_data_end);  
        this_stance_foot_ankle_angle_data_free_end = this_stance_foot_ankle_angle_data_end - ankle_angle_data_end_means(i_stance);
                
        fx_means(i_stance) = mean(this_stance_fx_data);
        this_stance_fx_data_mean_free = this_stance_fx_data - fx_means(i_stance);

        
        % fit step regression model
        [fit_object, fit_stats] = fit([this_stance_com_from_stance_ankle_data_mean_free_midstance, this_stance_com_x_vel_midstance_data_mean_free_midstance], this_stance_step_placement_x_data_mean_free, 'poly11');
        step_linear_models{i_stance} = [fit_object.p10, fit_object.p01];

         [fit_object, fit_stats] = fit([this_stance_com_x_from_stance_ankle_data_mean_free_end, this_stance_foot_ankle_angle_data_free_end], this_stance_fx_data_mean_free, 'poly11');
        fx_linear_models{i_stance} = [fit_object.p10, fit_object.p01];
        
        [fit_object, fit_stats] = fit([this_stance_com_from_stance_ankle_data_mean_free_initial, this_stance_com_x_vel_data_mean_free_initial], this_stance_cop_com_x_int_data_mean_free, 'poly11');
        cop_linear_models{i_stance} = [fit_object.p10, fit_object.p01];
        
        if visualize
            figure; 
            plot(fit_object, [this_stance_com_from_stance_ankle_data_mean_free_midstance, this_stance_com_x_vel_midstance_data_mean_free_midstance], this_stance_step_placement_x_data_mean_free)
            xlabel('\Delta CoM from stance ankle'); ylabel('\Delta CoM vel'); zlabel('\Delta foot placement')
            title(['slopes = [' num2str(linear_model_slopes(1)) ', ' num2str(linear_model_slopes(2)) '], r^2 = ' num2str(fit_stats.rsquare)]);
        end    
        
        if visualize
            figure; 
            plot(fit_object, [this_stance_com_from_stance_ankle_data_mean_free_midstance, this_stance_com_x_vel_midstance_data_mean_free_midstance], this_stance_fx_data_mean_free)
            xlabel('\Delta CoM from stance ankle'); ylabel('\Delta CoM vel'); zlabel('\Delta foot placement')
            title(['slopes = [' num2str(linear_model_slopes(1)) ', ' num2str(linear_model_slopes(2)) '], r^2 = ' num2str(fit_stats.rsquare)]);
        end    
    end

    %% for each data point, calculate difference from prediction
    foot_placement_adjusted_x_data = zeros(size(step_placement_x_data)) * NaN;
    foot_placement_adjusted_x_directions  = step_placement_x_directions;
    fx_adjusted_data = zeros(size(step_placement_x_data)) * NaN;
    fx_adjusted_directions = fx_directions;
    cop_com_int_adjusted_directions = cop_com_int_directions;
    
    for i_stretch = 1 : number_of_stretches
        this_stretch_stance_foot_data = conditions_session.stance_foot_data(i_stretch, :);
        this_condition_midstance_index_data = midstance_index_data(:, i_stretch);
        this_stretch_com_x_pos_data = com_x_data(:, i_stretch);
        this_stretch_com_x_vel_data = com_x_vel_data(:, i_stretch);
        this_stretch_lankle_x_data = lankle_x_data(:, i_stretch);
        this_stretch_rankle_x_data = rankle_x_data(:, i_stretch);
        this_stretch_fx_data = fx_data(:, i_stretch);
        
        this_stretch_lankle_angle_data = lankle_angle_data(:, i_stretch);
        this_stretch_rankle_angle_data = rankle_angle_data(:, i_stretch);
        this_stretch_cop_com_int_data = cop_com_int_data(:, i_stretch);
        this_stretch_step_placement_x_data = step_placement_x_data(:, i_stretch);
        
        % go through bands
        for i_band = 1 : number_of_bands_per_stretch
            if ~isnan(this_stretch_step_placement_x_data(i_band))
                % get data
                foot_index = find(strcmp(stance_feet, this_stretch_stance_foot_data{i_band}));
                com_x_pos_midstance_data = this_stretch_com_x_pos_data(this_condition_midstance_index_data(i_band));
                com_x_vel_midstance_data = this_stretch_com_x_vel_data(this_condition_midstance_index_data(i_band));
                com_x_vel_initial_data = this_stretch_com_x_vel_data(100*i_band - 99);
                com_x_pos_end_data = this_stretch_com_x_pos_data(100*i_band - 1);
                com_x_pos_initial_data = this_stretch_com_x_pos_data(100*i_band - 99);
                lankle_x_midstance_data = this_stretch_lankle_x_data(this_condition_midstance_index_data(i_band));
                rankle_x_midstance_data = this_stretch_rankle_x_data(this_condition_midstance_index_data(i_band));
                lankle_x_end_data = this_stretch_lankle_x_data(100*i_band-1);
                rankle_x_end_data = this_stretch_rankle_x_data(100*i_band-1);
                lankle_x_initial_data = this_stretch_lankle_x_data(100*i_band-99);
                rankle_x_initial_data = this_stretch_rankle_x_data(100*i_band-99);
                lankle_angle_end_data = this_stretch_lankle_angle_data(100*i_band-5);
                rankle_angle_end_data = this_stretch_rankle_angle_data(100*i_band-5);
                cop_com_int_data_here = this_stretch_cop_com_int_data(i_band);
                fx_data_here = this_stretch_fx_data(end);
                step_placement_x_here = this_stretch_step_placement_x_data(i_band);
                
                if strcmp(this_stretch_stance_foot_data{i_band}, 'STANCE_LEFT')
                    this_stance_foot_ankle_x_data_midstance = lankle_x_midstance_data;
                end
                if strcmp(this_stretch_stance_foot_data{i_band}, 'STANCE_RIGHT')
                    this_stance_foot_ankle_x_data_midstance = rankle_x_midstance_data;
                end
                com_from_stance_ankle_data_midstance = com_x_pos_midstance_data - this_stance_foot_ankle_x_data_midstance;
                
                if strcmp(this_stretch_stance_foot_data{i_band}, 'STANCE_LEFT')
                    this_stance_foot_ankle_x_data_end = lankle_x_end_data;
                end
                if strcmp(this_stretch_stance_foot_data{i_band}, 'STANCE_RIGHT')
                    this_stance_foot_ankle_x_data_end = rankle_x_end_data;
                end
                com_from_stance_ankle_data_end = com_x_pos_end_data - this_stance_foot_ankle_x_data_end;
                
                if strcmp(this_stretch_stance_foot_data{i_band}, 'STANCE_LEFT')
                    this_stance_foot_ankle_x_data_initial = lankle_x_initial_data;
                end
                if strcmp(this_stretch_stance_foot_data{i_band}, 'STANCE_RIGHT')
                    this_stance_foot_ankle_x_data_initial = rankle_x_initial_data;
                end
                com_from_stance_ankle_data_initial = com_x_pos_initial_data - this_stance_foot_ankle_x_data_initial;
                
                
                if strcmp(stance_feet{i_band}, 'STANCE_LEFT')
                    this_stance_foot_ankle_angle_data_end = lankle_angle_end_data;
                end
                if strcmp(stance_feet{i_band}, 'STANCE_RIGHT')
                    this_stance_foot_ankle_angle_data_end = rankle_angle_end_data;
                end
                
                this_stance_foot_ankle_angle_end_delta = this_stance_foot_ankle_angle_data_end - ankle_angle_data_end_means(foot_index);
                
                % calculated step deltas
                com_from_stance_ankle_midstance_delta = com_from_stance_ankle_data_midstance - com_from_ankle_midstance_means(foot_index);
                com_from_stance_ankle_data_end_delta = com_from_stance_ankle_data_end - com_from_ankle_end_means(foot_index); 
                com_from_stance_ankle_initial_delta = com_from_stance_ankle_data_initial - com_from_ankle_initial_means(foot_index);
                com_vel_midstance_delta = com_x_vel_midstance_data - com_vel_midstance_means(foot_index);
                com_vel_initial_delta = com_x_vel_initial_data - com_x_vel_initial_means(foot_index);
                foot_placement_delta = step_placement_x_here - foot_placement_means(foot_index);
                fx_delta = fx_data_here - fx_means(foot_index);
                cop_com_int_delta = cop_com_int_data_here - cop_com_int_data_means(foot_index);
                
                step_model_slopes = step_linear_models{foot_index};
                fx_model_slopes = fx_linear_models{foot_index};
                cop_model_slopes = cop_linear_models{foot_index};
                
                predicted_foot_placement_change = com_from_stance_ankle_midstance_delta * step_model_slopes(1) + com_vel_midstance_delta * step_model_slopes(2);
                actual_foot_placement_change = foot_placement_delta;
                foot_placement_adjusted = actual_foot_placement_change - predicted_foot_placement_change;
                foot_placement_adjusted_x_data(i_band, i_stretch) = foot_placement_adjusted;
                
                predicted_fx_change = com_from_stance_ankle_data_end_delta * fx_model_slopes(1) + this_stance_foot_ankle_angle_end_delta * fx_model_slopes(2);
                actual_fx_change = fx_delta;
                fx_adjusted = actual_fx_change - predicted_fx_change;
                fx_adjusted_data(i_band, i_stretch) = fx_adjusted;
                
                predicted_cop_com_int_change = com_from_stance_ankle_initial_delta * cop_model_slopes(1) + com_vel_initial_delta * cop_model_slopes(2);
                actual_cop_com_int_change = cop_com_int_delta;
                cop_com_int_adjusted = actual_cop_com_int_change - predicted_cop_com_int_change;
                cop_com_int_adjusted_data(i_band, i_stretch) = cop_com_int_adjusted;
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
            foot_placement_adjusted_x_data, ...
            'foot_placement_adjusted_x',...
            foot_placement_adjusted_x_directions ...
          );
     [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
        addOrReplaceResultsData ...
          ( ...
            analysis_data_session, ...
            analysis_names_session, ...
            analysis_directions_session, ...
            fx_adjusted_data, ...
            'fx_adjusted',...
            fx_adjusted_directions ...
          );  
      [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
        addOrReplaceResultsData ...
          ( ...
            analysis_data_session, ...
            analysis_names_session, ...
            analysis_directions_session, ...
            cop_com_int_adjusted_data, ...
            'cop_com_int_adjusted',...
            cop_com_int_adjusted_directions ...
          );
%       
%       
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












