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
    study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
    
    results_file_name = ['results' filesep makeFileName(collection_date, subject_id, 'results')];
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
    lankle_x_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'marker:LANK_x')};
    rankle_x_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'marker:RANK_x')};
    midstance_index_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, 'midstance_index')};
    number_of_bands_per_stretch = size(step_placement_x_data, 1);
    
    % make condition data tables
    condition_data_all = cell(number_of_stretches, number_of_condition_labels);
    for i_condition = 1 : number_of_condition_labels
        condition_data_all(:, i_condition) = conditions_session.(condition_source_variables{i_condition});
    end
    labels_to_ignore = study_settings.get('conditions_to_ignore', 1);
    levels_to_remove = study_settings.get('levels_to_remove', 1);
    [condition_combination_labels, condition_combinations_control] = determineConditionCombinations(condition_data_all, conditions_settings, labels_to_ignore, levels_to_remove, 'control');
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
    
    if visualize
        figure; axes; hold on; set(gca, 'ylim', [-0.1, 0.1])
        plot(com_x_pos_midstance_control_data(strcmp(stance_foot_control_data, 'STANCE_RIGHT')), 'rx')
        plot(com_x_pos_midstance_control_data(strcmp(stance_foot_control_data, 'STANCE_LEFT')), 'gx')
        title('CoM position at mid-stance')
        xlabel('step'); ylabel('CoM (m)')

        figure; axes; hold on; set(gca, 'ylim', [-0.1, 0.1])
        plot(com_x_vel_midstance_control_data(strcmp(stance_foot_control_data, 'STANCE_RIGHT')), 'rx')
        plot(com_x_vel_midstance_control_data(strcmp(stance_foot_control_data, 'STANCE_LEFT')), 'gx')
        title('CoM velocity at mid-stance')
        xlabel('step'); ylabel('CoM vel (m/s)')

        figure; axes; hold on; set(gca, 'ylim', [-0.25, 0.25])
        plot(step_placement_x_control_data(strcmp(stance_foot_control_data, 'STANCE_RIGHT')), 'rx')
        plot(step_placement_x_control_data(strcmp(stance_foot_control_data, 'STANCE_LEFT')), 'gx')
        title('Foot placement')
        xlabel('step'); ylabel('foot placement (m)')

%         figure; axes; hold on; set(gca, 'ylim', [-0.25, 0.25])
%         plot(lankle_x_midstance_control_data(strcmp(stance_foot_control_data, 'STANCE_RIGHT')), 'rx')
%         plot(lankle_x_midstance_control_data(strcmp(stance_foot_control_data, 'STANCE_LEFT')), 'gx')
% 
%         figure; axes; hold on; set(gca, 'ylim', [-0.25, 0.25])
%         plot(rankle_x_midstance_control_data(strcmp(stance_foot_control_data, 'STANCE_RIGHT')), 'rx')
%         plot(rankle_x_midstance_control_data(strcmp(stance_foot_control_data, 'STANCE_LEFT')), 'gx')
    end
    
    %% fit linear model
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
        
        % calculate means and mean-free data
        com_from_ankle_means(i_stance) = mean(this_stance_com_from_stance_ankle_data);
        com_vel_means(i_stance) = mean(this_stance_com_x_vel_midstance_data);
        foot_placement_means(i_stance) = mean(this_stance_step_placement_x_data);
        this_stance_com_from_stance_ankle_data_mean_free = this_stance_com_from_stance_ankle_data - mean(this_stance_com_from_stance_ankle_data);
        this_stance_com_x_vel_midstance_data_mean_free = this_stance_com_x_vel_midstance_data - mean(this_stance_com_x_vel_midstance_data);
        this_stance_step_placement_x_data_mean_free = this_stance_step_placement_x_data - mean(this_stance_step_placement_x_data);
        
        % fit regression model
        [fit_object, fit_stats] = fit([this_stance_com_from_stance_ankle_data_mean_free, this_stance_com_x_vel_midstance_data_mean_free], this_stance_step_placement_x_data_mean_free, 'poly11'); %#ok<ASGLU>
        linear_models{i_stance} = [fit_object.p10, fit_object.p01];

        if visualize
            stance_styles = {'gx', 'rx'};
            figure; axes; hold on; set(gca, 'ylim', [-0.1, 0.1])
            plot(this_stance_com_from_stance_ankle_data, 'bx')
            plot(this_stance_com_from_stance_ankle_data_mean_free, stance_styles{i_stance})

            figure; axes; hold on; set(gca, 'ylim', [-0.1, 0.1])
            plot(this_stance_com_x_vel_midstance_data_mean_free, 'bx')
            plot(this_stance_com_x_vel_midstance_data, stance_styles{i_stance})

            figure; axes; hold on; set(gca, 'ylim', [-0.25, 0.25])
            plot(this_stance_step_placement_x_data_mean_free, 'bx')
            plot(this_stance_step_placement_x_data, stance_styles{i_stance})

            figure; axes; hold on; 
            plot(this_stance_com_from_stance_ankle_data, 'bx')
            plot(this_stance_com_x_pos_midstance_data, 'mo')
            plot(this_stance_foot_ankle_x_data, 'co')

            
            
            figure; 
            plot(fit_object, [this_stance_com_from_stance_ankle_data_mean_free, this_stance_com_x_vel_midstance_data_mean_free], this_stance_step_placement_x_data_mean_free)
            xlabel('\Delta CoM from stance ankle'); ylabel('\Delta CoM vel'); zlabel('\Delta foot placement')
%             title(['slopes = [' num2str(linear_model_slopes(1)) ', ' num2str(linear_model_slopes(2)) '], r^2 = ' num2str(fit_stats.rsquare)]);
            
            figure; hold on; 
            plot(this_stance_com_from_stance_ankle_data_mean_free, this_stance_step_placement_x_data_mean_free, 'o')
            plot([-0.04, 0.04], [0, 0], 'color', [1 1 1]*0.5);
            plot([0, 0], [-0.15, 0.15], 'color', [1 1 1]*0.5);
            xlabel('\Delta CoM at midstance')
            ylabel('\Delta foot placement')
%             title(['J = ' num2str(Jacobian(1)) ', c = ' num2str(correlation_c(1)) ', p = ' num2str(correlation_p(1))]);
            title('com from ankle - pos');
            set(gca, 'xlim', [-0.04, 0.04], 'ylim', [-0.15, 0.15])
%             
            figure; plot(this_stance_com_x_vel_midstance_data_mean_free, this_stance_step_placement_x_data_mean_free, 'x')
            title('com from ankle - vel');
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
%     if isfield(loaded_data, 'analysis_data_session')
%         analysis_data_session = loaded_data.analysis_data_session;
%         analysis_names_session = loaded_data.analysis_names_session;
%         analysis_directions_session = loaded_data.analysis_directions_session;
%     else
%         analysis_data_session = {};
%         analysis_names_session = {};
%         analysis_directions_session = {};
%     end
%     [analysis_data_session, analysis_names_session, analysis_directions_session] = ...
%         addOrReplaceResultsData ...
%           ( ...
%             analysis_data_session, ...
%             analysis_names_session, ...
%             analysis_directions_session, ...
%             stimulus_response_x_data, ...
%             'stimulus_response_x',...
%             stimulus_response_x_directions ...
%           );
%     variables_to_save = loaded_data;
%     variables_to_save.analysis_data_session = analysis_data_session;
%     variables_to_save.analysis_names_session = analysis_names_session;
%     variables_to_save.analysis_directions_session = analysis_directions_session;
      
    new_data = struct;
    new_data.data = stimulus_response_x_data;
    new_data.directions = stimulus_response_x_directions;
    new_data.name = 'stimulus_response_x';
    data = addOrReplaceResultsData(loaded_data, new_data, 'analysis');
      
    
    save(results_file_name, '-struct', 'data');
end












