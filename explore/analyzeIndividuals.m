%% %% %% analyzeIndividuals
function analyzeIndividuals(varargin)
    load('C:\Users\Tyler Fettrow\Documents\Vision_HY\results.mat')
    test_settings = SettingsCustodian('statIndividuals.txt');
    
    % initialize
    conditions_to_test = test_settings.get('conditions_to_test');
    
    subject_group_assignment = test_settings.get('subject_group_assignment');
    
    conditions_control = test_settings.get('conditions_control');
    number_of_conditions_to_test = size(conditions_to_test, 1);
    factor_to_analyze = test_settings.get('factor_to_analyze');
    variable_info = test_settings.get('variables_to_test');
    variables_to_test = variable_info(:,1);
    number_of_variables_to_test = size(variables_to_test, 1);
    stance_foot = {'STANCE_LEFT','STANCE_RIGHT'};
    stance_foot_checks = {'left', 'right'};
    
    data_stimulus = cell(number_of_conditions_to_test, number_of_variables_to_test);
    
    for i_foot = 1:2
        for i_condition = 1 : number_of_conditions_to_test
            stimulus_condition = conditions_to_test(i_condition, :);
            control_condition = conditions_control(strcmp(conditions_control(:, 1), stimulus_condition{1}), :);

            % get condition indicators
            if strcmp(stimulus_condition{4},'LEFT')
                foot_assignment = test_settings.get('foot_assignment_sides');
                foot_assignment_column = 2;
            elseif strcmp(stimulus_condition{4},'RIGHT')
                foot_assignment = test_settings.get('foot_assignment_sides');
                foot_assignment_column = 3;
            elseif strcmp(stimulus_condition{4},'EARLY')
                foot_assignment = test_settings.get('foot_assignment_onset');
                foot_assignment_column = 2;
            elseif strcmp(stimulus_condition{4},'LATE')
                
                foot_assignment_column = 3;
            elseif strcmp(stimulus_condition{4},'NO')
                foot_assignment_column = 4;
            end
            foot_assignment_indicator = strcmp(foot_assignment(:,foot_assignment_column), 'both');
            foot_assignment_left_indicator = strcmp(foot_assignment(:,foot_assignment_column), 'left');
            foot_assignment_right_indicator = strcmp(foot_assignment(:,foot_assignment_column), 'right');

            subject_indices = find(strcmp(foot_assignment(:,foot_assignment_column),stance_foot_checks(i_foot)) | strcmp(foot_assignment(:,foot_assignment_column),'both'));
            subjects = foot_assignment(subject_indices,1);
            subject_indicator = ismember(subject_list,subjects);
            
            stance_foot_indicator_stimulus =  strcmp(condition_stance_foot_list, stance_foot(i_foot));
            stance_foot_indicator_control = strcmp(condition_stance_foot_list, stance_foot(i_foot));
      
            perturbation_indicator_stim_Right = strcmp(condition_perturbation_list, 'ILLUSION_RIGHT');
            perturbation_indicator_stim_Left = strcmp(condition_perturbation_list, 'ILLUSION_LEFT');
            perturbation_indicator_control = strcmp(condition_perturbation_list, 'CONTROL');
            index_indicator_stim = strcmp(condition_index_list, stimulus_condition{1}); % step to analyze 1-4

            index_indicator_control = strcmp(condition_index_list, 'CONTROL');



            indicator_stimulus_right = stance_foot_indicator_stimulus & perturbation_indicator_stim_Right & index_indicator_stim & subject_indicator;
            indicator_stimulus_left = stance_foot_indicator_stimulus & perturbation_indicator_stim_Left & index_indicator_stim & subject_indicator;
            indicator_control = stance_foot_indicator_control & perturbation_indicator_control & index_indicator_control & subject_indicator;

                % get data and perform t-tests
            for i_variable = 1 : number_of_variables_to_test
                this_variable = variables_to_test{i_variable};

                if strcmp(this_variable,'left_hip_abduction_angle')
                    step_time_here = variable_data{strcmp(variable_names, 'step_time')}(:, indicator_stimulus_left);
                    data_here = response_data{strcmp(variable_names, this_variable)}(:, indicator_stimulus_left);
                elseif strcmp(this_variable,'right_hip_abduction_angle')
                    step_time_here = variable_data{strcmp(variable_names, 'step_time')}(:, indicator_stimulus_right);
                    data_here = response_data{strcmp(variable_names, this_variable)}(:, indicator_stimulus_right);
                else    
                    step_time_here_right = variable_data{strcmp(variable_names, 'step_time')}(:, indicator_stimulus_right);
                    step_time_here_left = variable_data{strcmp(variable_names, 'step_time')}(:, indicator_stimulus_left);
                    step_time_here = [step_time_here_right step_time_here_left];

                    data_here_right = response_data{strcmp(variable_names, this_variable)}(:, indicator_stimulus_right);
                    data_here_left = response_data{strcmp(variable_names, this_variable)}(:, indicator_stimulus_left);
                    data_here_left_inverted = -data_here_left; 
                    data_here = [data_here_right data_here_left_inverted];
%                     [val, ind] = max(abs(data_here));

                end
            end
            data_to_plot(i_condition) = {data_here};
            mean_step_time = mean(step_time_here);
            number_of_time_steps_normalized = size(data_here, 1);
            time_normalized(:,i_condition) = linspace(0, mean_step_time, number_of_time_steps_normalized);
        end  
        direction_test = mean(data_to_plot{1},2);
        [val, ind] = max(abs(direction_test));
        if direction_test(ind) < 0
%             data_to_plot = -data_to_plot;
            first_condition = -mean(data_to_plot{1},2);
            second_condition = -mean(data_to_plot{2},2);
            third_condition = -mean(data_to_plot{3},2);
        else
            first_condition = mean(data_to_plot{1},2);
            second_condition = mean(data_to_plot{2},2);
            third_condition = mean(data_to_plot{3},2);
        end
        
        abscissa_unscaled = linspace(0, 100, test_settings.get('number_of_time_steps_normalized'));

        % scale abscissae
        conditions_this_comparison = comparison_indices{i_comparison};
        step_time_means_this_comparison = zeros(size(conditions_this_comparison));
        for i_condition = 1 : length(conditions_this_comparison)
             % find correct condition indicator
             condition_identifier = conditions_to_plot(conditions_this_comparison(i_condition), :);
             stance_foot_indicator = strcmp(condition_stance_foot_list_all, condition_identifier{1});
             perturbation_indicator = strcmp(condition_perturbation_list_all, condition_identifier{2});
             delay_indicator = strcmp(condition_delay_list_all, condition_identifier{3});
             index_indicator = strcmp(condition_index_list_all, condition_identifier{4});
             experimental_indicator = strcmp(condition_experimental_list_all, condition_identifier{5});
             stimulus_indicator = strcmp(condition_stimulus_list_all, condition_identifier{6});
             day_indicator = strcmp(condition_day_list_all, condition_identifier{7});
             this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
             step_time_data_this_condition = step_time_data(:, this_condition_indicator);
             step_time_means_this_comparison(i_condition) = mean(step_time_data_this_condition);
        end
        for i_condition = 1 : length(conditions_this_comparison)
             if strcmp(plot_settings.get('time_plot_style'), 'scaled_to_comparison_mean')
                  abscissa_scaled = abscissa_unscaled * mean(step_time_means_this_comparison) / 100;
                  abscissae_cell{i_comparison, i_variable}(i_condition, :) = abscissa_scaled;
             elseif strcmp(plot_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                  abscissa_scaled = abscissa_unscaled * step_time_means_this_comparison(i_condition) / 100;
                  abscissae_cell{i_comparison, i_variable}(i_condition, :) = abscissa_scaled;
             else
                  abscissae_cell{i_comparison, i_variable}(i_condition, :) = abscissa_unscaled;
             end
        end
             
        for i_comparison = 1 : number_of_comparisons            
            target_abscissa = abscissae_cell{i_comparison, i_variable}(1, :);

            % set x-limits accordingly
            set(axes_handles(i_comparison, i_variable), 'xlim', [target_abscissa(1) target_abscissa(end)]);
        end            
        abscissae_cell_unscaled = cell(size(abscissae_cell));
        % make one figure per episode and variable
        figure_handles = zeros(number_of_episodes, number_of_variables_to_plot);
        axes_handles = zeros(number_of_episodes, number_of_variables_to_plot);
        pos_text_handles = zeros(number_of_episodes, number_of_variables_to_plot);
        neg_text_handles = zeros(number_of_episodes, number_of_variables_to_plot);
        step_start_times_cell = cell(number_of_episodes, number_of_variables_to_plot);
        step_end_times_cell = cell(number_of_episodes, number_of_variables_to_plot);
        step_pushoff_times_cell = cell(number_of_episodes, number_of_variables_to_plot);
        step_stance_foot_cell = cell(number_of_episodes, number_of_variables_to_plot);

        for i_variable = 1 : number_of_variables_to_plot
            for i_episode = 1 : number_of_episodes
                % make figure and axes and store handles
                new_figure = figure; new_axes = axes; hold on;
                figure_handles(i_episode, i_variable) = new_figure;
                axes_handles(i_episode, i_variable) = new_axes;
                this_episode = episode_indices{i_episode};

                % store handles and determine abscissa data for all comparisons in this episode
                xtick = [];
                for i_comparison = 1 : length(this_episode)
                    comparison_variable_to_axes_index_map(this_episode(i_comparison)) = i_episode;
                    
                    % determine which step this is
                    this_comparison = this_episode(i_comparison);
                    conditions_this_comparison = comparison_indices{this_comparison};
                    example_condition_index = conditions_this_comparison(1);
                    condition_identifier = conditions_to_plot(example_condition_index, :);
                    gap_between_steps = 1;
                    if strcmp(condition_identifier{4}, 'ONE')
                        step_index = 1;
                    elseif strcmp(condition_identifier{4}, 'TWO')
                        step_index = 2;
                    elseif strcmp(condition_identifier{4}, 'THREE')
                        step_index = 3;
                    elseif strcmp(condition_identifier{4}, 'FOUR')
                        step_index = 4;
                    end
                    
                     abscissae_cell_unscaled{this_episode(i_comparison), i_variable} = (linspace(0, 100, study_settings.get('number_of_time_steps_normalized')));

                end
                
                 % calculate average step times and scale abscissa
        for i_variable = 1 : number_of_variables_to_plot
            for i_episode = 1 : number_of_episodes
                this_episode = episode_indices{i_episode};
                for i_comparison = 1 : length(this_episode)
                    
                    % determine which step this is
                    this_comparison = this_episode(i_comparison);
                    conditions_this_comparison = comparison_indices{this_comparison};
                    step_time_means_this_comparison = zeros(size(conditions_this_comparison));
                    for i_condition = 1 : length(conditions_this_comparison)
                        % find correct condition indicator
                        condition_identifier = conditions_to_plot(conditions_this_comparison(i_condition), :);
                        stance_foot_indicator = strcmp(condition_stance_foot_list_all, condition_identifier{1});
                        perturbation_indicator = strcmp(condition_perturbation_list_all, condition_identifier{2});
                        delay_indicator = strcmp(condition_delay_list_all, condition_identifier{3});
                        index_indicator = strcmp(condition_index_list_all, condition_identifier{4});
                        experimental_indicator = strcmp(condition_experimental_list_all, condition_identifier{5});
                        stimulus_indicator = strcmp(condition_stimulus_list_all, condition_identifier{6});
                        day_indicator = strcmp(condition_day_list_all, condition_identifier{7});
                        this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
                        if strcmp(plot_settings.get('time_plot_style'), 'scaled_to_comparison_mean')
                            % calculate average step time
                            step_time_data_this_condition = step_time_data(:, this_condition_indicator);
                            step_time_means_this_comparison(i_condition) = mean(step_time_data_this_condition);
                        end
                        
                    end
                    
                    % scale abscissa
                    if isContinuousVariable(i_variable, variable_data_all)
                        abscissa_unscaled = abscissae_cell_unscaled{this_episode(i_comparison), i_variable};
                        for i_condition = 1 : length(conditions_this_comparison)
                            if strcmp(plot_settings.get('time_plot_style'), 'scaled_to_comparison_mean')
                                abscissa_scaled = abscissa_unscaled * mean(step_time_means_this_comparison) / 100;
                                abscissae_cell{this_episode(i_comparison), i_variable}(i_condition, :) = abscissa_scaled;
                            elseif strcmp(plot_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                                abscissa_scaled = abscissa_unscaled * step_time_means_this_comparison(i_condition) / 100;
                                abscissae_cell{this_episode(i_comparison), i_variable}(i_condition, :) = abscissa_scaled;
                            else
                                abscissae_cell{this_episode(i_comparison), i_variable}(i_condition, :) = abscissa_unscaled;
                            end
                        end                    
                    end                    
                end
            end
        end
        
        % determine abscissa offsets
        for i_step = 2 : 4
            for i_variable = 1 : number_of_variables_to_plot
                if isContinuousVariable(i_variable, variable_data_all)
                    for i_episode = 1 : number_of_episodes
                        this_episode = episode_indices{i_episode};
                        for i_comparison = 1 : length(this_episode)

                            % determine which step this is
                            this_comparison = this_episode(i_comparison);
                            conditions_this_comparison = comparison_indices{this_comparison};
                            for i_condition = 1 : length(conditions_this_comparison)
                                % determine step index
                                condition_identifier = conditions_to_plot(conditions_this_comparison(i_condition), :);
                                if strcmp(condition_identifier{4}, 'ONE')
                                    step_index = 1;
                                elseif strcmp(condition_identifier{4}, 'TWO')
                                    step_index = 2;
                                    previous_step_label = 'ONE';
                                elseif strcmp(condition_identifier{4}, 'THREE')
                                    step_index = 3;
                                    previous_step_label = 'TWO';
                                elseif strcmp(condition_identifier{4}, 'FOUR')
                                    step_index = 4;
                                    previous_step_label = 'THREE';
                                end

                                if step_index == i_step
                                    % find condition index for previous step in same condition
                                    previous_step_condition_index = [];
                                    previous_step_comparison_index = [];
                                    for j_comparison = 1 : length(this_episode)
                                        candidate_comparison = this_episode(j_comparison);
                                        conditions_candidate_comparison = comparison_indices{candidate_comparison};

                                        for j_condition = 1 : length(conditions_candidate_comparison)
                                            candidate_condition_identifier = conditions_to_plot(conditions_candidate_comparison(j_condition), :);


                                            if strcmp(condition_identifier{2}, candidate_condition_identifier{2}) ...
                                            && strcmp(condition_identifier{3}, candidate_condition_identifier{3}) ...
                                            && strcmp(previous_step_label, candidate_condition_identifier{4}) ...
                                            && strcmp(condition_identifier{5}, candidate_condition_identifier{5}) ...
                                            && strcmp(condition_identifier{6}, candidate_condition_identifier{6}) ...
                                            && strcmp(condition_identifier{7}, candidate_condition_identifier{7})
                                                previous_step_comparison_index = j_comparison;
                                                previous_step_condition_index = j_condition;
                                            end
                                        end
                                    end

                                    % find out where abscissa for previous step ends
                                    previous_step_last_data_point = abscissae_cell{this_episode(previous_step_comparison_index), i_variable}(previous_step_condition_index, end);
                                    abscissae_cell{this_episode(i_comparison), i_variable}(i_condition, :) = abscissae_cell{this_episode(i_comparison), i_variable}(i_condition, :) + previous_step_last_data_point;
                                end
                            end
                        end
                    end
                end
            end
        end
        % set x-limits
        for i_variable = 1 : number_of_variables_to_plot
            if isContinuousVariable(i_variable, variable_data_all)
                for i_episode = 1 : number_of_episodes
                    % determine time window to show
                    this_episode = episode_indices{i_episode};
                    first_comparison_in_episode_index = this_episode(1);
                    first_step_abscissae = abscissae_cell{first_comparison_in_episode_index, i_variable};
                    episode_start_time = first_step_abscissae(1, 1);
                    last_comparison_in_episode_index = this_episode(end);
                    last_step_abscissae = abscissae_cell{last_comparison_in_episode_index, i_variable};
                    episode_end_time = last_step_abscissae(1, end);

                    % set x-limits accordingly
                    set(axes_handles(i_episode, i_variable), 'xlim', [episode_start_time episode_end_time]);
                end
            end
        end            
                
        
        target_abscissa = abscissae_cell{i_comparison, i_variable}(i_condition, :);
        
        % determine which control condition applies here
        representant_condition_index = conditions_this_comparison(1);
        applicable_control_condition_index = findApplicableControlConditionIndex(conditions_to_plot(representant_condition_index, :), conditions_control);
        applicable_control_condition_labels = conditions_control(applicable_control_condition_index, :);
                
        % extract data for control condition
        stance_foot_indicator = strcmp(condition_stance_foot_list_all, applicable_control_condition_labels{1});
        perturbation_indicator = strcmp(condition_perturbation_list_all, applicable_control_condition_labels{2});
        delay_indicator = strcmp(condition_delay_list_all, applicable_control_condition_labels{3});
        index_indicator = strcmp(condition_index_list_all, applicable_control_condition_labels{4});
        experimental_indicator = strcmp(condition_experimental_list_all, applicable_control_condition_labels{5});
        stimulus_indicator = strcmp(condition_stimulus_list_all, applicable_control_condition_labels{6});
        day_indicator = strcmp(condition_day_list_all, applicable_control_condition_labels{7});
        this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
        data_to_plot_this_condition = data_to_plot(:, this_condition_indicator);
        origin_indices = find(this_condition_indicator);
                
          plot_handles = shadedErrorBar ...
                          ( ...
                            target_abscissa, ...
                            mean(data_to_plot_this_condition, 2), ...
                            spread(data_to_plot_this_condition, spread_method), ...
                            { ...
                              'color', colors_comparison(i_condition, :), ...
                              'linewidth', 6 ...
                            }, ...
                            1, ...
                            target_axes_handle ...
                          ); 
                      
        for i_variable = 1 : number_of_variables_to_plot
            for i_episode = 1 : number_of_episodes
                these_axes = axes_handles(i_episode, i_variable);
                ylimits = get(these_axes, 'ylim');

                step_start_times = step_start_times_cell{i_episode, i_variable};
                step_end_times = step_end_times_cell{i_episode, i_variable};
                step_pushoff_times = step_pushoff_times_cell{i_episode, i_variable};
                step_stance_foot = step_stance_foot_cell{i_episode, i_variable};
                
                for i_step = 1 : length(step_start_times)
                    % double stance patch
                    double_stance_patch_color = plot_settings.get('stance_double_color');
                    stretch_start = step_start_times(i_step);
                    stretch_end = step_pushoff_times(i_step);
                    patch_x = [stretch_start stretch_end stretch_end stretch_start];
                    patch_y = [ylimits(1) ylimits(1) ylimits(2) ylimits(2)];
                    patch_handle = ...
                        patch ...
                          ( ...
                            patch_x, ...
                            patch_y, ...
                            double_stance_patch_color, ...
                            'parent', these_axes, ...
                            'EdgeColor', 'none', ...
                            'FaceAlpha', plot_settings.get('stance_alpha'), ...
                            'HandleVisibility', 'off' ...
                          ); 
                    uistack(patch_handle, 'bottom')
                    
                    % single stance patch
                    single_stance_patch_color = [1 1 1] * 0.8;
                    if step_stance_foot(i_step) == 0
                        single_stance_patch_color = plot_settings.get('stance_both_color');
                    end
                    if step_stance_foot(i_step) == 1
                        single_stance_patch_color = plot_settings.get('stance_left_color');
                    end
                    if step_stance_foot(i_step) == 2
                        single_stance_patch_color = plot_settings.get('stance_right_color');
                    end
                    stretch_start = step_pushoff_times(i_step);
                    stretch_end = step_end_times(i_step);
                    patch_x = [stretch_start stretch_end stretch_end stretch_start];
                    patch_y = [ylimits(1) ylimits(1) ylimits(2) ylimits(2)];
                    patch_handle = ...
                        patch ...
                          ( ...
                            patch_x, ...
                            patch_y, ...
                            single_stance_patch_color, ...
                            'parent', these_axes, ...
                            'EdgeColor', 'none', ...
                            'FaceAlpha', plot_settings.get('stance_alpha'), ...
                            'HandleVisibility', 'off' ...
                          ); 
                    uistack(patch_handle, 'bottom')
                end
            end
        end
        
        
        figure
        hold on
        this_title = strcat(this_variable, '  ---  ', stance_foot(i_foot));
        title(this_title,'Interpreter', 'none')
        ylim([str2num(variable_info{2}),str2num(variable_info{3})])
        xlabel('normalized time (s)')
        ylabel('Total Response from control  (m)')
        h1 = shadedErrorBar(time_normalized(:,1),first_condition,cinv(data_to_plot{1},2),{'color', 'm','linewidth', 6},1 )
        h2 = shadedErrorBar(time_normalized(:,2),second_condition,cinv(data_to_plot{2},2),{'color', 'c','linewidth', 6},1 )
        h3 = shadedErrorBar(time_normalized(:,3),third_condition,cinv(data_to_plot{3},2),{'color', 'b','linewidth', 6},1 )
        if strcmp(conditions_to_test(1, 4),'LEFT') | strcmp(conditions_to_test(1, 4),'RIGHT')
            legend([h1.mainLine h2.mainLine h3.mainLine], {'Left Dominant','Right Dominant','No Response'})
        else strcmp(conditions_to_test(1, 4),'EARLY') | strcmp(conditions_to_test(1, 4),'LATE')
            legend([h1.mainLine h2.mainLine h3.mainLine], {'EARLY CoP Manipulation','Late/No CoP Manipulation','No Response'})
        end
    end  
        end
    end
end