% run from groupFigs folder
clear,clc

loaded_data = load('G:\My Drive\CoBaL_Test_Data\Cadence\GVS\results.mat');
figure_settings = SettingsCustodian('figureSettings_CAD.txt');

% initialize
variables_to_plot = figure_settings.get('variables_to_plot');
number_of_variables_to_plot = size(variables_to_plot, 1);
double_stance_patch_color = figure_settings.get('stance_double_color');
towards_color_80 = figure_settings.get('towards_color_80');
towards_color_110 = figure_settings.get('towards_color_110');

stance_alpha = figure_settings.get('stance_alpha');
show_outliers = figure_settings.get('show_outliers');
box_width = figure_settings.get('box_width');

%% CAD

% count number of data points
% N_triggerRight_stimTowards_stepOne = sum(strcmp(loaded_data_PD.conditions.condition_stance_foot_list, 'STANCE_RIGHT') & strcmp(loaded_data_PD.conditions.condition_perturbation_list, 'ILLUSION_RIGHT') & strcmp(loaded_data_PD.conditions.condition_index_list, 'ONE'));
% N_triggerRight_stimTowards_stepTwo = sum(strcmp(loaded_data_PD.conditions.condition_stance_foot_list, 'STANCE_LEFT') & strcmp(loaded_data_PD.conditions.condition_perturbation_list, 'ILLUSION_RIGHT') & strcmp(loaded_data_PD.conditions.condition_index_list, 'TWO'));
% N_triggerRight_stimTowards_stepThree = sum(strcmp(loaded_data_PD.conditions.condition_stance_foot_list, 'STANCE_RIGHT') & strcmp(loaded_data_PD.conditions.condition_perturbation_list, 'ILLUSION_RIGHT') & strcmp(loaded_data_PD.conditions.condition_index_list, 'THREE'));
% N_triggerRight_stimTowards_stepFour = sum(strcmp(loaded_data_PD.conditions.condition_stance_foot_list, 'STANCE_LEFT') & strcmp(loaded_data_PD.conditions.condition_perturbation_list, 'ILLUSION_RIGHT') & strcmp(loaded_data_PD.conditions.condition_index_list, 'FOUR'));
% 
% N_table = ...
%   [ ...
%     N_triggerRight_stimTowards_stepOne N_triggerRight_stimTowards_stepTwo N_triggerRight_stimTowards_stepThree N_triggerRight_stimTowards_stepFour; ...
%   ];

% collect variables
variable_data_80 = cell(number_of_variables_to_plot, 4);
variable_data_towards_means_80 = cell(number_of_variables_to_plot, 1);
variable_data_towards_cinvs_80 = cell(number_of_variables_to_plot, 1);
variable_data_towards_stdvs = cell(number_of_variables_to_plot, 1);

for i_variable = 1 : number_of_variables_to_plot
    this_variable_name = variables_to_plot{i_variable, 1};
    this_variable_collected_data_firstStanceRight_80 = cell(1, 4);
    
        % change variable names
    if strcmp(this_variable_name, 'right_ankle_dorsiflexion')
       this_variable_name = 'joint_angle:right_ankle_dorsiflexion';     
    else
       this_variable_name = this_variable_name;
    end
    
    % will need to adjust for two PD variable names and concat them..
    
    this_condition_index_towards = strcmp(loaded_data.conditions.stimulus_list, 'STIM_RIGHT');
    this_cadence_index_80 =  strcmp(loaded_data.conditions.cadence_list, '80BPM');
    this_cadence_index_110 =  strcmp(loaded_data.conditions.cadence_list, '110BPM');
    
    this_variable_data_80 = loaded_data.variable_data{strcmp(loaded_data.variable_names, this_variable_name)}(:, this_condition_index_towards & this_cadence_index_80);
    this_variable_data_110 = loaded_data.variable_data{strcmp(loaded_data.variable_names, this_variable_name)}(:, this_condition_index_towards & this_cadence_index_110);

    % may need to concat CADGVS02
%     this_condition_index_towards = strcmp(loaded_data_PD02.conditions.trigger_foot_list, 'TRIGGER_RIGHT') & strcmp(loaded_data_PD02.conditions.stimulus_list, 'STIM_RIGHT');
%     this_variable_data_towards = [this_variable_data_towards loaded_data_PD02.variable_data{strcmp(loaded_data_PD02.variable_names, this_variable_name)}(:, this_condition_index_towards)];

    
    if size(this_variable_data_80, 1) > 4
        step_one_indices = 1 : 100;
        step_two_indices = 100 : 199;
        step_three_indices = 199 : 298;
        step_four_indices = 298 : 397;
    else
        step_one_indices = 1;
        step_two_indices = 2;
        step_three_indices = 3;
        step_four_indices = 4;
    end
    
    % add data
    this_variable_collected_data_firstStanceRight_80{1, 1} = [this_variable_data_80(step_one_indices, :)];
    this_variable_collected_data_firstStanceRight_80{1, 2} = [this_variable_data_80(step_two_indices, :)];
    this_variable_collected_data_firstStanceRight_80{1, 3} = [this_variable_data_80(step_three_indices, :)];
    this_variable_collected_data_firstStanceRight_80{1, 4} = [this_variable_data_80(step_four_indices, :)];
    
    % add data
    this_variable_collected_data_firstStanceRight_110{1, 1} = [this_variable_data_110(step_one_indices, :)];
    this_variable_collected_data_firstStanceRight_110{1, 2} = [this_variable_data_110(step_two_indices, :)];
    this_variable_collected_data_firstStanceRight_110{1, 3} = [this_variable_data_110(step_three_indices, :)];
    this_variable_collected_data_firstStanceRight_110{1, 4} = [this_variable_data_110(step_four_indices, :)];


    % merge and store
    if size(this_variable_collected_data_firstStanceRight_80{1, 1}, 1) == 1
        variable_data_towards_means_80{i_variable} = zeros(1, 4);
        variable_data_towards_cinvs_80{i_variable} = zeros(1, 4);
        for i_col = 1 : 4
            variable_data_80{i_variable, i_col} = [this_variable_collected_data_firstStanceRight_80{1, i_col}];
            variable_data_towards_means_80{i_variable}(i_col) = mean(variable_data_80{i_variable, i_col});
            variable_data_towards_cinvs_80{i_variable}(i_col) = cinv(variable_data_80{i_variable, i_col}');
            variable_data_towards_stdvs_80{i_variable}(i_col) = std(variable_data_80{i_variable, i_col}');
        end
    end
    
    
    % merge and store
    if size(this_variable_collected_data_firstStanceRight_110{1, 1}, 1) == 1
        variable_data_towards_means_110{i_variable} = zeros(1, 4);
        variable_data_towards_cinvs_110{i_variable} = zeros(1, 4);
        for i_col = 1 : 4
            variable_data_110{i_variable, i_col} = [this_variable_collected_data_firstStanceRight_110{1, i_col}];
            variable_data_towards_means_110{i_variable}(i_col) = mean(variable_data_110{i_variable, i_col});
            variable_data_towards_cinvs_110{i_variable}(i_col) = cinv(variable_data_110{i_variable, i_col}');
            variable_data_towards_stdvs_110{i_variable}(i_col) = std(variable_data_110{i_variable, i_col}');
        end
    end
    
    if size(this_variable_collected_data_firstStanceRight_80{1, 1}, 1) > 1
        % merge sources
        this_variable_collected_data_80 = cell(1, 4);
        this_variable_means_by_step_80 = cell(1, 4);
        this_variable_cinvs_by_step_80 = cell(1, 4);
        this_variable_stdvs_by_step_80 = cell(1, 4);
        this_variable_grand_means_by_step_80 = cell(1, 4);
        for i_col = 1 : 4
            this_variable_collected_data_80{1, i_col} = [this_variable_collected_data_firstStanceRight_80{1, i_col}];
            this_variable_means_by_step_80{1, i_col} = interp1(1:100, nanmean(this_variable_collected_data_80{1, i_col}, 2), linspace(1, 100, 101));
            this_variable_cinvs_by_step_80{1, i_col} = interp1(1:100, cinv(this_variable_collected_data_80{1, i_col}, 2), linspace(1, 100, 101));
            this_variable_stdvs_by_step_80{1, i_col} = interp1(1:100, nanstd(this_variable_collected_data_80{1, i_col}, 1, 2), linspace(1, 100, 101));
            
            % what is happening here?
            this_variable_grand_means_by_step_80{i_col} = interp1(1:100, this_variable_collected_data_80{1, i_col}, linspace(1, 100, 101));
            
            variable_data_80{i_variable, i_col} = this_variable_collected_data_80{1, i_col};

        end

        % merge steps
        this_variable_towards_means_80 = zeros(1, 401);
        this_variable_towards_cinvs_80 = zeros(1, 401);
        this_variable_towards_stdvs_80 = zeros(1, 401);
%         this_variable_grand_means = zeros(1, 401);
        for i_step = 1 : 4
            target_range = (i_step-1)*100 + (1:101);
            this_variable_towards_means_80(target_range) = this_variable_means_by_step_80{1, i_step};
            this_variable_towards_cinvs_80(target_range) = this_variable_cinvs_by_step_80{1, i_step};
            this_variable_towards_stdvs_80(target_range) = this_variable_stdvs_by_step_80{1, i_step};
%             this_variable_grand_means(target_range) = this_variable_grand_means_by_step{i_step};

            if i_step > 1
                this_variable_towards_means_80(target_range(1)) = mean([this_variable_means_by_step_80{1, i_step-1}(end), this_variable_means_by_step_80{1, i_step}(1)]);
                this_variable_towards_cinvs_80(target_range(1)) = mean([this_variable_cinvs_by_step_80{1, i_step-1}(end), this_variable_cinvs_by_step_80{1, i_step}(1)]);
                this_variable_towards_stdvs_80(target_range(1)) = mean([this_variable_stdvs_by_step_80{1, i_step-1}(end), this_variable_stdvs_by_step_80{1, i_step}(1)]);
%                 this_variable_grand_means(target_range(1)) = mean([this_variable_grand_means_by_step{i_step-1}(end), this_variable_grand_means_by_step{i_step}(1)]);
            end
        end
        
        % smooth over connection points
        if figure_settings.get('smooth_connection_points')
            this_variable_towards_means_filtered = this_variable_towards_means_80;
            this_variable_towards_means_replaced = this_variable_towards_means_80;
            this_variable_towards_cinvs_filtered = this_variable_towards_cinvs_80;
            this_variable_towards_cinvs_replaced = this_variable_towards_cinvs_80;
            this_variable_towards_stdvs_filtered = this_variable_towards_stdvs_80;
            this_variable_towards_stdvs_replaced = this_variable_towards_stdvs_80;
 
            filter_order = 4;
            cutoff_frequency = 10; % in Hz
            sampling_rate = (0.6 / 100)^(-1);
            [b_filter, a_filter] = butter(filter_order, cutoff_frequency/(sampling_rate/2), 'low');
            smooth_range_width = 6;
            sinusoid_original = (cos(linspace(0, 2*pi, smooth_range_width*2+1)) + 1) * 0.5;
            sinusoid_filtered = 1 - sinusoid_original;
            for i_step = 1 : 3
                merge_index = 100*i_step + 1;
                merge_range = merge_index-smooth_range_width : merge_index+smooth_range_width;

                this_variable_towards_means_range = this_variable_towards_means_80(merge_range);
                this_variable_towards_means_range_filtered = filtfilt(b_filter, a_filter, this_variable_towards_means_range);
                this_variable_towards_means_range_replace = this_variable_towards_means_range .* sinusoid_original + this_variable_towards_means_range_filtered .* sinusoid_filtered;
                this_variable_towards_means_filtered(merge_range) = this_variable_towards_means_range_filtered;
                this_variable_towards_means_replaced(merge_range) = this_variable_towards_means_range_replace;

                this_variable_towards_cinvs_range = this_variable_towards_cinvs_80(merge_range);
                this_variable_towards_cinvs_range_filtered = filtfilt(b_filter, a_filter, this_variable_towards_cinvs_range);
                this_variable_towards_cinvs_range_replace = this_variable_towards_cinvs_range .* sinusoid_original + this_variable_towards_cinvs_range_filtered .* sinusoid_filtered;
                this_variable_towards_cinvs_filtered(merge_range) = this_variable_towards_cinvs_range_filtered;
                this_variable_towards_cinvs_replaced(merge_range) = this_variable_towards_cinvs_range_replace;

                this_variable_towards_stdvs_range = this_variable_towards_stdvs_80(merge_range);
                this_variable_towards_stdvs_range_filtered = filtfilt(b_filter, a_filter, this_variable_towards_stdvs_range);
                this_variable_towards_stdvs_range_replace = this_variable_towards_stdvs_range .* sinusoid_original + this_variable_towards_stdvs_range_filtered .* sinusoid_filtered;
                this_variable_towards_stdvs_filtered(merge_range) = this_variable_towards_stdvs_range_filtered;
                this_variable_towards_stdvs_replaced(merge_range) = this_variable_towards_stdvs_range_replace;

            end

            % store
            variable_data_towards_means_80{i_variable} = this_variable_towards_means_replaced;
            variable_data_towards_cinvs_80{i_variable} = this_variable_towards_cinvs_replaced;
            variable_data_towards_stdvs{i_variable} = this_variable_towards_stdvs_replaced;
        else
            variable_data_towards_means_80{i_variable} = this_variable_towards_means_80;
            variable_data_towards_cinvs_80{i_variable} = this_variable_towards_cinvs_80;
            variable_data_towards_stdvs{i_variable} = this_variable_towards_stdvs_80;
        end
    end
    
    if size(this_variable_collected_data_firstStanceRight_110{1, 1}, 1) > 1
        % merge sources
        this_variable_collected_data_110 = cell(1, 4);
        this_variable_means_by_step_110 = cell(1, 4);
        this_variable_cinvs_by_step_110 = cell(1, 4);
        this_variable_stdvs_by_step_110 = cell(1, 4);
        this_variable_grand_means_by_step_110 = cell(1, 4);
        for i_col = 1 : 4
            this_variable_collected_data_110{1, i_col} = [this_variable_collected_data_firstStanceRight_110{1, i_col}];
            this_variable_means_by_step_110{1, i_col} = interp1(1:100, nanmean(this_variable_collected_data_110{1, i_col}, 2), linspace(1, 100, 101));
            this_variable_cinvs_by_step_110{1, i_col} = interp1(1:100, cinv(this_variable_collected_data_110{1, i_col}, 2), linspace(1, 100, 101));
            this_variable_stdvs_by_step_110{1, i_col} = interp1(1:100, nanstd(this_variable_collected_data_110{1, i_col}, 1, 2), linspace(1, 100, 101));
            
            % what is happening here?
            this_variable_grand_means_by_step_110{i_col} = interp1(1:100, this_variable_collected_data_110{1, i_col}, linspace(1, 100, 101));
            
            variable_data_110{i_variable, i_col} = this_variable_collected_data_110{1, i_col};

        end

        % merge steps
        this_variable_towards_means_110 = zeros(1, 401);
        this_variable_towards_cinvs_110 = zeros(1, 401);
        this_variable_towards_stdvs_110 = zeros(1, 401);
%         this_variable_grand_means = zeros(1, 401);
        for i_step = 1 : 4
            target_range = (i_step-1)*100 + (1:101);
            this_variable_towards_means_110(target_range) = this_variable_means_by_step_110{1, i_step};
            this_variable_towards_cinvs_110(target_range) = this_variable_cinvs_by_step_110{1, i_step};
            this_variable_towards_stdvs_110(target_range) = this_variable_stdvs_by_step_110{1, i_step};
%             this_variable_grand_means(target_range) = this_variable_grand_means_by_step{i_step};

            if i_step > 1
                this_variable_towards_means_110(target_range(1)) = mean([this_variable_means_by_step_110{1, i_step-1}(end), this_variable_means_by_step_110{1, i_step}(1)]);
                this_variable_towards_cinvs_110(target_range(1)) = mean([this_variable_cinvs_by_step_110{1, i_step-1}(end), this_variable_cinvs_by_step_110{1, i_step}(1)]);
                this_variable_towards_stdvs_110(target_range(1)) = mean([this_variable_stdvs_by_step_110{1, i_step-1}(end), this_variable_stdvs_by_step_110{1, i_step}(1)]);
%                 this_variable_grand_means(target_range(1)) = mean([this_variable_grand_means_by_step{i_step-1}(end), this_variable_grand_means_by_step{i_step}(1)]);
            end
        end
        
        % smooth over connection points
        if figure_settings.get('smooth_connection_points')
            this_variable_towards_means_filtered = this_variable_towards_means_110;
            this_variable_towards_means_replaced = this_variable_towards_means_110;
            this_variable_towards_cinvs_filtered = this_variable_towards_cinvs_110;
            this_variable_towards_cinvs_replaced = this_variable_towards_cinvs_110;
            this_variable_towards_stdvs_filtered = this_variable_towards_stdvs_110;
            this_variable_towards_stdvs_replaced = this_variable_towards_stdvs_110;
 
            filter_order = 4;
            cutoff_frequency = 10; % in Hz
            sampling_rate = (0.6 / 100)^(-1);
            [b_filter, a_filter] = butter(filter_order, cutoff_frequency/(sampling_rate/2), 'low');
            smooth_range_width = 6;
            sinusoid_original = (cos(linspace(0, 2*pi, smooth_range_width*2+1)) + 1) * 0.5;
            sinusoid_filtered = 1 - sinusoid_original;
            for i_step = 1 : 3
                merge_index = 100*i_step + 1;
                merge_range = merge_index-smooth_range_width : merge_index+smooth_range_width;

                this_variable_towards_means_range = this_variable_towards_means_110(merge_range);
                this_variable_towards_means_range_filtered = filtfilt(b_filter, a_filter, this_variable_towards_means_range);
                this_variable_towards_means_range_replace = this_variable_towards_means_range .* sinusoid_original + this_variable_towards_means_range_filtered .* sinusoid_filtered;
                this_variable_towards_means_filtered(merge_range) = this_variable_towards_means_range_filtered;
                this_variable_towards_means_replaced(merge_range) = this_variable_towards_means_range_replace;

                this_variable_towards_cinvs_range = this_variable_towards_cinvs_110(merge_range);
                this_variable_towards_cinvs_range_filtered = filtfilt(b_filter, a_filter, this_variable_towards_cinvs_range);
                this_variable_towards_cinvs_range_replace = this_variable_towards_cinvs_range .* sinusoid_original + this_variable_towards_cinvs_range_filtered .* sinusoid_filtered;
                this_variable_towards_cinvs_filtered(merge_range) = this_variable_towards_cinvs_range_filtered;
                this_variable_towards_cinvs_replaced(merge_range) = this_variable_towards_cinvs_range_replace;

                this_variable_towards_stdvs_range = this_variable_towards_stdvs_110(merge_range);
                this_variable_towards_stdvs_range_filtered = filtfilt(b_filter, a_filter, this_variable_towards_stdvs_range);
                this_variable_towards_stdvs_range_replace = this_variable_towards_stdvs_range .* sinusoid_original + this_variable_towards_stdvs_range_filtered .* sinusoid_filtered;
                this_variable_towards_stdvs_filtered(merge_range) = this_variable_towards_stdvs_range_filtered;
                this_variable_towards_stdvs_replaced(merge_range) = this_variable_towards_stdvs_range_replace;

            end

            % store
            variable_data_towards_means_110{i_variable} = this_variable_towards_means_replaced;
            variable_data_towards_cinvs_110{i_variable} = this_variable_towards_cinvs_replaced;
            variable_data_towards_stdvs{i_variable} = this_variable_towards_stdvs_replaced;
        else
            variable_data_towards_means_110{i_variable} = this_variable_towards_means_110;
            variable_data_towards_cinvs_110{i_variable} = this_variable_towards_cinvs_110;
            variable_data_towards_stdvs{i_variable} = this_variable_towards_stdvs_110;
        end
    end
end

%% plot

% calculate pushoff times - use HY data only for now
step_times = loaded_data.step_time_data;
% pushoff_times = loaded_data_HY.pushoff_time_data;
% pushoff_time_ratios = pushoff_times ./ step_times;
% pushoff_time_ratio_mean = mean(pushoff_time_ratios);

% plot 
for i_variable = 1 : number_of_variables_to_plot
    % create figure
    figure; axes; hold on
    title(variables_to_plot{i_variable, 2});
    if figure_settings.get('dictate_axes')
        set(gca, 'ylim', [str2double(variables_to_plot{i_variable, 5}), str2double(variables_to_plot{i_variable, 6})]);
    end
    
    % plots
    if size(variable_data_80{i_variable, 1}, 1) == 1
        set(gca, 'xtick', [])
%         for i_step = 1 : 4
            if strcmp(figure_settings.get('discrete_variable_plot_style'), 'box')
                singleBoxPlot ...
                  ( ...
                    gca, ...
                    i_step - box_width*0.55, ...
                    variable_data_80{i_variable, i_step}, ...
                    towards_color, ...
                    'TOWARDS', ...
                    show_outliers, ...
                    box_width ...
                  )
                singleBoxPlot ...
                  ( ...
                    gca, ...
                    i_step + box_width*0.55, ...
                    variable_data_away{i_variable, i_step}, ...
                    away_color, ...
                    'AWAY', ...
                    show_outliers, ...
                    box_width ...
                  )
                % plot zero
                zero_plot = plot([0 400], [0 0], 'color', [0.7 0.7 0.7]);
                uistack(zero_plot, 'bottom')
            end
            if strcmp(figure_settings.get('discrete_variable_plot_style'), 'bar')
                bar_80 = bar ...
                  ( ...
                    [1 - box_width*1], ...
                    [variable_data_towards_means_80{i_variable}(1)], ...
                    box_width, ...
                    'edgecolor', 'none', ...
                    'facecolor', towards_color_80  ...
                  )
                bar_110 = bar ...
                  ( ...
                    [1 + box_width*1], ...
                    [variable_data_towards_means_110{i_variable}(1)], ...
                    box_width, ...
                    'edgecolor', 'none', ...
                    'facecolor', towards_color_110  ...
                  )
              
                errorbar ...
                  ( ...
                    [1 - box_width*1, 1 + box_width*1 ], ...
                    [variable_data_towards_means_80{i_variable}(1), variable_data_towards_means_110{i_variable}(1)], ...
                    [variable_data_towards_cinvs_80{i_variable}(1), variable_data_towards_cinvs_110{i_variable}(1)], ...
                    'color', [1 1 1]*0.6, ...
                    'LineStyle', 'none' ...
                  )

            end
%         end
%         set(gca, 'xlim', [0.5, 1.5]); % first step
%         set(gca, 'xlim', [0.5, 2.5]); % first two steps
    end
    
    if size(variable_data_80{i_variable, 1}, 1) > 1
        plot_time = 0:400;
        plot_80 = shadedErrorBar(0:400, variable_data_towards_means_80{i_variable}, variable_data_towards_cinvs_80{i_variable}, {'color', towards_color_80, 'linewidth', 6}, 1);
        plot_110 = shadedErrorBar(0:400, variable_data_towards_means_110{i_variable}, variable_data_towards_cinvs_110{i_variable}, {'color', towards_color_110, 'linewidth', 6}, 1);
%         legend([plot_80.mainLine, plot_110.mainLine], '80BPM', '110BPM')
        
        xlim([0 130])
        
        % plot zero
        zero_plot = plot([0 400], [0 0], 'color', [0.7 0.7 0.7]);
        uistack(zero_plot, 'bottom')
        
        % mark double stance
        double_stance_patch_color = figure_settings.get('stance_double_color');
        ylimits = get(gca, 'ylim');
        % shade double stance
        for i_step = [1  2]
            patch_x1 = (i_step - 1) * 100;
            patch_x2 = (i_step - 1) * 100 + 30;
            patch_x = [patch_x1 patch_x2 patch_x2 patch_x1];
            patch_y = [ylimits(1) ylimits(1) ylimits(2) ylimits(2)];
            patch_handle = patch(patch_x, patch_y, double_stance_patch_color, 'EdgeColor', 'none', 'FaceAlpha', stance_alpha); 
            uistack(patch_handle, 'bottom')
        end
    end
    
    if figure_settings.get('save')
        if figure_settings.get('remove_labels')
            set(get(gca, 'xaxis'), 'visible', 'off');
            set(get(gca, 'yaxis'), 'visible', 'off');
            set(get(gca, 'xlabel'), 'visible', 'off');
            set(get(gca, 'ylabel'), 'visible', 'off');
            set(get(gca, 'title'), 'visible', 'off');
            set(gca, 'xticklabel', '');
            set(gca, 'yticklabel', '');
            set(gca, 'position', [0 0 1 1]);
            legend(gca, 'hide');
        end
        save_folder = pwd;
        filename = [save_folder filesep variables_to_plot{i_variable, 9} filesep variables_to_plot{i_variable, 4}];
%         filename = ['rawWithLabels' filesep variables_to_plot{i_variable, 9} filesep variables_to_plot{i_variable, 4}];
        print(gcf, filename, '-dtiff', '-r300')
        close(gcf)
    end
    
end

return









