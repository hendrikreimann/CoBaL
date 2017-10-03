

% load data and settings
load('E:\DataStorage\Vision_HY\results.mat')
test_settings = SettingsCustodian('analyzeGroups.txt');


subjects = {'DXT','MTB','ZKY'}; % L dom

%% Contains all subjects with a slightly earlier or larger CoP deviation on R 
% subjects = {'DJB','EFU','FNA','IDA','NGY','RON','RRB','SLL','SPA','XDY','YMU'}; % R dom

                % symmetric 

%% Contains only subjects who were earlier or ealier and larger on the R
% subjects = {'DJB','EFU','IDA','RRB','SPA','XDY','YMU'} % possibly symmetric removed?

 % symmetric -- Includes those Rs with only larger

% subjects = {'GHJ','UJD','WHO'}; % No cop

 subject_indicator = ismember(subject_list,subjects);


% initialize
conditions_to_test = test_settings.get('conditions_to_test');
conditions_control = test_settings.get('conditions_control');
number_of_conditions_to_test = size(conditions_to_test, 1);
factor_to_analyze = test_settings.get('factor_to_analyze');
variables_to_test = test_settings.get('variables_to_test');
number_of_variables_to_test = length(variables_to_test);

% test
data_stimulus = cell(number_of_conditions_to_test, number_of_variables_to_test);
h_results = cell(number_of_conditions_to_test, number_of_variables_to_test);
p_results = cell(number_of_conditions_to_test, number_of_variables_to_test);
onset_times = ones(number_of_conditions_to_test, number_of_variables_to_test) * -1;
for i_condition = 1 : number_of_conditions_to_test
    % get condition descriptors
    stimulus_condition = conditions_to_test(i_condition, :);
    control_condition = conditions_control(strcmp(conditions_control(:, 1), stimulus_condition{1}), :);
    
    % get condition indicators
    stance_foot_indicator_stimulus = strcmp(condition_stance_foot_list, stimulus_condition{1});
    stance_foot_indicator_control = strcmp(condition_stance_foot_list, control_condition{1});
    perturbation_indicator_stim_Right = strcmp(condition_perturbation_list, 'ILLUSION_RIGHT');
    perturbation_indicator_stim_Left = strcmp(condition_perturbation_list, 'ILLUSION_LEFT');
    perturbation_indicator_control = strcmp(condition_perturbation_list, control_condition{2});
    index_indicator_stim = strcmp(condition_index_list, stimulus_condition{4});
    index_indicator_control = strcmp(condition_index_list, control_condition{4});
    indicator_stimulus_right = stance_foot_indicator_stimulus & perturbation_indicator_stim_Right & index_indicator_stim & subject_indicator;
    indicator_stimulus_left = stance_foot_indicator_stimulus & perturbation_indicator_stim_Left & index_indicator_stim & subject_indicator;
    indicator_control = stance_foot_indicator_control & perturbation_indicator_control & index_indicator_control & subject_indicator;
    
    % get data and perform t-tests
    for i_variable = 1 : number_of_variables_to_test
        this_variable = variables_to_test{i_variable};
        this_condition = conditions_to_test(i_condition, :);

        tail = 'right';

        step_time_here_right = variable_data{strcmp(variable_names, 'step_time')}(:, indicator_stimulus_right);
        step_time_here_left = variable_data{strcmp(variable_names, 'step_time')}(:, indicator_stimulus_left);
        step_time_here = [step_time_here_right step_time_here_left];
        
        data_here_right = response_data{strcmp(variable_names, this_variable)}(:, indicator_stimulus_right);
        data_here_left_inverted = -response_data{strcmp(variable_names, this_variable)}(:, indicator_stimulus_left);
        data_here = [data_here_right data_here_left_inverted];
        
        data_stimulus{i_condition, i_variable} = data_here;
        [h, p] = bhTest(data_here, 'tail', tail);

        h_results{i_condition, i_variable} = h;
        p_results{i_condition, i_variable} = p;

        % estimate onset time - find absolute peak (but exclude start and end)
        exclude_range_start = 0;
        exclude_range_end = 0;
        filter_order = 2;
        cutoff_frequency = 10; % cutoff frequency, in Hz
        mean_step_time = mean(step_time_here);
        number_of_time_steps_normalized = size(data_here, 1);
        time_normalized = linspace(0, mean_step_time, number_of_time_steps_normalized);
        sampling_rate = number_of_time_steps_normalized / mean_step_time; % sampling rate of the data, in Hz
        [filter_b, filter_a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;

        data_here_mean = mean(data_here, 2);
        data_here_filtered = filtfilt(filter_b, filter_a, data_here_mean);
        data_here_dot = deriveByTime(data_here_filtered, 1/sampling_rate);
        data_here_dot_pruned = [zeros(exclude_range_start, 1); data_here_dot(exclude_range_start+1 : end-exclude_range_end); zeros(exclude_range_end, 1)];
        [~, peak_index_absolutePeak] = max(abs(data_here_dot_pruned));
        data_here_slope_at_absolutePeak = data_here_dot(peak_index_absolutePeak);
        data_here_slope_intersect_absolutePeak =  - time_normalized(peak_index_absolutePeak)*data_here_slope_at_absolutePeak + data_here_filtered(peak_index_absolutePeak);
        onset_time_here_absolutePeak = - data_here_slope_intersect_absolutePeak / data_here_slope_at_absolutePeak;
        
        
%         figure; hold on;
%         plot([1 number_of_time_steps_normalized], [0 0], 'k');
%         plot(data_here_mean);
%         plot(data_here_filtered);
%         plot(data_here_dot);
%         plot(data_here_dot_max_index, data_here_dot(data_here_dot_max_index), 'o');
        
        % estimate onset time - find first peak after 95% confidence interval excludes zero
        data_here_cinv_width = cinv(data_here, 2);
        cinv_lower = data_here_mean - data_here_cinv_width;
        cinv_upper = data_here_mean + data_here_cinv_width;
        cinv_excludes_zero = (cinv_lower > 0);   % just need because of inverted signal
        critical_index = find(cinv_excludes_zero, 1, 'first');
%         while critical_index <= 25
%             critical_index = critical_index + 1;
%         end
        
        data_of_interest = data_here_dot(critical_index : end);
        if max(data_of_interest) > abs(min(data_of_interest))
            [~, peak_index_local] = findpeaks(data_of_interest);
            if isempty(peak_index_local)
                [~, peak_index_local] = max(data_of_interest);
            end
        else
            [~, peak_index_local] = findpeaks(-data_of_interest);
            if isempty(peak_index_local)
                [~, peak_index_local] = max(-data_of_interest);
            end
        end
        
        peak_index_local(peak_index_local < 12) = [];
        peak_index_localPeak = critical_index - 1 + peak_index_local(1);
        data_here_slope_at_localPeak = data_here_dot(peak_index_localPeak);
        data_here_slope_intersect_localPeak =  - time_normalized(peak_index_localPeak)*data_here_slope_at_localPeak + data_here_filtered(peak_index_localPeak);
        onset_time_here_localPeak = - data_here_slope_intersect_localPeak / data_here_slope_at_localPeak;
        
        figure; hold on;
        plot(time_normalized([1 number_of_time_steps_normalized]), [0 0], 'k');
%         plot(time_normalized, data_here_mean);
        shadedErrorBar(time_normalized,data_here_mean,data_here_cinv_width)
%         plot(time_normalized, data_here_filtered);
        plot(time_normalized, data_here_dot);
%         plot(time_normalized(peak_index_absolutePeak), data_here_mean(peak_index_absolutePeak), 'o');
%         plot(time_normalized(peak_index_absolutePeak), data_here_dot(peak_index_absolutePeak), 'o');
%         plot(time_normalized, time_normalized*data_here_slope_at_absolutePeak + data_here_slope_intersect_absolutePeak)
%         plot(onset_time_here_absolutePeak, 0, 'o')
        plot(time_normalized(peak_index_localPeak), data_here_mean(peak_index_localPeak), 'o');
        plot(time_normalized(peak_index_localPeak), data_here_dot(peak_index_localPeak), 'o');
        plot(time_normalized, time_normalized*data_here_slope_at_localPeak + data_here_slope_intersect_localPeak)
        plot(onset_time_here_localPeak, 0, 'o')
        
        
        onset_times(i_condition, i_variable);
        disp(['Variable "' this_variable '", stance ' this_condition{1} ', Averaged stim direction'  ', step ' this_condition{4} ' - onset time: ' num2str(onset_time_here_localPeak)]);
        
    end
end
return

if ~isempty(strcmp(variables_to_test, 'step_placement_x'))
    data_step_placement_x_stim = [];
    data_step_placement_x_control = [];

    % STIM_RIGHT
    stance_foot_indicator = strcmp(condition_stance_foot_list, 'STANCE_RIGHT');
    perturbation_indicator_stimRight = strcmp(condition_perturbation_list, 'ILLUSION_RIGHT');
    perturbation_indicator_stimLeft = strcmp(condition_perturbation_list, 'ILLUSION_LEFT');
    perturbation_indicator_control = strcmp(condition_perturbation_list, 'CONTROL');
    index_indicator_stim = strcmp(condition_index_list, index_to_analyze);
    index_indicator_control = strcmp(condition_index_list, 'CONTROL');
    indicator_stimRight = stance_foot_indicator & perturbation_indicator_stimRight & index_indicator_stim;
    indicator_stimLeft = stance_foot_indicator & perturbation_indicator_stimLeft & index_indicator_stim;
    indicator_control = stance_foot_indicator & perturbation_indicator_control & index_indicator_control;

    data_step_placement_x = variable_data{strcmp(variable_names, 'step_placement_x')};
    data_step_placement_x_stimRight = data_step_placement_x(indicator_stimRight);
    data_step_placement_x_stimLeft = data_step_placement_x(indicator_stimLeft);
    data_step_placement_x_control = data_step_placement_x(indicator_control);
    data_step_placement_x_stimRight_response = data_step_placement_x_stimRight - mean(data_step_placement_x_control);
    data_step_placement_x_stimLeft_response = data_step_placement_x_stimLeft - mean(data_step_placement_x_control);
    data_step_placement_x_control_response = data_step_placement_x_control - mean(data_step_placement_x_control);

    data_step_placement_x_stimLeft_response_inverted = - data_step_placement_x_stimLeft_response;
    data_step_placement_x_stim_response = [data_step_placement_x_stimRight_response data_step_placement_x_stimLeft_response_inverted];
%     [h, p] = ttest2(data_step_placement_x_stim_response, data_step_placement_x_control_response, 'tail', 'right');
%     [h, p] = ttest2(data_step_placement_x_stim_response, data_step_placement_x_control_response);
    [h, p] = ttest(data_step_placement_x_stim_response, 0, 'tail', 'right');
    disp('One-tailed t-test for step placement response on step ONE, stance RIGHT, pooled stim directions')
    disp(['p = ' num2str(p) ', N = ' num2str(length(data_step_placement_x_stim_response)) '\Delta = ' num2str(mean(data_step_placement_x_stim_response))])
end



% visualize
color_left = [0.2, 0.7, 0.07];
color_right = [0.7, 0.2, 0.07];
color_control = [0.2, 0.07, 0.7];
for i_variable = 1 : number_of_variables_to_test
    for i_condition = 1 : number_of_conditions_to_test
        if strcmp(conditions_to_test{i_condition, 2}, 'ILLUSION_LEFT')
            color_stim = color_left;
        end
        if strcmp(conditions_to_test{i_condition, 2}, 'ILLUSION_RIGHT')
            color_stim = color_right;
        end
        figure; my_axes = axes; hold on
        shadedErrorBar ...
          ( ...
            1:100, ...
            mean(data_stimulus{i_condition, i_variable}, 2), ...
            cinv(data_stimulus{i_condition, i_variable}, 2), ...
            { ...
              'color', color_stim, ...
              'linewidth', 6 ...
            }, ...
            1 ...
          );
        
        title(strrep([variables_to_test{i_variable} ' - ' conditions_to_test{i_condition, 1} ' - ' conditions_to_test{i_condition, 2} ' - ' conditions_to_test{i_condition, 4}], '_', ' '));
        
%         % plot significance marker
%         significance_marker_x = [];
%         p = p_results(i_condition, i_variable);
%         if p < 0.05
%             significance_marker_x = 1.5;
%         end
%         if p < 0.01
%             significance_marker_x = [1.45 1.55];
%         end
%         if p < 0.001
%             significance_marker_x = [1.4 1.5 1.6];
%         end
%         ylimits = get(gca, 'ylim');
%         if ~isempty(significance_marker_x)
%             significance_marker_y = ylimits(1) + 0.9*(ylimits(2) - ylimits(1));
%             plot(significance_marker_x, ones(size(significance_marker_x))*significance_marker_y, 'kh', 'markersize', 15, 'linewidth', 1, 'MarkerFaceColor', [0 0 0]);
%         end
%         text(1.5, ylimits(1) + 0.85*(ylimits(2) - ylimits(1)), ['p = ' num2str(p)], 'FontSize', 14, 'HorizontalAlignment', 'center')
%         delta_mean = mean(data_stimulus{i_condition, i_variable}) - mean(data_control{i_condition, i_variable});
%         text(1.5, ylimits(1) + 0.15*(ylimits(2) - ylimits(1)), ['\Delta = ' num2str(delta_mean)], 'FontSize', 14, 'HorizontalAlignment', 'center')
        
    end
end








