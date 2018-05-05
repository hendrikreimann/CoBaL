
plot_step_response          = 1;
plot_raw_trajectories       = 0;
plot_response_trajectories  = 1;

data_folder = '/Users/reimajbi/Neuro/walking/BalanceBeam/model/data/10_gainSearch';
data_folder = '/Users/reimajbi/Neuro/walking/BalanceBeam/model/data/13_withCopModulation';
% data_folder = '/Users/reimajbi/Neuro/walking/BalanceBeam/model/data/12_noCopModulation';

load_file_index = '01';
load_file_label = 'reference';
% load_file_label = 'gvs';
% load_file_label = 'alphaXiDot';
% load_file_label = 'gamma_a';
load_file_label = 'posGvs';

load_file = [load_file_index '_' load_file_label];
reference_file = [load_file_index '_reference'];

variable_label = load_file(4:end);
variable_title = variable_label;
variable_title(variable_title == '_') = '-';

%% prepare
load([data_folder filesep reference_file]);
x_trajectory_grid_squeezed_ref = squeeze(x_trajectory_grid);
p_trajectory_grid_squeezed_ref = squeeze(p_trajectory_grid);
u_trajectory_grid_squeezed_ref = squeeze(u_trajectory_grid);
p_trajectory_reference_cell = squeeze(p_trajectory_grid);
p_trajectory_reference = mean(p_trajectory_reference_cell{1}, 2);
p_trajectory_reference_normalized = p_trajectory_reference - p_trajectory_reference(1);

load([data_folder filesep load_file]);
x_trajectory_grid_squeezed_stm = squeeze(x_trajectory_grid);
p_trajectory_grid_squeezed_stm = squeeze(p_trajectory_grid);
u_trajectory_grid_squeezed_stm = squeeze(u_trajectory_grid);

title_label = load_file(4:end);
if strcmp(variable_label, 'alphaXiDot')
    parameter_values = alphaXiDot_values;
elseif strcmp(variable_label, 'gamma_p')
    parameter_values = gamma_p_values;
elseif strcmp(variable_label, 'gamma_v')
    parameter_values = gamma_v_values;
elseif strcmp(variable_label, 'gamma_a')
    parameter_values = gamma_a_values;
elseif strcmp(variable_label, 'sigma_motor')
    parameter_values = sigma_motor_values;
elseif strcmp(variable_label, 'sigma_head_pos')
    parameter_values = sigma_head_pos_values;
elseif strcmp(variable_label, 'sigma_head_vel')
    parameter_values = sigma_head_vel_values;
elseif strcmp(variable_label, 'sigma_head_acc')
    parameter_values = sigma_head_acc_values;
else
    parameter_values = 0;
end

load '/Users/reimajbi/Dropbox/BalanceBeamPaper/responseData.mat';




%% calculate means and confidence intervals
time_start_index = 506;
time_end_index = 1013;
time_start = time(time_start_index);
time_end = time(time_end_index);
number_of_trajectories = size(p_trajectory_grid_squeezed_stm{1}, 2);
number_of_parameters = length(p_trajectory_grid_squeezed_stm);

p_means = zeros(length(time), size(p_trajectory_grid_squeezed_stm, 1));
p_civs = zeros(length(time), size(p_trajectory_grid_squeezed_stm, 1));
p_response_means = zeros(length(time), size(p_trajectory_grid_squeezed_stm, 1));
p_response_civs = zeros(length(time), size(p_trajectory_grid_squeezed_stm, 1));
for i_parameter = 1 : length(p_trajectory_grid_squeezed_stm);
    % raw trajectories
    p_trajectories = p_trajectory_grid_squeezed_stm{i_parameter};
    p_means(:, i_parameter) = mean(p_trajectories, 2);
    p_sem = std(p_trajectories, 0, 2) * 1/(size(p_trajectories, 2))^(0.5);
    p_civs(:, i_parameter) = tinv(0.025, size(p_trajectories, 2)-1) * p_sem;
    
    % responses
%     p_response_trajectories = p_trajectory_grid_squeezed{i_parameter} - p_means(1, i_parameter) - repmat(p_trajectory_reference, 1, size(p_trajectories, 2));
    p_trajectories_normalized = p_trajectory_grid_squeezed_stm{i_parameter} - p_means(1, i_parameter);
    p_response_trajectories = p_trajectories_normalized - repmat(p_trajectory_reference_normalized, 1, size(p_trajectories, 2));
    
    p_response_means(:, i_parameter) = mean(p_response_trajectories, 2);
    p_response_sem = std(p_response_trajectories, 0, 2) * 1/(size(p_response_trajectories, 2))^(0.5);
    p_response_civs(:, i_parameter) = tinv(0.025, size(p_response_trajectories, 2)-1) * p_response_sem;
end

%% calculate responses
step_responses_model = zeros(number_of_trajectories, number_of_parameters);
stimulus_responses_model = zeros(number_of_trajectories, number_of_parameters);

p_trajectories_ref = p_trajectory_grid_squeezed_ref{1};
x_trajectories_ref = x_trajectory_grid_squeezed_ref{1};
foot_trajectories_ref = p_trajectories_ref(time_start_index : time_end_index, :);
com_position_trajectories_ref = squeeze(x_trajectories_ref(1, time_start_index : time_end_index, :));
com_velocity_trajectories_ref = squeeze(x_trajectories_ref(2, time_start_index : time_end_index, :));
for i_parameter = 1 : number_of_parameters
    p_trajectories_stm = p_trajectory_grid_squeezed_stm{i_parameter};
    x_trajectories_stm = x_trajectory_grid_squeezed_stm{i_parameter};
    foot_trajectories_stm = p_trajectories_stm(time_start_index : time_end_index, :);
    com_position_trajectories_stm = squeeze(x_trajectories_stm(1, time_start_index : time_end_index, :));
    com_velocity_trajectories_stm = squeeze(x_trajectories_stm(2, time_start_index : time_end_index, :));
    
    [ ...
      Jacobian, ...
      correlation_c, ...
      correlation_p, ...
      step_response_model, ...
      stimulus_response_model ...
    ] ...
    = ...
    calculateStepResponse ...
      ( ...
        foot_trajectories_stm, ...
        foot_trajectories_ref, ...
        com_position_trajectories_stm, ...
        com_position_trajectories_ref, ...
        com_velocity_trajectories_stm, ...
        com_velocity_trajectories_ref ...
      );
  Jacobian
  correlation_c
  step_responses_model(:, i_parameter) = step_response_model;
  stimulus_responses_model(:, i_parameter) = stimulus_response_model;
end








%% define stuff for visualizations
color_human_pos = [1 0.3 .1] * 0.7;
color_human_neg = [0.3 1 0.1] * 0.7;
color_human_no = [0.3 0.1 1];

colors = ...
  [ ...
    0 0.4470 0.7410; ...
    0.8500 0.3250 0.0980; ...
    0.9290 0.6940 0.1250; ...
    0.4940 0.1840 0.5560; ...
    0.4660 0.6740 0.1880; ...
    0.3010 0.7450 0.9330; ...
    0.6350 0.0780 0.1840; ...
    0 0.4470 0.7410; ...
    0.8500 0.3250 0.0980; ...
    0.9290 0.6940 0.1250; ...
  ];




%% plot step response
if plot_step_response
    edges = linspace(-0.07, 0.07, 29);
    
%     % step response
%     step_response_figure = figure; axes('fontsize', 12); hold on; title(['step response - ' load_file_label], 'fontsize', 16);
%     human_histogram_pos = histogram(step_responses_pos, edges, 'displaystyle', 'stairs', 'edgecolor', color_human_pos, 'Normalization', 'probability');
%     human_histogram_neg = histogram(step_responses_neg, edges, 'displaystyle', 'stairs', 'edgecolor', color_human_neg, 'Normalization', 'probability');
%     model_histogram = histogram(step_response_model, edges, 'displaystyle', 'stairs', 'edgecolor', colors(1, :), 'Normalization', 'probability');
%     
%     fit_data_x = edges(1 : end-1) + mean(diff(edges))/2;
%     plot_data_x = linspace(-0.08, 0.08, 100);
%     gaussian_fit_human_pos = fit(fit_data_x.', human_histogram_pos.Values.', 'gauss1');
%     gaussian_fit_human_neg = fit(fit_data_x.', human_histogram_neg.Values.', 'gauss1');
%     gaussian_fit_model = fit(fit_data_x.',  model_histogram.Values.', 'gauss1', 'Robust', 'Bisquare');
%     legend('human pos', 'human neg', 'model')
%     
%     plot(plot_data_x, gaussian_fit_human_pos.a1*exp(-((plot_data_x-gaussian_fit_human_pos.b1)/gaussian_fit_human_pos.c1).^2), 'linewidth', 3, 'color', color_human_pos)
%     plot(plot_data_x, gaussian_fit_human_neg.a1*exp(-((plot_data_x-gaussian_fit_human_neg.b1)/gaussian_fit_human_neg.c1).^2), 'linewidth', 3, 'color', color_human_neg)
%     plot(plot_data_x, gaussian_fit_model.a1*exp(-((plot_data_x-gaussian_fit_model.b1)/gaussian_fit_model.c1).^2), 'linewidth', 3, 'color', colors(1, :))

    % stimulus response
    stimulus_response_figure = figure; axes('fontsize', 12); hold on; title(['stimulus response - ' load_file_label], 'fontsize', 16);
    human_histogram_pos = histogram(stimulus_responses_pos, edges, 'displaystyle', 'stairs', 'edgecolor', color_human_pos, 'Normalization', 'probability');
    human_histogram_neg = histogram(stimulus_responses_neg, edges, 'displaystyle', 'stairs', 'edgecolor', color_human_neg, 'Normalization', 'probability');
    model_histograms = cell(number_of_parameters, 1);
    for i_parameter = 1 : number_of_parameters
        model_histograms{i_parameter} = histogram(stimulus_responses_model(:, i_parameter), edges, 'displaystyle', 'stairs', 'edgecolor', colors(i_parameter, :), 'Normalization', 'probability');
    end
    
    fit_data_x = edges(1 : end-1) + mean(diff(edges))/2;
    plot_data_x = linspace(-0.08, 0.08, 100);
    gaussian_fit_human_pos = fit(fit_data_x.', human_histogram_pos.Values.', 'gauss1');
    gaussian_fit_human_neg = fit(fit_data_x.', human_histogram_neg.Values.', 'gauss1');
    gaussian_fits_model = cell(number_of_parameters, 1);
    for i_parameter = 1 : number_of_parameters
        gaussian_fits_model{i_parameter} = fit(fit_data_x.',  model_histograms{i_parameter}.Values.', 'gauss1', 'Robust', 'Bisquare');
    end
%     legend('human pos', 'human neg', 'model')
    
    figure; hold on
    plot(plot_data_x, gaussian_fit_human_pos.a1*exp(-((plot_data_x-gaussian_fit_human_pos.b1)/gaussian_fit_human_pos.c1).^2), 'linewidth', 3, 'color', color_human_pos, 'displayname', 'human, pos')
    plot(plot_data_x, gaussian_fit_human_neg.a1*exp(-((plot_data_x-gaussian_fit_human_neg.b1)/gaussian_fit_human_neg.c1).^2), 'linewidth', 3, 'color', color_human_neg, 'displayname', 'human, neg')
    model_plots = zeros(number_of_parameters, 1);
    for i_parameter = 1 : number_of_parameters
        model_plots(i_parameter) = plot(plot_data_x, gaussian_fits_model{i_parameter}.a1*exp(-((plot_data_x-gaussian_fits_model{i_parameter}.b1)/gaussian_fits_model{i_parameter}.c1).^2), 'linewidth', 3, 'color', colors(i_parameter, :), 'displayname', num2str(parameter_values(i_parameter)));
    end
    legend('show')
%     legend(model_plots, parameter_values)
    
end
% return

%% plot raw trajectories
if plot_raw_trajectories
    figure; hold on; title(variable_title)
    p_plots = cell(length(p_trajectory_grid_squeezed_stm), 1);
    for i_parameter = 1 : length(p_trajectory_grid_squeezed_stm);
        p_plots{i_parameter} = shadedErrorBar(time, p_means(:, i_parameter), p_civs(:, i_parameter), {'color', colors(i_parameter, :), 'linewidth', 5}, 1);
    end
    xlim([time_start time_end])
end
    
%% plot CoP response
if plot_response_trajectories
    figure; hold on; title(variable_title)
    legend_plot_handles = zeros(length(p_trajectory_grid_squeezed_stm)+1, 1);
    legend_plot_labels = cell(length(p_trajectory_grid_squeezed_stm)+1, 1);
    plot_handles = shadedErrorBar(time_normalized + time_start, cop_responses_pos_mean, cop_responses_pos_civ, {'color', color_human_pos, 'linewidth', 5}, 1);
    legend_plot_handles(1) = plot_handles.mainLine;
    plot_handles = shadedErrorBar(time_normalized + time_start, cop_responses_neg_mean, cop_responses_neg_civ, {'color', color_human_neg, 'linewidth', 5}, 1);
    % plot_handles = shadedErrorBar(time_normalized + time_start, cop_responses_no_mean, cop_responses_no_civ, {'color', color_human_no, 'linewidth', 5}, 1);
    legend_plot_labels{1} = 'human';
    for i_parameter = 1 : length(p_trajectory_grid_squeezed_stm);
%         plot_handles = shadedErrorBar(time, p_response_means(:, i_parameter), p_response_civs(:, i_parameter), {'color', colors(i_parameter, :), 'linewidth', 5}, 1);
        plot_handles = plot(time, p_response_means(:, i_parameter), 'color', colors(i_parameter, :), 'linewidth', 5);
%         plot_handles = shadedErrorBar(time, -p_response_means(:, i_parameter), p_response_civs(:, i_parameter), {'color', colors(i_parameter, :), 'linewidth', 5}, 1);
        legend_plot_handles(i_parameter+1) = plot_handles;%.mainLine;
        legend_plot_labels{i_parameter+1} = num2str(parameter_values(i_parameter));
    end
    xlim([time_start time_end - 0.005])
    xlim([time_start time_end + 0.005])
    legend(legend_plot_handles, legend_plot_labels);
end
% plot(com_acc_noboard_axes, f_human, com_acc_psds_mean_noboard_human, 'color', [0 0 0], 'linewidth', 3, 'displayname', 'human');
% plot(com_acc_withboard_axes, f_human, com_acc_psds_mean_withboard_human, 'color', [0 0 0], 'linewidth', 3, 'displayname', 'human');
% plot(board_axes, f_human, board_psds_mean_withboard_human, 'color', [0 0 0], 'linewidth', 3, 'displayname', 'human');

% legend(com_acc_noboard_axes, 'show')





