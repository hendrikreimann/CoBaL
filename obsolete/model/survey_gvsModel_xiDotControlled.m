

% data_folder = '/Users/reimajbi/Neuro/walking/BalanceBeam/model/data/10_gainSearch';
% data_folder = '/Users/reimajbi/Neuro/walking/BalanceBeam/model/data/13_withCopModulation';
% data_folder = '/Users/reimajbi/Neuro/walking/BalanceBeam/model/data/13_noCopModulation';
% data_folder = '/Users/reimajbi/Neuro/walking/BalanceBeam/model/data/13_proportionalControl';

save_file_index = '01';
save_file_label = 'reference';
% save_file_label = 'withCopModulationNeg';
% save_file_label = 'noCopModulation';
save_file_label = 'posGvs';
% save_file_label = 'alphaXiDot';
% save_file_label = 'gamma_a';
% save_file_label = 'gvs';

save_file = [save_file_index '_' save_file_label];


simulate_data       = 1;
save_results        = 1;
visualize_results   = 0;

%% set parameters
first_seed = 12e4;
number_of_repetitions = 433;
% number_of_repetitions = 32;
% number_of_repetitions = 1e4;
% number_of_repetitions = 128;
total_time = 2;
time_step = 0.001;
repetition_seeds = (1 : number_of_repetitions) + first_seed - 1;

% standard values
alphaXiDot_values = 0;
gamma_p_values = 0.0;
gamma_v_values = 0.0;
gamma_a_values = 0.0;
sigma_motor_values = 2.2608e-06;
sigma_head_pos_values = 4.9e-5;
sigma_head_vel_values = 3.7e-5;
sigma_head_acc_values = 3.7e-5;

% test values

% gamma_v_values = linspace(0, 0.05, 5);
% gamma_a_values = linspace(0, 0.01, 5);

alphaXiDot_values = 5;
gamma_p_values = 1e-1;
gamma_v_values = 1e-2;
gamma_a_values = 1e-2;

% alphaXiDot_values = [0.5 1 2 3 5];

% gamma_p_values = [1e-2 5e-2 1e-1 5e-1 1];
% gamma_v_values = [1e-2 3e-2 5e-2 8e-2 1e-1];
% gamma_a_values = [1e-3 5e-3 1e-2 3e-2 5e-2];
% gamma_a_values = logspace(-3, -2, 5);
% 
% gamma_p_values = 1e-1;
% gamma_v_values = 1e-2;
% gamma_a_values = 5e-3;

% factor = 1e1;
% sigma_head_pos_values = 4.9e-5 * factor;
% sigma_head_vel_values = 3.7e-5 * factor;
% sigma_head_acc_values = 3.7e-5 * factor;

%% simulate
if simulate_data
    tic
    number_of_dofs = 3;
    time = time_step : time_step : total_time;
    number_of_time_steps = length(time);

    x_trajectory_grid = cell(length(alphaXiDot_values), size(gamma_p_values, 3), size(gamma_v_values, 3), length(gamma_a_values), length(sigma_motor_values), length(sigma_head_pos_values), length(sigma_head_vel_values), length(sigma_head_acc_values));
    p_trajectory_grid = cell(length(alphaXiDot_values), size(gamma_p_values, 3), size(gamma_v_values, 3), length(gamma_a_values), length(sigma_motor_values), length(sigma_head_pos_values), length(sigma_head_vel_values), length(sigma_head_acc_values));
    u_trajectory_grid = cell(length(alphaXiDot_values), size(gamma_p_values, 3), size(gamma_v_values, 3), length(gamma_a_values), length(sigma_motor_values), length(sigma_head_pos_values), length(sigma_head_vel_values), length(sigma_head_acc_values));
    for i_alpha_xi_dot = 1 : length(alphaXiDot_values)
        for i_gamma_p = 1 : length(gamma_p_values)
            for i_gamma_v = 1 : length(gamma_v_values)
                for i_gamma_a = 1 : length(gamma_a_values)
                    for i_sigma_motor = 1 : length(sigma_motor_values)
                        for i_sigma_head_pos = 1 : length(sigma_head_pos_values)
                            for i_sigma_head_vel = 1 : length(sigma_head_vel_values)
                                for i_sigma_head_acc = 1 : length(sigma_head_acc_values)
                                    alpha_xi_dot_value = alphaXiDot_values(i_alpha_xi_dot);
                                    gamma_p_value = gamma_p_values(i_gamma_p);
                                    gamma_v_value = gamma_v_values(i_gamma_v);
                                    gamma_a_value = gamma_a_values(i_gamma_a);
                                    sigma_motor_value = sigma_motor_values(i_sigma_motor);
                                    sigma_head_pos_value = sigma_head_pos_values(i_sigma_head_pos);
                                    sigma_head_vel_value = sigma_head_vel_values(i_sigma_head_vel);
                                    sigma_head_acc_value = sigma_head_acc_values(i_sigma_head_acc);

                                    % simulate
                                    x_trajectories = zeros(number_of_dofs, number_of_time_steps, number_of_repetitions);
                                    p_trajectory = zeros(number_of_time_steps, number_of_repetitions);
                                    u_trajectories = zeros(number_of_time_steps, number_of_repetitions);
                                    parfor i_repetition = 1 : number_of_repetitions
                                        [ ...
                                          x_trajectories(:, :, i_repetition), ...
                                          p_trajectory(:, i_repetition), ...
                                          u_trajectories(:, i_repetition) ...
                                        ] ...
                                        = gvsModel_xiDotControlled ...
                                          ( ...
                                            alpha_xi_dot_value, ...
                                            gamma_p_value, ...
                                            gamma_v_value, ...
                                            gamma_a_value, ...
                                            sigma_motor_value, ...
                                            sigma_head_pos_value, ...
                                            sigma_head_vel_value, ...
                                            sigma_head_acc_value, ...
                                            time_step, ...
                                            total_time, ...
                                            repetition_seeds(i_repetition) ...
                                          );
                                    end
                                    x_trajectory_grid{i_alpha_xi_dot, i_gamma_p, i_gamma_v, i_gamma_a, i_sigma_motor, i_sigma_head_pos, i_sigma_head_vel, i_sigma_head_acc} = x_trajectories;
                                    p_trajectory_grid{i_alpha_xi_dot, i_gamma_p, i_gamma_v, i_gamma_a, i_sigma_motor, i_sigma_head_pos, i_sigma_head_vel, i_sigma_head_acc} = p_trajectory;
                                    u_trajectory_grid{i_alpha_xi_dot, i_gamma_p, i_gamma_v, i_gamma_a, i_sigma_motor, i_sigma_head_pos, i_sigma_head_vel, i_sigma_head_acc} = u_trajectories;
                                    
%                                     step = 100;
%                                     if floor(i_repetition/step) == (i_repetition/step)
%                                         disp(['Finished repetition ' num2str(i_repetition) ' of ' num2str(number_of_repetitions) '.']);
%                                     end
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    toc
end


%% save results
if save_results
    [~, ~, ~] = mkdir(data_folder);
    save([data_folder filesep save_file], ...
          'alphaXiDot_values', ...
          'gamma_p_values', ...
          'gamma_v_values', ...
          'gamma_a_values', ...
          'sigma_motor_values', ...
          'sigma_head_pos_values', ...
          'sigma_head_vel_values', ...
          'sigma_head_acc_values', ...
          'x_trajectory_grid', ...
          'p_trajectory_grid', ...
          'u_trajectory_grid', ...
          'time', ...
          'repetition_seeds' ...
        );
end

% visualizeSurvey_gvsModel_xiDotControlled

%% visualize_results
if visualize_results
    load '/Users/reimajbi/Dropbox/BalanceBeamPaper/responseData.mat';
    
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
    time_start = 0.506;
    time_end = 1.012;
    time_start_index = round(time_start*1000);
    
    p_trajectory_grid_squeezed = squeeze(p_trajectory_grid);
    p_means = zeros(length(time), size(p_trajectory_grid_squeezed, 1));
    p_civs = zeros(length(time), size(p_trajectory_grid_squeezed, 1));
    p_response_means = zeros(length(time), size(p_trajectory_grid_squeezed, 1));
    p_response_civs = zeros(length(time), size(p_trajectory_grid_squeezed, 1));

    p_trajectory = p_trajectory_grid_squeezed{1};
    p_mean = mean(p_trajectory, 2);
    p_std = std(p_trajectory, 0, 2);
    p_sem = std(p_trajectory, 0, 2) * 1/(size(p_trajectory, 2))^(0.5);
    p_civ = tinv(0.025, size(p_trajectory, 2)-1) * p_sem;
    
    figure; axes; hold on
    plot_human = shadedErrorBar(time_forceplate_normalized + time_start, cop_response_ml_no_mean, cop_response_ml_no_std, {'color', colors(1, :), 'linewidth', 5}, 1);
    plot_model = shadedErrorBar(time, p_mean - p_mean(time_start_index), p_std, {'color', colors(2, :), 'linewidth', 5}, 1);
    
    xlim([time_start time_end])
    legend([plot_human.mainLine plot_model.mainLine], {'human', 'model'});
    title(save_file)
end
