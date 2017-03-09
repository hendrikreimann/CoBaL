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


function analyzeAlpha(varargin)

    % parse arguments
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    addParameter(parser, 'alpha', 1 : 5)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;
    
    load('subjectInfo.mat', 'date', 'subject_id');

    % calculate RMS
    alpha_labels = 'undefined';
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        root_mean_square_errors_left = zeros(1, length(trials_to_process));
        root_mean_square_errors_right = zeros(1, length(trials_to_process));
        coefficients_of_multiple_correlation_left = zeros(1, length(trials_to_process));
        coefficients_of_multiple_correlation_right = zeros(1, length(trials_to_process));
        for i_trial = trials_to_process
            %% prepare
            % load data
            condition = condition_list{i_condition};
            [inclination_angle_mocap_left_trajectory, time_marker, sampling_rate_marker, marker_labels] = loadData(date, subject_id, condition, i_trial, 'inclination_angle_mocap_left_trajectory');
            [inclination_angle_mocap_right_trajectory, time_marker, sampling_rate_marker, marker_labels] = loadData(date, subject_id, condition, i_trial, 'inclination_angle_mocap_right_trajectory');
            [inclination_angle_armsense_left_trajectories, time_marker, sampling_rate_marker, alpha_labels_left] = loadData(date, subject_id, condition, i_trial, 'inclination_angle_armsense_left_trajectories');
            [inclination_angle_armsense_right_trajectories, time_marker, sampling_rate_marker, alpha_labels_right] = loadData(date, subject_id, condition, i_trial, 'inclination_angle_armsense_right_trajectories');
            if strcmp(alpha_labels, 'undefined')
                alpha_labels = alpha_labels_left;
            end
            if ~isequal(alpha_labels, alpha_labels_left) || ~isequal(alpha_labels, alpha_labels_right)
                error('Alpha values are different between trials, aborting.')
            end
            
            for i_alpha = 1 : length(alpha_labels)
                % calculate error
                error_left = inclination_angle_mocap_left_trajectory - inclination_angle_armsense_left_trajectories(:, i_alpha);
                error_right = inclination_angle_mocap_right_trajectory - inclination_angle_armsense_right_trajectories(:, i_alpha);
                
                % root mean square error
                root_mean_square_errors_left(i_alpha, i_trial) = mean(error_left.^2).^0.5;
                root_mean_square_errors_right(i_alpha, i_trial) = mean(error_right.^2).^0.5;
                
                % coefficient of correlation left
                number_of_time_steps_left = length(inclination_angle_mocap_left_trajectory);
                number_of_signals_left = 2;
                y_mc_left = inclination_angle_mocap_left_trajectory;
                y_as_left = inclination_angle_armsense_left_trajectories(:, i_alpha);
                y_left_mean_trajectory = (y_mc_left + y_as_left) / 2;
                y_left_grand_mean = mean(y_left_mean_trajectory);
                
                y_mc_left_mean_free_trajectory = y_mc_left - y_left_mean_trajectory;
                y_mc_left_deviation_squared_sum = sum(y_mc_left_mean_free_trajectory.^2);
                y_as_left_mean_free_trajectory = y_as_left - y_left_mean_trajectory;
                y_as_left_deviation_squared_sum = sum(y_as_left_mean_free_trajectory.^2);
                numerator_left = (y_mc_left_deviation_squared_sum + y_as_left_deviation_squared_sum) / (number_of_time_steps_left) * (number_of_signals_left - 1);
                
                y_mc_left_grand_mean_free_trajectory = y_mc_left - y_left_grand_mean;
                y_mc_left_grand_deviation_squared_sum = sum(y_mc_left_grand_mean_free_trajectory.^2);
                y_as_left_grand_mean_free_trajectory = y_as_left - y_left_mean_trajectory;
                y_as_left_grand_deviation_squared_sum = sum(y_as_left_grand_mean_free_trajectory.^2);
                denominator_left = (y_mc_left_grand_deviation_squared_sum + y_as_left_grand_deviation_squared_sum) / (number_of_time_steps_left * number_of_signals_left - 1);
                
                coefficients_of_multiple_correlation_left(i_alpha, i_trial) = sqrt(1 - numerator_left/denominator_left);

                % coefficient of correlation right
                number_of_time_steps_right = length(inclination_angle_mocap_right_trajectory);
                number_of_signals_right = 2;
                y_mc_right = inclination_angle_mocap_right_trajectory;
                y_as_right = inclination_angle_armsense_right_trajectories(:, i_alpha);
                y_right_mean_trajectory = (y_mc_right + y_as_right) / 2;
                y_right_grand_mean = mean(y_right_mean_trajectory);
                
                y_mc_right_mean_free_trajectory = y_mc_right - y_right_mean_trajectory;
                y_mc_right_deviation_squared_sum = sum(y_mc_right_mean_free_trajectory.^2);
                y_as_right_mean_free_trajectory = y_as_right - y_right_mean_trajectory;
                y_as_right_deviation_squared_sum = sum(y_as_right_mean_free_trajectory.^2);
                numerator_right = (y_mc_right_deviation_squared_sum + y_as_right_deviation_squared_sum) / (number_of_time_steps_right) * (number_of_signals_right - 1);
                
                y_mc_right_grand_mean_free_trajectory = y_mc_right - y_right_grand_mean;
                y_mc_right_grand_deviation_squared_sum = sum(y_mc_right_grand_mean_free_trajectory.^2);
                y_as_right_grand_mean_free_trajectory = y_as_right - y_right_mean_trajectory;
                y_as_right_grand_deviation_squared_sum = sum(y_as_right_grand_mean_free_trajectory.^2);
                denominator_right = (y_mc_right_grand_deviation_squared_sum + y_as_right_grand_deviation_squared_sum) / (number_of_time_steps_right * number_of_signals_right - 1);
                
                coefficients_of_multiple_correlation_right(i_alpha, i_trial) = sqrt(1 - numerator_right/denominator_right);
            end
        end
    end
    
    % average
    mean_root_mean_square_error_left = mean(root_mean_square_errors_left, 2);
    mean_root_mean_square_error_right = mean(root_mean_square_errors_right, 2);
    mean_coefficient_of_multiple_correlation_left = mean(coefficients_of_multiple_correlation_left, 2);
    mean_coefficient_of_multiple_correlation_right = mean(coefficients_of_multiple_correlation_right, 2);
    
    % visualize
    figure; axes; hold on; title('RMS')
    plot(mean_root_mean_square_error_left, 'x-', 'displayname', 'left')
    plot(mean_root_mean_square_error_right, 'x-', 'displayname', 'right')
    legend('show')
    
    % visualize
    figure; axes; hold on; title('CMC')
    plot(mean_coefficient_of_multiple_correlation_left, 'x-', 'displayname', 'left')
    plot(mean_coefficient_of_multiple_correlation_right, 'x-', 'displayname', 'right')
    legend('show')
    
    
end













