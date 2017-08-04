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
    addParameter(parser, 'alpha', 1 : 8)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;
    
    load('subjectInfo.mat', 'date', 'subject_id');

    % calculate RMS
    alpha_labels = 'undefined';
    
    root_mean_square_errors_left = cell(length(condition_list),1);
    root_mean_square_errors_right = cell(length(condition_list),1);
    coefficients_of_multiple_correlation_left = cell(length(condition_list),1);
    coefficients_of_multiple_correlation_right = cell(length(condition_list),1);
    trialaveraged_mocap_peak_amplitude_left = cell(length(condition_list),1);
    trialaveraged_mocap_peak_amplitude_right = cell(length(condition_list),1);
    angular_velocity_left = cell(length(condition_list),1);
    angular_velocity_right = cell(length(condition_list),1);
        
    adjust_i_condition = 0;
    for i_condition = 1: length(condition_list)
        condition = condition_list{i_condition};
            
        % remove the adaptation trials from armsense analysis
        if strcmp(condition,"adaptation")
            length_of_conditions = length(condition_list) - 1;
            adjust_i_condition = 1;
            
            root_mean_square_errors_left = cell(length_of_conditions,1);
            root_mean_square_errors_right = cell(length_of_conditions,1);
            coefficients_of_multiple_correlation_left = cell(length_of_conditions,1);
            coefficients_of_multiple_correlation_right = cell(length_of_conditions,1);
            trialaveraged_mocap_peak_amplitude_left = cell(length_of_conditions,1);
            trialaveraged_mocap_peak_amplitude_right = cell(length_of_conditions,1);
            continue
        end
        adjusted_i_condition = i_condition;
        if adjust_i_condition
            adjusted_i_condition = i_condition - 1;
        end
       
        trials_to_process = trial_number_list{i_condition};
        root_mean_square_errors_left{adjusted_i_condition} = zeros(1, length(trials_to_process));
        root_mean_square_errors_right{adjusted_i_condition} = zeros(1, length(trials_to_process));
        coefficients_of_multiple_correlation_left{adjusted_i_condition} = zeros(1, length(trials_to_process));
        coefficients_of_multiple_correlation_right{adjusted_i_condition} = zeros(1, length(trials_to_process));
        angular_velocity_left{adjusted_i_condition} = zeros(1,length(trials_to_process));
        angular_velocity_right{adjusted_i_condition} = zeros(1,length(trials_to_process));
        
        for i_trial = trials_to_process
            %% prepare
            % load data
            [inclination_angle_mocap_left_trajectory, time_marker, sampling_rate_marker, marker_labels] = loadData(date, subject_id, condition, i_trial, 'inclination_angle_mocap_left_trajectory');
            [inclination_angle_mocap_right_trajectory, time_marker, sampling_rate_marker, marker_labels] = loadData(date, subject_id, condition, i_trial, 'inclination_angle_mocap_right_trajectory');
            [inclination_angle_armsense_left_trajectories, time_marker, sampling_rate_marker, alpha_labels_left] = loadData(date, subject_id, condition, i_trial, 'inclination_angle_armsense_left_trajectories');
            [inclination_angle_armsense_right_trajectories, time_marker, sampling_rate_marker, alpha_labels_right] = loadData(date, subject_id, condition, i_trial, 'inclination_angle_armsense_right_trajectories');
            
            % clip data due to missing marker data in beginning and end of
            % trials
            inclination_angle_mocap_left_trajectory = inclination_angle_mocap_left_trajectory(500:29500,:);
            inclination_angle_mocap_right_trajectory = inclination_angle_mocap_right_trajectory(500:29500,:);
            inclination_angle_armsense_left_trajectories = inclination_angle_armsense_left_trajectories(500:29500,:);
            inclination_angle_armsense_right_trajectories = inclination_angle_armsense_right_trajectories(500:29500,:);
            time_marker = time_marker(500:29500,:);
            % find velocity of arm movement via mocap
            angular_velocity_left{adjusted_i_condition}(i_trial) = mean(abs(diff(inclination_angle_mocap_left_trajectory)./diff(time_marker)));
            angular_velocity_right{adjusted_i_condition}(i_trial) = mean(abs(diff(inclination_angle_mocap_right_trajectory)./diff(time_marker)));
            
            if strcmp(alpha_labels, 'undefined')
                alpha_labels = alpha_labels_left;
            end
            if ~isequal(alpha_labels, alpha_labels_left) || ~isequal(alpha_labels, alpha_labels_right)
                error('Alpha values are different between trials, aborting.')
            end
            
            % find peaks for average 
            [maxPeaks_left, maxPeaks_left_indices] = findpeaks(inclination_angle_mocap_left_trajectory,'MinPeakHeight',30);
            [maxPeaks_right, maxPeaks_right_indices] = findpeaks(inclination_angle_mocap_right_trajectory,'MinPeakHeight',30);
            [minPeaks_left, minPeaks_left_indices] = findpeaks(-inclination_angle_mocap_left_trajectory);
            [minPeaks_right, minPeaks_right_indices] = findpeaks(-inclination_angle_mocap_right_trajectory);
            minPeaks_left = abs(minPeaks_left);
            minPeaks_right = abs(minPeaks_right);
            trialaveraged_mocap_peak_amplitude_left{adjusted_i_condition}(i_trial) = mean(maxPeaks_left) - mean(minPeaks_left);
            trialaveraged_mocap_peak_amplitude_right{adjusted_i_condition}(i_trial) = mean(maxPeaks_right) - mean(minPeaks_right);
            
%             figure; hold on;
%             plot(inclination_angle_mocap_left_trajectory,'-')
%             plot(maxPeaks_left_indices,maxPeaks_left,'*')
%             plot(minPeaks_left_indices,minPeaks_left,'*')
%             plot(inclination_angle_mocap_right_trajectory,'-')
%             plot(maxPeaks_right_indices,maxPeaks_right,'*')
%             plot(minPeaks_right_indices,minPeaks_right,'*')
%             hold off;
%             close all;
            
            for i_alpha = 1 : length(alpha_labels)
                % calculate error
                error_left = inclination_angle_mocap_left_trajectory - inclination_angle_armsense_left_trajectories(:, i_alpha);
                error_right = inclination_angle_mocap_right_trajectory - inclination_angle_armsense_right_trajectories(:, i_alpha);
                
                % root mean square error
                root_mean_square_errors_left{adjusted_i_condition}(i_alpha, i_trial) = mean(error_left.^2).^0.5;
                root_mean_square_errors_right{adjusted_i_condition}(i_alpha, i_trial) = mean(error_right.^2).^0.5;
                
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
                numerator_left = (y_mc_left_deviation_squared_sum + y_as_left_deviation_squared_sum) / (number_of_time_steps_left * (number_of_signals_left - 1));
                
                y_mc_left_grand_mean_free_trajectory = y_mc_left - y_left_grand_mean;
                y_mc_left_grand_deviation_squared_sum = sum(y_mc_left_grand_mean_free_trajectory.^2);
                y_as_left_grand_mean_free_trajectory = y_as_left - y_left_grand_mean;
                y_as_left_grand_deviation_squared_sum = sum(y_as_left_grand_mean_free_trajectory.^2);
                denominator_left = (y_mc_left_grand_deviation_squared_sum + y_as_left_grand_deviation_squared_sum) / (number_of_time_steps_left * number_of_signals_left - 1);
                
                coefficients_of_multiple_correlation_left{adjusted_i_condition}(i_alpha, i_trial) = sqrt(1 - numerator_left/denominator_left);

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
                numerator_right = (y_mc_right_deviation_squared_sum + y_as_right_deviation_squared_sum) / (number_of_time_steps_right * (number_of_signals_right - 1));
                
                y_mc_right_grand_mean_free_trajectory = y_mc_right - y_right_grand_mean;
                y_mc_right_grand_deviation_squared_sum = sum(y_mc_right_grand_mean_free_trajectory.^2);
                y_as_right_grand_mean_free_trajectory = y_as_right - y_right_grand_mean;
                y_as_right_grand_deviation_squared_sum = sum(y_as_right_grand_mean_free_trajectory.^2);
                denominator_right = (y_mc_right_grand_deviation_squared_sum + y_as_right_grand_deviation_squared_sum) / (number_of_time_steps_right * number_of_signals_right - 1);
                
                coefficients_of_multiple_correlation_right{adjusted_i_condition}(i_alpha, i_trial) = sqrt(1 - numerator_right/denominator_right);
                
                % alternative calculation
                y_left_jt = [inclination_angle_mocap_left_trajectory'; inclination_angle_armsense_left_trajectories(:, i_alpha)'];
                N = size(y_left_jt, 1);
                T = size(y_left_jt, 2);
                numerator_left = sum(sum((y_left_jt - repmat(mean(y_left_jt, 1), N, 1)).^2, 2), 1) / (T * (N - 1));
                denominator_left = sum(sum((y_left_jt - mean(mean(y_left_jt, 1))).^2, 2), 1) / (T * N - 1);
                coefficients_of_multiple_correlation_left{adjusted_i_condition}(i_alpha, i_trial) = sqrt(1 - numerator_left/denominator_left);
                
                y_right_jt = [inclination_angle_mocap_right_trajectory'; inclination_angle_armsense_right_trajectories(:, i_alpha)'];
                N = size(y_right_jt, 1);
                T = size(y_right_jt, 2);
                numerator_right = sum(sum((y_right_jt - repmat(mean(y_right_jt, 1), N, 1)).^2, 2), 1) / (T * (N - 1));
                denominator_right = sum(sum((y_right_jt - mean(mean(y_right_jt, 1))).^2, 2), 1) / (T * N - 1);
                coefficients_of_multiple_correlation_right{adjusted_i_condition}(i_alpha, i_trial) = sqrt(1 - numerator_right/denominator_right);
                
            end
        end
    end
    
    % cell2mat
    root_mean_square_errors_left = cell2mat(root_mean_square_errors_left');
    root_mean_square_errors_right = cell2mat(root_mean_square_errors_right');
    coefficients_of_multiple_correlation_left = cell2mat(coefficients_of_multiple_correlation_left');
    coefficients_of_multiple_correlation_right = cell2mat(coefficients_of_multiple_correlation_right');
    trialaveraged_mocap_peak_amplitude_left = cell2mat(trialaveraged_mocap_peak_amplitude_left');
    trialaveraged_mocap_peak_amplitude_right = cell2mat(trialaveraged_mocap_peak_amplitude_right');
    angular_velocity_left = cell2mat(angular_velocity_left');
    angular_velocity_right = cell2mat(angular_velocity_right');
    
    % average
    mean_root_mean_square_error_left = mean(root_mean_square_errors_left, 2);
    mean_root_mean_square_error_right = mean(root_mean_square_errors_right, 2);
    mean_mocap_peak_amplitude_left = mean(trialaveraged_mocap_peak_amplitude_left);
    mean_mocap_peak_amplitude_right = mean(trialaveraged_mocap_peak_amplitude_right);
    percent_error_left = (mean_root_mean_square_error_left/mean_mocap_peak_amplitude_left)*100;
    percent_error_right = (mean_root_mean_square_error_right/mean_mocap_peak_amplitude_right)*100;
    std_root_mean_square_error_left = std(root_mean_square_errors_left,0, 2);
    std_root_mean_square_error_right = std(root_mean_square_errors_right,0, 2);
    mean_coefficient_of_multiple_correlation_left = mean(coefficients_of_multiple_correlation_left, 2);
    mean_coefficient_of_multiple_correlation_right = mean(coefficients_of_multiple_correlation_right, 2);
    std_coefficient_of_multiple_correlation_left = std(coefficients_of_multiple_correlation_left,0,2);
    std_coefficient_of_multiple_correlation_right = std(coefficients_of_multiple_correlation_right,0,2);
    mean_angular_velocity_left = mean(angular_velocity_left);
    mean_angular_velocity_right = mean(angular_velocity_right);
    std_angular_velocity_left = std(angular_velocity_left);
    std_angular_velocity_right = std(angular_velocity_right);
  
    
    % visualize
    alpha_values = zeros(size(alpha_labels));
    for i_alpha = 1 : length(alpha_values)
        alpha_values(i_alpha) = str2double(alpha_labels{i_alpha});
    end
    
    figure; axes; hold on; title('RMS')
    errorbar(alpha_values, mean_root_mean_square_error_left, std_root_mean_square_error_left, 'o-', 'displayname', 'left')
    errorbar(alpha_values, mean_root_mean_square_error_right, std_root_mean_square_error_right, 'o-', 'displayname', 'right')
    set(gca, 'xtick', alpha_values)
    set(gca, 'xticklabels', alpha_labels)
    xlabel('alpha')
    ylabel('RMS (deg)')
    legend('show')
    
    figure; axes; hold on; title('% Error')
    plot(alpha_values, percent_error_left, 'x-', 'displayname', 'left')
    plot(alpha_values, percent_error_right, 'x-', 'displayname', 'right')
    xlabel('Alpha')
    ylabel('% Error (Relative to Avg Mocap Amplitude)')
    legend('show')


    % visualize
    figure; axes; hold on; title('CMC')
    errorbar(alpha_values, mean_coefficient_of_multiple_correlation_left, std_coefficient_of_multiple_correlation_left, 'o-', 'displayname', 'left')
    errorbar(alpha_values, mean_coefficient_of_multiple_correlation_right, std_coefficient_of_multiple_correlation_right, 'o-', 'displayname', 'right')
    set(gca, 'xtick', alpha_values)
    set(gca, 'xticklabels', alpha_labels)
    xlabel('alpha')
    ylabel('CMC')
    legend('show')
    
    results_file_name = ['analysis' filesep makeFileName(date, subject_id, 'validation')];
    save ...
      ( ...
        results_file_name, ...
        'mean_coefficient_of_multiple_correlation_left',...
        'mean_coefficient_of_multiple_correlation_right',...
        'std_coefficient_of_multiple_correlation_left', ...
        'std_coefficient_of_multiple_correlation_right', ...
        'mean_root_mean_square_error_left', ...
        'mean_root_mean_square_error_right', ...
        'std_root_mean_square_error_left', ...
        'std_root_mean_square_error_right', ...
        'mean_angular_velocity_left', ...
        'mean_angular_velocity_right', ...
        'std_angular_velocity_left', ...
        'std_angular_velocity_right' ...
      )
    
    
end













