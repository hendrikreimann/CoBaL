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


function mocapArmsenseComparison_AlphaOnly_AllTrials(varargin)
    % parse arguments
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    addParameter(parser, 'alpha', 1 : 8)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;
    
    load('subjectInfo.mat', 'date', 'subject_id');

    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            %% prepare
            % load data
            condition = condition_list{i_condition};
            [marker_trajectories, time_marker, sampling_rate_marker, marker_labels] = loadData(date, subject_id, condition, i_trial, 'marker_trajectories');
            [left_acc_x_trajectories, time_armsense_left, sampling_rate_armsense_left, labels] = loadData(date, subject_id, condition, i_trial, 'acc_x_left_trajectory');
            [left_acc_y_trajectories, time_armsense_left, sampling_rate_armsense_left, labels] = loadData(date, subject_id, condition, i_trial, 'acc_y_left_trajectory');
            [left_acc_z_trajectories, time_armsense_left, sampling_rate_armsense_left, labels] = loadData(date, subject_id, condition, i_trial, 'acc_z_left_trajectory');
            [left_gyro_x_trajectories, time_armsense_left, sampling_rate_armsense_left, labels] = loadData(date, subject_id, condition, i_trial, 'gyro_x_left_trajectory');
            [left_gyro_y_trajectories, time_armsense_left, sampling_rate_armsense_left, labels] = loadData(date, subject_id, condition, i_trial, 'gyro_y_left_trajectory');
            [left_gyro_z_trajectories, time_armsense_left, sampling_rate_armsense_left, labels] = loadData(date, subject_id, condition, i_trial, 'gyro_z_left_trajectory');
            [right_acc_x_trajectories, time_armsense_right, sampling_rate_armsense_right, labels] = loadData(date, subject_id, condition, i_trial, 'acc_x_right_trajectory');
            [right_acc_y_trajectories, time_armsense_right, sampling_rate_armsense_right, labels] = loadData(date, subject_id, condition, i_trial, 'acc_y_right_trajectory');
            [right_acc_z_trajectories, time_armsense_right, sampling_rate_armsense_right, labels] = loadData(date, subject_id, condition, i_trial, 'acc_z_right_trajectory');
            [right_gyro_x_trajectories, time_armsense_right, sampling_rate_armsense_right, labels] = loadData(date, subject_id, condition, i_trial, 'gyro_x_right_trajectory');
            [right_gyro_y_trajectories, time_armsense_right, sampling_rate_armsense_right, labels] = loadData(date, subject_id, condition, i_trial, 'gyro_y_right_trajectory');
            [right_gyro_z_trajectories, time_armsense_right, sampling_rate_armsense_right, labels] = loadData(date, subject_id, condition, i_trial, 'gyro_z_right_trajectory');
            
            % change time to start at 0 and scale appropriately
            time_armsense_left = (time_armsense_left - time_armsense_left(1)) * (time_marker(end) - time_marker(1)) / (time_armsense_left(end) - time_armsense_left(1)) + time_marker(1);
            time_armsense_right = (time_armsense_right - time_armsense_right(1)) * (time_marker(end) - time_marker(1)) / (time_armsense_right(end) - time_armsense_right(1)) + time_marker(1);
            
            % remove outliers
            outliers = find(diff(time_armsense_left) - median(diff(time_armsense_left)) > 2*median(diff(time_armsense_left))) + 1;
            time_armsense_left(outliers) = [];
            left_acc_x_trajectories(outliers) = [];
            left_acc_y_trajectories(outliers) = [];
            left_acc_z_trajectories(outliers) = [];
            left_gyro_x_trajectories(outliers) = [];
            left_gyro_y_trajectories(outliers) = [];
            left_gyro_z_trajectories(outliers) = [];
            
            outliers = find(diff(time_armsense_right) - median(diff(time_armsense_right)) > 2*median(diff(time_armsense_right))) + 1;
            time_armsense_right(outliers) = [];
            right_acc_x_trajectories(outliers) = [];
            right_acc_y_trajectories(outliers) = [];
            right_acc_z_trajectories(outliers) = [];
            right_gyro_x_trajectories(outliers) = [];
            right_gyro_y_trajectories(outliers) = [];
            right_gyro_z_trajectories(outliers) = [];
            
            % resample to mocap time
            left_acc_x_trajectories = spline(time_armsense_left, left_acc_x_trajectories, time_marker);
            left_acc_y_trajectories = spline(time_armsense_left, left_acc_y_trajectories, time_marker);
            left_acc_z_trajectories = spline(time_armsense_left, left_acc_z_trajectories, time_marker);
            left_gyro_x_trajectories = spline(time_armsense_left, left_gyro_x_trajectories, time_marker);
            left_gyro_y_trajectories = spline(time_armsense_left, left_gyro_y_trajectories, time_marker);
            left_gyro_z_trajectories = spline(time_armsense_left, left_gyro_z_trajectories, time_marker);
            right_acc_x_trajectories = spline(time_armsense_right, right_acc_x_trajectories, time_marker);
            right_acc_y_trajectories = spline(time_armsense_right, right_acc_y_trajectories, time_marker);
            right_acc_z_trajectories = spline(time_armsense_right, right_acc_z_trajectories, time_marker);
            right_gyro_x_trajectories = spline(time_armsense_right, right_gyro_x_trajectories, time_marker);
            right_gyro_y_trajectories = spline(time_armsense_right, right_gyro_y_trajectories, time_marker);
            right_gyro_z_trajectories = spline(time_armsense_right, right_gyro_z_trajectories, time_marker);
            
            % integrate sensor axis and calculate inclination angles
            delta_t = sampling_rate_marker^(-1);
            vertical_sensor_left_init = -normVector([left_acc_x_trajectories(1); left_acc_y_trajectories(1); left_acc_z_trajectories(1)]);
            inclination_angle_armsense_left_trajectories = zeros(length(time_marker), length(parser.Results.alpha));
            vertical_sensor_right_init = -normVector([right_acc_x_trajectories(1); right_acc_y_trajectories(1); right_acc_z_trajectories(1)]);
            inclination_angle_armsense_right_trajectories = zeros(length(time_marker), length(parser.Results.alpha));
            for i_alpha = 1 : length(parser.Results.alpha)
                vertical_sensor_left_dot = zeros(length(time_marker), 3);
                vertical_sensor_left_trajectory = zeros(length(time_marker), 3);
                vertical_sensor_left_trajectory(1, :) = vertical_sensor_left_init;
                vertical_sensor_right_dot = zeros(length(time_marker), 3);
                vertical_sensor_right_trajectory = zeros(length(time_marker), 3);
                vertical_sensor_right_trajectory(1, :) = vertical_sensor_right_init;
                for i_time = 2 : length(time_marker)
                    % calculate rate of change left
                    vertical_sensor_left_accel = -normVector([left_acc_x_trajectories(1); left_acc_y_trajectories(1); left_acc_z_trajectories(1)]);
                    angular_body_velocity_gyro = [left_gyro_x_trajectories(i_time); left_gyro_y_trajectories(i_time); left_gyro_z_trajectories(i_time);];
                    angular_velocity_body_matrix = wedgeAxis(angular_body_velocity_gyro);
                    f_g = - angular_velocity_body_matrix * vertical_sensor_left_trajectory(i_time-1, :)';
                    f_a = - parser.Results.alpha(i_alpha) * (vertical_sensor_left_trajectory(i_time-1, :) - vertical_sensor_left_accel')';
                    vertical_sensor_left_dot(i_time, :) = f_g + f_a;
                    
                    % calculate rate of change right
                    vertical_sensor_right_accel = -normVector([right_acc_x_trajectories(1); right_acc_y_trajectories(1); right_acc_z_trajectories(1)]);
                    angular_body_velocity_gyro = [right_gyro_x_trajectories(i_time); right_gyro_y_trajectories(i_time); right_gyro_z_trajectories(i_time);];
                    angular_velocity_body_matrix = wedgeAxis(angular_body_velocity_gyro);
                    f_g = - angular_velocity_body_matrix * vertical_sensor_right_trajectory(i_time-1, :)';
                    f_a = - parser.Results.alpha(i_alpha) * (vertical_sensor_right_trajectory(i_time-1, :) - vertical_sensor_right_accel')';
                    vertical_sensor_right_dot(i_time, :) = f_g + f_a;

                    % integrate
                    vertical_sensor_left_trajectory(i_time, :) = normVector(vertical_sensor_left_trajectory(i_time-1, :)' + delta_t * vertical_sensor_left_dot(i_time, :)');
                    vertical_sensor_right_trajectory(i_time, :) = normVector(vertical_sensor_right_trajectory(i_time-1, :)' + delta_t * vertical_sensor_right_dot(i_time, :)');
                end
                
                % calculate inclination angles
                inclination_angle_armsense_left_trajectories(:, i_alpha) = rad2deg(asin(sum(vertical_sensor_left_trajectory(:, [1 3]).^2, 2).^(0.5)));
                inclination_angle_armsense_right_trajectories(:, i_alpha) = rad2deg(asin(sum(vertical_sensor_right_trajectory(:, [1 3]).^2, 2).^(0.5)));
            end
            
            % calculate arm trajectories from mocap data
            LELB_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'LELB');
            LWRB_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'LWRB');
            left_arm_vector_trajectory = normVector(LWRB_trajectory' - LELB_trajectory')';
            RELB_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'RELB');
            RWRB_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'RWRB');
            right_arm_vector_trajectory = normVector(RWRB_trajectory' - RELB_trajectory')';
            
            % calculate inclination angles
            inclination_angle_mocap_left_trajectory = rad2deg(asin(sum(left_arm_vector_trajectory(:, 1:2).^2, 2).^(0.5)));
            inclination_angle_mocap_right_trajectory = rad2deg(asin(sum(right_arm_vector_trajectory(:, 1:2).^2, 2).^(0.5)));
            
            % visualize
            if visualize
                figure; axes; hold on
                plot(inclination_angle_mocap_left_trajectory, 'linewidth', 3, 'displayname', 'mocap');
                plot(inclination_angle_armsense_left_trajectories, 'linewidth', 1, 'displayname', 'mocap');
                
                figure; axes; hold on
                plot(inclination_angle_mocap_right_trajectory, 'linewidth', 3, 'displayname', 'mocap');
                plot(inclination_angle_armsense_right_trajectories, 'linewidth', 1, 'displayname', 'mocap');
            end
            
            % save
            alpha_labels = cell(size(parser.Results.alpha));
            for i_alpha = 1 : length(parser.Results.alpha)
                alpha_labels{i_alpha} = num2str(parser.Results.alpha(i_alpha));
            end
            
            save_folder = 'processed';
            save_file_name = makeFileName(date, subject_id, condition, i_trial, 'inclinationAngleTrajectories.mat');
            save ...
              ( ...
                [save_folder filesep save_file_name], ...
                'inclination_angle_mocap_left_trajectory', ...
                'inclination_angle_mocap_right_trajectory', ...
                'inclination_angle_armsense_left_trajectories', ...
                'inclination_angle_armsense_right_trajectories', ...
                'alpha_labels', ...
                'time_marker', ...
                'sampling_rate_marker' ...
              );
            addAvailableData('inclination_angle_mocap_left_trajectory', 'time_marker', 'sampling_rate_marker', '', save_folder, save_file_name);
            addAvailableData('inclination_angle_mocap_right_trajectory', 'time_marker', 'sampling_rate_marker', '', save_folder, save_file_name);
            addAvailableData('inclination_angle_armsense_left_trajectories', 'time_marker', 'sampling_rate_marker', 'alpha_labels', save_folder, save_file_name);
            addAvailableData('inclination_angle_armsense_right_trajectories', 'time_marker', 'sampling_rate_marker', 'alpha_labels', save_folder, save_file_name);
            

            disp(['Calculating inclination angles: condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' save_file_name]);


        end



    end
end
    
    





