% optimizeAlpha
% need to determine 


% calculateInclinationAngles;
OptimalAlpha = fminunc(@myfun,3);

function RMS_error_ratio = myfun(this_alpha,varargin)

% for this alpha access a specific subject, load the data (mocap, raw
% sensor data, calculate inclination angle, 
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'subjects', [])
    addParameter(parser, 'save', false)
    addParameter(parser, 'settings', 'plotSettings.txt')
    parse(parser, varargin{:})
    subjects = parser.Results.subjects;
    
    subject_dir = cd;
    data_folder_list = determineDataStructure(subjects);
    
    root_mean_square_errors_left = cell(length(data_folder_list),3);
    root_mean_square_errors_right = cell(length(data_folder_list),3);
    trialaveraged_mocap_peak_amplitude_left = cell(length(data_folder_list),3);
    trialaveraged_mocap_peak_amplitude_right = cell(length(data_folder_list),3);
    
    for i_folder = 1 : length(data_folder_list)
        % load data
        data_path = data_folder_list{i_folder};
        cd(data_path);
        load([data_path filesep 'subjectInfo.mat'], 'date', 'subject_id', 'condition_list','trial_number_list');
%          load('subjectInfo.mat', 'condition_list','trial_number_list')

        remove_conditions = {'adaptation', 'calibration', 'emg'};
        remove_conditions_indices = find(ismember(condition_list, remove_conditions));
        condition_list(remove_conditions_indices) = [];
        trial_number_list(remove_conditions_indices) = [];
        alpha_labels = 'undefined';
        
        for i_condition = 1: length(condition_list)
            condition = condition_list{i_condition};
            trials_to_process = trial_number_list{i_condition};
            root_mean_square_errors_left{i_folder,i_condition} = zeros(1, length(trials_to_process));
            root_mean_square_errors_right{i_folder,i_condition} = zeros(1, length(trials_to_process));
            trialaveraged_mocap_peak_amplitude_left{i_folder,i_condition} = zeros(1, length(trials_to_process));
            trialaveraged_mocap_peak_amplitude_right{i_folder,i_condition} = zeros(1, length(trials_to_process));
               
            for i_trial = trials_to_process
                [inclination_angle_mocap_left_trajectory, time_marker, sampling_rate_marker, marker_labels] = loadData(date, subject_id, condition, i_trial, 'inclination_angle_mocap_left_trajectory');
                [inclination_angle_mocap_right_trajectory, time_marker, sampling_rate_marker, marker_labels] = loadData(date, subject_id, condition, i_trial, 'inclination_angle_mocap_right_trajectory');
                [inclination_angle_armsense_left_trajectories, time_marker, sampling_rate_marker, alpha_labels_left] = loadData(date, subject_id, condition, i_trial, 'inclination_angle_armsense_left_trajectories');
                [inclination_angle_armsense_right_trajectories, time_marker, sampling_rate_marker, alpha_labels_right] = loadData(date, subject_id, condition, i_trial, 'inclination_angle_armsense_right_trajectories');
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
                vertical_sensor_right_init = -normVector([right_acc_x_trajectories(1); right_acc_y_trajectories(1); right_acc_z_trajectories(1)]);
                inclination_angle_armsense_left_trajectories = zeros(length(time_marker),1);
                % change the acos to angle/atan2 when get a chance
                inclination_angle_left_rad = acos(vertical_sensor_left_init' * [0; 1; 0]);
                inclination_angle_right_rad = acos(vertical_sensor_right_init' * [0; 1; 0]);
                vertical_sensor_right_init = -normVector([right_acc_x_trajectories(1); right_acc_y_trajectories(1); right_acc_z_trajectories(1)]);
                inclination_angle_armsense_right_trajectories = zeros(length(time_marker), 1);
                
                % calculate inclination angle armsense with this_alpha
                vertical_sensor_left_dot = zeros(length(time_marker), 3);
                vertical_sensor_left_trajectory = zeros(length(time_marker), 3);
                vertical_sensor_left_trajectory(1, :) = vertical_sensor_left_init;
                vertical_sensor_right_dot = zeros(length(time_marker), 3);
                vertical_sensor_right_trajectory = zeros(length(time_marker), 3);
                vertical_sensor_right_trajectory(1, :) = vertical_sensor_right_init;
                for i_time = 2 : length(time_marker)
                    % calculate rate of change left
                    vertical_sensor_left_accel = -normVector([left_acc_x_trajectories(i_time); left_acc_y_trajectories(i_time); left_acc_z_trajectories(i_time)]);
                    angular_body_velocity_gyro = [left_gyro_x_trajectories(i_time); left_gyro_y_trajectories(i_time); left_gyro_z_trajectories(i_time);];
                    angular_velocity_body_matrix = wedgeAxis(angular_body_velocity_gyro);
                    f_g = - angular_velocity_body_matrix * vertical_sensor_left_trajectory(i_time-1, :)';
                    f_a = -this_alpha * (vertical_sensor_left_trajectory(i_time-1, :) - vertical_sensor_left_accel')';
                    vertical_sensor_left_dot(i_time, :) = f_g + f_a;
                    
                
                    % calculate rate of change right
                    vertical_sensor_right_accel = -normVector([right_acc_x_trajectories(i_time); right_acc_y_trajectories(i_time); right_acc_z_trajectories(i_time)]);
                    angular_body_velocity_gyro = [right_gyro_x_trajectories(i_time); right_gyro_y_trajectories(i_time); right_gyro_z_trajectories(i_time);];
                    angular_velocity_body_matrix = wedgeAxis(angular_body_velocity_gyro);
                    f_g = - angular_velocity_body_matrix * vertical_sensor_right_trajectory(i_time-1, :)';
                    f_a = - this_alpha * (vertical_sensor_right_trajectory(i_time-1, :) - vertical_sensor_right_accel')';
                    vertical_sensor_right_dot(i_time, :) = f_g + f_a;
                
                    % integrate
                    vertical_sensor_left_trajectory(i_time, :) = normVector(vertical_sensor_left_trajectory(i_time-1, :)' + delta_t * vertical_sensor_left_dot(i_time, :)');
                    vertical_sensor_right_trajectory(i_time, :) = normVector(vertical_sensor_right_trajectory(i_time-1, :)' + delta_t * vertical_sensor_right_dot(i_time, :)');
                end
                % calculate inclination angles
                inclination_angle_armsense_left_trajectories = rad2deg(asin(sum(vertical_sensor_left_trajectory(:, [1 3]).^2, 2).^(0.5)));
                inclination_angle_armsense_right_trajectories = rad2deg(asin(sum(vertical_sensor_right_trajectory(:, [1 3]).^2, 2).^(0.5)));
                
                % begin error calculation
                begin_time = 2;
                end_time = 118;
                [begin_marker_clip begin_marker_clip] = min(abs(time_marker-begin_time));
                [end_marker_clip end_marker_clip] = min(abs(time_marker-end_time));
                inclination_angle_mocap_left_trajectory = inclination_angle_mocap_left_trajectory(begin_marker_clip:end_marker_clip,:);
                inclination_angle_mocap_right_trajectory = inclination_angle_mocap_right_trajectory(begin_marker_clip:end_marker_clip,:);
                inclination_angle_armsense_left_trajectories = inclination_angle_armsense_left_trajectories(begin_marker_clip:end_marker_clip,:);
                inclination_angle_armsense_right_trajectories = inclination_angle_armsense_right_trajectories(begin_marker_clip:end_marker_clip,:);
                time_marker = time_marker(begin_marker_clip:end_marker_clip,:);
                
                % find peaks for average 
                [maxPeaks_left, maxPeaks_left_indices] = findpeaks(inclination_angle_mocap_left_trajectory,'MinPeakHeight',30);
                [maxPeaks_right, maxPeaks_right_indices] = findpeaks(inclination_angle_mocap_right_trajectory,'MinPeakHeight',30);
                [minPeaks_left, minPeaks_left_indices] = findpeaks(-inclination_angle_mocap_left_trajectory);
                [minPeaks_right, minPeaks_right_indices] = findpeaks(-inclination_angle_mocap_right_trajectory);
                minPeaks_left = abs(minPeaks_left);
                minPeaks_right = abs(minPeaks_right);
                trialaveraged_mocap_peak_amplitude_left{i_folder,i_condition}(i_trial) = mean(maxPeaks_left) - mean(minPeaks_left);
                trialaveraged_mocap_peak_amplitude_right{i_folder,i_condition}(i_trial) = mean(maxPeaks_right) - mean(minPeaks_right);
                
                % calculate error
                error_left = inclination_angle_mocap_left_trajectory - inclination_angle_armsense_left_trajectories;
                error_right = inclination_angle_mocap_right_trajectory - inclination_angle_armsense_right_trajectories;
                
                % root mean square error
                root_mean_square_errors_left{i_folder,i_condition}(i_trial) = mean(error_left.^2).^0.5;
                root_mean_square_errors_right{i_folder,i_condition}(i_trial) = mean(error_right.^2).^0.5;
            end
        end
        cd(subject_dir)
    end
        root_mean_square_errors_left = cell2mat(root_mean_square_errors_left);
        root_mean_square_errors_right = cell2mat(root_mean_square_errors_right);
        trialaveraged_mocap_peak_amplitude_left = cell2mat(trialaveraged_mocap_peak_amplitude_left);
        trialaveraged_mocap_peak_amplitude_right = cell2mat(trialaveraged_mocap_peak_amplitude_right);
    
        
        % do we find a % error for every subject or calculate a single %
        % error for all subject data?
        % average
        mean_root_mean_square_error_left = mean(root_mean_square_errors_left, 2);
        mean_root_mean_square_error_right = mean(root_mean_square_errors_right, 2);
        mean_mocap_peak_amplitude_left = mean(trialaveraged_mocap_peak_amplitude_left, 2);
        mean_mocap_peak_amplitude_right = mean(trialaveraged_mocap_peak_amplitude_right, 2);
        percent_error_left = (mean_root_mean_square_error_left./mean_mocap_peak_amplitude_left);
        percent_error_right = (mean_root_mean_square_error_right./mean_mocap_peak_amplitude_right);
        
        RMS_error_ratio = mean(mean([percent_error_left, percent_error_right]));
end



