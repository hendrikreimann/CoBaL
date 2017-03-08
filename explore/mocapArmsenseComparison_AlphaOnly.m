%% Alpha Testing 4 Armsense
%% Choose Parameters and Methods
trial_number = 2;
%% Choose Sides
    Sides2Analyze = {'Left','Right'};
%     Sides2Analyze = {'Left'};
%     Sides2Analyze = {'Right'};

%% Choose Anatomy
    AnatomicalMethods = {'Forearm'};
%     AnatomicalMethods = {'Arm'};
%% Choose Alpha Values to Investigate
    % alpha = [1:1:6];
%     alpha = [1:1:5];
%     alpha = [2:.5:5];
%     alpha = [3:.2:4];
      alpha = [3:1:4];








%% LOAD AND INITIALIZE
% load and extract marker trajectories
close all
mocap_filename = sprintf('20170227_IFS_walking_00%d_markerTrajectoriesRaw.mat',trial_number);
load(mocap_filename);
time_step_mocap = mean(diff(time_mocap));
for i_side = 1:length(Sides2Analyze) 
    as_filename_right = sprintf('20170227_IFS_walking_00%d_armSense%sData.mat',trial_number,Sides2Analyze{i_side});
    load(as_filename_right);
    time_armsense = [0;time]*1000';

    % time_armsense_right = time * 1000;
    % time_armsense_right = linspace(0,time_mocap(end),length(gyro_x_trajectory));
    gyroscope_x_raw = [0;gyro_x_trajectory]';
    gyroscope_y_raw = [0;gyro_y_trajectory]';
    gyroscope_z_raw = [0;gyro_z_trajectory]';
    accelerometer_x_raw = [0;acc_x_trajectory]';
    accelerometer_y_raw = [0;acc_y_trajectory]';
    accelerometer_z_raw = [0;acc_z_trajectory]';

    outliers = find(diff(time_armsense) - median(diff(time_armsense)) > 2*median(diff(time_armsense))) + 1;
    time_armsense(outliers) = [];

    gyroscope_x_raw(outliers) = [];
    gyroscope_y_raw(outliers) = [];
    gyroscope_z_raw(outliers) = [];
    accelerometer_x_raw(outliers) = [];
    accelerometer_y_raw(outliers) = [];
    accelerometer_z_raw(outliers) = [];

    gyroscope_x = interp1(time_armsense, gyroscope_x_raw, time_mocap);
    gyroscope_y = interp1(time_armsense, gyroscope_y_raw, time_mocap);
    gyroscope_z = interp1(time_armsense, gyroscope_z_raw, time_mocap);
    accelerometer_x = interp1(time_armsense, accelerometer_x_raw, time_mocap);
    accelerometer_y = interp1(time_armsense, accelerometer_y_raw, time_mocap);
    accelerometer_z = interp1(time_armsense, accelerometer_z_raw, time_mocap);

%% ARMSENSE INCLINATION ANGLE    
    for i_alpha = 1:length(alpha)
        % integrate direction of vertical from gyroscope readings

        vertical_sensor_init = normVector([-accelerometer_x(1); -accelerometer_y(1); -accelerometer_z(1)]);
        vertical_sensor_dot = zeros(length(time_mocap), 3);
        vertical_sensor_trajectory = zeros(length(time_mocap), 3);
        vertical_sensor_trajectory(1, :) = vertical_sensor_init;
        vertical_sensor_trajectory_accel = zeros(length(time_mocap), 3);
        vertical_sensor_trajectory_accel(1, :) = normVector([accelerometer_x(1), accelerometer_y(1), accelerometer_z(1)]);
        for i_time = 2 : length(time_mocap)
            vertical_sensor_trajectory_accel(i_time, :) = normVector([accelerometer_x(i_time), accelerometer_y(i_time), accelerometer_z(i_time)]);

            angular_body_velocity_gyro = [gyroscope_x(i_time); gyroscope_y(i_time); gyroscope_z(i_time);];
            angular_velocity_body_matrix = wedgeAxis(angular_body_velocity_gyro);
            f_g = - angular_velocity_body_matrix * vertical_sensor_trajectory(i_time-1, :)';
            f_a = - alpha(i_alpha) * (vertical_sensor_trajectory(i_time-1, :) - vertical_sensor_trajectory_accel(i_time-1, :))';

            vertical_sensor_dot(i_time, :) = f_g + f_a;

            % integrate
            vertical_sensor_trajectory(i_time, :) = normVector(vertical_sensor_trajectory(i_time-1, :) + time_step_mocap * vertical_sensor_dot(i_time, :));
        end

        % Calculate Inclination Angle
        inclination_angle_gyro_trajectory = zeros(length(time_mocap), 1);
        inclination_angle_accel_trajectory = zeros(length(time_mocap), 1);
        inclination_angle_combined_trajectory = zeros(length(time_mocap), 1);
        inclination_angle_combined_trajectory_right = zeros(length(time_mocap), 1);
        inclination_angle_combined_trajectory_left = zeros(length(time_mocap), 1);
        sensor_axis_sensor = [0; 1; 0];
        for i_time = 1 : length(time_mocap)
            inclination_angle_accel_trajectory(i_time) = rad2deg(subspace(vertical_sensor_trajectory_accel(i_time, :)', sensor_axis_sensor));
            inclination_angle_combined_trajectory(i_time) = rad2deg(subspace(vertical_sensor_trajectory(i_time, :)', sensor_axis_sensor));
        end

%% MOCAP INCLINATION ANGLE

        if strcmp('Right',Sides2Analyze{i_side})
            marker_index = find(strcmp(marker_labels, 'RSHO'), 1, 'last');
            coordinate_indices = [1 2 3] + (marker_index-1)*3;
            shoulder_marker_trajectory = marker_trajectories_raw(:, coordinate_indices);
            
            marker_index = find(strcmp(marker_labels, 'RELB'), 1, 'last');
            coordinate_indices = [1 2 3] + (marker_index-1)*3;
            elbow_marker_trajectory = marker_trajectories_raw(:, coordinate_indices);
            
            marker_index = find(strcmp(marker_labels, 'RWRB'), 1, 'last');
            coordinate_indices = [1 2 3] + (marker_index-1)*3;
            medialWrist_marker_trajectory = marker_trajectories_raw(:, coordinate_indices);
            
            marker_index = find(strcmp(marker_labels, 'RWRA'), 1, 'last');
            coordinate_indices = [1 2 3] + (marker_index-1)*3;
            lateralWrist_marker_trajectory = marker_trajectories_raw(:, coordinate_indices);
 
        elseif strcmp('Left',Sides2Analyze{i_side})
            marker_index = find(strcmp(marker_labels, 'LSHO'), 1, 'last');
            coordinate_indices = [1 2 3] + (marker_index-1)*3;
            shoulder_marker_trajectory = marker_trajectories_raw(:, coordinate_indices);
            
            marker_index = find(strcmp(marker_labels, 'LELB'), 1, 'last');
            coordinate_indices = [1 2 3] + (marker_index-1)*3;
            elbow_marker_trajectory = marker_trajectories_raw(:, coordinate_indices);
            
            marker_index = find(strcmp(marker_labels, 'LWRB'), 1, 'last');
            coordinate_indices = [1 2 3] + (marker_index-1)*3;
            medialWrist_marker_trajectory = marker_trajectories_raw(:, coordinate_indices);
            
            marker_index = find(strcmp(marker_labels, 'LWRA'), 1, 'last');
            coordinate_indices = [1 2 3] + (marker_index-1)*3;
            lateralWrist_marker_trajectory = marker_trajectories_raw(:, coordinate_indices);
        end

        vertical_vector = [0; 0; 1];
        inclination_angle_mocap_trajectory = [];
%         proximal_marker_trajectory;
        for i_time = 1 : length(time_mocap) 
            
            if strcmp('Forearm',AnatomicalMethods)     
%                 arm_vector = normVector(((medialWrist_marker_trajectory(i_time,:) + lateralWrist_marker_trajectory(i_time,:)) / 2) - elbow_marker_trajectory(i_time,:));
                arm_vector = normVector(medialWrist_marker_trajectory(i_time,:) - elbow_marker_trajectory(i_time,:));
            elseif strcmp('Arm',AnatomicalMethods)
                arm_vector = normVector(((medialWrist_marker_trajectory(i_time,:) + lateralWrist_marker_trajectory(i_time,:)) / 2) - shoulder_marker_trajectory(i_time,:));
            end
            inclination_angle_mocap_trajectory(i_time) = rad2deg(subspace(arm_vector',vertical_vector));                
        end
        
        rms_check = rms(inclination_angle_mocap_trajectory(10*250:110*250)' - inclination_angle_combined_trajectory(10*250:110*250));
%         mean_signals = (inclination_angle_mocap_trajectory' + inclination_angle_combined_trajectory) / 2;
%         signals = [inclination_angle_mocap_trajectory' inclination_angle_combined_trajectory];
%         
%         % still need to set frame cutoffs...
%         for i_signal = 1:2
%            for i_time = 1:length(time_mocap)
%                CMC_Numerator(i_time,i_signal) = (signals(i_time,i_signal) - mean_signals(i_time))^2 / length(time_mocap);
%                CMC_Denominator(i_time,i_signal) = (signals(i_time,i_signal) - mean_signals)^2 / (2*length(time_mocap) - 1); 
%            end
%         end
%         CMC_Numerator = sum(reshape(CMC_Numerator,[1,2*length(time_mocap)]));
%         CMC_Denominator = sum(reshape(CMC_Denominator,[1,2*length(time_mocap)]));
%         CMC = sqrt(1 - (CMC_Numerator/CMC_Denominator));
        
%% PLOT VARIABLES
        figure; axes; hold on
        plot(time_mocap, inclination_angle_combined_trajectory, 'linewidth', 2, 'displayname', 'angle from AS')
        plot(time_mocap, inclination_angle_mocap_trajectory, 'linewidth', 2, 'displayname', 'angle from mocap')
        plot(time_mocap, inclination_angle_accel_trajectory, 'linewidth', 2, 'displayname', 'angle from accel')
        legend('toggle')
        xlim([10 20]);
        xlabel('Degrees');
        

        Arm_title = ['Arm Side = ' num2str(Sides2Analyze{i_side})]; 
        Alpha_title = ['Alpha Value = ' num2str(alpha(i_alpha))];
        RMS_title = ['RMS Error(deg) = ' num2str(rms_errors)];
        RMS_check = ['RMS Error(deg) = ' num2str(rms_check)];
%         CMC_title = ['CMC = ' num2str(CMC)];
        title({Arm_title, Alpha_title, RMS_title, CMC_title});
        
    end
end
handles=findall(0,'type','axes');
linkaxes([handles],'xy');
distFig('rows', 2)
