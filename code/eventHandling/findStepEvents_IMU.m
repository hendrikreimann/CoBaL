% this function finds turning events

% input: 
% subjectInfo.mat
% Gyroscope Trajectories 
% turn start and end in subjectSettings
%
% output:
% file stepEvents.mat, containing
% - turn start times
% - turn end times


function findStepEvents_IMU(varargin)

    % parse arguments
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'APDM', true)
    addParameter(parser, 'plots', false)
    parse(parser, varargin{:})

    % figure out folders
    if ~exist('analysis', 'dir')
        mkdir('analysis')
    end
    
    % load settings
    subject_settings = loadSettingsFromFile('subject');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
    
    
    % load turn identification settings
    min_peak_distance_threshold = subject_settings.get('midturn_peak_distance_threshold'); 
    min_peak_height_threshold = subject_settings.get('midturn_peak_amplitude_threshold');
    start_turn_modifier = subject_settings.get('start_turn_modifier');
    end_turn_modifier = subject_settings.get('end_turn_modifier');
    
    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            %% prepare
            % load data
            condition = condition_list{i_condition};
            [lumbar_x_gyro, sensor_time, sampling_rate, sensor_label] = loadData(collection_date, subject_id, condition, i_trial, 'lumbar_x_gyro');
%             [lumbar_y_gyro, sensor_time, sampling_rate, sensor_label] = loadData(collection_date, subject_id, condition, i_trial, 'lumbar_y_gyro');
%             [lumbar_z_gyro, sensor_time, sampling_rate, sensor_label] = loadData(collection_date, subject_id, condition, i_trial, 'lumbar_z_gyro');
%             [left_foot_x_accel, sensor_time, sampling_rate, sensor_label] = loadData(collection_date, subject_id, condition, i_trial, 'left_foot_x_accel');
%             [left_foot_y_accel, sensor_time, sampling_rate, sensor_label] = loadData(collection_date, subject_id, condition, i_trial, 'left_foot_y_accel');
%             [left_foot_z_accel, sensor_time, sampling_rate, sensor_label] = loadData(collection_date, subject_id, condition, i_trial, 'left_foot_z_accel');
%             [right_foot_x_accel, sensor_time, sampling_rate, sensor_label] = loadData(collection_date, subject_id, condition, i_trial, 'right_foot_x_accel');
%             [right_foot_y_accel, sensor_time, sampling_rate, sensor_label] = loadData(collection_date, subject_id, condition, i_trial, 'right_foot_y_accel');
%             [right_foot_z_accel, sensor_time, sampling_rate, sensor_label] = loadData(collection_date, subject_id, condition, i_trial, 'right_foot_z_accel');
%             
%             left_foot_total_accel = (left_foot_x_accel.^2 + left_foot_y_accel.^2 + left_foot_y_accel.^2).^(0.5);
%             right_foot_total_accel = (right_foot_x_accel.^2 + right_foot_y_accel.^2 + right_foot_y_accel.^2).^(0.5);
%             
%             figure; axes_left = axes; hold on
%             plot(sensor_time, left_foot_x_accel, 'linewidth', 2, 'displayname', 'left x')
%             plot(sensor_time, left_foot_y_accel, 'linewidth', 2, 'displayname', 'left y')
%             plot(sensor_time, left_foot_z_accel, 'linewidth', 2, 'displayname', 'left z')
%             plot(sensor_time, left_foot_total_accel, 'linewidth', 2, 'displayname', 'left total')
%             legend('show')
%             
%             figure; axes_right = axes; hold on
%             plot(sensor_time, right_foot_x_accel, 'linewidth', 2, 'displayname', 'right x')
%             plot(sensor_time, right_foot_y_accel, 'linewidth', 2, 'displayname', 'right y')
%             plot(sensor_time, right_foot_z_accel, 'linewidth', 2, 'displayname', 'right z')
%             plot(sensor_time, right_foot_total_accel, 'linewidth', 2, 'displayname', 'right total')
%             legend('show')
%             
%             figure; axes_lumbar = axes; hold on
%             plot(sensor_time, lumbar_x_gyro, 'linewidth', 2, 'displayname', 'lumbar gyro x')
%             plot(sensor_time, lumbar_y_gyro, 'linewidth', 2, 'displayname', 'lumbar gyro y')
%             plot(sensor_time, lumbar_z_gyro, 'linewidth', 2, 'displayname', 'lumbar gyro z')
%             legend('show')
%             
%             figure; axes_both = axes; hold on
%             plot(sensor_time, lumbar_x_gyro * 10, 'linewidth', 2, 'displayname', 'lumbar gyro')
%             plot(sensor_time, left_foot_total_accel, 'linewidth', 2, 'displayname', 'left total')
%             plot(sensor_time, right_foot_total_accel, 'linewidth', 2, 'displayname', 'right total')
%             legend('show')
%             
%             linkaxes([axes_left, axes_right, axes_both axes_lumbar], 'x')
            
            %% find turn events
            turn_start_times = [];
            turn_end_times = [];
            
            if any(strcmp(subject_settings.get('midturn_method', 1), 'lumbar_x_peak'))
                % get pushoff and touchdown indices
               this_absolute_data = abs(lumbar_x_gyro);
               
               % Plots the turns to spot check the thresholds are correct.
                   % can comment this out or make it as a parser...
               figure; grid off;
               set(gcf, 'Position',[ 250 10 1000 1000]); 
               findpeaks(this_absolute_data, 128, 'MinPeakDistance', min_peak_distance_threshold, 'MinPeakHeight', min_peak_height_threshold)
               title(['Turns - ' condition])
               xlabel('Time')
               grid off; 
               ylim([0,9]);
               
               [peakValues, turn_times] = findpeaks(this_absolute_data, 128, 'MinPeakDistance',min_peak_distance_threshold, 'MinPeakHeight', min_peak_height_threshold);
               

                % go through and assign each stant and end to each turn
            
                for i_turns = 1 : length(turn_times)
                    turn_start_times(i_turns) = turn_times(i_turns) - start_turn_modifier;
                    turn_end_times(i_turns) = turn_times(i_turns) + end_turn_modifier;
                    
                    
                    %%%% Plots vert lines at stat and end of turn
           
                    gcf;        
                    xline(turn_start_times(i_turns), '--', 'color', '#0B0');
                    xline(turn_end_times(i_turns), '--ro');
                  
               end
            end
           
            
            % find heel-strikes
            
            
            
            
            
            
            
            
        end
    

        
           

        %% save
        % struct for saving
        variables_to_save = struct;

        event_data = ...
          { ...
            turn_times; ...
            turn_start_times; ...
            turn_end_times; ...
          };
        event_labels = ...
          { ...
            'turn_times'; ...
            'turn_start_time'; ...
            'turn_end_time'; ...
          };

        % add new variables to be saved
        variables_to_save.event_data = event_data;
        variables_to_save.event_labels = event_labels;

        step_events_file_name = ['analysis' filesep makeFileName(collection_date, subject_id, condition, i_trial, 'events.mat')];
        saveDataToFile(step_events_file_name, variables_to_save);

        disp(['Finding Turn Events: condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' step_events_file_name]);


    end
end
   


   














