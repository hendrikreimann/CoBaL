function preprocessRawData()

    visualize                   = 0;
    process_emg                 = 1;
    process_forceplate          = 0;
    process_marker              = 0;
    map_emg_sensors_manually    = 1;
    transform_to_belt_space     = 1;
    
    conditions_to_transform_to_belt_space = {'baselineTM', 'feedback', 'postTM'};
    
    load('subjectInfo.mat');

    %% emg
    if process_emg
        data_dir = dir(['raw' filesep '*_emgTrajectoriesRaw.mat']);
        clear file_name_list;
        [file_name_list{1:length(data_dir)}] = deal(data_dir.name);

        number_of_files = length(file_name_list);
        for i_trial = 1 : number_of_files
            % load data
            raw_emg_file_name = file_name_list{i_trial};
            [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(raw_emg_file_name);
            load(['raw' filesep raw_emg_file_name]);

            % define filters
            % initial low pass filter
            filter_order_low = 4;
            cutoff_frequency_low = 500; % in Hz
            [b_low, a_low] = butter(filter_order_low, cutoff_frequency_low/(sampling_rate_emg/2), 'low');

            % high pass filter at 20 hz to get rid of DC offset
            filter_order_high = 4;
            cutoff_frequency_high = 20; % in Hz
            [b_high, a_high] = butter(filter_order_high, cutoff_frequency_high/(sampling_rate_emg/2), 'high');

            % low pass filter below 10 Hz -- aggressive smoothing after rectification
            filter_order_final = 4;
            cutoff_frequency_final = 10; % in Hz
            [b_final, a_final] = butter(filter_order_final, cutoff_frequency_high/(sampling_rate_emg/2), 'low');

            % filter, then rectify
            emg_trajectories_filtered_lowpass = filtfilt(b_low, a_low, emg_trajectories_raw);
            emg_trajectories_filtered_highpass = filtfilt(b_high, a_high, emg_trajectories_filtered_lowpass);
            emg_trajectories_rectified = abs(emg_trajectories_filtered_highpass);
            emg_trajectories = filtfilt(b_final, a_final, emg_trajectories_rectified);
            
            % rectify, then filter
        %     emg_trajectories_rectified = abs(emg_trajectories_raw);
        %     emg_trajectories_filtered_lowpass = filtfilt(b_low, a_low, emg_trajectories_rectified);
        %     emg_trajectories_filtered_highpass = filtfilt(b_high, a_high, emg_trajectories_filtered_lowpass);
        %     emg_trajectories = filtfilt(b_final, a_final, emg_trajectories_filtered_highpass);

            % rectify, then smooth
        %     emg_trajectories_rectified = abs(emg_trajectories_raw);
        %     rms_smooth_window_length_indices = rms_smooth_window_length * sampling_rate_emg;
        %     time_smoothed = rms_gbiomech(time_emg, rms_smooth_window_length_indices, 0, 0);
        %     rms_smoothed = rms_gbiomech(emg_trajectories_rectified(:, 1), rms_smooth_window_length_indices, 0, 0);
            
            % Determine max EMG reading for each muscle
            emg_max_trajectories = max(emg_trajectories, [],1);
            
            if map_emg_sensors_manually
                lglutmed_sensor_index = 1;
                ltibiant_sensor_index = 2;
                lgastroc_sensor_index = 3;
                lperolng_sensor_index = 4;
                rglutmed_sensor_index = 5;
                rtibiant_sensor_index = 6;
                rgastroc_sensor_index = 7;
                rperolng_sensor_index = 8;

                emg_headers{lglutmed_sensor_index} = 'LGLUTMED';
                emg_headers{ltibiant_sensor_index} = 'LTIBIANT';
                emg_headers{lgastroc_sensor_index} = 'LGASTROC';
                emg_headers{lperolng_sensor_index} = 'LPEROLNG';
                emg_headers{rglutmed_sensor_index} = 'RGLUTMED';
                emg_headers{rtibiant_sensor_index} = 'RTIBIANT';
                emg_headers{rgastroc_sensor_index} = 'RGASTROC';
                emg_headers{rperolng_sensor_index} = 'RPEROLNG';
            end

            % save
            save_file_name = ['processed' filesep makeFileName(date, subject_id, trial_type, trial_number, 'emgTrajectories.mat')];
            save ...
              ( ...
                save_file_name, ...
                'emg_trajectories', ...
                'emg_max_trajectories', ...
                'time_emg', ...
                'sampling_rate_emg', ...
                'emg_headers' ...
              );
            disp(['filtered and saved as ' save_file_name])

            % visualize
            if visualize
                i_channel = 1;
                figure; axes; hold on
        %         plot(time_emg, emg_trajectories_raw(:, i_channel), 'DisplayName', 'raw');
                plot(time_emg, emg_trajectories_rectified(:, i_channel), 'DisplayName', 'rectified');
        %         plot(time_smoothed, rms_smoothed, 'DisplayName', 'rms smoothed', 'linewidth', 2);
        %         plot(time_emg, emg_trajectories_filtered_lowpass(:, i_channel), 'DisplayName', 'lowpass');
        %         plot(time_emg, emg_trajectories_filtered_highpass(:, i_channel), 'DisplayName', 'highpass');
                plot(time_emg, emg_trajectories(:, i_channel), 'linewidth', 2, 'DisplayName', 'final');
                plot([time_emg(1) time_emg(end)], [emg_max_trial(i_channel) emg_max_trial(i_channel)], 'r');
                legend('toggle');
                keyboard
                distFig
            end    




        end
    end
    
    %% forceplate data
    if process_forceplate
        data_dir = dir(['raw' filesep '*_forceplateTrajectoriesRaw.mat']);
        clear file_name_list;
        [file_name_list{1:length(data_dir)}] = deal(data_dir.name);


        number_of_files = length(file_name_list);
        for i_trial = 1 : number_of_files

            % load data
            raw_forceplate_file_name = file_name_list{i_trial};
            [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(raw_forceplate_file_name);
            load(['raw' filesep raw_forceplate_file_name]);
            
            % filter
            filter_order_low = 4;
            cutoff_frequency_low = 50; % in Hz
            [b_lowpass, a_lowpass] = butter(filter_order_low, cutoff_frequency_low/(sampling_rate_forceplate/2), 'low');
            forceplate_trajectories_filtered = filtfilt(b_lowpass, a_lowpass, forceplate_trajectories_raw);


            % extract
            fxl_trajectory = forceplate_trajectories_filtered(:, 1);
            fyl_trajectory = forceplate_trajectories_filtered(:, 2);
            fzl_trajectory = forceplate_trajectories_filtered(:, 3);
            mxl_trajectory = forceplate_trajectories_filtered(:, 4);
            myl_trajectory = forceplate_trajectories_filtered(:, 5);
            mzl_trajectory = forceplate_trajectories_filtered(:, 6);
            fxr_trajectory = forceplate_trajectories_filtered(:, 7);
            fyr_trajectory = forceplate_trajectories_filtered(:, 8);
            fzr_trajectory = forceplate_trajectories_filtered(:, 9);
            mxr_trajectory = forceplate_trajectories_filtered(:, 10);
            myr_trajectory = forceplate_trajectories_filtered(:, 11);
            mzr_trajectory = forceplate_trajectories_filtered(:, 12);

            % calculate CoP
            copxl_trajectory = - myl_trajectory ./ fzl_trajectory;
            copyl_trajectory = mxl_trajectory ./ fzl_trajectory;
            copxr_trajectory = - myr_trajectory ./ fzr_trajectory;
            copyr_trajectory = mxr_trajectory ./ fzr_trajectory;

            % zero CoP for no contact times
            fz_threshold = 60;
            copxl_trajectory(fzl_trajectory < fz_threshold) = 0;
            copyl_trajectory(fzl_trajectory < fz_threshold) = 0;
            copxr_trajectory(fzr_trajectory < fz_threshold) = 0;
            copyr_trajectory(fzr_trajectory < fz_threshold) = 0;

            % transform forceplate data to CoBaL world frame A_cw
            left_forceplate_wrench_Acl = [fxl_trajectory fyl_trajectory fzl_trajectory mxl_trajectory myl_trajectory mzl_trajectory];
            left_forceplate_cop_Acl = [copxl_trajectory copyl_trajectory zeros(size(copxl_trajectory))];
            right_forceplate_wrench_Acr = [fxr_trajectory fyr_trajectory fzr_trajectory mxr_trajectory myr_trajectory mzr_trajectory];
            right_forceplate_cop_Acr = [copxr_trajectory copyr_trajectory zeros(size(copxr_trajectory))];

            % define forceplate rotation and translation
            Acr_to_Acw_rotation = [-1 0 0; 0 1 0; 0 0 -1];
            Acr_to_Acw_translation = [0.5588; 0; 0];
            Acr_to_Acw_trafo = [Acr_to_Acw_rotation Acr_to_Acw_translation; 0 0 0 1];
            Acl_to_Acw_rotation = [-1 0 0; 0 1 0; 0 0 -1];
            Acl_to_Acw_translation = [-0.5588; 0; 0];
            Acl_to_Acw_trafo = [Acl_to_Acw_rotation Acl_to_Acw_translation; 0 0 0 1];
            Acr_to_Acw_adjoint = rigidToAdjointTransformation(Acr_to_Acw_trafo);
            Acl_to_Acw_adjoint = rigidToAdjointTransformation(Acl_to_Acw_trafo);

            % transform
            left_forceplate_wrench_Acw = (Acl_to_Acw_adjoint' * left_forceplate_wrench_Acl')';
            left_forceplate_cop_Acw = (eye(2, 4) * Acl_to_Acw_trafo * [left_forceplate_cop_Acl ones(size(left_forceplate_cop_Acl, 1), 1)]')';
            right_forceplate_wrench_Acw = (Acr_to_Acw_adjoint' * right_forceplate_wrench_Acr')';
            right_forceplate_cop_Acw = (eye(2, 4) * Acr_to_Acw_trafo * [right_forceplate_cop_Acr ones(size(right_forceplate_cop_Acr, 1), 1)]')';

            % calculate wrenches and CoP for complete plate
            total_forceplate_wrench_Acw = left_forceplate_wrench_Acw + right_forceplate_wrench_Acw;
            copx_trajectory = - total_forceplate_wrench_Acw(:, 5) ./ total_forceplate_wrench_Acw(:, 3);
            copy_trajectory = total_forceplate_wrench_Acw(:, 4) ./ total_forceplate_wrench_Acw(:, 3);
            total_forceplate_cop_Acw = [copx_trajectory copy_trajectory];

            % re-zero CoP for low loads
            left_forceplate_low_load_indicator = copxl_trajectory == 0;
            left_forceplate_cop_Acw(left_forceplate_low_load_indicator, :) = 0;
            right_forceplate_low_load_indicator = copxr_trajectory == 0;
            right_forceplate_cop_Acw(right_forceplate_low_load_indicator, :) = 0;
            total_forceplate_low_load_indicator = copx_trajectory == 0;
            total_forceplate_cop_Acw(total_forceplate_low_load_indicator, :) = 0;


        %     % visualize
        %     figure; axes; hold on
        %     plot(time_forceplate, total_forceplate_cop_Acw(:, 1));
        %     plot(time_forceplate, left_forceplate_cop_Acw(:, 1));
        %     plot(time_forceplate, right_forceplate_cop_Acw(:, 1));
        %     legend('total', 'left', 'right')
        % 
        %     figure; axes; hold on
        %     plot(time_forceplate, total_forceplate_cop_Acw(:, 2));
        %     plot(time_forceplate, left_forceplate_cop_Acw(:, 2));
        %     plot(time_forceplate, right_forceplate_cop_Acw(:, 2));
        %     legend('total', 'left', 'right')    

            % re-zero CoP for low loads
            left_forceplate_low_load_indicator = copxl_trajectory == 0;
            left_forceplate_cop_Acw(left_forceplate_low_load_indicator, :) = 0;
            right_forceplate_low_load_indicator = copxr_trajectory == 0;
            right_forceplate_cop_Acw(right_forceplate_low_load_indicator, :) = 0;

            % save
            save_file_name = ['processed' filesep makeFileName(date, subject_id, trial_type, trial_number, 'forceplateTrajectories.mat')];
            save ...
              ( ...
                save_file_name, ...
                'left_forceplate_wrench_Acw', ...
                'left_forceplate_cop_Acw', ...
                'right_forceplate_wrench_Acw', ...
                'right_forceplate_cop_Acw', ...
                'total_forceplate_wrench_Acw', ...
                'total_forceplate_cop_Acw', ...
                'fxl_trajectory', ...
                'fyl_trajectory', ...
                'fzl_trajectory', ...
                'mxl_trajectory', ...
                'myl_trajectory', ...
                'mzl_trajectory', ...
                'copxl_trajectory', ...
                'copyl_trajectory', ...
                'fxr_trajectory', ...
                'fyr_trajectory', ...
                'fzr_trajectory', ...
                'mxr_trajectory', ...
                'myr_trajectory', ...
                'mzr_trajectory', ...
                'copxr_trajectory', ...
                'copyr_trajectory', ...
                'time_forceplate', ...
                'sampling_rate_forceplate' ...
              );
            disp(['processed ' raw_forceplate_file_name ' and saved as ' save_file_name])        

        end
    end

    %% marker data
    if process_marker
        data_dir = dir(['raw' filesep '*_markerTrajectoriesRaw.mat']);
        clear file_name_list;
        [file_name_list{1:length(data_dir)}] = deal(data_dir.name);
        number_of_files = length(file_name_list);
        for i_trial = 1 : number_of_files

            % load data
            raw_marker_file_name = file_name_list{i_trial};
            [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(raw_marker_file_name);
            load(['raw' filesep raw_marker_file_name]);
            
            % save
            save_file_name = ['processed' filesep makeFileName(date, subject_id, trial_type, trial_number, 'markerTrajectories.mat')];
            save ...
              ( ...
                save_file_name, ...
                'marker_trajectories', ...
                'time_mocap', ...
                'sampling_rate_mocap', ...
                'marker_headers' ...
              );
            disp(['processed ' raw_marker_file_name ' and saved as ' save_file_name])        
        end
    end

    
    %% if transform to belt space
    if transform_to_belt_space
        [condition_list, trial_number_list] = parseTrialArguments();
    
        for i_condition = 1 : length(condition_list)
            trials_to_process = trial_number_list{i_condition};
            for i_trial = trials_to_process
                % load data
                condition = condition_list{i_condition};
                if any(strcmp(conditions_to_transform_to_belt_space, condition))
                    % extract data for new structure
                    load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'plcData')]);
                    time_belts = time_plcData - time_plcData(1);

                    % calculate shift
                    belt_speed_trajectory = mean([belt_speed_left_trajectory belt_speed_right_trajectory], 2);
                    delta_t = diff(time_belts);
                    belt_position_trajectory_plcData = zeros(size(belt_speed_trajectory));
                    for i_time = 2 : length(belt_speed_trajectory)
                        belt_position_trajectory_plcData(i_time) = belt_position_trajectory_plcData(i_time-1) + delta_t(i_time-1) * belt_speed_trajectory(i_time-1);
                    end

                    
                    % apply shift to marker trajectories
                    file_name_raw = ['raw' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectoriesRaw.mat')];
                    load(file_name_raw);
%                     marker_trajectories_original = marker_trajectories;
                    belt_position_trajectory_mocap = spline(time_belts, belt_position_trajectory_plcData, time_mocap)';
                    for i_marker = 1 : size(marker_headers, 2)
                        marker_trajectories(:, (i_marker-1)*3+2) = marker_trajectories(:, (i_marker-1)*3+2) + belt_position_trajectory_mocap;
                    end
                    
                    
%                     % do the same thing, but move to mocap time first
%                     belt_speed_trajectory = mean([belt_speed_left_trajectory belt_speed_right_trajectory], 2);
%                     belt_speed_trajectory_mocap = spline(time_belts, belt_speed_trajectory, time_mocap)';
%                     delta_t = 1/sampling_rate_mocap;
%                     belt_position_trajectory_mocap = zeros(size(belt_speed_trajectory_mocap));
%                     for i_time = 2 : length(belt_speed_trajectory_mocap)
%                         belt_position_trajectory_mocap(i_time) = belt_position_trajectory_mocap(i_time-1) + delta_t * belt_speed_trajectory_mocap(i_time-1);
%                     end
%                     marker_trajectories_alt = marker_trajectories_original;
%                     for i_marker = 1 : size(marker_headers, 2)
%                         marker_trajectories_alt(:, (i_marker-1)*3+2) = marker_trajectories_alt(:, (i_marker-1)*3+2) + belt_position_trajectory_mocap;
%                     end
                    
%                     figure; axes; hold on
                    plot(marker_trajectories(:, 101), 'linewidth', 2)
%                     plot(marker_trajectories_alt(:, 101), 'linewidth', 1)

%                     % things do not make sense, play around
%                     left_heel_y_marker_trajectory_raw = marker_trajectories_original(:, 101);
%                     left_heel_y_marker_speed_raw = deriveByTime(left_heel_y_marker_trajectory_raw, delta_t);
%                     
%                     figure; axes; hold on
% %                     plot(left_heel_y_marker_trajectory_raw, 'linewidth', 2)
%                     plot(time_mocap, left_heel_y_marker_speed_raw, 'linewidth', 2)
%                     plot(time_belts, -belt_speed_trajectory, 'linewidth', 2)



                    file_name_shifted = ['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories.mat')];
                    save ...
                      ( ...
                        file_name_shifted, ...
                        'marker_trajectories', ...
                        'time_mocap', ...
                        'sampling_rate_mocap', ...
                        'marker_headers' ...
                      );
                    disp(['Transformed marker data in ' file_name_raw ' to belt space and saved to ' file_name_shifted])                    
                end
                

    



    %             % apply shift to forceplate trajectories
    %             load(makeFileName(date, subject_id, condition, i_trial, 'forceplateTrajectories'));
    %             belt_position_trajectory_forceplate = spline(time_belts, belt_position_trajectory_plcData, time_forceplate)';
    %             
    %             for i_time = 1 : length(time_forceplate)
    %                 % define forceplate rotation and translation
    %                 Acw_to_Acb_rotation = [1 0 0; 0 1 0; 0 0 1];
    %                 Acw_to_Acb_translation = [0.5588; 0; 0];
    %                 Acw_to_Acb_trafo = [Acw_to_Acb_rotation Acw_to_Acb_translation; 0 0 0 1];
    %                 Acw_to_Acb_adjoint = rigidToAdjointTransformation(Acw_to_Acb_trafo);
    % 
    %                 % transform
    %                 left_forceplate_wrench_Acb = (Acw_to_Acb_adjoint' * left_forceplate_wrench_Acw')';
    %                 right_forceplate_wrench_Acb = (Acw_to_Acb_adjoint' * right_forceplate_wrench_Acw')';
    % 
    %             end
    % 
    %             % calculate wrenches and CoP for complete plate
    %             total_forceplate_wrench_Acb = left_forceplate_wrench_Acb + right_forceplate_wrench_Acb;
    %             copx_trajectory = - total_forceplate_wrench_Acb(:, 5) ./ total_forceplate_wrench_Acb(:, 3);
    %             copy_trajectory = total_forceplate_wrench_Acb(:, 4) ./ total_forceplate_wrench_Acb(:, 3);
    %             total_forceplate_cop_Acb = [copx_trajectory copy_trajectory];
    %                 
    %             % re-zero CoP for low loads
    %             left_forceplate_low_load_indicator = copxl_trajectory == 0;
    %             left_forceplate_cop_Acw(left_forceplate_low_load_indicator, :) = 0;
    %             right_forceplate_low_load_indicator = copxr_trajectory == 0;
    %             right_forceplate_cop_Acw(right_forceplate_low_load_indicator, :) = 0;
    %             total_forceplate_low_load_indicator = copx_trajectory == 0;
    %             total_forceplate_cop_Acw(total_forceplate_low_load_indicator, :) = 0;            









            end
        end

    end
    
    


    
    
    
    
    
    
    
    
end