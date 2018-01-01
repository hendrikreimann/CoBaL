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

% this function applies several basic processing steps to experimental data, e.g. filtering

% input: 
% Experimental data files generated by importAscii.m, in the subfolder "raw"
%
% output: 
% multiple files with processed data for each trial, in the subfolder "processed"


function preprocessRawData(varargin)
    load('subjectInfo.mat');
    
    % parse arguments
    [condition_list, trial_number_list, calibration_trials, emg_trials] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;

    % add calibration and EMG trials if no specific condition was specified
    if ~isempty(calibration_trials)
        condition_list = [condition_list; 'calibration'];
        trial_number_list = [trial_number_list; calibration_trials];
    end
    if ~isempty(emg_trials)
        condition_list = [condition_list; 'emg'];
        trial_number_list = [trial_number_list; emg_trials];
    end
    
    % load settings
    study_settings_file = '';
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    subject_settings = SettingsCustodian('subjectSettings.txt');
    
    %% emg
    data_dir = dir(['raw' filesep '*_emgTrajectoriesRaw.mat']);
    clear file_name_list;
    [file_name_list{1:length(data_dir)}] = deal(data_dir.name);
    number_of_files = length(file_name_list);
    for i_trial = 1 : number_of_files
        raw_emg_file_name = file_name_list{i_trial};
        [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(raw_emg_file_name);

        % does the caller want to process this file?
        if any(strcmp(trial_type, condition_list))
            % condition is set to be processed, now check trial number
            trial_number_list_this_condition = trial_number_list{strcmp(trial_type, condition_list)};
            if ismember(trial_number, trial_number_list_this_condition)
                % process file
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
                cutoff_frequency_final = 6; % in Hz
                [b_final, a_final] = butter(filter_order_final, cutoff_frequency_final/(sampling_rate_emg/2), 'low');

                % filter, then rectify
%                 emg_trajectories_filtered_lowpass = filtfilt(b_low, a_low, emg_trajectories_raw);
%                 emg_trajectories_filtered_highpass = filtfilt(b_high, a_high, emg_trajectories_filtered_lowpass);
%                 emg_trajectories_rectified = abs(emg_trajectories_filtered_highpass);
%                 emg_trajectories = filtfilt(b_final, a_final, emg_trajectories_rectified);
                
                % rectify, then filter
                emg_trajectories_rectified = abs(emg_trajectories_raw);
                emg_trajectories = filtfilt(b_final, a_final, emg_trajectories_rectified);
                
                
%                 emg_labels_from_source = emg_labels;
%                 emg_labels = cell(size(emg_labels_from_source));
%                 for i_label = 1 : size(emg_labels_from_source, 2)
%                     % find entry in emg_sensor_map
%                     match_column = find(strcmp(emg_labels_from_source{i_label}, emg_sensor_map(1, :)));
% 
%                     if ~isempty(match_column)
%                         emg_labels(i_label) = emg_sensor_map(2, match_column);
%                     end
%                 end
                emg_labels = subject_settings.get('emg_labels');
                
                % apply RMS smoothing
%                 window_size_time = 0.2; % time in seconds
%                 window_radius_steps = ceil(sampling_rate_emg * window_size_time / 2);
%                 emg_rms_smoothed = zeros(size(emg_trajectories_raw)) * NaN;
%                 emg_rms_rectified = abs(emg_trajectories_raw);
%                 for i_time = 1 + window_radius_steps : size(emg_trajectories_raw, 1) - window_radius_steps
%                     window_steps = i_time - window_radius_steps : i_time + window_radius_steps;
%                     window_data = emg_rms_rectified(window_steps, :);
%                     
%                     rms = mean(window_data.^2).^(0.5);
%                     emg_rms_smoothed(i_time, :) = mean(window_data.^2).^(0.5);
%                     
%                 end
                
                % save
                save_folder = 'processed';
                save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'emgTrajectories.mat');
                save ...
                  ( ...
                    [save_folder filesep save_file_name], ...
                    'emg_trajectories', ...
                    'time_emg', ...
                    'sampling_rate_emg', ...
                    'emg_labels' ...
                  );
                addAvailableData('emg_trajectories', 'time_emg', 'sampling_rate_emg', 'emg_labels', save_folder, save_file_name);
                disp(['filtered and saved as ' save_file_name])

                % visualize
                if visualize
                    i_channel = 3;
                    figure; axes; hold on; title(['EMG, condition ' trial_type ', trial ' num2str(trial_number)])
                    plot(time_emg, emg_trajectories_raw(:, i_channel), 'DisplayName', 'raw');
                    plot(time_emg, emg_trajectories_rectified(:, i_channel), 'DisplayName', 'rectified');
%                     plot(time_emg, emg_rms_rectified(:, i_channel), 'DisplayName', 'rms rectified', 'linewidth', 2);
%                     plot(time_emg, emg_rms_smoothed(:, i_channel), 'DisplayName', 'rms smoothed', 'linewidth', 2);
            %         plot(time_emg, emg_trajectories_filtered_lowpass(:, i_channel), 'DisplayName', 'lowpass');
            %         plot(time_emg, emg_trajectories_filtered_highpass(:, i_channel), 'DisplayName', 'highpass');
                    plot(time_emg, emg_trajectories(:, i_channel), 'linewidth', 2, 'DisplayName', 'final');
                    legend('toggle');
                end                        


            end
        end
    end
 
    %% forceplate data
    data_dir = dir(['raw' filesep '*_forceplateTrajectoriesRaw.mat']);
    clear file_name_list;
    [file_name_list{1:length(data_dir)}] = deal(data_dir.name);
    number_of_files = length(file_name_list);
    for i_trial = 1 : number_of_files
        raw_forceplate_file_name = file_name_list{i_trial};
        [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(raw_forceplate_file_name);
        
        % does the caller want to process this file?
        if any(strcmp(trial_type, condition_list))
            % condition is set to be processed, now check trial number
            trial_number_list_this_condition = trial_number_list{strcmp(trial_type, condition_list)};
            if ismember(trial_number, trial_number_list_this_condition)
                % load raw data
                load(['raw' filesep raw_forceplate_file_name]);

                % define filter
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
                
                % apply offset for cases where forceplate wasn't set to zero
                if subject_settings.get('apply_forceplate_offset')
                    fxl_trajectory = fxl_trajectory - subject_settings.get('offset_fxl');
                    fyl_trajectory = fyl_trajectory - subject_settings.get('offset_fyl');
                    fzl_trajectory = fzl_trajectory - subject_settings.get('offset_fzl');
                    mxl_trajectory = mxl_trajectory - subject_settings.get('offset_mxl');
                    myl_trajectory = myl_trajectory - subject_settings.get('offset_myl');
                    mzl_trajectory = mzl_trajectory - subject_settings.get('offset_mzl');
                    fxr_trajectory = fxr_trajectory - subject_settings.get('offset_fxr');
                    fyr_trajectory = fyr_trajectory - subject_settings.get('offset_fyr');
                    fzr_trajectory = fzr_trajectory - subject_settings.get('offset_fzr');
                    mxr_trajectory = mxr_trajectory - subject_settings.get('offset_mxr');
                    myr_trajectory = myr_trajectory - subject_settings.get('offset_myr');
                    mzr_trajectory = mzr_trajectory - subject_settings.get('offset_mzr');
                end
                

%                 % calculate CoP
%                 copxl_trajectory = - myl_trajectory ./ fzl_trajectory;
%                 copyl_trajectory = mxl_trajectory ./ fzl_trajectory;
%                 copxr_trajectory = - myr_trajectory ./ fzr_trajectory;
%                 copyr_trajectory = mxr_trajectory ./ fzr_trajectory;
% 
%                 % zero CoP for no contact times
%                 fz_threshold = 60;
                fz_threshold = subject_settings.get('forceplate_load_threshold');
%                 copxl_trajectory(fzl_trajectory < fz_threshold) = 0;
%                 copyl_trajectory(fzl_trajectory < fz_threshold) = 0;
%                 copxr_trajectory(fzr_trajectory < fz_threshold) = 0;
%                 copyr_trajectory(fzr_trajectory < fz_threshold) = 0;

                if strcmp(data_source, 'nexus')
                    % transform forceplate data to CoBaL world frame A_cw
                    left_forceplate_wrench_Acl = [fxl_trajectory fyl_trajectory fzl_trajectory mxl_trajectory myl_trajectory mzl_trajectory];
%                     left_forceplate_cop_Acl = [copxl_trajectory copyl_trajectory zeros(size(copxl_trajectory))];
                    right_forceplate_wrench_Acr = [fxr_trajectory fyr_trajectory fzr_trajectory mxr_trajectory myr_trajectory mzr_trajectory];
%                     right_forceplate_cop_Acr = [copxr_trajectory copyr_trajectory zeros(size(copxr_trajectory))];

                    % define forceplate rotation and translation
                    Acr_to_world_rotation = [-1 0 0; 0 1 0; 0 0 -1];
                    Acr_to_world_translation = [0.5588; 0; 0]; % origin of Acw in Acr frame
                    Acr_to_world_trafo = [Acr_to_world_rotation Acr_to_world_translation; 0 0 0 1];
                    Acl_to_world_rotation = [-1 0 0; 0 1 0; 0 0 -1];
                    Acl_to_world_translation = [-0.5588; 0; 0]; % origin of Acw in Acl frame
                    Acl_to_world_trafo = [Acl_to_world_rotation Acl_to_world_translation; 0 0 0 1];
                    Acr_to_world_adjoint = rigidToAdjointTransformation(Acr_to_world_trafo);
                    Acl_to_world_adjoint = rigidToAdjointTransformation(Acl_to_world_trafo);

                    % transform
                    left_forceplate_wrench_world = (Acl_to_world_adjoint' * left_forceplate_wrench_Acl')';
%                     left_forceplate_cop_world = (eye(2, 4) * Acl_to_world_trafo * [left_forceplate_cop_Acl ones(size(left_forceplate_cop_Acl, 1), 1)]')';
                    right_forceplate_wrench_world = (Acr_to_world_adjoint' * right_forceplate_wrench_Acr')';
%                     right_forceplate_cop_world = (eye(2, 4) * Acr_to_world_trafo * [right_forceplate_cop_Acr ones(size(right_forceplate_cop_Acr, 1), 1)]')';
                elseif strcmp(data_source, 'qtm')
                    % transform forceplate data to CoBaL world frame A_cw
                    left_forceplate_wrench_Acl = [fxl_trajectory fyl_trajectory fzl_trajectory mxl_trajectory myl_trajectory mzl_trajectory];
                    right_forceplate_wrench_Acr = [fxr_trajectory fyr_trajectory fzr_trajectory mxr_trajectory myr_trajectory mzr_trajectory];
                    % define forceplate rotation and translation
                    Acr_to_world_rotation = [-1 0 0; 0 1 0; 0 0 -1];
                    Acr_to_world_translation = [0.5588; 0; 0]; % origin of Acw in Acr frame
                    Acr_to_world_trafo = [Acr_to_world_rotation Acr_to_world_translation; 0 0 0 1];
                    Acl_to_world_rotation = [-1 0 0; 0 1 0; 0 0 -1];
                    Acl_to_world_translation = [-0.5588; 0; 0]; % origin of Acw in Acl frame
                    Acl_to_world_trafo = [Acl_to_world_rotation Acl_to_world_translation; 0 0 0 1];
                    Acr_to_world_adjoint = rigidToAdjointTransformation(Acr_to_world_trafo);
                    Acl_to_world_adjoint = rigidToAdjointTransformation(Acl_to_world_trafo);
                    
                    % transform
                    left_forceplate_wrench_world = (Acl_to_world_adjoint' * left_forceplate_wrench_Acl')';
                    right_forceplate_wrench_world = (Acr_to_world_adjoint' * right_forceplate_wrench_Acr')';
                elseif strcmp(data_source, 'neurocom')
                    % transform forceplate data to VEPO world frame A_vw
                    left_forceplate_wrench_Anl = [fxl_trajectory fyl_trajectory fzl_trajectory mxl_trajectory myl_trajectory mzl_trajectory];
%                     left_forceplate_cop_Anl = [copxl_trajectory copyl_trajectory zeros(size(copxl_trajectory))];
                    right_forceplate_wrench_Anr = [fxr_trajectory fyr_trajectory fzr_trajectory mxr_trajectory myr_trajectory mzr_trajectory];
%                     right_forceplate_cop_Anr = [copxr_trajectory copyr_trajectory zeros(size(copxr_trajectory))];
                    
                    % define forceplate rotation and translation
                    inchToMeter = 0.0254;
                    
                    Anl_to_world_rotation = [-1 0 0; 0 1 0; 0 0 -1];
                    Anl_to_world_translation = [5; -10; 0] * inchToMeter;
%                     Anl_to_world_translation = [5; 10; 0] * 0;
                    Anl_to_world_trafo = [Anl_to_world_rotation Anl_to_world_translation; 0 0 0 1];
                    Anl_to_world_adjoint = rigidToAdjointTransformation(Anl_to_world_trafo);
                    
                    Anr_to_world_rotation = [-1 0 0; 0 1 0; 0 0 -1];
                    Anr_to_world_translation = [-5; -10; 0] * inchToMeter;
                    Anr_to_world_trafo = [Anr_to_world_rotation Anr_to_world_translation; 0 0 0 1];
                    Anr_to_world_adjoint = rigidToAdjointTransformation(Anr_to_world_trafo);
                    
                    % transform
                    left_forceplate_wrench_world = (Anl_to_world_adjoint' * left_forceplate_wrench_Anl')';
%                     left_forceplate_cop_world = (eye(2, 4) * Anl_to_world_trafo * [left_forceplate_cop_Anl ones(size(left_forceplate_cop_Anl, 1), 1)]')';
                    right_forceplate_wrench_world = (Anr_to_world_adjoint' * right_forceplate_wrench_Anr')';
%                     right_forceplate_cop_world = (eye(2, 4) * Anr_to_world_trafo * [right_forceplate_cop_Anr ones(size(right_forceplate_cop_Anr, 1), 1)]')';
                    
                else
                    error(['data source "' data_source '" not recognized'])
                end

                % calculate wrench for complete plate
                total_forceplate_wrench_world = left_forceplate_wrench_world + right_forceplate_wrench_world;
                
                % calculate CoP
                copxl_raw = - myl_trajectory ./ fzl_trajectory;
                copyl_raw = mxl_trajectory ./ fzl_trajectory;
                copxr_raw = - myr_trajectory ./ fzr_trajectory;
                copyr_raw = mxr_trajectory ./ fzr_trajectory;
                
                
                copxl_world = - left_forceplate_wrench_world(:, 5) ./ left_forceplate_wrench_world(:, 3);
                copyl_world = left_forceplate_wrench_world(:, 4) ./ left_forceplate_wrench_world(:, 3);
                copxr_world = - right_forceplate_wrench_world(:, 5) ./ right_forceplate_wrench_world(:, 3);
                copyr_world = right_forceplate_wrench_world(:, 4) ./ right_forceplate_wrench_world(:, 3);
                copx_world = - total_forceplate_wrench_world(:, 5) ./ total_forceplate_wrench_world(:, 3);
                copy_world = total_forceplate_wrench_world(:, 4) ./ total_forceplate_wrench_world(:, 3);
%                 copxr_world = - left_forceplate_wrench_world(:, 5) ./ left_forceplate_wrench_world(:, 3);
%                 copyr_world = left_forceplate_wrench_world(:, 4) ./ left_forceplate_wrench_world(:, 3);
%                 copx_world = - left_forceplate_wrench_world(:, 5) ./ left_forceplate_wrench_world(:, 3);
%                 copy_world = left_forceplate_wrench_world(:, 4) ./ left_forceplate_wrench_world(:, 3);
                
                % re-zero CoP for low loads
                left_forceplate_low_load_indicator = (fzl_trajectory < fz_threshold);
                copxl_raw(left_forceplate_low_load_indicator, :) = 0;
                copyl_raw(left_forceplate_low_load_indicator, :) = 0;
                copxl_world(left_forceplate_low_load_indicator, :) = 0;
                copyl_world(left_forceplate_low_load_indicator, :) = 0;
                
                right_forceplate_low_load_indicator = (fzr_trajectory < fz_threshold);
                copxr_raw(right_forceplate_low_load_indicator, :) = 0;
                copyr_raw(right_forceplate_low_load_indicator, :) = 0;
                copxr_world(right_forceplate_low_load_indicator, :) = 0;
                copyr_world(right_forceplate_low_load_indicator, :) = 0;
                
                total_forceplate_low_load_indicator = (left_forceplate_low_load_indicator & right_forceplate_low_load_indicator);
                copx_world(total_forceplate_low_load_indicator, :) = 0;
                copy_world(total_forceplate_low_load_indicator, :) = 0;
                
                left_forceplate_cop_world = [copxl_world copyl_world];
                right_forceplate_cop_world = [copxr_world copyr_world];
                total_forceplate_cop_world = [copx_world copy_world];
                
                % define wrench and CoP trajectories for feet instead of forceplate sides
                if strcmp(data_source, 'nexus')
                    left_foot_cop_world = left_forceplate_cop_world;
                    left_foot_wrench_world = left_forceplate_wrench_world;
                    right_foot_cop_world = right_forceplate_cop_world;
                    right_foot_wrench_world = right_forceplate_wrench_world;
                end
                if strcmp(data_source, 'qtm')
                    left_foot_cop_world = left_forceplate_cop_world;
                    left_foot_wrench_world = left_forceplate_wrench_world;
                    right_foot_cop_world = right_forceplate_cop_world;
                    right_foot_wrench_world = right_forceplate_wrench_world;
                end
                if strcmp(data_source, 'neurocom')
                    left_foot_cop_world = right_forceplate_cop_world;
                    left_foot_wrench_world = right_forceplate_wrench_world;
                    right_foot_cop_world = left_forceplate_cop_world;
                    right_foot_wrench_world = left_forceplate_wrench_world;
                end

                % save
                save_folder = 'processed';
                save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'forceplateTrajectories.mat');
                save ...
                  ( ...
                    [save_folder filesep save_file_name], ...
                    'left_forceplate_wrench_world', ...
                    'left_forceplate_cop_world', ...
                    'right_forceplate_wrench_world', ...
                    'right_forceplate_cop_world', ...
                    'left_foot_wrench_world', ...
                    'left_foot_cop_world', ...
                    'right_foot_wrench_world', ...
                    'right_foot_cop_world', ...
                    'total_forceplate_wrench_world', ...
                    'total_forceplate_cop_world', ...
                    'fxl_trajectory', ...
                    'fyl_trajectory', ...
                    'fzl_trajectory', ...
                    'mxl_trajectory', ...
                    'myl_trajectory', ...
                    'mzl_trajectory', ...
                    'fxr_trajectory', ...
                    'fyr_trajectory', ...
                    'fzr_trajectory', ...
                    'mxr_trajectory', ...
                    'myr_trajectory', ...
                    'mzr_trajectory', ...
                    'time_forceplate', ...
                    'sampling_rate_forceplate' ...
                  );
                addAvailableData('left_forceplate_wrench_world', 'time_forceplate', 'sampling_rate_forceplate', 'left_forceplate_wrench_world', save_folder, save_file_name);
                addAvailableData('left_forceplate_cop_world', 'time_forceplate', 'sampling_rate_forceplate', 'left_forceplate_cop_world', save_folder, save_file_name);
                addAvailableData('right_forceplate_wrench_world', 'time_forceplate', 'sampling_rate_forceplate', 'right_forceplate_wrench_world', save_folder, save_file_name);
                addAvailableData('right_forceplate_cop_world', 'time_forceplate', 'sampling_rate_forceplate', 'right_forceplate_cop_world', save_folder, save_file_name);

                addAvailableData('left_foot_wrench_world', 'time_forceplate', 'sampling_rate_forceplate', 'left_foot_wrench_world', save_folder, save_file_name);
                addAvailableData('left_foot_cop_world', 'time_forceplate', 'sampling_rate_forceplate', 'left_foot_cop_world', save_folder, save_file_name);
                addAvailableData('right_foot_wrench_world', 'time_forceplate', 'sampling_rate_forceplate', 'right_foot_wrench_world', save_folder, save_file_name);
                addAvailableData('right_foot_cop_world', 'time_forceplate', 'sampling_rate_forceplate', 'right_foot_cop_world', save_folder, save_file_name);

                addAvailableData('total_forceplate_wrench_world', 'time_forceplate', 'sampling_rate_forceplate', 'total_forceplate_wrench_world', save_folder, save_file_name);
                addAvailableData('total_forceplate_cop_world', 'time_forceplate', 'sampling_rate_forceplate', 'total_forceplate_cop_world', save_folder, save_file_name);
                
                addAvailableData('fxl_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'fxl_trajectory', save_folder, save_file_name);
                addAvailableData('fyl_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'fyl_trajectory', save_folder, save_file_name);
                addAvailableData('fzl_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'fzl_trajectory', save_folder, save_file_name);
                addAvailableData('mxl_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'mxl_trajectory', save_folder, save_file_name);
                addAvailableData('myl_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'myl_trajectory', save_folder, save_file_name);
                addAvailableData('mzl_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'mzl_trajectory', save_folder, save_file_name);
                
                addAvailableData('fxr_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'fxr_trajectory', save_folder, save_file_name);
                addAvailableData('fyr_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'fyr_trajectory', save_folder, save_file_name);
                addAvailableData('fzr_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'fzr_trajectory', save_folder, save_file_name);
                addAvailableData('mxr_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'mxr_trajectory', save_folder, save_file_name);
                addAvailableData('myr_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'myr_trajectory', save_folder, save_file_name);
                addAvailableData('mzr_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'mzr_trajectory', save_folder, save_file_name);
                disp(['processed ' raw_forceplate_file_name ' and saved as ' [save_folder filesep save_file_name]])        
            end
        end
    end

    %% marker data
    data_dir = dir(['raw' filesep '*_markerTrajectoriesRaw.mat']);
    clear file_name_list;
    [file_name_list{1:length(data_dir)}] = deal(data_dir.name);
    number_of_files = length(file_name_list);
    for i_trial = 1 : number_of_files
        % load data
        raw_marker_file_name = file_name_list{i_trial};
        [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(raw_marker_file_name);
        % does the caller want to process this file?
        if any(strcmp(trial_type, condition_list))
            % condition is set to be processed, now check trial number
            trial_number_list_this_condition = trial_number_list{strcmp(trial_type, condition_list)};
            if ismember(trial_number, trial_number_list_this_condition)
                load(['raw' filesep raw_marker_file_name]);
                
                if study_settings.get('filter_marker_data')
                    filter_order = 4;
                    cutoff_frequency = study_settings.get('marker_data_cutoff_frequency'); % in Hz
                    [b_marker, a_marker] = butter(filter_order, cutoff_frequency/(sampling_rate_mocap/2));
                    marker_trajectories = nanfiltfilt(b_marker, a_marker, marker_trajectories_raw);
                else
                    marker_trajectories = marker_trajectories_raw;
                end
                
                marker_trajectories = marker_trajectories';
                
                % save
                save_folder = 'processed';
                save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'markerTrajectories.mat');
                save ...
                  ( ...
                    [save_folder filesep save_file_name], ...
                    'marker_trajectories', ...
                    'time_mocap', ...
                    'sampling_rate_mocap', ...
                    'marker_labels' ...
                  );
                addAvailableData('marker_trajectories', 'time_mocap', 'sampling_rate_mocap', 'marker_labels', save_folder, save_file_name);
                disp(['processed ' raw_marker_file_name ' and saved as ' save_file_name])
            end
        end
    end

    %% transform to belt space
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            % load data
            this_condition = condition_list{i_condition};
            if any(strcmp(study_settings.get('conditions_to_transform_to_belt_space'), this_condition))
                % extract data for new structure
                if exist(['processed' filesep makeFileName(date, subject_id, trial_type, i_trial, 'plcData')], 'file')
                    load(['processed' filesep makeFileName(date, subject_id, trial_type, i_trial, 'plcData')])
                else
                    error(['Failed to load PLC data file for condition ' trial_type ', trial ' num2str(i_trial)]);
                end
                time_belts = time_plcData - time_plcData(1);

                % calculate shift
                belt_speed_trajectory = mean([belt_speed_left_trajectory belt_speed_right_trajectory], 2);
                delta_t = diff(time_belts);
                belt_position_trajectory_plcData = zeros(size(belt_speed_trajectory));
                for i_time = 2 : length(belt_speed_trajectory)
                    belt_position_trajectory_plcData(i_time) = belt_position_trajectory_plcData(i_time-1) + delta_t(i_time-1) * belt_speed_trajectory(i_time-1);
                end

                % apply shift to marker trajectories
                file_name_raw = ['raw' filesep makeFileName(date, subject_id, this_condition, i_trial, 'markerTrajectoriesRaw.mat')];
                load(file_name_raw);
                marker_trajectories = marker_trajectories_raw;
                belt_position_trajectory_mocap = spline(time_belts, belt_position_trajectory_plcData, time_mocap)';
                for i_marker = 1 : size(marker_headers, 2)
                    marker_trajectories(:, (i_marker-1)*3+2) = marker_trajectories(:, (i_marker-1)*3+2) + belt_position_trajectory_mocap;
                end
              
                save_folder = 'processed';
                save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'markerTrajectories.mat');
                save ...
                  ( ...
                    [save_folder filesep save_file_name], ...
                    'marker_trajectories', ...
                    'time_mocap', ...
                    'sampling_rate_mocap', ...
                    'marker_labels' ...
                  );
                addAvailableData('marker_trajectories', 'time_mocap', 'marker_labels', save_folder, save_file_name);
                disp(['Transformed marker data in ' file_name_raw ' to belt space and saved to ' file_name_shifted])                    
            end






%             % apply shift to forceplate trajectories
%             load(makeFileName(date, subject_id, condition, i_trial, 'forceplateTrajectories'));
%             belt_position_trajectory_forceplate = spline(time_belts, belt_position_trajectory_plcData, time_forceplate)';
%             
%             for i_time = 1 : length(time_forceplate)
%                 % define forceplate rotation and translation
%                 world_to_Acb_rotation = [1 0 0; 0 1 0; 0 0 1];
%                 world_to_Acb_translation = [0.5588; 0; 0];
%                 world_to_Acb_trafo = [world_to_Acb_rotation world_to_Acb_translation; 0 0 0 1];
%                 world_to_Acb_adjoint = rigidToAdjointTransformation(world_to_Acb_trafo);
% 
%                 % transform
%                 left_forceplate_wrench_Acb = (world_to_Acb_adjoint' * left_forceplate_wrench_world')';
%                 right_forceplate_wrench_Acb = (world_to_Acb_adjoint' * right_forceplate_wrench_world')';
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
%             left_forceplate_cop_world(left_forceplate_low_load_indicator, :) = 0;
%             right_forceplate_low_load_indicator = copxr_trajectory == 0;
%             right_forceplate_cop_world(right_forceplate_low_load_indicator, :) = 0;
%             total_forceplate_low_load_indicator = copx_trajectory == 0;
%             total_forceplate_cop_world(total_forceplate_low_load_indicator, :) = 0;            

        end
    end

    
    


    
    
    
    
    
    
    
    
end