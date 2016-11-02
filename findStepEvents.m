function findStepEvents(varargin)
    
    
    % figure out folders
    if ~exist('analysis', 'dir')
        mkdir('analysis')
    end

    
    [condition_list, trial_number_list] = parseTrialArguments(varargin);
    load('subjectInfo.mat', 'date', 'subject_id');


    % findStepEvents

    visualize_left          = 1;
    visualize_right         = 1;
    visualize_position      = 0;
    visualize_derivatives   = 0;

    show_forceplate         = 0;



%     trials_to_process = 1 : 20;


    % Choose Identification method for each event and foot
    left_method_touchdown = 'left_heel_position_minima';
    % left_method_touchdown = 'left_toe_velocity_minima';
    % left_method_touchdown = 'left_toe_position_minima';
    % left_method_touchdown = 'left_first_acceleration_peak';

    left_method_pushoff = 'left_first_velocity_peak';

    right_method_touchdown = 'right_heel_position_minima';
    % right_method_touchdown = 'right_toe_position_minima';
    % right_method_touchdown = 'right_toe_velocity_minima';
    % right_method_touchdown = 'right_first_acceleration_peak';
    right_method_pushoff = 'right_first_velocity_peak';

    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process

            %% prepare
            % load data
            condition = condition_list{i_condition};
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);
%             load(makeFileName(date, subject_id, condition, i_trial, 'labviewTrajectories'));
%             load(makeFileName(date, subject_id, condition, i_trial, 'forceplateTrajectories'));


            % extract data
            left_heel_marker = find(strcmp(marker_headers, 'LHEE'));
            left_toes_marker = find(strcmp(marker_headers, 'LTOE'));
            right_heel_marker = find(strcmp(marker_headers, 'RHEE'));
            right_toes_marker = find(strcmp(marker_headers, 'RTOE'));
            left_heel_marker_indices = reshape([(left_heel_marker - 1) * 3 + 1; (left_heel_marker - 1) * 3 + 2; (left_heel_marker - 1) * 3 + 3], 1, length(left_heel_marker)*3);
            left_toes_marker_indices = reshape([(left_toes_marker - 1) * 3 + 1; (left_toes_marker - 1) * 3 + 2; (left_toes_marker - 1) * 3 + 3], 1, length(left_toes_marker)*3);
            right_heel_marker_indices = reshape([(right_heel_marker - 1) * 3 + 1; (right_heel_marker - 1) * 3 + 2; (right_heel_marker - 1) * 3 + 3], 1, length(right_heel_marker)*3);
            right_toes_marker_indices = reshape([(right_toes_marker - 1) * 3 + 1; (right_toes_marker - 1) * 3 + 2; (right_toes_marker - 1) * 3 + 3], 1, length(right_toes_marker)*3);
            left_heel_marker_z_trajectory = marker_trajectories(:, left_heel_marker_indices(3));
            left_toes_marker_z_trajectory = marker_trajectories(:, left_toes_marker_indices(3));
            right_heel_marker_z_trajectory = marker_trajectories(:, right_heel_marker_indices(3));
            right_toes_marker_z_trajectory = marker_trajectories(:, right_toes_marker_indices(3));

            % calculate derivatives
            filter_order = 2;
            cutoff_frequency = 20; % cutoff frequency, in Hz
            [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate_mocap/2));	% set filter parameters for butterworth filter: 2=order of filter;
            left_heel_marker_z_vel_trajectory = deriveByTime(filtfilt(b, a, left_heel_marker_z_trajectory), 1/sampling_rate_mocap);
            right_heel_marker_z_vel_trajectory = deriveByTime(filtfilt(b, a, right_heel_marker_z_trajectory), 1/sampling_rate_mocap);
            left_heel_marker_z_acc_trajectory = deriveByTime(filtfilt(b, a, left_heel_marker_z_vel_trajectory), 1/sampling_rate_mocap);
            right_heel_marker_z_acc_trajectory = deriveByTime(filtfilt(b, a, right_heel_marker_z_vel_trajectory), 1/sampling_rate_mocap);
            left_toes_marker_z_vel_trajectory = deriveByTime(filtfilt(b, a, spline(time_mocap, left_toes_marker_z_trajectory, time_mocap)'), 1/sampling_rate_mocap);
            right_toes_marker_z_vel_trajectory = deriveByTime(filtfilt(b, a, spline(time_mocap, right_toes_marker_z_trajectory, time_mocap)'), 1/sampling_rate_mocap);
            left_toes_marker_z_acc_trajectory = deriveByTime(filtfilt(b, a, spline(time_mocap, left_toes_marker_z_vel_trajectory, time_mocap)'), 1/sampling_rate_mocap);
            right_toes_marker_z_acc_trajectory = deriveByTime(filtfilt(b, a, spline(time_mocap, right_toes_marker_z_vel_trajectory, time_mocap)'), 1/sampling_rate_mocap);



            %% find events for left foot
%             peak_width_threshold = sampling_rate_mocap * 0.25;
%             peak_prominence_threshold = 0.5;
            peak_width_threshold = sampling_rate_mocap * 0.25;
            peak_prominence_threshold = 0.25;
            % For Left Foot   
            if strcmp(left_method_touchdown, 'left_heel_position_minima')
                % find touch down indices as negative peaks of the heel marker z-position


                % left
                [~, left_heel_peak_locations, left_heel_peak_widths] = findpeaks(-left_heel_marker_z_trajectory);
                left_heel_peak_locations = left_heel_peak_locations';
                left_touchdown_indices_mocap = left_heel_peak_locations(left_heel_peak_widths > peak_width_threshold);


            elseif strcmp(left_method_touchdown, 'left_toe_position_minima')
                % left
                [~, left_toe_peak_locations, left_toe_peak_widths] = findpeaks(-left_toes_marker_z_trajectory);
                left_toe_peak_locations = left_toe_peak_locations';
                left_touchdown_indices_mocap = left_toe_peak_locations(left_toe_peak_widths > peak_width_threshold);

            elseif strcmp(left_method_touchdown, 'left_toe_velocity_minima')
                min_peak_prominence = 1.3;
                [~, left_toe_peak_locations, left_toe_peak_widths, p] = findpeaks(-left_toes_marker_z_vel_trajectory, 'MinPeakProminence', min_peak_prominence);
                left_toe_peak_locations = left_toe_peak_locations';
                left_touchdown_indices_mocap = left_toe_peak_locations;

            elseif strcmp(left_method_touchdown, 'left_first_acceleration_peak')
                % left

                % find mid-swing as peaks of heel position
                min_peak_prominence = 0.05;
                [~, left_heel_midswing_locations] = findpeaks(left_heel_marker_z_trajectory, 'MinPeakProminence', min_peak_prominence);
                left_heel_midswing_locations = left_heel_midswing_locations';
                % find acceleration peaks
                min_peak_prominence = 5;
                [~, left_heel_acc_peak_locations] = findpeaks(left_heel_marker_z_acc_trajectory, 'MinPeakProminence', min_peak_prominence);
                left_heel_acc_peak_locations = left_heel_acc_peak_locations';
                % identify acceleration peaks as touchdowns
                left_touchdown_indices_mocap = zeros(size(left_heel_midswing_locations));
                for i_step = 1 : length(left_heel_midswing_locations)
                    touchdown_index_index = left_heel_acc_peak_locations(find(left_heel_acc_peak_locations > left_heel_midswing_locations(i_step), 1, 'first'));
                    if ~isempty(touchdown_index_index)
                        left_touchdown_indices_mocap(i_step) = touchdown_index_index;
                    end
                end
                left_touchdown_indices_mocap(left_touchdown_indices_mocap==0) = [];

            end

            if strcmp(left_method_pushoff, 'left_first_velocity_peak')
                % for pushoff, find the first significant toes z-velocity peak after each heelstrike
                min_peak_prominence = 0.35;
                [~, left_toes_vel_peak_locations] = findpeaks(left_toes_marker_z_vel_trajectory, 'MinPeakProminence', min_peak_prominence);
                left_toes_vel_peak_locations = left_toes_vel_peak_locations';
                left_pushoff_indices_mocap = zeros(size(left_touchdown_indices_mocap));
                for i_touchdown = 1 : length(left_touchdown_indices_mocap)
                    pushoff_index_index = find(left_toes_vel_peak_locations > left_touchdown_indices_mocap(i_touchdown), 1, 'first');
                    if ~isempty(pushoff_index_index)
                        left_pushoff_indices_mocap(i_touchdown) = left_toes_vel_peak_locations(pushoff_index_index);
                    end
                end
                left_pushoff_indices_mocap(left_pushoff_indices_mocap==0) = [];
            end

            %% find events for right foot

            if strcmp(right_method_touchdown, 'right_heel_position_minima')
                [~, right_heel_peak_locations, right_heel_peak_widths] = findpeaks(-right_heel_marker_z_trajectory);
                right_heel_peak_locations = right_heel_peak_locations';
                right_touchdown_indices_mocap = right_heel_peak_locations(right_heel_peak_widths > peak_width_threshold);

            elseif strcmp(right_method_touchdown, 'right_toe_position_minima')
                [~, right_toe_peak_locations, right_toe_peak_widths] = findpeaks(-right_toes_marker_z_trajectory);
                right_toe_peak_locations = right_toe_peak_locations';
                right_touchdown_indices_mocap = right_toe_peak_locations(right_toe_peak_widths > peak_width_threshold);    

            elseif strcmp(right_method_touchdown, 'right_toe_velocity_minima')
                min_peak_prominence = 1.2;
                [~, right_toe_peak_locations, right_toe_peak_widths, p] = findpeaks(-right_toes_marker_z_vel_trajectory, 'MinPeakProminence', min_peak_prominence);
                right_toe_peak_locations = right_toe_peak_locations';
                right_touchdown_indices_mocap = right_toe_peak_locations;

            elseif strcmp(right_method_touchdown, 'right_first_acceleration_peak')
                % find mid-swing as peaks of heel position
                min_peak_prominence = 0.05;
                [~, right_heel_midswing_locations] = findpeaks(right_heel_marker_z_trajectory, 'MinPeakProminence', min_peak_prominence);
                right_heel_midswing_locations = right_heel_midswing_locations';
                % find acceleration peaks
                min_peak_prominence = 5;
                [~, right_heel_acc_peak_locations] = findpeaks(right_heel_marker_z_acc_trajectory, 'MinPeakProminence', min_peak_prominence);
                right_heel_acc_peak_locations = right_heel_acc_peak_locations';
                % identify acceleration peaks as touchdowns
                right_touchdown_indices_mocap = zeros(size(right_heel_midswing_locations));
                for i_step = 1 : length(right_heel_midswing_locations)
                    touchdown_index_index = right_heel_acc_peak_locations(find(right_heel_acc_peak_locations > right_heel_midswing_locations(i_step), 1, 'first'));
                    if ~isempty(touchdown_index_index)
                        right_touchdown_indices_mocap(i_step) = touchdown_index_index;
                    end
                end
                right_touchdown_indices_mocap(right_touchdown_indices_mocap==0) = [];
            end

            if strcmp(right_method_pushoff, 'right_first_velocity_peak')
                min_peak_prominence = 0.35;
                [~, right_toes_vel_peak_locations] = findpeaks(right_toes_marker_z_vel_trajectory, 'MinPeakProminence', min_peak_prominence);
                right_toes_vel_peak_locations = right_toes_vel_peak_locations';
                right_pushoff_indices_mocap = zeros(size(right_touchdown_indices_mocap));
                for i_touchdown = 1 : length(right_touchdown_indices_mocap)
                    pushoff_index_index = find(right_toes_vel_peak_locations > right_touchdown_indices_mocap(i_touchdown), 1, 'first');
                    if ~isempty(pushoff_index_index)
                        right_pushoff_indices_mocap(i_touchdown) = right_toes_vel_peak_locations(pushoff_index_index);
                    end
                end
                right_pushoff_indices_mocap(right_pushoff_indices_mocap==0) = [];
            end

            % remove duplicates
            % Achtung, untested!
            left_pushoff_indices_mocap = unique(left_pushoff_indices_mocap);
            left_touchdown_indices_mocap = unique(left_touchdown_indices_mocap);
            right_pushoff_indices_mocap = unique(right_pushoff_indices_mocap);
            right_touchdown_indices_mocap = unique(right_touchdown_indices_mocap);


%             % transform to forceplate time
%             left_pushoff_indices_forceplate = zeros(size(left_pushoff_indices_mocap));
%             for i_index = 1 : length(left_pushoff_indices_mocap)
%                 [~, index_forceplate] = min(abs(time_forceplate - time_mocap(left_pushoff_indices_mocap(i_index))));
%                 left_pushoff_indices_forceplate(i_index) = index_forceplate;
%             end
%             left_touchdown_indices_forceplate = zeros(size(left_touchdown_indices_mocap));
%             for i_index = 1 : length(left_touchdown_indices_mocap)
%                 [~, index_forceplate] = min(abs(time_forceplate - time_mocap(left_touchdown_indices_mocap(i_index))));
%                 left_touchdown_indices_forceplate(i_index) = index_forceplate;
%             end
%             right_pushoff_indices_forceplate = zeros(size(right_pushoff_indices_mocap));
%             for i_index = 1 : length(right_pushoff_indices_mocap)
%                 [~, index_forceplate] = min(abs(time_forceplate - time_mocap(right_pushoff_indices_mocap(i_index))));
%                 right_pushoff_indices_forceplate(i_index) = index_forceplate;
%             end
%             right_touchdown_indices_forceplate = zeros(size(right_touchdown_indices_mocap));
%             for i_index = 1 : length(right_touchdown_indices_mocap)
%                 [~, index_forceplate] = min(abs(time_forceplate - time_mocap(right_touchdown_indices_mocap(i_index))));
%                 right_touchdown_indices_forceplate(i_index) = index_forceplate;
%             end
% 
%             % transform to labview time
%             left_pushoff_indices_labview = zeros(size(left_pushoff_indices_mocap));
%             for i_index = 1 : length(left_pushoff_indices_mocap)
%                 [~, index_labview] = min(abs(time_labview - time_mocap(left_pushoff_indices_mocap(i_index))));
%                 left_pushoff_indices_labview(i_index) = index_labview;
%             end
%             left_touchdown_indices_labview = zeros(size(left_touchdown_indices_mocap));
%             for i_index = 1 : length(left_touchdown_indices_mocap)
%                 [~, index_labview] = min(abs(time_labview - time_mocap(left_touchdown_indices_mocap(i_index))));
%                 left_touchdown_indices_labview(i_index) = index_labview;
%             end
%             right_pushoff_indices_labview = zeros(size(right_pushoff_indices_mocap));
%             for i_index = 1 : length(right_pushoff_indices_mocap)
%                 [~, index_labview] = min(abs(time_labview - time_mocap(right_pushoff_indices_mocap(i_index))));
%                 right_pushoff_indices_labview(i_index) = index_labview;
%             end
%             right_touchdown_indices_labview = zeros(size(right_touchdown_indices_mocap));
%             for i_index = 1 : length(right_touchdown_indices_mocap)
%                 [~, index_labview] = min(abs(time_labview - time_mocap(right_touchdown_indices_mocap(i_index))));
%                 right_touchdown_indices_labview(i_index) = index_labview;
%             end

            % calculate times
            left_pushoff_times = time_mocap(left_pushoff_indices_mocap);
            left_touchdown_times = time_mocap(left_touchdown_indices_mocap);
            right_pushoff_times = time_mocap(right_pushoff_indices_mocap);
            right_touchdown_times = time_mocap(right_touchdown_indices_mocap);

            %% form contact indicators
            number_of_time_steps_mocap = length(time_mocap);
%             number_of_time_steps_forceplate = length(time_forceplate);
%             number_of_time_steps_labview = length(time_labview);

            left_contact_indicators_mocap = formContactIndicatorTrajectory(left_pushoff_indices_mocap, left_touchdown_indices_mocap, number_of_time_steps_mocap);
            right_contact_indicators_mocap = formContactIndicatorTrajectory(right_pushoff_indices_mocap, right_touchdown_indices_mocap, number_of_time_steps_mocap);
%             left_contact_indicators_forceplate = formContactIndicatorTrajectory(left_pushoff_indices_forceplate, left_touchdown_indices_forceplate, number_of_time_steps_forceplate);
%             right_contact_indicators_forceplate = formContactIndicatorTrajectory(right_pushoff_indices_forceplate, right_touchdown_indices_forceplate, number_of_time_steps_forceplate);
%             left_contact_indicators_labview = formContactIndicatorTrajectory(left_pushoff_indices_labview, left_touchdown_indices_labview, number_of_time_steps_labview);
%             right_contact_indicators_labview = formContactIndicatorTrajectory(right_pushoff_indices_labview, right_touchdown_indices_labview, number_of_time_steps_labview);

            % form contact trajectories
            left_heel_contact_trajectories = left_heel_marker_z_trajectory; left_heel_contact_trajectories(~left_contact_indicators_mocap, :) = NaN;
            left_toes_contact_trajectories = left_toes_marker_z_trajectory; left_toes_contact_trajectories(~left_contact_indicators_mocap, :) = NaN;
            left_toes_velocity_contact_trajectories = left_toes_marker_z_vel_trajectory; left_toes_velocity_contact_trajectories(~left_contact_indicators_mocap, :) = NaN;

            right_heel_contact_trajectories = right_heel_marker_z_trajectory; right_heel_contact_trajectories(~right_contact_indicators_mocap, :) = NaN;
            right_toes_contact_trajectories = right_toes_marker_z_trajectory; right_toes_contact_trajectories(~right_contact_indicators_mocap, :) = NaN;
            right_toes_velocity_contact_trajectories = right_toes_marker_z_vel_trajectory; right_toes_velocity_contact_trajectories(~right_contact_indicators_mocap, :) = NaN;





            %% visualize
            color_heelstrike = [1 0 0];
            color_pushoff = [0 1 0];

            force_scaler = 2e-4;
            vel_scaler = 1;
            acc_scaler = .2;

%             fzl_trajectory_mocap = spline(time_forceplate, fzl_trajectory, time_mocap);
%             fzr_trajectory_mocap = spline(time_forceplate, fzr_trajectory, time_mocap);


            if visualize_left
                if visualize_position
                    % left events
                    figure; axes_left = axes; hold on; title('left foot marker positions')
                    plot(time_mocap, left_heel_marker_z_trajectory, 'linewidth', 1);
                    plot(time_mocap, left_toes_marker_z_trajectory, 'linewidth', 1);
                    plot(time_mocap(left_touchdown_indices_mocap), left_heel_marker_z_trajectory(left_touchdown_indices_mocap), 'o', 'linewidth', 2, 'color', color_heelstrike);
                    plot(time_mocap(left_pushoff_indices_mocap), left_toes_marker_z_trajectory(left_pushoff_indices_mocap), 'o', 'linewidth', 2, 'color', color_pushoff);
                    if show_forceplate
                        plot(time_forceplate, fzl_trajectory*force_scaler)
                        legend('heel', 'toes', 'touchdown', 'pushoff', 'fzl');
                        plot(time_mocap(left_touchdown_indices_mocap), fzl_trajectory_mocap(left_touchdown_indices_mocap)*force_scaler, 'o', 'linewidth', 2, 'color', color_heelstrike);
                        plot(time_mocap(left_pushoff_indices_mocap), fzl_trajectory_mocap(left_pushoff_indices_mocap)*force_scaler, 'o', 'linewidth', 2, 'color', color_pushoff);
                    else
                        legend('heel', 'toes', 'touchdown', 'pushoff');
                    end
                end

                if visualize_derivatives
                    figure; axes_left_derivatives = axes; hold on;  title('left foot marker derivatives')
            %         plot(time_mocap, left_heel_marker_z_vel_trajectory, 'linewidth', 1);
                    plot(time_mocap, left_heel_marker_z_acc_trajectory*acc_scaler, 'linewidth', 1);
                    plot(time_mocap, left_toes_marker_z_vel_trajectory*vel_scaler, 'linewidth', 1);
            %         plot(time_mocap, left_toes_marker_z_acc_trajectory, 'linewidth', 1);
            %         plot(time_mocap(left_touchdown_indices_mocap), left_heel_marker_z_acc_trajectory(left_touchdown_indices_mocap), 'o', 'linewidth', 2);
                    legend('heel acc', 'toes vel');
                    plot(time_mocap(left_touchdown_indices_mocap), left_heel_marker_z_acc_trajectory(left_touchdown_indices_mocap)*acc_scaler, 'o', 'linewidth', 2, 'color', color_heelstrike);
                    plot(time_mocap(left_pushoff_indices_mocap), left_toes_marker_z_vel_trajectory(left_pushoff_indices_mocap)*vel_scaler, 'o', 'linewidth', 2, 'color', color_pushoff);

                    linkaxes([axes_left, axes_left_derivatives], 'x')
            %         distFig('rows', 2)
                end
            end
            if visualize_right
                if visualize_position
                    % right events
                    figure; axes_right = axes; hold on; title('right foot marker positions')
                    plot(time_mocap, right_heel_marker_z_trajectory, 'linewidth', 1);
                    plot(time_mocap, right_toes_marker_z_trajectory, 'linewidth', 1);
                    plot(time_mocap(right_touchdown_indices_mocap), right_heel_marker_z_trajectory(right_touchdown_indices_mocap), 'o', 'linewidth', 2, 'color', color_heelstrike);
                    plot(time_mocap(right_pushoff_indices_mocap), right_toes_marker_z_trajectory(right_pushoff_indices_mocap), 'o', 'linewidth', 2, 'color', color_pushoff);
                    if show_forceplate
                        plot(time_forceplate, fzr_trajectory*force_scaler)
                        legend('heel', 'toes', 'touchdown', 'pushoff', 'fzl');
                        plot(time_mocap(right_touchdown_indices_mocap), fzr_trajectory_mocap(right_touchdown_indices_mocap)*force_scaler, 'o', 'linewidth', 2, 'color', color_heelstrike);
                        plot(time_mocap(right_pushoff_indices_mocap), fzr_trajectory_mocap(right_pushoff_indices_mocap)*force_scaler, 'o', 'linewidth', 2, 'color', color_pushoff);
                    else
                        legend('heel', 'toes', 'touchdown', 'pushoff');
                    end
                end

                if visualize_derivatives
                    figure; axes_right_derivatives = axes; hold on; title('right foot marker derivatives')
            %         plot(time_mocap, right_heel_marker_z_vel_trajectory, 'linewidth', 1);
                    plot(time_mocap, right_heel_marker_z_acc_trajectory*acc_scaler, 'linewidth', 1);
                    plot(time_mocap, right_toes_marker_z_vel_trajectory*vel_scaler, 'linewidth', 1);
            %         plot(time_mocap, right_toes_marker_z_acc_trajectory, 'linewidth', 1);
                    legend('heel acc', 'toes vel');
                    plot(time_mocap(right_touchdown_indices_mocap), right_heel_marker_z_acc_trajectory(right_touchdown_indices_mocap)*acc_scaler, 'o', 'linewidth', 2, 'color', color_heelstrike);
                    plot(time_mocap(right_pushoff_indices_mocap), right_toes_marker_z_vel_trajectory(right_pushoff_indices_mocap)*vel_scaler, 'o', 'linewidth', 2, 'color', color_pushoff);

                    linkaxes([axes_right, axes_right_derivatives], 'x')
                end
            end


        %     linkaxes(getAllAxes, 'x')
            if isrow(left_pushoff_times)
                left_pushoff_times = left_pushoff_times';
            end
            if isrow(left_touchdown_times)
                left_touchdown_times = left_touchdown_times';
            end
            if isrow(right_pushoff_times)
                right_pushoff_times = right_pushoff_times';
            end
            if isrow(right_touchdown_times)
                right_touchdown_times = right_touchdown_times';
            end

            %% save
            step_events_file_name = ['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'stepEvents')];
            save ...
              ( ...
                step_events_file_name, ...
                'left_pushoff_times', ...
                'left_touchdown_times', ...
                'right_pushoff_times', ...
                'right_touchdown_times' ...
              );
%                 'left_touchdown_indices_mocap', ...
%                 'right_touchdown_indices_mocap', ...
%                 'left_pushoff_indices_mocap', ...
%                 'right_pushoff_indices_mocap', ...
%                 'left_contact_indicators_mocap', ...
%                 'right_contact_indicators_mocap', ...
%                 'left_touchdown_indices_forceplate', ...
%                 'right_touchdown_indices_forceplate', ...
%                 'left_pushoff_indices_forceplate', ...
%                 'right_pushoff_indices_forceplate', ...
%                 'left_contact_indicators_forceplate', ...
%                 'right_contact_indicators_forceplate', ...
%                 'left_touchdown_indices_labview', ...
%                 'right_touchdown_indices_labview', ...
%                 'left_pushoff_indices_labview', ...
%                 'right_pushoff_indices_labview', ...
%                 'left_contact_indicators_labview', ...
%                 'right_contact_indicators_labview', ...

            disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' step_events_file_name]);
        end
    end
    if ~isempty(findobj('type','figure'))
        distFig('rows', 2)
    end
end















