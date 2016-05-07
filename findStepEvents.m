% findStepEvents

visualize_left = 1;
visualize_right = 1;

trials_to_process = 0 : 21;
% trials_to_process = [0:7 9:23];
trials_to_process = 4;


% method_touchdown = 'position_minima';
method_touchdown = 'first_acceleration_peak';

method_pushoff = 'first_velocity_peak';



for i_trial = trials_to_process

    %% prepare

    % load data
    load subjectInfo.mat;
    load(makeFileName(date, subject_id, 'walking', i_trial, 'markerTrajectories'));
    load(makeFileName(date, subject_id, 'walking', i_trial, 'forcePlateData'));

    % extract data
    left_heel_marker = 34;
    left_toes_marker = 35;
    right_heel_marker = 42;
    right_toes_marker = 43;
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
    [b, a] = butter(filter_order, cutoff_frequency/(samplingRate/2));	% set filter parameters for butterworth filter: 2=order of filter;
    left_heel_marker_z_vel_trajectory = deriveByTime(filtfilt(b, a, left_heel_marker_z_trajectory), 1/samplingRate);
    right_heel_marker_z_vel_trajectory = deriveByTime(filtfilt(b, a, right_heel_marker_z_trajectory), 1/samplingRate);
    left_heel_marker_z_acc_trajectory = deriveByTime(filtfilt(b, a, left_heel_marker_z_vel_trajectory), 1/samplingRate);
    right_heel_marker_z_acc_trajectory = deriveByTime(filtfilt(b, a, right_heel_marker_z_vel_trajectory), 1/samplingRate);
    left_toes_marker_z_vel_trajectory = deriveByTime(filtfilt(b, a, left_toes_marker_z_trajectory), 1/samplingRate);
    right_toes_marker_z_vel_trajectory = deriveByTime(filtfilt(b, a, right_toes_marker_z_trajectory), 1/samplingRate);
    left_toes_marker_z_acc_trajectory = deriveByTime(filtfilt(b, a, left_toes_marker_z_vel_trajectory), 1/samplingRate);
    right_toes_marker_z_acc_trajectory = deriveByTime(filtfilt(b, a, right_toes_marker_z_vel_trajectory), 1/samplingRate);
    
    

    %% find events
    
    if strcmp(method_touchdown, 'position_minima')
        % find touch down indices as negative peaks of the heel marker z-position
        peak_width_threshold = samplingRate * 0.25;
        peak_prominence_threshold = 0.5;
        peak_width_threshold = samplingRate * 0.25;
        peak_prominence_threshold = 0.25;
        
        % left
        [~, left_heel_peak_locations, left_heel_peak_widths] = findpeaks(-left_heel_marker_z_trajectory);
        left_touchdown_indices_mocap = left_heel_peak_locations(left_heel_peak_widths > peak_width_threshold);
        
        % right
        [~, right_heel_peak_locations, right_heel_peak_widths] = findpeaks(-right_heel_marker_z_trajectory);
        right_touchdown_indices_mocap = right_heel_peak_locations(right_heel_peak_widths > peak_width_threshold);
    elseif strcmp(method_touchdown, 'first_acceleration_peak')
        % left
        
        % find mid-swing as peaks of heel position
        min_peak_prominence = 0.05;
        [~, left_heel_midswing_locations] = findpeaks(left_heel_marker_z_trajectory, 'MinPeakProminence', min_peak_prominence);
        % find acceleration peaks
        min_peak_prominence = 5;
        [~, left_heel_acc_peak_locations] = findpeaks(left_heel_marker_z_acc_trajectory, 'MinPeakProminence', min_peak_prominence);
        % identify acceleration peaks as touchdowns
        left_touchdown_indices_mocap = zeros(size(left_heel_midswing_locations));
        for i_step = 1 : length(left_heel_midswing_locations)
            touchdown_index_index = left_heel_acc_peak_locations(find(left_heel_acc_peak_locations > left_heel_midswing_locations(i_step), 1, 'first'));
            if ~isempty(touchdown_index_index)
                left_touchdown_indices_mocap(i_step) = touchdown_index_index;
            end
        end
        left_touchdown_indices_mocap(left_touchdown_indices_mocap==0) = [];
        
        % right
        
        % find mid-swing as peaks of heel position
        min_peak_prominence = 0.05;
        [~, right_heel_midswing_locations] = findpeaks(right_heel_marker_z_trajectory, 'MinPeakProminence', min_peak_prominence);
        % find acceleration peaks
        min_peak_prominence = 5;
        [~, right_heel_acc_peak_locations] = findpeaks(right_heel_marker_z_acc_trajectory, 'MinPeakProminence', min_peak_prominence);
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
    
    if strcmp(method_pushoff, 'first_velocity_peak')
        % for pushoff, find the first significant toes z-velocity peak after each heelstrike
        min_peak_prominence = 0.1;
        [~, left_toes_vel_peak_locations] = findpeaks(left_toes_marker_z_vel_trajectory, 'MinPeakProminence', min_peak_prominence);
        left_pushoff_indices_mocap = zeros(size(left_touchdown_indices_mocap));
        for i_touchdown = 1 : length(left_touchdown_indices_mocap)
            pushoff_index_index = find(left_toes_vel_peak_locations > left_touchdown_indices_mocap(i_touchdown), 1, 'first');
            if ~isempty(pushoff_index_index)
                left_pushoff_indices_mocap(i_touchdown) = left_toes_vel_peak_locations(pushoff_index_index);
            end
        end
        left_pushoff_indices_mocap(left_pushoff_indices_mocap==0) = [];

        min_peak_prominence = 0.1;
        [~, right_toes_vel_peak_locations] = findpeaks(right_toes_marker_z_vel_trajectory, 'MinPeakProminence', min_peak_prominence);
        right_pushoff_indices_mocap = zeros(size(right_touchdown_indices_mocap));
        for i_touchdown = 1 : length(right_touchdown_indices_mocap)
            pushoff_index_index = find(right_toes_vel_peak_locations > right_touchdown_indices_mocap(i_touchdown), 1, 'first');
            if ~isempty(pushoff_index_index)
                right_pushoff_indices_mocap(i_touchdown) = right_toes_vel_peak_locations(pushoff_index_index);
            end
        end
        right_pushoff_indices_mocap(right_pushoff_indices_mocap==0) = [];
    end
%     right_pushoff_indices_mocap = right_toes_vel_peak_locations;
    
    % transform to force plate time
    left_pushoff_indices_force_plate = zeros(size(left_pushoff_indices_mocap));
    for i_index = 1 : length(left_pushoff_indices_mocap)
        [~, index_force_plate] = min(abs(time_force_plate - time_mocap(left_pushoff_indices_mocap(i_index))));
        left_pushoff_indices_force_plate(i_index) = index_force_plate;
    end
    left_touchdown_indices_force_plate = zeros(size(left_touchdown_indices_mocap));
    for i_index = 1 : length(left_touchdown_indices_mocap)
        [~, index_force_plate] = min(abs(time_force_plate - time_mocap(left_touchdown_indices_mocap(i_index))));
        left_touchdown_indices_force_plate(i_index) = index_force_plate;
    end
    right_pushoff_indices_force_plate = zeros(size(right_pushoff_indices_mocap));
    for i_index = 1 : length(right_pushoff_indices_mocap)
        [~, index_force_plate] = min(abs(time_force_plate - time_mocap(right_pushoff_indices_mocap(i_index))));
        right_pushoff_indices_force_plate(i_index) = index_force_plate;
    end
    right_touchdown_indices_force_plate = zeros(size(right_touchdown_indices_mocap));
    for i_index = 1 : length(right_touchdown_indices_mocap)
        [~, index_force_plate] = min(abs(time_force_plate - time_mocap(right_touchdown_indices_mocap(i_index))));
        right_touchdown_indices_force_plate(i_index) = index_force_plate;
    end



    %% form contact indicators
    number_of_time_steps = length(time_mocap);

    left_contact_indicators_mocap = formContactIndicatorTrajectory(left_pushoff_indices_mocap, left_touchdown_indices_mocap, number_of_time_steps);
    right_contact_indicators_mocap = formContactIndicatorTrajectory(right_pushoff_indices_mocap, right_touchdown_indices_mocap, number_of_time_steps);
    left_contact_indicators_force_plate = formContactIndicatorTrajectory(left_pushoff_indices_force_plate, left_touchdown_indices_force_plate, number_of_time_steps);
    right_contact_indicators_force_plate = formContactIndicatorTrajectory(right_pushoff_indices_force_plate, right_touchdown_indices_force_plate, number_of_time_steps);

    % form contact trajectories
    left_heel_contact_trajectories = left_heel_marker_z_trajectory; left_heel_contact_trajectories(~left_contact_indicators_mocap, :) = NaN;
    left_toes_contact_trajectories = left_toes_marker_z_trajectory; left_toes_contact_trajectories(~left_contact_indicators_mocap, :) = NaN;
    left_toes_velocity_contact_trajectories = left_toes_marker_z_vel_trajectory; left_toes_velocity_contact_trajectories(~left_contact_indicators_mocap, :) = NaN;

    right_heel_contact_trajectories = right_heel_marker_z_trajectory; right_heel_contact_trajectories(~right_contact_indicators_mocap, :) = NaN;
    right_toes_contact_trajectories = right_toes_marker_z_trajectory; right_toes_contact_trajectories(~right_contact_indicators_mocap, :) = NaN;
    right_toes_velocity_contact_trajectories = right_toes_marker_z_vel_trajectory; right_toes_velocity_contact_trajectories(~right_contact_indicators_mocap, :) = NaN;





    %% visualize
    force_scaler = 2e-4;
    vel_scaler = 5e-2;
    acc_scaler = 1e-2;


    if visualize_left
        % left events
        figure; axes_left = axes; hold on
        plot(time_force_plate, fzl_trajectory*force_scaler)
        plot(time_mocap, left_heel_marker_z_trajectory, 'linewidth', 1);
        plot(time_mocap, left_toes_marker_z_trajectory, 'linewidth', 1);
        plot(time_mocap(left_touchdown_indices_mocap), left_heel_marker_z_trajectory(left_touchdown_indices_mocap)*0, 'o', 'linewidth', 2);
        plot(time_mocap(left_pushoff_indices_mocap), left_toes_marker_z_trajectory(left_pushoff_indices_mocap)*0, 'o', 'linewidth', 2);
        legend('fzl', 'heel', 'toes', 'touchdown', 'pushoff');
        
        figure; axes_left_derivatives = axes; hold on
        plot(time_mocap, left_heel_marker_z_vel_trajectory, 'linewidth', 1);
        plot(time_mocap, left_toes_marker_z_vel_trajectory, 'linewidth', 1);
%         plot(time_mocap, left_heel_marker_z_acc_trajectory, 'linewidth', 1);
%         plot(time_mocap, left_toes_marker_z_acc_trajectory, 'linewidth', 1);
%         plot(time_mocap(left_touchdown_indices_mocap), left_heel_marker_z_acc_trajectory(left_touchdown_indices_mocap), 'o', 'linewidth', 2);
        legend('heel vel', 'toes vel', 'heel acc', 'toes acc');
        plot(time_mocap(left_pushoff_indices_mocap), left_toes_marker_z_vel_trajectory(left_pushoff_indices_mocap), 'o', 'linewidth', 2);
        
        linkaxes([axes_left, axes_left_derivatives], 'x')
%         distFig('rows', 2)
    end
    if visualize_right
        % right events
        figure; axes_right = axes; hold on
        plot(time_force_plate, fzr_trajectory*force_scaler)
        plot(time_mocap, right_heel_marker_z_trajectory, 'linewidth', 1);
        plot(time_mocap, right_toes_marker_z_trajectory, 'linewidth', 1);
        plot(time_mocap(right_touchdown_indices_mocap), right_toes_marker_z_trajectory(right_touchdown_indices_mocap)*0, 'o', 'linewidth', 2);
        plot(time_mocap(right_pushoff_indices_mocap), right_toes_marker_z_trajectory(right_pushoff_indices_mocap)*0, 'o', 'linewidth', 2);
        legend('fzl', 'heel', 'toes', 'touchdown', 'pushoff');
        
        figure; axes_right_derivatives = axes; hold on
        plot(time_mocap, right_heel_marker_z_vel_trajectory, 'linewidth', 1);
        plot(time_mocap, right_toes_marker_z_vel_trajectory, 'linewidth', 1);
%         plot(time_mocap, right_heel_marker_z_acc_trajectory, 'linewidth', 1);
%         plot(time_mocap, right_toes_marker_z_acc_trajectory, 'linewidth', 1);
        legend('heel vel', 'toes vel', 'heel acc', 'toes acc');
        plot(time_mocap(right_pushoff_indices_mocap), right_toes_marker_z_vel_trajectory(right_pushoff_indices_mocap), 'o', 'linewidth', 2);

        linkaxes([axes_right, axes_right_derivatives], 'x')
        distFig('rows', 2)
    end






    step_events_file_name = makeFileName(date, subject_id, 'walking', i_trial, 'stepEvents');
    save ...
      ( ...
        step_events_file_name, ...
        'left_touchdown_indices_mocap', ...
        'right_touchdown_indices_mocap', ...
        'left_pushoff_indices_mocap', ...
        'right_pushoff_indices_mocap', ...
        'left_contact_indicators_mocap', ...
        'right_contact_indicators_mocap', ...
        'left_touchdown_indices_force_plate', ...
        'right_touchdown_indices_force_plate', ...
        'left_pushoff_indices_force_plate', ...
        'right_pushoff_indices_force_plate', ...
        'left_contact_indicators_force_plate', ...
        'right_contact_indicators_force_plate' ...
      );
    
    disp(['Trial ' num2str(i_trial) ' completed']);


end















