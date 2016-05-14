% findStepEvents

visualize_left = 1;
visualize_right = 1;
visualize_position = 1;
visualize_derivatives = 1;
show_forceplate = 0;

trials_to_process = 1 : 1 : 21;
% trials_to_process = [0:7 9:23];
trials_to_process = 1;


% Choose Identification method for each event and foot
% left_method_touchdown = 'left_heel_position_minima';
left_method_touchdown = 'left_toe_position_minima';
% method_touchdown = 'zpos_threshold';
left_method_pushoff = 'left_first_velocity_peak';

right_method_touchdown = 'right_heel_position_minima';
% right_method_touchdown = 'right_toe_position_minima';
right_method_pushoff = 'right_first_velocity_peak';

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
    %         peak_width_threshold = samplingRate * 0.25;
%         peak_prominence_threshold = 0.5;
        peak_width_threshold = samplingRate * 0.25;
        peak_prominence_threshold = 0.25;
    % For Left Foot   
    if strcmp(left_method_touchdown, 'left_heel_position_minima')
        % find touch down indices as negative peaks of the heel marker z-position

        
        % left
        [~, left_heel_peak_locations, left_heel_peak_widths] = findpeaks(-left_heel_marker_z_trajectory);
        left_touchdown_indices_mocap = left_heel_peak_locations(left_heel_peak_widths > peak_width_threshold);
        
       
    elseif strcmp(left_method_touchdown, 'left_toe_position_minima')
        % left
        [~, left_toe_peak_locations, left_toe_peak_widths] = findpeaks(-left_toes_marker_z_trajectory);
        left_touchdown_indices_mocap = left_toe_peak_locations(left_toe_peak_widths > peak_width_threshold);
        
    elseif strcmp(left_method_touchdown, 'left_first_acceleration_peak')
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
 
    end
    
    if strcmp(left_method_pushoff, 'left_first_velocity_peak')
        % for pushoff, find the first significant toes z-velocity peak after each heelstrike
        min_peak_prominence = 0.25;
        [~, left_toes_vel_peak_locations] = findpeaks(left_toes_marker_z_vel_trajectory, 'MinPeakProminence', min_peak_prominence);
        left_pushoff_indices_mocap = zeros(size(left_touchdown_indices_mocap));
        for i_touchdown = 1 : length(left_touchdown_indices_mocap)
            pushoff_index_index = find(left_toes_vel_peak_locations > left_touchdown_indices_mocap(i_touchdown), 1, 'first');
            if ~isempty(pushoff_index_index)
                left_pushoff_indices_mocap(i_touchdown) = left_toes_vel_peak_locations(pushoff_index_index);
            end
        end
        left_pushoff_indices_mocap(left_pushoff_indices_mocap==0) = [];
    end
     
        %% RIGHT FOOT IDENTIFY
     
    if strcmp(right_method_touchdown, 'right_heel_position_minima')
        [~, right_heel_peak_locations, right_heel_peak_widths] = findpeaks(-right_heel_marker_z_trajectory);
        right_touchdown_indices_mocap = right_heel_peak_locations(right_heel_peak_widths > peak_width_threshold);
    
    elseif strcmp(left_method_touchdown, 'left_toe_position_minima')
        [~, right_toe_peak_locations, right_toe_peak_widths] = findpeaks(-right_toes_marker_z_trajectory);
        right_touchdown_indices_mocap = right_toe_peak_locations(right_toe_peak_widths > peak_width_threshold);    
         
    elseif strcmp(left_method_touchdown, 'left_first_acceleration_peak')       
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
        
 if strcmp(right_method_pushoff, 'right_first_velocity_peak')
 min_peak_prominence = 0.25;
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
    number_of_time_steps_mocap = length(time_mocap);
    number_of_time_steps_force_plate = length(time_force_plate);

    left_contact_indicators_mocap = formContactIndicatorTrajectory(left_pushoff_indices_mocap, left_touchdown_indices_mocap, number_of_time_steps_mocap);
    right_contact_indicators_mocap = formContactIndicatorTrajectory(right_pushoff_indices_mocap, right_touchdown_indices_mocap, number_of_time_steps_mocap);
    left_contact_indicators_force_plate = formContactIndicatorTrajectory(left_pushoff_indices_force_plate, left_touchdown_indices_force_plate, number_of_time_steps_force_plate);
    right_contact_indicators_force_plate = formContactIndicatorTrajectory(right_pushoff_indices_force_plate, right_touchdown_indices_force_plate, number_of_time_steps_force_plate);

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

    fzl_trajectory_mocap = spline(time_force_plate, fzl_trajectory, time_mocap);
    fzr_trajectory_mocap = spline(time_force_plate, fzr_trajectory, time_mocap);
    

    if visualize_left
        if visualize_position
            % left events
            figure; axes_left = axes; hold on; title('left foot marker positions')
            plot(time_mocap, left_heel_marker_z_trajectory, 'linewidth', 1);
            plot(time_mocap, left_toes_marker_z_trajectory, 'linewidth', 1);
            plot(time_mocap(left_touchdown_indices_mocap), left_heel_marker_z_trajectory(left_touchdown_indices_mocap), 'o', 'linewidth', 2, 'color', color_heelstrike);
            plot(time_mocap(left_pushoff_indices_mocap), left_toes_marker_z_trajectory(left_pushoff_indices_mocap), 'o', 'linewidth', 2, 'color', color_pushoff);
            if show_forceplate
                plot(time_force_plate, fzl_trajectory*force_scaler)
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
                plot(time_force_plate, fzr_trajectory*force_scaler)
                legend('heel', 'toes', 'touchdown', 'pushoff', 'fzl');
                plot(time_mocap(right_touchdown_indices_mocap), fzr_trajectory_mocap(right_touchdown_indices_mocap)*force_scaler, 'o', 'linewidth', 2, 'color', color_heelstrike);
                plot(time_mocap(right_pushoff_indices_mocap), fzr_trajectory_mocap(right_pushoff_indices_mocap)*force_scaler, 'o', 'linewidth', 2, 'color', color_pushoff);
            else
                legend('heel', 'toes', 'touchdown', 'pushoff');
            end
        end
        
        if visualize_derivatives
            figure; axes_right_derivatives = axes; hold on; title('left foot marker derivatives')
    %         plot(time_mocap, right_heel_marker_z_vel_trajectory, 'linewidth', 1);
            plot(time_mocap, right_heel_marker_z_acc_trajectory*acc_scaler, 'linewidth', 1);
            plot(time_mocap, right_toes_marker_z_vel_trajectory*vel_scaler, 'linewidth', 1);
    %         plot(time_mocap, right_toes_marker_z_acc_trajectory, 'linewidth', 1);
            legend('heel acc', 'toes vel');
            plot(time_mocap(right_touchdown_indices_mocap), right_heel_marker_z_acc_trajectory(right_touchdown_indices_mocap)*acc_scaler, 'o', 'linewidth', 2, 'color', color_heelstrike);
            plot(time_mocap(right_pushoff_indices_mocap), right_toes_marker_z_vel_trajectory(right_pushoff_indices_mocap)*vel_scaler, 'o', 'linewidth', 2, 'color', color_pushoff);

            linkaxes([axes_right, axes_right_derivatives], 'x')
            distFig('rows', 2)
        end
    end


%     linkaxes(getAllAxes, 'x')



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
    
    if visualize_left || visualize_right
        keyboard
        distFig
    end
end















