% data source flags
% data_source = 'current_folder';
data_source = 'sim';
save_data       = 0;

plot_source     = 0;

save_figures    = 1; save_label = 'stepResponse_fitSms';

% update flags
update_linear   = 0;
update_splined  = 0;
update_pid      = 0;
update_pidAvg   = 1;
update_sms      = 1;
update_pidf     = 0;
update_mmj      = 0;

update_all      = 0;

% plot flags - controllers
plot_linear   = 0;
plot_splined  = 0;
plot_pid      = 0;
plot_pidAvg   = 1;
plot_sms      = 1;
plot_pidf     = 0;
plot_mmj      = 0;

% plot flags - data
plot_pos_data = 1;
plot_vel_data = 1;
plot_acc_data = 0;
plot_jrk_data = 0;
plot_rms_data = 0;


t_min = 10;
t_max = 120;

xlim_1 = 0;
xlim_2 = 30;

% these values are what we usually use, corresponding to K_p = 0.3, K_d = 0.03 min
k_p_pid = 0.3;
k_v_pid = 0.54;

% hand-fitted values to be similar to the SMS output for 0.1, 0.25
% k_p_pid = 0.2;
% k_v_pid = 0.4;

% optimized values to be similar to the SMS output for 0.1, 0.25
% k_p_pid = 0.215380825420512;
% k_v_pid = 0.446725659501051;

% values from SMS paper
k_p_sms = 0.1;
k_v_sms = 0.25;

% these values are fit to have position error output close to what our standard settings for PD do
k_p_sms = 0.1370;
k_v_sms = 0.3071;

% preparing data, depending on source
% data needed:
% pos - body position overground
% vel - body velocity, d pos / d t
% time - time data
% heelstrike_indices - indices in time array at which heel-strikes occurred

if strcmp(data_source, 'current_folder')
    % define data trial
%     trial_type = 'default'; % for SP01
%     trial_type = 'pdNormal'; % for SP02
    trial_type = 'pdAveraged'; % for SP02
%     trial_type = 'seungmoon'; % for SP02
    trial_number = 1;

    % load data
    subject_settings = loadSettingsFromFile('subject');
    marker_data_belt = load(['processed' filesep makeFileName(subject_settings.get('collection_date'), subject_settings.get('subject_id'), trial_type, trial_number, 'markerBeltTrajectories.mat')]);
    marker_data_world = load(['processed' filesep makeFileName(subject_settings.get('collection_date'), subject_settings.get('subject_id'), trial_type, trial_number, 'markerTrajectories.mat')]);
    time = marker_data_belt.time_mocap;
    plc_data = load(['processed' filesep makeFileName(subject_settings.get('collection_date'), subject_settings.get('subject_id'), trial_type, trial_number, 'PLCData.mat')]);
    event_data = load(['analysis' filesep makeFileName(subject_settings.get('collection_date'), subject_settings.get('subject_id'), trial_type, trial_number, 'events.mat')]);

    % extract and pre-process data
    belt_speed = (plc_data.belt_speed_left_trajectory + plc_data.belt_speed_right_trajectory) * 0.5;
    belt_time = plc_data.time;
    belt_speed_resampled = spline(belt_time, belt_speed, time);
    LPSI_pos = extractMarkerData(marker_data_belt.marker_trajectories, marker_data_belt.marker_labels, 'LPSI');
    RPSI_pos = extractMarkerData(marker_data_belt.marker_trajectories, marker_data_belt.marker_labels, 'RPSI');
    pos = (LPSI_pos(:, 2) + RPSI_pos(:, 2)) * 0.5;
    v_init = belt_speed_resampled(1);

    LPSI_pos_world = extractMarkerData(marker_data_world.marker_trajectories, marker_data_world.marker_labels, 'LPSI');
    RPSI_pos_world = extractMarkerData(marker_data_world.marker_trajectories, marker_data_world.marker_labels, 'RPSI');
    pos_actual = (LPSI_pos_world(:, 2) + RPSI_pos_world(:, 2)) * 0.5;

    heelstrike_left_times = event_data.event_data{strcmp(event_data.event_labels, 'left_touchdown')};
    heelstrike_left_indices = findClosestIndex(heelstrike_left_times, time);
    heelstrike_right_times = event_data.event_data{strcmp(event_data.event_labels, 'right_touchdown')};
    heelstrike_right_indices = findClosestIndex(heelstrike_right_times, time);
    heelstrike_indices = sort([heelstrike_left_indices; heelstrike_right_indices]);
    
    vel = deriveByTime(pos, time);
    
    time_indices_to_analyze = (time > t_min & time < t_max);
end
if strcmp(data_source, 'sim')
    T = 60;
    dt = 0.001;
    time = (dt : dt : T)';
    
    step_frequency = 2;
    oscillation_amplitude = 0;
%     oscillation_amplitude = 0.25;
    
    vel_base = ones(size(time));
    
    t_acc_start = 5;
    t_acc_duration = 2;
    
    speed_base = 1;
    speed_delta = 1;
    
    vel_base = speed_base + speed_delta * sinusSigmoid(time, t_acc_start, t_acc_start + t_acc_duration);
    
    
    vel_oscillation = oscillation_amplitude * sin(2 * pi * time * step_frequency);
    
    heelstrike_oscillation = sin(2 * pi * time * step_frequency);
    [~, heelstrike_indices] = findpeaks(heelstrike_oscillation);
    
    vel = vel_base + vel_oscillation;
    pos = cumtrapz(time, vel);
    v_init = vel(1);
    
    time_indices_to_analyze = (time > t_min & time < t_max);
    
end
if save_data
    heelstrike_indicator = zeros(size(time));
    heelstrike_indicator(heelstrike_indices) = 1;
    data_to_save_header = {'time', 'pos', 'vel', 'pos_actual', 'heelstrike_indicator', 'belt_speed_resampled'};
    data_to_save_table = table(time, pos, vel, pos_actual, heelstrike_indicator, belt_speed_resampled, 'VariableNames', data_to_save_header);
    filename = ['test_data_' data_source '.csv'];
    writetable(data_to_save_table, filename);
    
    
end
    
if plot_source
    % plot belt space
    if plot_pos_data
        pos_fig = figure; axes('fontsize', 16); hold on; xlim([0, 20])
        plot(time, pos, 'linewidth', 5)
        xlabel('time (s)'); ylabel('pelvis position (m)')
    % plot(time(heelstrike_indices), pos(heelstrike_indices), 'v');
        print(pos_fig, 'posBeltSpace', '-djpeg', '-r150')
    end
    
    if plot_vel_data
        vel_fig = figure; axes('fontsize', 16); hold on; xlim([0, 20])
        xlabel('time (s)'); ylabel('pelvis velocity (m/s)')
        plot(time, vel, 'linewidth', 5)
        print(vel_fig, 'velBeltSpace', '-djpeg', '-r150')
    end
end



% pace - linear
if update_linear || update_all
    [pos_paced_linear, pace_linear] = selfPacer_linear(pos, time);
    pace_accel_linear = deriveByTime(pace_linear, time);
    pace_jerk_linear = deriveByTime(pace_accel_linear, time);
    rms_pos_linear = rms(pos_paced_linear(time_indices_to_analyze));
    rms_accel_linear = rms(pace_accel_linear(time_indices_to_analyze));
    rms_jerk_linear = rms(pace_jerk_linear(time_indices_to_analyze));
end

% pace - ideal (splined)
if update_splined || update_all
    [pos_paced_splined, pace_splined] = selfPacer_spline(pos, time, heelstrike_indices);
    pace_accel_splined = deriveByTime(pace_splined, time);
    pace_jerk_splined = deriveByTime(pace_accel_splined, time);
    rms_pos_splined = rms(pos_paced_splined(time_indices_to_analyze));
    rms_accel_splined = rms(pace_accel_splined(time_indices_to_analyze));
    rms_jerk_splined = rms(pace_jerk_splined(time_indices_to_analyze));

end

% pace - pid
if update_pid || update_all
    [pos_paced_pid, pace_pid] = selfPacer_pid(pos, time, k_p_pid, k_v_pid, v_init);
    pace_accel_pid = deriveByTime(pace_pid, time);
    pace_jerk_pid = deriveByTime(pace_accel_pid, time);
    rms_pos_pid = rms(pos_paced_pid(time_indices_to_analyze));
    rms_accel_pid = rms(pace_accel_pid(time_indices_to_analyze));
    rms_jerk_pid = rms(pace_jerk_pid(time_indices_to_analyze));

end

% pace - pid with position averaged across previous step
if update_pidAvg || update_all
    [pos_paced_pidAvg, pace_pidAvg] = selfPacer_pidWithAveraging(pos, time, k_p_pid, k_v_pid, v_init, heelstrike_indices);
    pace_accel_pidAvg = deriveByTime(pace_pidAvg, time);
    pace_jerk_pidAvg = deriveByTime(pace_accel_pidAvg, time);
    rms_pos_pidAvg = rms(pos_paced_pidAvg(time_indices_to_analyze));
    rms_accel_pidAvg = rms(pace_accel_pidAvg(time_indices_to_analyze));
    rms_jerk_pidAvg = rms(pace_jerk_pidAvg(time_indices_to_analyze));

end

% pace - Seungmoon Song controller
if update_sms || update_all
    [pos_paced_sms, pace_sms] = selfPacer_sms(pos, time, v_init, k_p_sms, k_v_sms, heelstrike_indices);
%     [pos_paced_sms, pace_sms] = selfPacer_sms(pos, time, v_init, 0.1, 0.25, heelstrike_indices);
    pace_accel_sms = deriveByTime(pace_sms, time);
    pace_jerk_sms = deriveByTime(pace_accel_sms, time);
    rms_pos_sms = rms(pos_paced_sms(time_indices_to_analyze));
    rms_accel_sms = rms(pace_accel_sms(time_indices_to_analyze));
    rms_jerk_sms = rms(pace_jerk_sms(time_indices_to_analyze));
end

% pace - minimum jerk
if update_mmj || update_all
    [pos_paced_mmj, pace_mmj] = selfPacer_mmj(pos, time, v_init, heelstrike_indices);
    pace_accel_mmj = deriveByTime(pace_mmj, time);
    pace_jerk_mmj = deriveByTime(pace_accel_mmj, time);
    rms_pos_mmj = rms(pos_paced_mmj(time_indices_to_analyze));
    rms_accel_mmj = rms(pace_accel_mmj(time_indices_to_analyze));
    rms_jerk_mmj = rms(pace_jerk_mmj(time_indices_to_analyze));

end

% pace - pid filtered
if update_pidf || update_all
    [pos_paced_pidf, pace_pidf] = selfPacer_pidWithFilter(pos, time, k_p_pid, k_v_pid, v_init);
    pace_accel_pidf = deriveByTime(pace_pidf, time);
    pace_jerk_pidf = deriveByTime(pace_accel_pidf, time);
    rms_pos_pidf = rms(pos_paced_pidf(time_indices_to_analyze));
    rms_accel_pidf = rms(pace_accel_pidf(time_indices_to_analyze));
    rms_jerk_pidf = rms(pace_jerk_pidf(time_indices_to_analyze));

end

% create figures and axes
colors = lines(8);
colors = colors([3 4], :); % choose these lines to be consistent with earlier plots

axes_to_link = [];
if plot_pos_data
    fig_pos = figure; axes_pos = axes('fontsize', 16); hold on; title('position - world space'); legend('show'); xlabel('time (s)'); ylabel('$p_w (m)$', 'interpreter', 'latex');
%     plot(axes_pos, time, pos_actual, 'displayname', 'actual', 'color', [0.4 0.4 0.4]);
    axes_to_link = [axes_to_link, axes_pos];
end
if plot_vel_data
    fig_vel = figure; axes_vel = axes('fontsize', 16); hold on; title('belt speed'); legend('show'); xlabel('time (s)'); ylabel('$s~(m s^{-1})$', 'interpreter', 'latex');
    plot(axes_vel, time, vel, 'displayname', 'pelvis velocity', 'color', [0.4 0.4 0.4], 'linewidth', 2);
    axes_to_link = [axes_to_link, axes_vel];
end
if plot_acc_data
    fig_acc = figure; axes_acc = axes('fontsize', 16); hold on; title('belt acceleration'); legend('show'); xlabel('time (s)'); ylabel('$\dot s~(m s^{-2})$', 'interpreter', 'latex');
    axes_to_link = [axes_to_link, axes_acc];
end
if plot_jrk_data
    fig_jrk = figure; axes_jrk = axes('fontsize', 16); hold on; title('belt jerk'); legend('show'); xlabel('time (s)'); ylabel('$\ddot s~(m s^{-3})$', 'interpreter', 'latex');
    axes_to_link = [axes_to_link, axes_jrk];
end
if plot_rms_data
    fig_rms = figure;
    axes_rms_pos = subplot(3, 1, 1); hold on; title('RMS pelvis position'); set(axes_rms_pos, 'fontsize', 16);
    axes_rms_acc = subplot(3, 1, 2); hold on; title('RMS belt acceleration'); set(axes_rms_acc, 'fontsize', 16);
    axes_rms_jrk = subplot(3, 1, 3); hold on; title('RMS belt jerk'); set(axes_rms_jrk, 'fontsize', 16);
    bar_labels = {};
    set(axes_rms_pos, 'xtick', []);
    set(axes_rms_acc, 'xtick', []);
    ylabel(axes_rms_pos, '$\mathrm{rms}(p_w)$', 'interpreter', 'latex')
    ylabel(axes_rms_acc, '$\mathrm{rms}(\dot s)$', 'interpreter', 'latex')
    ylabel(axes_rms_jrk, '$\mathrm{rms}(\ddot s)$', 'interpreter', 'latex')
end
count = 0;


linewidth = 5;
if plot_linear
    count = count + 1;

    if plot_pos_data
        plot(axes_pos, time, pos_paced_linear, 'displayname', 'linear', 'color', colors(count, :));
    end
    if plot_vel_data
        plot(axes_vel, time, pace_linear, 'displayname', 'linear', 'color', colors(count, :));
    end
    if plot_acc_data
        plot(axes_acc, time, pace_accel_linear, 'displayname', 'linear', 'color', colors(count, :));
    end
    if plot_jrk_data
        plot(axes_jrk, time, pace_jerk_linear, 'displayname', 'linear', 'color', colors(count, :));
    end    
end

if plot_splined
    count = count + 1;

    if plot_pos_data
        plot(axes_pos, time, pos_paced_splined, 'displayname', 'splined', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_vel_data
        plot(axes_vel, time, pace_splined, 'displayname', 'splined', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_acc_data
        plot(axes_acc, time, pace_accel_splined, 'displayname', 'splined', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_jrk_data
        plot(axes_jrk, time, pace_jerk_splined, 'displayname', 'splined', 'linewidth', linewidth, 'color', colors(count, :));
    end    
    if plot_rms_data
        bar_labels = [bar_labels, 'splined'];
        bar(axes_rms_pos, count, rms_pos_splined, 'facecolor', colors(count, :))
        bar(axes_rms_acc, count, rms_accel_splined, 'facecolor', colors(count, :))
        bar(axes_rms_jrk, count, rms_jerk_splined, 'facecolor', colors(count, :))
    end
end

if plot_pid
    count = count + 1;

    if plot_pos_data
        plot(axes_pos, time, pos_paced_pid, 'displayname', 'PD', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_vel_data
        plot(axes_vel, time, pace_pid, 'displayname', 'PD', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_acc_data
        plot(axes_acc, time, pace_accel_pid, 'displayname', 'PD', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_jrk_data
        plot(axes_jrk, time, pace_jerk_pid, 'displayname', 'PD', 'linewidth', linewidth, 'color', colors(count, :));
    end    
    if plot_rms_data
        bar_labels = [bar_labels, 'PD'];
        bar(axes_rms_pos, count, rms_pos_pid, 'facecolor', colors(count, :))
        bar(axes_rms_acc, count, rms_accel_pid, 'facecolor', colors(count, :))
        bar(axes_rms_jrk, count, rms_jerk_pid, 'facecolor', colors(count, :))
    end    
end

if plot_pidAvg
    count = count + 1;

    if plot_pos_data
        plot(axes_pos, time, pos_paced_pidAvg, 'displayname', 'PD + averaging', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_vel_data
        plot(axes_vel, time, pace_pidAvg, 'displayname', 'PD + averaging', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_acc_data
        plot(axes_acc, time, pace_accel_pidAvg, 'displayname', 'PD + averaging', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_jrk_data
        plot(axes_jrk, time, pace_jerk_pidAvg, 'displayname', 'PD + averaging', 'linewidth', linewidth, 'color', colors(count, :));
    end    
    if plot_rms_data
        bar_labels = [bar_labels, 'PD+avg'];
        bar(axes_rms_pos, count, rms_pos_pidAvg, 'facecolor', colors(count, :))
        bar(axes_rms_acc, count, rms_accel_pidAvg, 'facecolor', colors(count, :))
        bar(axes_rms_jrk, count, rms_jerk_pidAvg, 'facecolor', colors(count, :))
    end    
end

if plot_sms
    count = count + 1;

    if plot_pos_data
        plot(axes_pos, time, pos_paced_sms, 'displayname', 'Song2020', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_vel_data
        plot(axes_vel, time, pace_sms, 'displayname', 'Song2020', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_acc_data
        plot(axes_acc, time, pace_accel_sms, 'displayname', 'Song2020', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_jrk_data
        plot(axes_jrk, time, pace_jerk_sms, 'displayname', 'Song2020', 'linewidth', linewidth, 'color', colors(count, :));
    end    
    if plot_rms_data
        bar_labels = [bar_labels, 'Song2020'];
        bar(axes_rms_pos, count, rms_pos_sms, 'facecolor', colors(count, :))
        bar(axes_rms_acc, count, rms_accel_sms, 'facecolor', colors(count, :))
        bar(axes_rms_jrk, count, rms_jerk_sms, 'facecolor', colors(count, :))
    end    
end

if plot_mmj
    count = count + 1;

    if plot_pos_data
        plot(axes_pos, time, pos_paced_mmj, 'displayname', 'minimum jerk', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_vel_data
        plot(axes_vel, time, pace_mmj, 'displayname', 'minimum jerk', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_acc_data
        plot(axes_acc, time, pace_accel_mmj, 'displayname', 'minimum jerk', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_jrk_data
        plot(axes_jrk, time, pace_jerk_mmj, 'displayname', 'minimum jerk', 'linewidth', linewidth, 'color', colors(count, :));
    end    
    if plot_rms_data
        bar_labels = [bar_labels, 'mmj'];
        bar(axes_rms_pos, count, rms_pos_mmj, 'facecolor', colors(count, :))
        bar(axes_rms_acc, count, rms_accel_mmj, 'facecolor', colors(count, :))
        bar(axes_rms_jrk, count, rms_jerk_mmj, 'facecolor', colors(count, :))
    end
end

if plot_pidf
    count = count + 1;

    if plot_pos_data
        plot(axes_pos, time, pos_paced_pidf, 'displayname', 'PD+filter', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_vel_data
        plot(axes_vel, time, pace_pidf, 'displayname', 'PD+filter', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_acc_data
        plot(axes_acc, time, pace_accel_pidf, 'displayname', 'PD+filter', 'linewidth', linewidth, 'color', colors(count, :));
    end
    if plot_jrk_data
        plot(axes_jrk, time, pace_jerk_pidf, 'displayname', 'PD+filter', 'linewidth', linewidth, 'color', colors(count, :));
    end    
    if plot_rms_data
        bar_labels = [bar_labels, 'PD+filter'];
        bar(axes_rms_pos, count, rms_pos_pidf, 'facecolor', colors(count, :))
        bar(axes_rms_acc, count, rms_accel_pidf, 'facecolor', colors(count, :))
        bar(axes_rms_jrk, count, rms_jerk_pidf, 'facecolor', colors(count, :))
    end
end

if plot_rms_data
    set(axes_rms_jrk, 'xtick', 1:count, 'xticklabels', bar_labels);
end
linkaxes(axes_to_link, 'x')
set(axes_to_link(1), 'xlim', [xlim_1, xlim_2]);

if save_figures
    if plot_pos_data
        set(axes_pos, 'ylim', [-0.4 1.4])
        print(fig_pos, [save_label '_pos'], '-djpeg', '-r150')
    end
    if plot_vel_data
        set(axes_vel, 'ylim', [1 2.5])
        print(fig_vel, [save_label '_vel'], '-djpeg', '-r150')
    end
    if plot_acc_data
        print(fig_acc, [save_label '_acc'], '-djpeg', '-r150')
    end
    if plot_jrk_data
        print(fig_jrk, [save_label '_jrk'], '-djpeg', '-r150')
    end
    if plot_rms_data
        print(fig_rms, [save_label '_rms'], '-djpeg', '-r150')
    end
end















