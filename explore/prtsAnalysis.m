create_stimulus = 1;
analyze_data = 1;

plot_results_1 = 0;
plot_results_2 = 1;

save_figure_1 = 0;
save_figure_2 = 1;

trial_type = 'standing'; type_label = 'standing';
% trial_type = 'stepping'; type_label = 'stepping';
% trial_type = 'stimulus'; type_label = 'walking';

% response_variable = 'com_position'; response_label = 'CoM'; response_unit = 'm'; response_filename_label = 'com'; 
% response_variable = 'trunk_angle'; response_label = 'trunk roll angle'; response_unit = 'deg'; response_filename_label = 'trunkAngle';
% response_variable = 'com_from_cop_position'; response_label = 'CoM - CoP'; response_unit = 'm'; response_filename_label = 'comFromCop';
response_variable = 'com_angle'; response_label = 'CoM angle'; response_unit = 'deg'; response_filename_label = 'comAngle';
% response_variable = 'com_velocity'; response_label = 'CoM vel'; response_unit = 'm/s'; filename = 'comVel';

trials_to_analyze = 1 : 3; trial_number_label = 'fake stimulus'; trial_number_filename_label = 'fake'; use_fake_gains = 1;
trials_to_analyze = 4 : 6; trial_number_label = 'real stimulus'; trial_number_filename_label = 'real'; use_fake_gains = 1;
% trials_to_analyze = 4 : 10; trial_number_label = 'larger stimuli'; trial_number_filename_label = 'largerStimuli'; use_fake_gains = 0;

% trials_to_analyze = 8 : 10; trial_number_label = 'testing'; trial_number_filename_label = 'testing'; use_fake_gains = 0;

filename = [type_label '_' response_filename_label '_' trial_number_filename_label];

stimulus_order = 'pos';
% stimulus_order = 'vel';

% set gain ylimits
dictate_axes = 1;
if strcmp(trial_type, 'standing')
    if strcmp(response_variable, 'com_position')
        gain_ylim = [1e-5 1e-2];
    end
    if strcmp(response_variable, 'trunk_angle')
        gain_ylim = [1e-3 1e-0];
    end
    if strcmp(response_variable, 'com_from_cop_position')
        gain_ylim = [1e-6 1e-2];
    end
    if strcmp(response_variable, 'com_angle')
        gain_ylim = [1e-3 1e-0];
        response_ylim = [-0.1 0.1];
        response_rms_ylim = [0 0.05];
    end
end
if strcmp(trial_type, 'stepping')
    if strcmp(response_variable, 'com_position')
        gain_ylim = [1e-4 1e-1];
    end
    if strcmp(response_variable, 'trunk_angle')
        gain_ylim = [1e-2 1e1];
    end
    if strcmp(response_variable, 'com_from_cop_position')
        gain_ylim = [1e-5 1e-1];
    end
    if strcmp(response_variable, 'com_angle')
        gain_ylim = [1e-2 1e1];
        response_ylim = [-0.5 0.5];
        response_rms_ylim = [0 0.3];
    end
end
if strcmp(trial_type, 'stimulus')
    if strcmp(response_variable, 'com_position')
        gain_ylim = [1e-4 1e-1];
    end
    if strcmp(response_variable, 'trunk_angle')
        gain_ylim = [1e-2 1e1];
    end
    if strcmp(response_variable, 'com_from_cop_position')
        gain_ylim = [1e-5 1e-0];
    end
    if strcmp(response_variable, 'com_angle')
        gain_ylim = [1e-3 1e1];
        response_ylim = [-0.5 0.5];
        response_rms_ylim = [0 0.2];
    end
end



% create stimulus
if create_stimulus
    % set prts parameters from memory -- this is not ideal, should store these instead
    ppdt=30;
    sreg=5;
    
    rate=200;       % sampling rate
    amp = 1;
    fmax=200;     % number of frequency points used in frequency analysis
    SampRate = rate; % sampling rate, this is implicit information
    
%     [p,p1]=pseudogen3(sreg,[0 0 1 -1 -1],[0 0 2 0 2],ppdt);               % symmetrical PRTS
    
    
    [s,x,xi] = pseudogen4(sreg,[0 0 1 -1 -1],[0 0 2 0 2],ppdt);
    x = x * amp;
    tt=(0:length(x)-1)/rate; % define sample time vector
    x_integrated = cumtrapz(tt, x);
    x_integrated_normalized = x_integrated * 1/range(x_integrated);
    x_i_normalized = xi * 1/range(xi);
    
    if strcmp(stimulus_order, 'pos')
        stimulus_base = xi;
    end
    if strcmp(stimulus_order, 'vel')
        stimulus_base = x;
    end
    
    

    ncyc=6;
    nprts=3^sreg-1;
    pts=ppdt*nprts; % number of samples per PRTS cycle
%     ss2=amp*p1/(max(p1)-min(p1)); % amp = p-p position stim
    ss2=amp*stimulus_base/(max(stimulus_base)-min(stimulus_base)); % amp = p-p position stim
    PRTSperiod=pts/rate;
    ss3=ss2;      % extend to specified number of cycles
    for i=1:(ncyc-1)
        ss3=[ss3 ss2];
    end

    ss=[zeros(1,pts/2) ss3 zeros(1,pts/2)];
    stimulus=ss;

    time = (1:max(size(stimulus)))/rate;

    points_per_cycle = pts;
    number_of_cycles = ncyc;


    % figure; hold on
    % plot(time_stim_loaded, stimulus_trajectory_loaded)
    % plot(t, VSstim)
end 


if analyze_data
    number_of_trials = length(trials_to_analyze);
    
    protocol_data = load('protocolInfo.mat');
    gains = zeros(size(trials_to_analyze));
    for i_trial = 1 : number_of_trials
        trial_number = trials_to_analyze(i_trial);
        trial_index = (protocol_data.trial_number == trial_number) & strcmp(protocol_data.trial_type, trial_type);
        this_trial_gain = protocol_data.stimulus_gain(trial_index);
        gains(i_trial) = this_trial_gain;
    end
    
    % fake gains
    if use_fake_gains
        gains = [1 3 6];
    end
    
    % define filter
    filter_order = 6;
    cutoff_frequency = 0.5; % in Hz
    [b_filter, a_filter] = butter(filter_order, cutoff_frequency/(SampRate/2), 'low');
    
    FRFd = cell(1, number_of_trials);
    Cohd = cell(1, number_of_trials);
    avg_stim = cell(1, number_of_trials);
    avg_resp = cell(1, number_of_trials);
    avg_filtered_resp = cell(1, number_of_trials);
    stimulus_all = cell(1, number_of_trials);
    response_all = cell(1, number_of_trials);
    FRF_gains_low = zeros(1, number_of_trials);
    coherence_average_low = zeros(1, number_of_trials);
    avg_filtered_resp_rms = zeros(1, number_of_trials);

    for i_trial = 1 : length(trials_to_analyze)
        subject_settings = loadSettingsFromFile('subject');
        collection_date = subject_settings.get('collection_date');
        subject_id = subject_settings.get('subject_id');
        trial_number = trials_to_analyze(i_trial);
        
        % load CoM
        [com_trajectories, time_com, com_sampling_rate, com_labels, com_directions] = loadData(collection_date, subject_id, trial_type, trial_number, 'com_position_trajectories');
        com_x = com_trajectories(:, strcmp(com_labels, 'center_of_mass_x'));
        com_z = com_trajectories(:, strcmp(com_labels, 'center_of_mass_z'));
        [com_vel_trajectories, time_com_vel, com_vel_sampling_rate, com_vel_labels, com_vel_directions] = loadData(collection_date, subject_id, trial_type, trial_number, 'com_velocity_trajectories');
        com_x_vel = com_vel_trajectories(:, strcmp(com_labels, 'center_of_mass_x'));

        % load CoP
        [forceplate_trajectories, time_forceplate, forceplate_sampling_rate, forceplate_labels, forceplate_directions] = loadData(collection_date, subject_id, trial_type, trial_number, 'forceplate_trajectories');
        cop_x = forceplate_trajectories(:, strcmp(forceplate_labels, 'copx'));
        
        % calculate trunk angle
        [marker_trajectories, time_marker, sampling_rate_marker, marker_labels] = loadData(collection_date, subject_id, trial_type, trial_number, 'marker_trajectories');
        LPSI_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'LPSI', 'trajectories');
        RPSI_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'RPSI', 'trajectories');
        LHEE_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'LHEE', 'trajectories');
        RHEE_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'RHEE', 'trajectories');
        C7_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'C7', 'trajectories');
        MPSI_trajectory = (LPSI_trajectory + RPSI_trajectory) * 0.5;
        foot_basis_trajectory = (LHEE_trajectory + RHEE_trajectory) * 0.5;

        C7_x_trajectory = C7_trajectory(:, 1);
        C7_z_trajectory = C7_trajectory(:, 3);
        MPSI_x_trajectory = MPSI_trajectory(:, 1);
        MPSI_z_trajectory = MPSI_trajectory(:, 3);
        foot_basis_x_trajectory = foot_basis_trajectory(:, 1);
        foot_basis_x_trajectory_splined = spline(time_marker, foot_basis_x_trajectory, time_com);

        trunk_vector_x = C7_x_trajectory - MPSI_x_trajectory;
        trunk_vector_z = C7_z_trajectory - MPSI_z_trajectory;
        trunk_angle_trajectory = rad2deg(atan2(trunk_vector_x, trunk_vector_z));

        
        
        
        com_vector_x = com_x - foot_basis_x_trajectory_splined;
        com_vector_z = com_z;
        com_angle_trajectory = rad2deg(atan2(com_vector_x, com_vector_z));

        [stimulus_trajectory_loaded, time_stim_loaded, stim_sampling_rate, stim_labels, stim_directions] = loadData(collection_date, subject_id, trial_type, trial_number, 'visual_rotation_angle_trajectory');

        stim_imported = importdata(['labview' filesep 'prts.csv']);
        time_stim_new = stim_imported.data(:, 1);
        stimulus_trajectory_new = stim_imported.data(:, 2);

    %     figure; hold on
    %     plot(time_stim_loaded, -stimulus_trajectory_loaded * 0.01)
    %     plot(time_stim_new, -stimulus_trajectory_new * 0.15)
    %     plot(time_com, com_x)
    %     plot(time_marker, trunk_angle_trajectory)



        % analyze
        if strcmp(response_variable, 'com_position')
            time_response = time_com;
            response = - com_x;
            response = spline(time_response, response, time);
        end
        if strcmp(response_variable, 'com_velocity')
            time_response = time_com_vel;
            response = - com_x_vel;
            response = spline(time_response, response, time);
        end
        if strcmp(response_variable, 'trunk_angle')
            time_response = time_marker;
            response = - trunk_angle_trajectory;
            response = spline(time_response, response, time);
        end
        if strcmp(response_variable, 'com_from_cop_position')
            com_splined = spline(time_com, com_x, time);
            cop_splined = spline(time_forceplate, cop_x, time);
            
            response = -(com_splined - cop_splined);
        end
        if strcmp(response_variable, 'com_angle')
            time_response = time_com;
            response = - com_angle_trajectory;
            response = spline(time_response, response, time);
        end
        
        this_gain = gains(i_trial);
        this_stimulus = stimulus * this_gain;
        
        response=detrend(response);
        

        response_filtered = filtfilt(b_filter, a_filter, response);

        % Frequency domain variables for accumulating results
        time_single_cycle = (0:(points_per_cycle-1))/SampRate; % define sample time vector
        f=(SampRate/points_per_cycle)*(1:fmax);   % frequency vector for frequency analysis

        yi_trial = zeros(1,fmax);       % stimulus dft  - fmax freqs
        yo_trial = zeros(1,fmax);		% response dft - fmax freqs
        yoi_trial = zeros(1,fmax);      % used for coherence calculations
        yii_trial = zeros(1,fmax);
        yoo_trial = zeros(1,fmax);
                                %
        % Time domain variables for accumulating results
        avg_stim_trial = zeros(1,points_per_cycle);
        avg_resp_trial = zeros(1,points_per_cycle);
        avg_filtered_resp_trial = zeros(1,points_per_cycle);
        psfactor_trial = 1/(2*SampRate*points_per_cycle);      % Factor for scaling power spectra such that integration across
                                                % all frequencies (deltaf*sum(power spectral values) gives the 
                                                % mean square value. The power spectra are calculated beginning
                                                % with 2*fft (2*fft gives the one sided spectrum).
        cyc = 0;      % cycle counter
        startidx = 1; % starting index of first stimulus cycle
        startidx = 3570; % starting index of first stimulus cycle
        startidx = points_per_cycle/2;
        for k = 2 : number_of_cycles          % skip 1st cycles to avoid transient response
            cyc=cyc+1;
            yi_raw=fft(this_stimulus((startidx+(k-1)*points_per_cycle):(startidx+k*points_per_cycle-1))');	    % stimulus dft
            yi_raw=conj(yi_raw);
            yo_raw=fft(response((startidx+(k-1)*points_per_cycle):(startidx+k*points_per_cycle-1))');     % response dft
            yo_raw=conj(yo_raw);

            yi_raw2=2*yi_raw(2:(fmax+1))';       % ignore DC values and multiply by 2 for one-sided dft spectra
            yo_raw2=2*yo_raw(2:(fmax+1))';

            yi_trial(cyc,:)=yi_raw2;                 % accumulate dft spectra - fmax points
            yo_trial(cyc,:)=yo_raw2;

            yoi_trial(cyc,:)=psfactor_trial*yo_raw2.*conj(yi_raw2);     % accumulate scaled power spectra and cross power spectra
            yii_trial(cyc,:)=psfactor_trial*abs(yi_raw2).*abs(yi_raw2);
            yoo_trial(cyc,:)=psfactor_trial*abs(yo_raw2).*abs(yo_raw2);

            avg_stim_trial=avg_stim_trial+this_stimulus((startidx+(k-1)*points_per_cycle):(startidx+k*points_per_cycle-1));  % time averages
            avg_resp_trial=avg_resp_trial+response((startidx+(k-1)*points_per_cycle):(startidx+k*points_per_cycle-1));                                          % resp time averages
            avg_filtered_resp_trial=avg_filtered_resp_trial+response_filtered((startidx+(k-1)*points_per_cycle):(startidx+k*points_per_cycle-1));                                          % resp time averages


        end	%for k=2:ncyc
        avg_stim_trial=avg_stim_trial/cyc;	% average stimulus
        avg_resp_trial=avg_resp_trial/cyc;	% average response
        avg_filtered_resp_trial=avg_filtered_resp_trial/cyc;	% average response

        yi_mean_trial = mean(yi_trial,1);
        yo_mean_trial = mean(yo_trial,1);
        yoi_mean_trial = mean(yoi_trial,1);
        yii_mean_trial = mean(yii_trial,1);
        yoo_mean_trial = mean(yoo_trial,1);

        yi_var_trial = zeros(1,fmax);  % calculate variance of spectra - Pintelon & Schoukens eq 2-31
        yo_var_trial = zeros(1,fmax);
        yoi_var_trial = zeros(1,fmax);
        for k=1:cyc
            yi_var_trial = yi_var_trial+abs(yi_trial(cyc,:)-yi_mean_trial).^2;
            yo_var_trial = yo_var_trial+abs(yo_trial(cyc,:)-yo_mean_trial).^2;
            yoi_var_trial = yoi_var_trial+(yo_trial(cyc,:)-yo_mean_trial).*conj(yi_trial(cyc,:)-yi_mean_trial);
        end
        yi_var_trial = psfactor_trial*yi_var_trial/(cyc-1);
        yo_var_trial = psfactor_trial*yo_var_trial/(cyc-1);
        yoi_var_trial = psfactor_trial*yoi_var_trial/(cyc-1);

        % Calculate variance of FRF (frequency response function) - Pintelon & Schoukens eq 2-32
        %   Calculate at all frequencies, then decimate to get final values at
        %   PRTS frequencies that have stimulus energy.
        FRF_trial=yo_mean_trial./yi_mean_trial;   % FRF calculation at all frequencies - Pintelon & Schoukens eq 2-30

        FRF_var_trial=(abs(FRF_trial).^2).*(yo_var_trial./(psfactor_trial*abs(yo_mean_trial).^2)+yi_var_trial./(psfactor_trial*abs(yi_mean_trial).^2)+2*real(yoi_var_trial./(psfactor_trial*yo_mean_trial.*conj(yi_mean_trial))));
        FRFd_var_trial=decimate2(FRF_var_trial,2);

        FRFd_trial=decimate2(FRF_trial,2);  % FRF at frequencies where stimulus had energy

        fd=decimate2(f,2);      % frequencies with stimulus PRTS energy
        yid_mean_trial=decimate2(yi_mean_trial,2);
        yod_mean_trial=decimate2(yo_mean_trial,2);

        yoid_mean_trial=decimate2(yoi_mean_trial,2);
        yiid_mean_trial=decimate2(yii_mean_trial,2);
        yood_mean_trial=decimate2(yoo_mean_trial,2);
        
        % smooth
%         fd = SiPsmooth(fd);
%         FRFd_trial = SiPsmooth(FRFd_trial);
%         yoid_mean_trial = SiPsmooth(yoid_mean_trial);
%         yiid_mean_trial = SiPsmooth(yiid_mean_trial);
%         yood_mean_trial = SiPsmooth(yood_mean_trial);
        
        fd = SiPsmooth_more1(fd);
        FRFd_trial = SiPsmooth_more1(FRFd_trial);
        yoid_mean_trial = SiPsmooth_more1(yoid_mean_trial);
        yiid_mean_trial = SiPsmooth_more1(yiid_mean_trial);
        yood_mean_trial = SiPsmooth_more1(yood_mean_trial);
        
%         fd = SiPsmooth_more2(fd);
%         FRFd_trial = SiPsmooth_more2(FRFd_trial);
%         yoid_mean_trial = SiPsmooth_more2(yoid_mean_trial);
%         yiid_mean_trial = SiPsmooth_more2(yiid_mean_trial);
%         yood_mean_trial = SiPsmooth_more2(yood_mean_trial);
        

        Cohd_trial=(abs(yoid_mean_trial).^2)./(yiid_mean_trial.*yood_mean_trial); % Coherence function calculation
        
        
        
        FRF_gains_low(i_trial) = mean(abs(FRFd_trial(1:2)));
        coherence_average_low(i_trial) = mean(Cohd_trial(1:2));
        avg_filtered_resp_rms(i_trial) = rms(avg_filtered_resp_trial);
        

        fdd=decimate2(f(2:fmax),2);		% frequencies without stimulus PRTS energy
        yipdd_trial=decimate2(yi_mean_trial(2:fmax),2);
        yopdd_trial=decimate2(yo_mean_trial(2:fmax),2);

        FRFd{i_trial} = FRFd_trial;
        Cohd{i_trial} = Cohd_trial;
        avg_stim{i_trial} = avg_stim_trial;
        avg_resp{i_trial} = avg_resp_trial;
        avg_filtered_resp{i_trial} = avg_filtered_resp_trial;
        response_all{i_trial} = response;
        stimulus_all{i_trial} = this_stimulus;
    end
end


colors = copper(length(trials_to_analyze));
[~, plot_order] = sort(gains);
if plot_results_1
    figure('position', [0 0 1200 1800])
    tiledlayout(3,2,'TileSpacing','Compact', 'padding', 'compact');
    
    axes_321 = nexttile;
    hold(axes_321, 'on')
    set(axes_321, 'XScale', 'log', 'YScale', 'log', 'fontsize', 12)
    xlim([0.02 4])
    if dictate_axes
        ylim(gain_ylim)
    end
    ylabel(['FRF gain (' response_unit '/deg)'])

%     axes_322 = nexttile;
%     hold(axes_322, 'on')
%     set(axes_322, 'fontsize', 12)
%     axis([0 max(time_single_cycle) -max(avg_stim_trial)*1.2 max(avg_stim_trial)*1.2])
%     ylabel('Average stimulus')
%     legend('show', 'Location', 'northeast')

    axes_322 = nexttile;
    hold(axes_322, 'on')
    set(axes_322, 'fontsize', 12)
    xlim([0.5 max(gains)+0.5])
    ylabel('Gain average below 0.1 Hz')
    legend('show', 'Location', 'northeast')

    axes_323 = nexttile;
    hold(axes_323, 'on')
    set(axes_323, 'XScale', 'log', 'fontsize', 12)
    xlim([0.02 4])
    ylabel('FRF phase')

    axes_324 = nexttile;
    hold(axes_324, 'on')
    set(axes_324, 'fontsize', 12)
    xlim([0 max(time_single_cycle)])
    ylabel(['Average ' response_label ' response (' response_unit ')'])

    axes_325 = nexttile;
    hold(axes_325, 'on')
    set(axes_325, 'XScale', 'log', 'fontsize', 12)
    xlim([0.02 4]); ylim([0 1])
    ylabel('Coherence')
    xlabel('Frequency (Hz)')

%     axes_326 = nexttile;
%     hold(axes_326, 'on')
%     set(axes_326, 'fontsize', 12)
%     xlim([0 max(time)])
%     ylabel('Response over time')
%     xlabel('Time (s)')

    axes_326 = nexttile;
    hold(axes_326, 'on')
    set(axes_326, 'fontsize', 12)
    xlim([0.5 max(gains)+0.5])
    ylim([0 1])
    ylabel('Coherence average below 0.1 Hz')

    linewidth = 1;
%     plot(axes_322, gains, FRF_gains_low, ':', 'color', [1 1 1]*0.7, 'linewidth', 3)
    for i_trial = 1 : number_of_trials
        this_trial_index = plot_order(i_trial);
        if gains(this_trial_index) ~= 0
            plot(axes_321, fd, abs(FRFd{this_trial_index}),'-o', 'color', colors(i_trial, :), 'linewidth', linewidth);

            plot(axes_323, fd, 180/pi*phase(FRFd{this_trial_index}),'-o', 'color', colors(i_trial, :), 'linewidth', linewidth);

            plot(axes_325, fd, Cohd{this_trial_index},'-o', 'color', colors(i_trial, :), 'linewidth', linewidth);

%             plot(axes_322, time_single_cycle, avg_stim{this_trial_index}, 'displayname', num2str(gains(this_trial_index)), 'color', colors(i_trial, :), 'linewidth', linewidth);
            plot(axes_322, gains(i_trial), FRF_gains_low(i_trial), 'o', 'displayname', num2str(gains(this_trial_index)), 'color', colors(i_trial, :), 'markersize', 18, 'MarkerFaceColor', colors(i_trial, :));

            plot(axes_324, time_single_cycle, avg_resp{this_trial_index}, 'color', colors(i_trial, :), 'linewidth', linewidth);

%             plot(axes_326, time, response_all{this_trial_index}, 'color', colors(i_trial, :), 'linewidth', linewidth);
%             plot(axes_326, time, stimulus_all{this_trial_index}, ':', 'color', colors(i_trial, :), 'linewidth', linewidth);
            plot(axes_326, gains(i_trial), coherence_average_low(i_trial), 'o', 'displayname', num2str(gains(this_trial_index)), 'color', colors(i_trial, :), 'markersize', 18, 'MarkerFaceColor', colors(i_trial, :));
        end
        
    end
    
    
    this_title = [response_label ' - ' type_label ' - ' trial_number_label];
    uicontrol('style', 'text', 'string', this_title, 'units', 'normalized', 'position', [0, 0.97, 1, 0.03], 'fontsize', 16, 'FontWeight', 'bold')
    
    % resize and print to pdf
    if save_figure_1
        saveas(gcf, ['figures' filesep 'fig' filesep filename])
        set(gcf, 'Position', [100 100 600 800])
        saveas(gcf, ['figures' filesep 'pdf' filesep filename], 'pdf')
    end    
%     print('test.pdf', '-bestfit', '-dpdf')
    
end

if plot_results_2
    figure('position', [0 0 600 900])
    tiledlayout(3, 1, 'TileSpacing', 'Compact', 'padding', 'compact');
    
    axes_21 = nexttile;
    hold(axes_21, 'on')
    set(axes_21, 'fontsize', 12)
    xlabel('Time (s)')
    ylabel('Stimulus')
    xlim([0, max(time_single_cycle)])

    axes_22 = nexttile;
    hold(axes_22, 'on')
    set(axes_22, 'fontsize', 12)
    xlim([0, max(time_single_cycle)])
    if dictate_axes
        ylim(response_ylim)
    end
    xlabel('Time (s)')
    ylabel('Average response, low-pass filtered')

    axes_23 = nexttile;
    hold(axes_23, 'on')
    set(axes_23, 'fontsize', 12)
%     ylabel('Average stimulus')
    ylabel('Average response, low-pass filtered, RMS')
%     legend('show', 'Location', 'northwest')
    if dictate_axes
        ylim(response_rms_ylim)
    end
    xlabel('Stimulus peak-to-peak amplitude (deg)')
    xlim([0.5 max(gains)+0.5])
    legend('show', 'Location', 'NW')
    
    linewidth = 1;
%     plot(axes_322, gains, FRF_gains_low, ':', 'color', [1 1 1]*0.7, 'linewidth', 3)
    for i_trial = 1 : number_of_trials
        this_trial_index = plot_order(i_trial);
        if gains(this_trial_index) ~= 0

            plot(axes_21, time_single_cycle, avg_stim{this_trial_index}, 'displayname', num2str(gains(this_trial_index)), 'color', colors(i_trial, :), 'linewidth', linewidth);

            plot(axes_22, time_single_cycle, avg_filtered_resp{this_trial_index}, 'color', colors(i_trial, :), 'linewidth', linewidth);

            plot(axes_23, gains(i_trial), avg_filtered_resp_rms(i_trial), 'o', 'displayname', num2str(gains(this_trial_index)), 'color', colors(i_trial, :), 'markersize', 18, 'MarkerFaceColor', colors(i_trial, :));
        end
        
    end
    
    
    this_title = [response_label ' - ' type_label ' - ' trial_number_label];
    uicontrol('style', 'text', 'string', this_title, 'units', 'normalized', 'position', [0, 0.97, 1, 0.03], 'fontsize', 16, 'FontWeight', 'bold')
    
    % resize and print to pdf
    if save_figure_2
        saveas(gcf, ['figures' filesep 'fig' filesep filename '_filtered'])
        set(gcf, 'Position', [100 100 600 800])
        saveas(gcf, ['figures' filesep 'pdf' filesep filename '_filtered'], 'pdf')
    end    
%     print('test.pdf', '-bestfit', '-dpdf')
    
end


