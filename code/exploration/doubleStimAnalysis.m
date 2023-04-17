% flags
analyze_data            = 1;
plot_results            = 0;
dictate_axes            = 0;
save_figure             = 0;

export_data             = 0;
export_results          = 1;
show_metronome_figure   = 0;

trial_type = 'stimulus'; type_label = 'walking';

% select response variable
response_variable = 'com_position'; response_label = 'CoM'; response_unit = 'm'; response_filename_label = 'com'; fit_model = 0;
% response_variable = 'com_angle'; response_label = 'CoM angle'; response_unit = 'deg'; response_filename_label = 'comAngle'; fit_model = 1;
% response_variable = 'foot_base'; response_label = 'foot base'; response_unit = 'm'; response_filename_label = 'footBase'; fit_model = 0;

% select trials to analyze
trials_to_analyze = 1:6; trial_label = 'low cadence'; trial_filename_label = 'lowCadence';
% trials_to_analyze = 7:12; trial_label = 'high cadence'; trial_filename_label = 'highCadence';

trials_to_analyze = 1:12; trial_label = 'all cadences'; trial_filename_label = 'allCadences';

% create filename and title label
filename = [response_filename_label '_' trial_filename_label];
title_label = [response_label '_' trial_label];

% set parameters for analysis
fmax                            = 200;  % number of frequency points used in frequency analysis
number_of_cycles_to_ignore_end  = 0;    % number of stimulus cycles ignored at the end

subject_settings = loadSettingsFromFile('subject');
subject_id = subject_settings.get('subject_id');

% analyze data
if analyze_data
    stimulus_data = load(['labview' filesep 'stimulus.mat']);
    number_of_trials = length(trials_to_analyze);

    protocol_data = load('protocolInfo.mat');
    vis_stimulus_amplitudes = zeros(size(trials_to_analyze));
    paces = zeros(size(trials_to_analyze));
    gvs_stimulus_amplitudes = zeros(size(trials_to_analyze));
    metronome_files = cell(number_of_trials, 1);
    vis_stimulus_files = cell(number_of_trials, 1);
    gvs_stimulus_files = cell(number_of_trials, 1);
    for i_trial = 1 : number_of_trials
        trial_number = trials_to_analyze(i_trial);
        trial_index = (protocol_data.trial_number == trial_number) & strcmp(protocol_data.trial_type, trial_type);
        
        % get pace
        metronome_files{i_trial} = protocol_data.metronome_file_name{trial_index};%this_trial_metronome_filename{1};
        this_trial_metronome_filename = metronome_files{i_trial};
        this_trial_metronome_strsplit = strsplit(this_trial_metronome_filename, '_');
        this_trial_metronome = str2double(this_trial_metronome_strsplit{2});
        metronome_frequency = this_trial_metronome / stimulus_data.single_cycle_duration;
        paces(i_trial) = metronome_frequency;
        
        % get visual stimulus amplitude
        vis_stimulus_files{i_trial} = protocol_data.UseVisualStimulus_FromFile{trial_index};
        this_trial_vis_stimulus_filename = vis_stimulus_files{i_trial};
        this_trial_vis_stimulus_strsplit = strsplit(this_trial_vis_stimulus_filename, '_');
        this_trial_vis_stimulus_strsplit = strsplit(this_trial_vis_stimulus_strsplit{2}, '.');
        this_trial_vis_stimulus_amplitude = str2double(this_trial_vis_stimulus_strsplit{1});
        vis_stimulus_amplitudes(i_trial) = this_trial_vis_stimulus_amplitude;
        
        % get gvs stimulus
        gvs_stimulus_files{i_trial} = protocol_data.UseGVSStimulus_FromFile{trial_index};
        this_trial_gvs_stimulus_filename = gvs_stimulus_files{i_trial};
        if isempty(this_trial_gvs_stimulus_filename)
            gvs_stimulus_amplitudes(i_trial) = 0;
        else
            gvs_stimulus_amplitudes(i_trial) = 0.25;
        end
        
    end

    % define filter
    signal_data = load(['labview' filesep 'stimulus.mat']);
    filter_order = 6;
    cutoff_frequency = 0.5; % in Hz
    [b_filter, a_filter] = butter(filter_order, cutoff_frequency/(signal_data.sampling_rate/2), 'low');
    
    % create containers for results
    FRF_decimated = cell(1, number_of_trials);
    Coh_decimated = cell(1, number_of_trials);
    average_vis_stimulus = cell(1, number_of_trials);
    average_gvs_stimulus = cell(1, number_of_trials);
    average_response = cell(1, number_of_trials);
    average_filtered_response = cell(1, number_of_trials);
    average_pulse_response = cell(1, number_of_trials);
    average_pulse_filtered_response = cell(1, number_of_trials);
    stimulus_all = cell(1, number_of_trials);
    response_all = cell(1, number_of_trials);
    FRF_gains_low = zeros(1, number_of_trials);
    coherence_average_low = zeros(1, number_of_trials);
    average_filtered_response_rms = zeros(1, number_of_trials);
    model_fits = cell(1, number_of_trials);
    
    % process
    subject_settings = loadSettingsFromFile('subjectSettings.txt');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
    for i_trial = 1 : length(trials_to_analyze)
        % load data
        trial_number = trials_to_analyze(i_trial);
        
        metronome_frequency = paces(i_trial);
        metronome_data = importdata(['labview' filesep metronome_files{i_trial}]);
        metronome_reset_indices = find(diff(metronome_data)) + 1;
        metronome_phase = wrapToPi(2 * pi * signal_data.time * metronome_frequency);
        
        if export_data
            data_for_export = struct;
            data_for_export.metronome_frequency = metronome_frequency;
            data_for_export.units.metronome_frequency = 'Hz';
        end
        
        % load heel-strikes and calculate phase
        events_file_name = makeFileName(collection_date, subject_id, trial_type, trial_number, 'events.mat');
        events_data = load(['analysis' filesep events_file_name]);
        left_heel_strikes = events_data.event_data{strcmp(events_data.event_labels, 'left_touchdown')};
        right_heel_strikes = events_data.event_data{strcmp(events_data.event_labels, 'right_touchdown')};
        left_heel_strike_indices = findClosestIndex(left_heel_strikes, signal_data.time);
        right_heel_strike_indices = findClosestIndex(right_heel_strikes, signal_data.time);
        left_heel_strike_phase = metronome_phase(left_heel_strike_indices);
        right_heel_strike_phase = metronome_phase(right_heel_strike_indices);
        
        if show_metronome_figure
            figure; hold on;
            plot(signal_data.time(left_heel_strike_indices), left_heel_strike_phase, 'x')
            plot(signal_data.time(right_heel_strike_indices), right_heel_strike_phase, 'x')
        end
        
        % load CoM
        if strcmp(response_variable, 'com_position') || strcmp(response_variable, 'com_angle') || export_data
            [com_trajectories, time_com, com_sampling_rate, com_labels, com_directions] = loadData(collection_date, subject_id, trial_type, trial_number, 'com_position_trajectories');
            com_x = com_trajectories(:, strcmp(com_labels, 'center_of_mass_x'));
            com_z = com_trajectories(:, strcmp(com_labels, 'center_of_mass_z'));
        end
        
        % calculate angles
        if strcmp(response_variable, 'trunk_angle') || strcmp(response_variable, 'com_angle') || strcmp(response_variable, 'foot_base') || export_data
            [marker_trajectories, time_marker, sampling_rate_marker, marker_labels] = loadData(collection_date, subject_id, trial_type, trial_number, 'marker_trajectories');
            LPSI_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'LPSI', 'trajectories');
            RPSI_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'RPSI', 'trajectories');
            LHEE_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'LHEE', 'trajectories');
            RHEE_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'RHEE', 'trajectories');
            C7_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'C7', 'trajectories');
            MPSI_trajectory = (LPSI_trajectory + RPSI_trajectory) * 0.5;
            
            events_file_name = makeFileName(collection_date, subject_id, trial_type, trial_number, 'events.mat');
            events_data = load(['analysis' filesep events_file_name]);
            
            foot_base_x_trajectory = estimateFootBasisTrajectory_boxcar ...
              ( ...
                time_marker, ...
                marker_trajectories, ...
                marker_labels, ...
                events_data, ...
                sampling_rate_marker, ...
                metronome_frequency ...
              );
          
            C7_x_trajectory = C7_trajectory(:, 1);
            C7_z_trajectory = C7_trajectory(:, 3);
            MPSI_x_trajectory = MPSI_trajectory(:, 1);
            MPSI_z_trajectory = MPSI_trajectory(:, 3);
            
            trunk_vector_x = C7_x_trajectory - MPSI_x_trajectory;
            trunk_vector_z = C7_z_trajectory - MPSI_z_trajectory;
            trunk_angle_trajectory = rad2deg(atan2(trunk_vector_x, trunk_vector_z));

            if strcmp(response_variable, 'com_angle')
                com_x_trajectory_splined = spline(time_com, com_x, time_marker);
                com_z_trajectory_splined = spline(time_com, com_z, time_marker);
                com_vector_x = com_x_trajectory_splined - foot_base_x_trajectory;
                com_vector_z = com_z_trajectory_splined;
                com_angle_trajectory = rad2deg(atan2(com_vector_x, com_vector_z));
            end            
            
        end        
        
        % identify response
        time = signal_data.time;
        if strcmp(response_variable, 'com_position')
            time_response = time_com;
            time(time>time_response(end)) = [];
            response = - com_x;
            response = spline(time_response, response, time);
        end
        if strcmp(response_variable, 'trunk_angle')
            time_response = time_marker;
            response = - trunk_angle_trajectory;
            time(time>time_response(end)) = [];
            response = spline(time_response, response, time);
        end
        if strcmp(response_variable, 'com_angle')
            time_response = time_marker;
            time(time>time_response(end)) = [];
            response = - com_angle_trajectory;
            response = spline(time_response, response, time);
        end
        if strcmp(response_variable, 'foot_base')
            time_response = time_marker;
            time(time>time_response(end)) = [];
            response = - foot_base_x_trajectory;
            response = spline(time_response, response, time);
        end
        if export_data
            % figure out time
            time_max = min([time_com(end), time_marker(end)]);
            data_for_export.time = time;
            data_for_export.time(data_for_export.time > time_max) = [];
            
            % com position
            data_for_export.com_x = spline(time_com, - com_x, data_for_export.time);
            data_for_export.com_z = spline(time_com, com_z, data_for_export.time);
            data_for_export.directions.com_x.positive = 'left';
            data_for_export.directions.com_x.negative = 'right';
            data_for_export.directions.com_z.positive = 'up';
            data_for_export.directions.com_z.negative = 'down';
            data_for_export.units.com_x = 'm';
            data_for_export.units.com_z = 'm';
            
            % com angle
            data_for_export.com_angle = spline(time_marker, - com_angle_trajectory, data_for_export.time);
            data_for_export.directions.com_angle.positive = 'left roll';
            data_for_export.directions.com_angle.negative = 'right roll';
            data_for_export.units.com_angle = 'deg';
            
            % foot base
            data_for_export.foot_base_x = spline(time_marker, - foot_base_x_trajectory, data_for_export.time);
            data_for_export.directions.foot_base_x.positive = 'left';
            data_for_export.directions.foot_base_x.negative = 'right';
            data_for_export.units.foot_base_x = 'm';
        end
        
        % clean up
        this_pace = paces(i_trial);
        this_vis_stim_amplitude = vis_stimulus_amplitudes(i_trial);
        this_vis_stimulus = signal_data.multisine_stimulus_normalized * this_vis_stim_amplitude;
        this_vis_stimulus(signal_data.time>time_response(end)) = [];
        
        this_gvs_stim_amplitude = gvs_stimulus_amplitudes(i_trial);
        this_gvs_stimulus = - signal_data.pulse_stimulus_normalized * this_gvs_stim_amplitude; % invert to make responses look nice, direction is arbitrary, depending on anode/cathode placement
        this_gvs_stimulus(signal_data.time>time_response(end)) = [];
        
        if export_data
            data_for_export.stimulus = this_vis_stimulus;
            data_for_export.units.stimulus = 'deg';
        end        
        
        response = detrend(response);
        response_filtered = filtfilt(b_filter, a_filter, response);
        
        % determine cycles to analyze
        cycles_to_analyze = signal_data.number_of_warmup_cycles + (1 : signal_data.number_of_stimulus_cycles);
        cycles_to_ignore_end = (length(cycles_to_analyze) + 1 - number_of_cycles_to_ignore_end) : length(cycles_to_analyze);
        cycles_to_analyze(cycles_to_ignore_end) = [];
        number_of_cycles_to_analyze = length(cycles_to_analyze);
        
        % prepare frequency domain analysis
        time_single_cycle = signal_data.single_cycle_time;
        f = (signal_data.sampling_rate/signal_data.points_per_cycle)*(1:fmax);   % frequency vector for frequency analysis
        yi_trial = zeros(number_of_cycles_to_analyze, fmax);    % stimulus dft  - fmax freqs
        yo_trial = zeros(number_of_cycles_to_analyze, fmax);    % response dft - fmax freqs
        yoi_trial = zeros(number_of_cycles_to_analyze, fmax);   % used for coherence calculations
        yii_trial = zeros(number_of_cycles_to_analyze, fmax);
        yoo_trial = zeros(number_of_cycles_to_analyze, fmax);
                                %
        % Time domain variables for accumulating results
        average_vis_stimulus_trial = zeros(1, signal_data.points_per_cycle);
        average_gvs_stimulus_trial = zeros(1, signal_data.points_per_cycle);
        average_response_trial = zeros(1, signal_data.points_per_cycle);
        average_pulse_response_trial = zeros(1, signal_data.points_per_pulse);
        average_filtered_response_trial = zeros(1, signal_data.points_per_cycle);
        psfactor_trial = 1/(2*signal_data.sampling_rate*signal_data.points_per_cycle);      % Factor for scaling power spectra such that integration across
                                                % all frequencies (deltaf*sum(power spectral values) gives the 
                                                % mean square value. The power spectra are calculated beginning
                                                % with 2*fft (2*fft gives the one sided spectrum).

        
        % analyze cycles
        ppc = signal_data.points_per_cycle;
        for i_cycle = 1 : number_of_cycles_to_analyze
            this_cycle_index = cycles_to_analyze(i_cycle);
            this_cycle_start_index = 1 + ppc * (this_cycle_index-1);
            this_cycle_end_index = ppc * this_cycle_index;
            this_cycle_indices = this_cycle_start_index : this_cycle_end_index;
            yi_raw = fft(this_vis_stimulus(this_cycle_indices)');	    % stimulus dft
            yi_raw = conj(yi_raw);
            yo_raw = fft(response(this_cycle_indices)');            % response dft
            yo_raw = conj(yo_raw);

            % ignore DC values and multiply by 2 for one-sided dft spectra
            yi_raw2 = 2*yi_raw(2:(fmax+1))';
            yo_raw2 = 2*yo_raw(2:(fmax+1))';

            % accumulate dft spectra - fmax points
            yi_trial(i_cycle, :) = yi_raw2;
            yo_trial(i_cycle, :) = yo_raw2;

            % accumulate scaled power spectra and cross power spectra
            yoi_trial(i_cycle, :) = psfactor_trial*yo_raw2.*conj(yi_raw2);
            yii_trial(i_cycle, :) = psfactor_trial*abs(yi_raw2).*abs(yi_raw2);
            yoo_trial(i_cycle, :) = psfactor_trial*abs(yo_raw2).*abs(yo_raw2);

            % add up time domain variables
            average_vis_stimulus_trial = average_vis_stimulus_trial + this_vis_stimulus(this_cycle_indices);
            average_gvs_stimulus_trial = average_gvs_stimulus_trial + this_gvs_stimulus(this_cycle_indices);
            average_response_trial = average_response_trial + response(this_cycle_indices);
            average_filtered_response_trial = average_filtered_response_trial + response_filtered(this_cycle_indices);
        end

        % calculate averages
        average_vis_stimulus_trial = average_vis_stimulus_trial/number_of_cycles_to_analyze;	% average stimulus
        average_gvs_stimulus_trial = average_gvs_stimulus_trial/number_of_cycles_to_analyze;	% average stimulus
        average_response_trial = average_response_trial/number_of_cycles_to_analyze;	% average response
        average_filtered_response_trial = average_filtered_response_trial/number_of_cycles_to_analyze;	% average filtered response
        
        % pulse average
        for i_pulse = 1 : signal_data.number_of_pulses_per_cycle
            this_pulse_start_index = 1 + signal_data.points_per_pulse * (i_pulse-1);
            this_pulse_end_index = signal_data.points_per_pulse * i_pulse;
            this_pulse_indices = this_pulse_start_index : this_pulse_end_index;
            this_pulse_response = average_response_trial(this_pulse_indices);
            this_pulse_filtered_response = average_filtered_response_trial(this_pulse_indices);
            this_pulse_stimulus = average_gvs_stimulus_trial(this_pulse_indices);
            % is this pulse positive or negative?
            this_sign = (-1)^i_pulse;
            
%             [~, extremum_index] = max(abs(this_pulse_stimulus));
%             this_sign = sign(this_pulse_stimulus(extremum_index));
            
            % store
            average_pulse_response_trial = average_pulse_response_trial + this_pulse_response * this_sign;
            average_pulse_filtered_response_trial = average_pulse_response_trial + this_pulse_filtered_response * this_sign;
        end
        average_pulse_response_trial = average_pulse_response_trial/signal_data.number_of_pulses_per_cycle;	% average response
        average_pulse_filtered_response_trial = average_pulse_filtered_response_trial/signal_data.number_of_pulses_per_cycle;	% average response
        time_single_pulse = signal_data.single_pulse_time;
        
        yi_mean_trial = mean(yi_trial, 1);
        yo_mean_trial = mean(yo_trial, 1);
        yoi_mean_trial = mean(yoi_trial, 1);
        yii_mean_trial = mean(yii_trial, 1);
        yoo_mean_trial = mean(yoo_trial, 1);
        
        % calculate variance of spectra - Pintelon & Schoukens eq 2-31
        yi_var_trial = zeros(1, fmax);  
        yo_var_trial = zeros(1, fmax);
        yoi_var_trial = zeros(1, fmax);
        % TODO: check this against the mentioned equations
        % the iterator isn't used in this loop, is this correct?
        for i_cycle = 1 : number_of_cycles_to_analyze
            yi_var_trial = yi_var_trial + abs(yi_trial(number_of_cycles_to_analyze, :) - yi_mean_trial).^2;
            yo_var_trial = yo_var_trial + abs(yo_trial(number_of_cycles_to_analyze, :) - yo_mean_trial).^2;
            yoi_var_trial = yoi_var_trial + (yo_trial(number_of_cycles_to_analyze, :) - yo_mean_trial).*conj(yi_trial(number_of_cycles_to_analyze, :) - yi_mean_trial);
        end
        yi_var_trial = psfactor_trial * yi_var_trial / (number_of_cycles_to_analyze - 1);
        yo_var_trial = psfactor_trial * yo_var_trial / (number_of_cycles_to_analyze - 1);
        yoi_var_trial = psfactor_trial * yoi_var_trial / (number_of_cycles_to_analyze - 1);
        
        % Calculate variance of FRF (frequency response function) - Pintelon & Schoukens eq 2-32
        %   Calculate at all frequencies, then decimate to get final values at
        %   PRTS frequencies that have stimulus energy.
        FRF_trial = yo_mean_trial./yi_mean_trial;   % FRF calculation at all frequencies - Pintelon & Schoukens eq 2-30
        FRF_var_trial = (abs(FRF_trial).^2).*(yo_var_trial./(psfactor_trial*abs(yo_mean_trial).^2)+yi_var_trial./(psfactor_trial*abs(yi_mean_trial).^2)+2*real(yoi_var_trial./(psfactor_trial*yo_mean_trial.*conj(yi_mean_trial))));
        FRFd_var_trial = decimate2(FRF_var_trial,2);

        % limit to frequencies where stimulus had energy
        FRFd_trial = decimate2(FRF_trial,2);
        fd = decimate2(f, 2);      % frequencies with stimulus PRTS energy
        yid_mean_trial = decimate2(yi_mean_trial, 2);
        yod_mean_trial = decimate2(yo_mean_trial, 2);
        yoid_mean_trial = decimate2(yoi_mean_trial, 2);
        yiid_mean_trial = decimate2(yii_mean_trial, 2);
        yood_mean_trial = decimate2(yoo_mean_trial, 2);
        
        FRFd_trial = FRFd_trial(1 : signal_data.number_of_excited_frequencies);
        fd = fd(1 : signal_data.number_of_excited_frequencies);
        yid_mean_trial = yid_mean_trial(1 : signal_data.number_of_excited_frequencies);
        yod_mean_trial = yod_mean_trial(1 : signal_data.number_of_excited_frequencies);
        yoid_mean_trial = yoid_mean_trial(1 : signal_data.number_of_excited_frequencies);
        yiid_mean_trial = yiid_mean_trial(1 : signal_data.number_of_excited_frequencies);
        yood_mean_trial = yood_mean_trial(1 : signal_data.number_of_excited_frequencies);
        
        % smooth frequencies
%         fd = frequencySmoothing_multisine(fd);
%         FRFd_trial = frequencySmoothing_multisine(FRFd_trial);
%         yoid_mean_trial = frequencySmoothing_multisine(yoid_mean_trial);
%         yiid_mean_trial = frequencySmoothing_multisine(yiid_mean_trial);
%         yood_mean_trial = frequencySmoothing_multisine(yood_mean_trial);
        
        % calculate coherence
        Cohd_trial = (abs(yoid_mean_trial).^2) ./ (yiid_mean_trial .* yood_mean_trial);
        
        % calculate average response to an individual GVS pulse
        
        
        
        
        % model fit
        if fit_model || export_data
            F_model = fd;
            FRF_model = FRFd_trial;
            Coh_model = Cohd_trial;

            Mass_kg = 85;
            height_m = 1.87;
            Len_L_m = height_m * (0.285 - 0.039);
            Len_T_m = height_m * (0.530 - 0.285);
            Len_HAT_m = height_m * (0.818 - 0.530);
            [J, COM] = BodyCalc(Mass_kg, Len_L_m, Len_T_m, Len_HAT_m);
            [Fit, mse] = FRF_fit_and_plot_Reimann ...
              ( ...
                FRF_model, ...
                F_model, ...
                Coh_model, ...
                average_vis_stimulus_trial, ...
                average_filtered_response_trial, ...
                ['metronome=' strrep(num2str(this_pace), '.', '_'), '_amplitude=' num2str(vis_stimulus_amplitudes(i_trial))], ...
                '', ...
                J, ...
                Mass_kg, ...
                COM, ...
                Mass_kg * 9.81 * COM, ...
                0, ...
                1, ...
                5, ...
                1 ...
              );
            model_fits{i_trial} = Fit;
            
            if export_data
                data_for_export.biomechanics.moment_of_inertia = J;
                data_for_export.biomechanics.mass = Mass_kg;
                data_for_export.biomechanics.center_of_mass_height = COM;
                data_for_export.units.biomechanics.moment_of_inertia = 'kg m^2';
                data_for_export.units.biomechanics.mass = 'kg';
                data_for_export.units.biomechanics.center_of_mass_height = 'm';
            end            
        end
        
        % store results from this trial
        FRF_decimated{i_trial} = FRFd_trial;
        Coh_decimated{i_trial} = Cohd_trial;
        average_vis_stimulus{i_trial} = average_vis_stimulus_trial;
        average_gvs_stimulus{i_trial} = average_gvs_stimulus_trial;
        average_response{i_trial} = average_response_trial;
        average_pulse_response{i_trial} = average_pulse_response_trial;
        average_pulse_filtered_response{i_trial} = average_pulse_filtered_response_trial;
        average_filtered_response{i_trial} = average_filtered_response_trial;
        stimulus_all{i_trial} = this_vis_stimulus;
        response_all{i_trial} = response;
        
        if export_data
            filename = makeFileName(collection_date, subject_id, trial_type, trial_number, 'dataForBob.mat');
            save(filename, '-struct', 'data_for_export');
        end
    end

end


if plot_results
    if dictate_axes
        gain_ylim = [0.0005 0.1];
        phase_ylim = [-200 500];
        com_response_ylim = [-0.15 0.15];
        gvs_response_ylim = [-0.16 0.04];
    end
    
    
    frequency_limits = [0.023 0.5];
    unique_paces = sort(unique(paces));

    number_of_vis_amplitudes = length(unique(vis_stimulus_amplitudes));
    number_of_paces = length(unique(paces));
    marker_pace_cell = ...
      { ...
        0.6, 'o'; ...
        1, 'd'; ...
      };
    marker_pace_table = cell2table(marker_pace_cell, 'VariableNames', {'pace', 'marker'});
    all_vis_stim_amplitudes = [0; 6; 12];
    colors_noGvs = copper(numel(all_vis_stim_amplitudes));
    colors_withGvs = winter(numel(all_vis_stim_amplitudes));
    color_gain_table = table(all_vis_stim_amplitudes, colors_noGvs, colors_withGvs, 'VariableNames', {'vis_stim_amplitude', 'color_noGvs', 'color_withGvs'});
    

    % create figure
    figure('position', [0 0 1200 1800])
    tiledlayout(3,2,'TileSpacing','Compact', 'padding', 'compact');
    
    axes_321 = nexttile;
    hold(axes_321, 'on')
    set(axes_321, 'XScale', 'log', 'YScale', 'log', 'fontsize', 12)
    xlim(frequency_limits)
    if dictate_axes
        ylim(gain_ylim)
    end
    ylabel(['FRF gain (' response_unit '/deg)'], 'fontsize', 18)

    axes_322 = nexttile;
    hold(axes_322, 'on')
    set(axes_322, 'fontsize', 12)
    axis([0 max(time_single_cycle) -max(average_vis_stimulus_trial)*1.2 max(average_vis_stimulus_trial)*1.2])
    ylabel('Average stimulus', 'fontsize', 18)

    axes_323 = nexttile;
    hold(axes_323, 'on')
    set(axes_323, 'XScale', 'log', 'fontsize', 12)
    xlim(frequency_limits)
    if dictate_axes
        ylim(phase_ylim)
    end
    ylabel('FRF phase', 'fontsize', 18)
    legend('show', 'Location', 'best')

    axes_324 = nexttile;
    hold(axes_324, 'on')
    set(axes_324, 'fontsize', 12)
    xlim([0 max(time_single_cycle)])
    if dictate_axes
        ylim(com_response_ylim)
    end
    ylabel(['Average ' response_label ' response (' response_unit ')'], 'fontsize', 18)

    axes_325 = nexttile;
    hold(axes_325, 'on')
    set(axes_325, 'XScale', 'log', 'fontsize', 12)
    xlim(frequency_limits); ylim([0 1])
    ylabel('Coherence', 'fontsize', 18)
    xlabel('Frequency (Hz)', 'fontsize', 18)

    axes_326 = nexttile;
    hold(axes_326, 'on')
    set(axes_326, 'fontsize', 12)
%     xlim([0.5 max(gains)+0.5])
%     ylabel('Sensory Weights from model fit', 'fontsize', 18)
    if dictate_axes
        ylim(gvs_response_ylim)
    end
    ylabel('Average response to GVS pulse', 'fontsize', 18)
%     legend('show', 'Location', 'best')
%     xlim([0, 16])
%     ylim([0, 1.1])

    
    linewidth = 2;
    for i_trial = 1 : number_of_trials
        this_vis_stim_amplitude = vis_stimulus_amplitudes(i_trial);
        this_gvs_stim_amplitude = gvs_stimulus_amplitudes(i_trial);
        this_pace = paces(i_trial);
        
        this_vis_stim_index = find(color_gain_table.vis_stim_amplitude == this_vis_stim_amplitude);
        this_pace_index = find(marker_pace_table.pace == this_pace);
        
        if this_gvs_stim_amplitude == 0
            this_color = color_gain_table.color_noGvs(this_vis_stim_index, :);
        else
            this_color = color_gain_table.color_withGvs(this_vis_stim_index, :);
        end
        this_marker = marker_pace_table.marker{this_pace_index};
        this_label = [num2str(vis_stimulus_amplitudes(i_trial)) ' deg, ' num2str(this_pace) ' Hz'];
        
        % plot FRF
        plot(axes_321, fd, abs(FRF_decimated{i_trial}), '-', 'marker', this_marker, 'color', this_color, 'linewidth', linewidth, 'MarkerFaceColor', this_color, 'MarkerEdgeColor', 'none', 'displayname', this_label);

        % plot phase
        plot(axes_323, fd, 180/pi*phase(FRF_decimated{i_trial}), '-', 'marker', this_marker, 'color', this_color, 'linewidth', linewidth, 'MarkerFaceColor', this_color, 'MarkerEdgeColor', 'none', 'displayname', this_label);

        % plot coherence
        plot(axes_325, fd, Coh_decimated{i_trial}, '-', 'marker', this_marker, 'color', this_color, 'linewidth', linewidth, 'MarkerFaceColor', this_color, 'MarkerEdgeColor', 'none');

        % plot single cycle stimulus
        plot(axes_322, time_single_cycle, average_vis_stimulus{i_trial}, 'displayname', this_label, 'color', this_color, 'linewidth', linewidth, 'MarkerFaceColor', this_color, 'MarkerEdgeColor', 'none');

        % plot single cycle filtered response
        plot(axes_324, time_single_cycle, average_filtered_response{i_trial}, 'color', this_color, 'linewidth', linewidth, 'MarkerFaceColor', this_color, 'MarkerEdgeColor', 'none');

%         if fit_model
%             plot(axes_326, vis_stimulus_amplitudes(i_trial), model_fits{i_trial}.gn, 'marker', this_marker, 'displayname', this_label, 'color', this_color, 'markersize', 18, 'MarkerFaceColor', this_color, 'MarkerEdgeColor', 'none');
%         end
        % plot GVS pulse response
        plot(axes_326, time_single_pulse, average_pulse_filtered_response{i_trial}, 'color', this_color, 'linewidth', linewidth, 'MarkerFaceColor', this_color);
        
    end
        
    yyaxis(axes_322, 'right')
    token_gvs_index = find(gvs_stimulus_amplitudes, 1);
    plot(axes_322, time_single_cycle, average_gvs_stimulus{token_gvs_index}, 'color', [0.5 0.5 0.5], 'linewidth', 3);
    axes_322.YColor = [0.5 0.5 0.5];

    yyaxis(axes_326, 'right')
    this_time = signal_data.single_pulse_time;
    
    this_pulse = signal_data.pulse_stimulus_normalized(signal_data.points_per_cycle + (1 : signal_data.points_per_pulse)) * 0.25;
    plot(axes_326, this_time, this_pulse, 'color', [0.5 0.5 0.5], 'linewidth', 3);
    axes_326.YColor = [0.5 0.5 0.5];
%     ylim(gvs_response_ylim * 2)
    
    % connect weights
    if fit_model
        weights = zeros(size(vis_stimulus_amplitudes));
        for i_trial = 1 : number_of_trials
            weights(i_trial) = model_fits{i_trial}.gn;
        end
        for i_pace = 1 : number_of_paces
            this_pace = unique_paces(i_pace);
            this_pace_indices = sort(find(paces == this_pace));
            this_pace_gains = vis_stimulus_amplitudes(this_pace_indices);
            this_pace_weights = weights(this_pace_indices);
            connector = plot(axes_326, this_pace_gains, this_pace_weights, '-', 'color', [0.6 0.6 0.6], 'linewidth', 2, 'HandleVisibility', 'off');
            uistack(connector, 'bottom')
        end
        
    end
    
    % make title
    this_title = title_label;
    uicontrol('style', 'text', 'string', this_title, 'units', 'normalized', 'position', [0, 0.97, 1, 0.03], 'fontsize', 14, 'FontWeight', 'bold')
    
    % resize and print to pdf
    if save_figure
        if ~directoryExists('figures')
            mkdir('figures')
        end
        set(gcf, 'Position', [100 100 600 800])
        saveas(gcf, ['figures' filesep filename], 'pdf')
        close all
    end    
    
    
end


% export results
if export_results
    data_to_export = struct;
    data_to_export.FRF_decimated = FRF_decimated;
    data_to_export.Coh_decimated = Coh_decimated;
    data_to_export.f_decimated = fd;
    data_to_export.time_single_cycle = time_single_cycle;
    data_to_export.average_stimulus = average_vis_stimulus;
    data_to_export.average_filtered_response = average_filtered_response;
    data_to_export.average_pulse_filtered_response = average_pulse_filtered_response;
    data_to_export.multisine_stimulus_amplitudes = vis_stimulus_amplitudes;
    data_to_export.pulse_stimulus_amplitudes = gvs_stimulus_amplitudes;
    data_to_export.paces = paces;
    
    save(['results' filesep filename '_' subject_id], '-struct', 'data_to_export');
end




