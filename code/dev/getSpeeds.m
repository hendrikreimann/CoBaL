% script to extract last 10 s of belt speed

function getSpeeds(collection_date, subject_id)

    % set parameters
    trial_type = 'adaptation';
    trials_to_analyze = 2 : 6;
    duration_to_analyze = 30; % duration of the interval at the end of the trial that will be analyzed, in seconds
    correction_factor = 1.175; % the commanded speed in labview is multiplied by this factor before sending to the PLC, so we need to correct for it
    reduction_factor = 0.9;
    decimals_to_show = 2;

    % load settings and protocol
    % subject_settings = loadSettingsFromFile('subjectSettings.txt');
    % collection_date = subject_settings.get('collection_date');
    % subject_id = subject_settings.get('subject_id');
    importProtocolData('source', pwd)
    protocol_data = load('protocolInfo.mat');
    importLabviewData('sampling_rate', 10, 'source', pwd)

    % extract data
    number_of_trials = length(trials_to_analyze);
    speeds = zeros(number_of_trials, 1);
    cadences = zeros(number_of_trials, 1);
    for i_trial = 1 : number_of_trials
        % get info for this file
        trial_number = trials_to_analyze(i_trial);
        filename = makeFileName(collection_date, subject_id, trial_type, trial_number, 'PLCData.mat');

        % load belt speed data
        [belt_speed_left, time] = loadData(collection_date, subject_id, trial_type, trial_number, 'belt_speed_left_trajectory');
        belt_speed_right = loadData(collection_date, subject_id, trial_type, trial_number, 'belt_speed_right_trajectory');
        belt_speed = (belt_speed_left + belt_speed_right) * 0.5;

        % extract last stretch median
        time_to_analyze_start = time(end) - duration_to_analyze;
        belt_speed_to_analyze = belt_speed(time > time_to_analyze_start);
        speed_this_trial = mean(belt_speed_to_analyze);
        speeds(i_trial) = abs(speed_this_trial);

        trial_index = (protocol_data.trial_number == trial_number) & strcmp(protocol_data.trial_type, trial_type);
        this_trial_cadence = protocol_data.metronome_cadence(trial_index);
        cadences(i_trial) = this_trial_cadence;    
    end

    % sort
    [cadences, sort_indices] = sort(cadences);
    speeds = speeds(sort_indices);

    % apply correction factor and reduction factor and round
    speeds = speeds * 1/correction_factor;
    speeds = speeds * reduction_factor;
    speeds = round(speeds, decimals_to_show);

    % write to csv
    speed_table = table(cadences, speeds, 'VariableNames', {'Cadence', 'Starting Speed'} );
    writetable(speed_table, 'speeds.csv')

    % plot
    figure; hold on
    color = lines(1);
    plot(cadences, speeds, 's-', 'linewidth', 3, 'markersize', 18, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none')
    xlim([min(cadences)-5, max(cadences)+5]);
    xlabel('Cadence'); ylabel('Speed')

end