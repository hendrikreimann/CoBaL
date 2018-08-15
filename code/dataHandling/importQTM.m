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

% this script transform raw data from ascii into matlab format

% input:
% Experimental data files generated by e.g. Nexus, QTM, labview etc.
% Files are expected to be in subfolders named <subject code>_<source type>, e.g. XYZ_labview
% any source type is viable, but if it is not in the default list, it must be specified as a name-value pair,
% e.g. importAscii('sources', 'someSource') would look for data in the subfolder XYZ_someSource
%
% output:
% Files containing the same data in .mat format, with some additional information about where they came from.
% Output files will be saved to folders "raw" and "processed".

function importQTM()

    %% prepare
    % set some parameters
    millimeter_to_meter = 1e-3;
    centimeter_to_meter = 1e-2;
    milliseconds_to_seconds = 1e-3;
    qtm_emg_scale = 1;

    subject_settings = SettingsCustodian('subjectSettings.txt');
    analog_to_protocol_mapping = subject_settings.get('analog_to_protocol_mapping');
    
    % initialize
    total_number_of_trials_extracted_this_subject = 0;

    % create folders if necessary
    if ~directoryExists('raw')
        mkdir('raw')
    end
    if ~directoryExists('processed')
        mkdir('processed')
    end
    if ~directoryExists('analysis')
        mkdir('analysis')
    end
    current_path = pwd;
    path_split = strsplit(current_path, filesep);
    subject_code = path_split{end};

    labview_source_dir = 'labview';
    qtm_source_dir = 'qtm';

    %% import labview data
    % get list of files to import from this directory
    clear file_name_list;
    data_dir_csv = dir([labview_source_dir filesep '*csv']);
    [file_name_list{1:length(data_dir_csv)}] = deal(data_dir_csv.name);

    % import protocol
    protocol_file_available = 0;
    if any(strcmp(file_name_list, 'protocol.csv'))
        protocol_file_available = 1;

        % we have a protocol file, so import that
        [imported_data, delimiter, nheaderlines] = importdata([labview_source_dir filesep 'protocol.csv'], ',', 1);
        headers = imported_data.textdata(1, :);
        textdata = imported_data.textdata(2:end, :);
        protocol_data = imported_data.data;
        protocol_headers = headers(2:end);
        protocol_trial_type = textdata(:, strcmp(headers, 'Trial Type'));
        if strcmp(protocol_trial_type{end}, '(stop)')
            protocol_trial_type(end) = [];
        end
        protocol_trial_number = protocol_data(:, strcmp(protocol_headers, 'Trial Number'));
        protocol_trial_duration = protocol_data(:, strcmp(protocol_headers, 'Duration (s)'));
        protocol_metronome_cadence = protocol_data(:, strcmp(protocol_headers, 'Metronome'));
        protocol_trial_saved = protocol_data(:, strcmp(protocol_headers, 'save data (0/1)'));
        protocol_count_left_step = protocol_data(:, strcmp(protocol_headers, 'Count left steps (0/1)'));
        protocol_count_right_step = protocol_data(:, strcmp(protocol_headers, 'Count right steps (0/1)'));
        protocol_stim_visual_intermittent = protocol_data(:, strcmp(protocol_headers, 'Use Visual Stimulus - intermittent'));
        protocol_stim_gvs_intermittent = protocol_data(:, strcmp(protocol_headers, 'GVS intermittent'));

        % save protocol data
        protocol_data = struct;
        protocol_data.trial_type = protocol_trial_type;
        protocol_data.trial_number = protocol_trial_number;
        protocol_data.trial_duration = protocol_trial_duration;
        protocol_data.metronome_cadence = protocol_metronome_cadence;
        protocol_data.trial_saved = protocol_trial_saved;
        protocol_data.count_left_step = protocol_count_left_step;
        protocol_data.count_right_step = protocol_count_right_step;
        protocol_data.stim_visual_intermittent = protocol_stim_visual_intermittent;
        protocol_data.stim_gvs_intermittent = protocol_stim_gvs_intermittent;
        save_file_name = 'protocolInfo.mat';
        save(save_file_name, '-struct', 'protocol_data');

        % remove the protocol file from the list
        file_name_list(strcmp(file_name_list, 'protocol.csv')) = [];
    end

    % import labview saved data
    number_of_files = length(file_name_list);
    for i_file = 1 : number_of_files
        % file name stuff
        data_file_name = file_name_list{i_file};
        [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(data_file_name);
        imported_data = importdata([labview_source_dir filesep data_file_name], ',', 2);
        labview_trajectories = imported_data.data;

        % extract headers
        column_name_string = imported_data.textdata{1, 1};
        labview_header = strsplit(column_name_string, ',');
        number_of_data_columns = size(imported_data.textdata, 2);

        % extract data into properly named variables
        variables_to_save = struct();
        variables_to_save_list = {};
        for i_column = 1 : number_of_data_columns
            variable_name = [strrep(labview_header{i_column}, ' ', '_'), '_trajectory'];
            extract_string = ['variables_to_save.' variable_name ' = labview_trajectories(:, i_column);'];
            eval(extract_string);
            variables_to_save_list = [variables_to_save_list; variable_name];
        end

        % take special care of time, transform to seconds and rename according to file type
        variables_to_save.time = variables_to_save.time_trajectory * milliseconds_to_seconds; % transform to seconds
        variables_to_save.time = variables_to_save.time - variables_to_save.time(1); % zero
        variables_to_save = rmfield(variables_to_save, 'time_trajectory'); % this is not a variable, so remove from list
        variables_to_save_list(strcmp(variables_to_save_list, 'time_trajectory')) = [];

        % add data source
        variables_to_save.data_source = 'labview';

        % add sampling rate
        variables_to_save.sampling_rate = NaN;

        % save
        save_folder = 'processed';
        save_file_name = makeFileName(date, subject_id, trial_type, trial_number, file_type);
        save([save_folder filesep save_file_name], '-struct', 'variables_to_save');

        for i_variable = 1 : length(variables_to_save_list)
            if ~checkDataAvailability(date, subject_id, trial_type, trial_number, variables_to_save_list{i_variable})
                addAvailableData ...
                  ( ...
                    variables_to_save_list{i_variable}, ...
                    'time', ...
                    'sampling_rate', ...
                    variables_to_save_list{i_variable}, ...
                    {'~', '~'}, ... % placeholder for direction
                    save_folder, ...
                    save_file_name ...
                  );
            end
        end

        disp(['imported ' labview_source_dir filesep data_file_name ' and saved as ' save_folder filesep save_file_name])



    end

    %% import qtm data
    clear file_name_list;
    data_dir_mat = dir([qtm_source_dir filesep '*mat']);
    [file_name_list{1:length(data_dir_mat)}] = deal(data_dir_mat.name);

    % go through files and import
    number_of_files = length(file_name_list);
    importing_trial_number = 0;
    for i_file = 1 : number_of_files
%     for i_file = 3 : number_of_files % XXX to speed things up during bug-fixing
        % file name stuff
        data_file_name = file_name_list{i_file};
        disp(['Importing ' qtm_source_dir filesep data_file_name])
        [date, subject_id, trial_type, trial_number] = getFileParameters(data_file_name);

        % this is marker data from QTM
        data_source = 'qtm';
        var_name = whos('-file', [qtm_source_dir, filesep, data_file_name]);
        temp_data = load([qtm_source_dir, filesep, data_file_name]);
        qtm_data = temp_data.(var_name.name);
        % TODO: this was a hack for having data missing, fix this!
        if strcmp(trial_type, 'calibration') || strcmp(trial_type, 'static')
            analog_fs = 2000;
            start_indices = 1;
            end_indices = 20000;
            this_trial_duration = 10;
        else
            analog_fs = qtm_data.Analog.Frequency;

            % this is walking data, so break up into chunks for trials if necessary
            trigger_mask = contains(qtm_data.Analog.Labels, 'labview_sync');
            trigger = qtm_data.Analog.Data(trigger_mask,:);
            % normalise to the range [0,1] and round
            trigger = round(trigger/max(trigger));

            % find edges
            trigger_edges = [diff(trigger), 0];
            analog_index = 1:length(qtm_data.Analog.Data);
            start_indices = analog_index(trigger_edges == -1); % trials start on negative edge
            end_indices = analog_index(trigger_edges == 1); % trials end on positive edge

            % figure out protocol steps
            protocol_step_mask = contains(qtm_data.Analog.Labels, 'currentStep');
            protocol_step_analog = qtm_data.Analog.Data(protocol_step_mask,:);
%             [~, protocol_step_up_indices] = findpeaks([diff(protocol_step_analog), 0], 'MinPeakProminence', 0.08);
            protocol_step_up_indices = find([diff(protocol_step_analog), 0] > 0.08); % TODO: changed this to a simple thresholding because of increased noise level, revisit this at some point!
            
            if isempty(start_indices) && isempty(end_indices)
                start_indices = 1;
                end_indices = length(analog_index);
            end

            % check whether each start index has a matching end index
            for i_start_index = 1 : length(start_indices)
                % remove end indices that come before the start index - HR: not tested yet
                while length(end_indices) >= i_start_index && end_indices(i_start_index) <= start_indices(i_start_index)
                    end_indices(i_start_index) = [];
                end
            end
            if length(end_indices) > length(start_indices)
                disp('More end indices than start indices, removing superfluous ones')
                end_indices(length(start_indices)+1 : end) = [];
            end
            if length(start_indices) > length(end_indices)
                disp('More start indices than end indices, removing superfluous ones')
                start_indices(length(end_indices)+1 : end) = [];
            end

            figure; hold on;
            time_analog = (1 : length(protocol_step_analog)) * 1/analog_fs;
            plot(time_analog, protocol_step_analog);
            plot(time_analog(protocol_step_up_indices), protocol_step_analog(protocol_step_up_indices), '^');
            plot(time_analog(start_indices), protocol_step_analog(start_indices), '>');
            plot(time_analog(end_indices), protocol_step_analog(end_indices), '<');
            drawnow
            
            if end_indices(end) > length(qtm_data.Analog.Data)
                end_indices(end) = length(qtm_data.Analog.Data);
            end
        end
        number_of_trials_in_this_qtm_file = length(start_indices);
        for i_trial_this_qtm_file = 1 : number_of_trials_in_this_qtm_file
            this_trial_start_index = start_indices(i_trial_this_qtm_file);
            this_trial_end_index = end_indices(i_trial_this_qtm_file);
            number_of_samples = this_trial_end_index - this_trial_start_index + 1;
            this_trial_length = number_of_samples / analog_fs;

            % figure out information from protocol, or define for calibration trials
            if strcmp(trial_type, 'calibration') || strcmp(trial_type, 'static')
                importing_trial_number = trial_number;
                importing_trial_type = trial_type;
                save_this_trial = 1;
                this_trial_duration = 10;
            else
                if ~protocol_file_available
                    error('Failed to locate protocol.csv in the labview folder, exiting.')
                end
                
                % check if a mapping was provided for this step
                use_analog_signal = true;
                trial_type_matches = strcmp(analog_to_protocol_mapping(:, 1), trial_type);
                trial_number_matches = strcmp(analog_to_protocol_mapping(:, 2), num2str(trial_number));
                index_number_matches = strcmp(analog_to_protocol_mapping(:, 3), num2str(i_trial_this_qtm_file));
                
                if any(trial_type_matches & trial_number_matches & index_number_matches)
                    use_analog_signal = false;
                    this_trial_protocol_index = str2num(analog_to_protocol_mapping{trial_type_matches & trial_number_matches & index_number_matches, 4});
                end
                    
                if use_analog_signal
                    % figure out protocol step from the analog signal
                    this_trial_protocol_step_analog = protocol_step_analog(this_trial_start_index:this_trial_end_index);

                    analog_to_step_range = (-10 : 0.1 : 10)';
                    offset = analogOffset();
                    this_trial_offset = interp1(analog_to_step_range, offset, this_trial_protocol_step_analog);
                    this_trial_protocol_step_analog_corrected = this_trial_protocol_step_analog + this_trial_offset;

                    % discrete map from [-10 10] --> [0 200]
                    this_trial_protocol_step_digital = round((this_trial_protocol_step_analog_corrected + 10) * 10);
                    this_trial_protocol_step = median(this_trial_protocol_step_digital);

                    if sum(this_trial_protocol_step_digital(2:end-1)~=median(this_trial_protocol_step_digital)) > 20
                        warning(['Ambiguous protocol step data for trial starting at time step ' num2str(this_trial_start_index)]);
                    end
                    this_trial_protocol_index = median(this_trial_protocol_step) + 1;
                end
                    
                importing_trial_type = protocol_trial_type{this_trial_protocol_index};
                importing_trial_number = protocol_trial_number(this_trial_protocol_index);
                this_trial_duration = protocol_trial_duration(this_trial_protocol_index);           
                save_this_trial = protocol_trial_saved(this_trial_protocol_index);

                
                % sanity check: recorded data should be within 0.1% of expected duration
                if this_trial_duration < this_trial_length*0.999 || this_trial_duration > this_trial_length*1.001
                    warning ...
                      ( ...
                        [ ...
                          'Problem with trial starting at time step ' ...
                          num2str(this_trial_start_index) ...
                          ', expected duration is ' ...
                          num2str(this_trial_duration) ...
                          's, but data is ' ...
                          num2str(this_trial_length) ...
                          's long.' ...
                        ] ...
                      );
                end

            end

            if ~strcmp(trial_type, 'calibration') && ~strcmp(trial_type, 'static')
                % EMG
                sampling_rate_emg = analog_fs;
                time_emg = (1 : number_of_samples)' / sampling_rate_emg;
                emg_labels = qtm_data.Analog.Labels(14:29);
                emg_trajectories_raw = qtm_data.Analog.Data(14:29, this_trial_start_index:this_trial_end_index)';
                data_type = 'emg';

                % make directions
                emg_directions = cell(2, length(emg_labels));
                [emg_directions{1, :}] = deal('positive');
                [emg_directions{2, :}] = deal('negative');

                % save emg data
                if save_this_trial
                    save_folder = 'raw';
                    save_file_name = makeFileName(date, subject_id, importing_trial_type, importing_trial_number, 'emgTrajectoriesRaw.mat');
                    save ...
                        ( ...
                        [save_folder filesep save_file_name], ...
                        'emg_trajectories_raw', ...
                        'time_emg', ...
                        'sampling_rate_emg', ...
                        'data_source', ...
                        'emg_labels', ...
                        'emg_directions' ...
                        );
                    addAvailableData('emg_trajectories_raw', 'time_emg', 'sampling_rate_emg', 'emg_labels', '_emg_directions', save_folder, save_file_name);
                end

                % Force data - data are in qtm_data.Force(n).Force
                data_type = 'forceplate';
                if any(any(qtm_data.Force(1).Force))
                    forceplate_tajectories_Left = ...
                        [ ...
                        qtm_data.Force(1).Force(:, this_trial_start_index : this_trial_end_index)', ...
                        qtm_data.Force(1).Moment(:, this_trial_start_index : this_trial_end_index)' ...
                        ];
                    forceplate_tajectories_Right = ...
                        [ ...
                        qtm_data.Force(2).Force(:, this_trial_start_index : this_trial_end_index)', ...
                        qtm_data.Force(2).Moment(:, this_trial_start_index : this_trial_end_index)' ...
                        ];
                    forceplate_trajectories_raw = [forceplate_tajectories_Left, forceplate_tajectories_Right];
                else % currently taking volts... need to scale accordinginly
                    Warning('No force data found, using analog data instead. This is currently not scaling correctly.')
                    forceplate_trajectories_raw = [qtm_data.Analog.Data(1:12, this_trial_start_index : this_trial_end_index)]';
                end
                forceplate_labels = qtm_data.Analog.Labels(1:12);
                %
                forceplate_location_left = mean(qtm_data.Force(1).ForcePlateLocation) * millimeter_to_meter; % mean of corner coordinates gives center
                forceplate_location_right = mean(qtm_data.Force(2).ForcePlateLocation) * millimeter_to_meter; % mean of corner coordinates gives center

                sampling_rate_forceplate = analog_fs;
                time_forceplate = (1 : number_of_samples)' / sampling_rate_forceplate;

                % make directions
                % NOTE: this defines directions and makes assumptions, make sure everything is right here
                forceplate_directions = cell(2, length(forceplate_labels));
                [forceplate_directions{1, [1 4 7 10]}] = deal('right');
                [forceplate_directions{2, [1 4 7 10]}] = deal('left');
                [forceplate_directions{1, [2 5 8 11]}] = deal('forward');
                [forceplate_directions{2, [2 5 8 11]}] = deal('backward');
                [forceplate_directions{1, [3 6 9 12]}] = deal('up');
                [forceplate_directions{2, [3 6 9 12]}] = deal('down');

                % save forceplate data
                if save_this_trial
                    save_folder = 'raw';
                    save_file_name = makeFileName(date, subject_id, importing_trial_type, importing_trial_number, 'forceplateTrajectoriesRaw.mat');
                    save ...
                        ( ...
                        [save_folder filesep save_file_name], ...
                        'forceplate_trajectories_raw', ...
                        'forceplate_labels', ...
                        'data_source', ...
                        'time_forceplate', ...
                        'sampling_rate_forceplate', ...
                        'forceplate_location_left', ...
                        'forceplate_location_right', ...
                        'forceplate_directions' ...
                        );
                    addAvailableData('forceplate_trajectories_raw', 'time_forceplate', 'sampling_rate_forceplate', 'forceplate_labels', '_forceplate_directions', save_folder, save_file_name);
                end
            end
            % Markers
            sampling_rate_mocap = qtm_data.FrameRate;
            % align indices
            start_indices_mocap = round(start_indices * sampling_rate_mocap/analog_fs);
            if start_indices_mocap == 0
                start_indices_mocap = 1;
            end
            end_indices_mocap = round(end_indices * sampling_rate_mocap/analog_fs);
            this_trial_start_index_mocap = start_indices_mocap(i_trial_this_qtm_file);
            this_trial_end_index_mocap = end_indices_mocap(i_trial_this_qtm_file);
            number_of_frames = this_trial_end_index_mocap - this_trial_start_index_mocap + 1;

            temp_markers = qtm_data.Trajectories.Labeled.Data(:, 1:3, this_trial_start_index_mocap:this_trial_end_index_mocap);

            data_type = 'markers';
            % deal with marker data

            marker_labels = qtm_data.Trajectories.Labeled.Labels;
            sampling_rate_mocap = qtm_data.FrameRate;
            time_mocap = (1 : number_of_frames)' / sampling_rate_mocap;

            marker_count = 1;
            marker_trajectories_raw = [];
            for i_marker = 1: size(temp_markers,1)
                this_marker = temp_markers(i_marker,:,:);
                marker_trajectories_raw(marker_count:marker_count+2,:) = reshape(this_marker, size(this_marker,2), size(this_marker,3)) * millimeter_to_meter; 
                marker_count = marker_count + 3;
            end
            marker_trajectories_raw = marker_trajectories_raw';

            
            % replace marker labels if necessary
            marker_label_replacement_map = subject_settings.get('marker_label_replacement_map');
            for i_label = 1 : size(marker_label_replacement_map, 1)
                old_label = marker_label_replacement_map{i_label, 1};
                new_label = marker_label_replacement_map{i_label, 2};
                marker_labels{strcmp(marker_labels, old_label)} = new_label;
            end
            
            % triplicate labels
            marker_labels_loaded = marker_labels;
            number_of_markers = length(marker_labels_loaded);
            marker_labels = cell(3, number_of_markers);
            for i_marker = 1 : length(marker_labels)
                marker_labels{1, i_marker} = [marker_labels_loaded{i_marker} '_x'];
                marker_labels{2, i_marker} = [marker_labels_loaded{i_marker} '_y'];
                marker_labels{3, i_marker} = [marker_labels_loaded{i_marker} '_z'];
            end
            marker_labels = reshape(marker_labels, 1, number_of_markers*3);

            % make directions
            % NOTE: this defines directions and makes assumptions, make sure everything is right here
            number_of_marker_trajectories = size(marker_trajectories_raw, 2);
            marker_directions = cell(2, number_of_marker_trajectories);
            [marker_directions{1, 1 : 3 : number_of_marker_trajectories}] = deal('right');
            [marker_directions{2, 1 : 3 : number_of_marker_trajectories}] = deal('left');
            [marker_directions{1, 2 : 3 : number_of_marker_trajectories}] = deal('forward');
            [marker_directions{2, 2 : 3 : number_of_marker_trajectories}] = deal('backward');
            [marker_directions{1, 3 : 3 : number_of_marker_trajectories}] = deal('up');
            [marker_directions{2, 3 : 3 : number_of_marker_trajectories}] = deal('down');


            % save
            if save_this_trial
                save_folder = 'raw';
                save_file_name = makeFileName(date, subject_id, importing_trial_type, importing_trial_number, 'markerTrajectoriesRaw.mat');
                save ...
                    ( ...
                    [save_folder filesep save_file_name], ...
                    'marker_trajectories_raw', ...
                    'time_mocap', ...
                    'data_source', ...
                    'sampling_rate_mocap', ...
                    'marker_labels', ...
                    'marker_directions' ...
                    );
                addAvailableData('marker_trajectories_raw', 'time_mocap', 'sampling_rate_mocap', 'marker_labels', '_marker_directions', save_folder, save_file_name);
            end

            if save_this_trial
%                 disp(['  Saved data files for type ' importing_trial_type ', trial ' num2str(importing_trial_number) ' with duration ' num2str(this_trial_duration) 's'])
                disp( ...
                      [ ...
                        '  Saved data files for type ' importing_trial_type ...
                        ', trial ' num2str(importing_trial_number) ...
                        ' with duration ' num2str(this_trial_duration) ...
                        's (' num2str(this_trial_start_index_mocap * 1/sampling_rate_mocap) ...
                        '-' num2str(this_trial_end_index_mocap * 1/sampling_rate_mocap) ...
                        ')' ...
                      ] ...
                    )
            end
        end
        total_number_of_trials_extracted_this_subject = total_number_of_trials_extracted_this_subject + number_of_trials_in_this_qtm_file;



    end


end




function offset = analogOffset()
    offset = ...
      [ ...
        0; ...
        -0.0317098627837868; ...
        -0.0355064291396836; ...
        -0.0383662432835923; ...
        -0.0407330287129284; ...
        -0.0417785032417743; ...
        -0.0414159429847061; ...
        -0.0417277364763837; ...
        -0.0412429925743254; ...
        -0.0388849311198740; ...
        -0.0363677844014454; ...
        -0.0350654941860018; ...
        -0.0371848832588526; ...
        -0.0339162265276105; ...
        -0.0308550321637160; ...
        -0.0258480019743690; ...
        -0.0227666814457557; ...
        -0.0199298398628880; ...
        -0.0173543278850836; ...
        -0.0132742999011271; ...
        -0.0106610228969846; ...
        -0.00936478599793489; ...
        -0.00473625185189075; ...
        -0.00257883532928460; ...
        -0.000572734957710708; ...
        0.00132398762374120; ...
        0.00347722371700066; ...
        0.00398356837392289; ...
        0.00599911144133358; ...
        0.00753927394536458; ...
        0.00753085203582149; ...
        0.00870865857276826; ...
        0.00853249776011822; ...
        0.00853628316823052; ...
        0.00788663034824477; ...
        0.00840332559286772; ...
        0.00815475542401511; ...
        0.00798220864381083; ...
        0.00709983160489891; ...
        0.00690247828238899; ...
        0.00643181355926981; ...
        0.00503683000999455; ...
        0.00413383944497703; ...
        0.00329203775458176; ...
        0.00242120802380175; ...
        0.00149387016846525; ...
        0.000821574193582642; ...
        0.000369138339894271; ...
        -0.000277877672783866; ...
        -0.00102545670460508; ...
        -0.00271989182455190; ...
        -0.00365981488830602; ...
        -0.00482093843100540; ...
        -0.00494539216657319; ...
        -0.00588487154153494; ...
        -0.00599837804270820; ...
        -0.00682649900498600; ...
        -0.00655203723933173; ...
        -0.00791718834077937; ...
        -0.00825049603177419; ...
        -0.00769935445543801; ...
        -0.00852675928681101; ...
        -0.00815403884460952; ...
        -0.00963015799901790; ...
        -0.00895956615075910; ...
        -0.00961831343282515; ...
        -0.00966309582918967; ...
        -0.00952247347543045; ...
        -0.00919973803649565; ...
        -0.00920595233365562; ...
        -0.00893583474366810; ...
        -0.00877341543027255; ...
        -0.00870756374502735; ...
        -0.00856975753089673; ...
        -0.00870519970119821; ...
        -0.00932097017889166; ...
        -0.00848747970298325; ...
        -0.00844678720233416; ...
        -0.00783984295034923; ...
        -0.00773505034545563; ...
        -0.00705187100005489; ...
        -0.00704416271013297; ...
        -0.00689469491530370; ...
        -0.00674065592462925; ...
        -0.00612953452561960; ...
        -0.00560228535979301; ...
        -0.00535692483822525; ...
        -0.00488826650123864; ...
        -0.00468637134357453; ...
        -0.00436269192422967; ...
        -0.00396357745781306; ...
        -0.00378863032705390; ...
        -0.00367707688493169; ...
        -0.00314717004448473; ...
        -0.00288107214429201; ...
        -0.00249950966782508; ...
        -0.00246150124316424; ...
        -0.00185530510153731; ...
        -0.00161144218362400; ...
        -0.00146637549603207; ...
        -0.00102024355158730; ...
        -0.000815409068261150; ...
        -0.000734256181997323; ...
        -0.000832852161098741; ...
        -0.000792806693304171; ...
        -0.000108319700746151; ...
        -0.000265702930942080; ...
        -0.000279534687805283; ...
        -0.000491790882055976; ...
        -0.000426932735417163; ...
        -0.00117144863524299; ...
        -0.00153105624076844; ...
        -0.000998698616601734; ...
        -0.000935137108814876; ...
        -0.000958296975685613; ...
        -0.00142491381875454; ...
        -0.00164261683167122; ...
        -0.00264091356187746; ...
        -0.00309681499754499; ...
        -0.00401206048587333; ...
        -0.00396441579734863; ...
        -0.00422148403188150; ...
        -0.00529833879784958; ...
        -0.00596774150663659; ...
        -0.00678126028756942; ...
        -0.00723620698253669; ...
        -0.00772259027777311; ...
        -0.00825504622265649; ...
        -0.00816044205792954; ...
        -0.00916009217047842; ...
        -0.00977628691772603; ...
        -0.0106469097669994; ...
        -0.0109333378042864; ...
        -0.0119922291978138; ...
        -0.0122519306587239; ...
        -0.0130286490528753; ...
        -0.0130401366798907; ...
        -0.0137444705007748; ...
        -0.0138314855721378; ...
        -0.0141161149084041; ...
        -0.0144106983101242; ...
        -0.0146721926696580; ...
        -0.0144676458229069; ...
        -0.0148205448241381; ...
        -0.0148786138861228; ...
        -0.0144537309394677; ...
        -0.0143732644422432; ...
        -0.0139970179820059; ...
        -0.0137701494310054; ...
        -0.0128238368055325; ...
        -0.0125910885262401; ...
        -0.0135005954160290; ...
        -0.0126732101628155; ...
        -0.0115996832421770; ...
        -0.0105618871208364; ...
        -0.00950531667495813; ...
        -0.00839970631782627; ...
        -0.00740346782180268; ...
        -0.00593091658393519; ...
        -0.00513883242166635; ...
        -0.00375714780028513; ...
        -0.00162831240730377; ...
        -0.000394209661337186; ...
        0.00243617559528087; ...
        0.00190042906742605; ...
        0.00386235131119417; ...
        0.00538588999007050; ...
        0.00703525654969717; ...
        0.00907773138043932; ...
        0.0108388384692422; ...
        0.0129475880893617; ...
        0.0138018623944065; ...
        0.0152044950641930; ...
        0.0162984264486976; ...
        0.0179685779542922; ...
        0.0188687694610845; ...
        0.0202168224764048; ...
        0.0208049101290397; ...
        0.0220546340257348; ...
        0.0223920981075532; ...
        0.0225321870361555; ...
        0.0229908094999018; ...
        0.0241814339901651; ...
        0.0240792613920782; ...
        0.0238237301350370; ...
        0.0233348590569662; ...
        0.0232842471463997; ...
        0.0227356714357168; ...
        0.0216668784253855; ...
        0.0215045626865980; ...
        0.0215450198313842; ...
        0.0203315176705363; ...
        0.0190577137879586; ...
        0.0179738315997273; ...
        0.0168931232604397; ...
        0.0156401403073261; ...
        0.0150204841270369; ...
        0.0147189566948640; ...
        0.0131844504728420; ...
        0.0123703604766217; ...
        0; ...
      ];
end


















  
