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

% analyze the data

% input
% relevantDataStretches.mat

function analyzeVarianceAcrossTime(varargin)
    [trial_type_list, trial_number_list] = parseTrialArguments(varargin{:});
    
    % load settings
    study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
    model_data = load('subjectModel.mat');
    
    across_time_conditions = study_settings.get('across_time_conditions', 1);
    variance_variables = study_settings.get('variance_variables');
    number_of_variance_variables = size(variance_variables, 1);
    group_label = subject_settings.get('group');
    
    % make containers to hold the data
    subject_list_session = {};
    condition_list_session = {};
    time_point_list_session = {};
    group_list_session = {};
    block_list_session = {};
    origin_trial_list_session = [];
    
    % define stretch variables from list of ucm and variance variables
    variance_data_session = cell(number_of_variance_variables, 1);
    variance_directions_session = cell(number_of_variance_variables, 2);
    [variance_directions_session{:, :}] = deal('~');
    variance_names_session = {};
    for i_variable = 1 : number_of_variance_variables
        this_variable_name = variance_variables{i_variable, 1};
        variance_names_session = [variance_names_session; this_variable_name]; %#ok<AGROW>
    end
    
    % calculate variance across time within a single trial
    for i_condition = 1 : length(across_time_conditions)
        % get list of trials for this condition
        this_condition_label = across_time_conditions{i_condition};
        this_condition_index = strcmp(trial_type_list, this_condition_label);
        if any(this_condition_index)
            this_condition_trial_numbers = trial_number_list{this_condition_index};
        else
            this_condition_trial_numbers = [];
        end
        
        % analyze across time
        number_of_trials_this_condition = length(this_condition_trial_numbers);
        expected_event_labels = {'trial_start_time';'trial_end_time'};
        number_of_events = length(expected_event_labels);

        for i_trial = 1 : number_of_trials_this_condition
            % load data
            loaded_data = load(['processed' filesep makeFileName(collection_date, subject_id, this_condition_label, this_condition_trial_numbers(i_trial), 'kinematicTrajectories.mat')]);
            joint_angle_trajectories = loaded_data.joint_angle_trajectories;
            joint_angle_labels = loaded_data.joint_labels; %#ok<NASGU>
            time_mocap = loaded_data.time_mocap;
            loaded_data = load(['processed' filesep makeFileName(collection_date, subject_id, this_condition_label, this_condition_trial_numbers(i_trial), 'comTrajectories.mat')]);
            com_trajectories = loaded_data.com_trajectories;
            com_labels = loaded_data.com_labels; %#ok<NASGU>
            loaded_data = load(['processed' filesep makeFileName(collection_date, subject_id, this_condition_label, this_condition_trial_numbers(i_trial), 'eefTrajectories.mat')]);
            eef_trajectories = loaded_data.eef_trajectories;
            eef_labels = loaded_data.eef_labels; %#ok<NASGU>

            % load and check events
            loaded_data = load(['analysis' filesep makeFileName(collection_date, subject_id, this_condition_label, this_condition_trial_numbers(i_trial), 'events.mat')]);
            if ~(length(loaded_data.event_data)==number_of_events)
                error(['Trial ' num2str(this_condition_trial_numbers(i_trial)) ' - expected ' num2str(number_of_events) ' events, but found ' num2str(length(event_indices_mocap))]);
            end

            % extract data data indices for events
            event_times = zeros(1, number_of_events);
            event_indices_mocap = zeros(1, number_of_events);
            for i_event = 1 : number_of_events
                if ~(any(strcmp(expected_event_labels{i_event}, loaded_data.event_labels)))
                    error(['Trial ' num2str(this_condition_trial_numbers(i_trial)) ' - event label "' expected_event_labels{i_event} '" not found']);
                end
                event_times(i_event) = loaded_data.event_data{strcmp(loaded_data.event_labels, expected_event_labels{i_event})};
                event_indices_mocap(i_event) = findClosestIndex(event_times(i_event), time_mocap);
            end
            
            % extract data at the events
            joint_data_to_analyze = joint_angle_trajectories(event_indices_mocap(1) : event_indices_mocap(2), :);
            com_data_to_analyze = com_trajectories(event_indices_mocap(1) : event_indices_mocap(2), :); %#ok<NASGU>
            eef_data_to_analyze = eef_trajectories(event_indices_mocap(1) : event_indices_mocap(2), :); %#ok<NASGU>
            
            % analyze variance of other variables
            for i_variable = 1 : number_of_variance_variables
                % get data
                this_variable_label = variance_variables{i_variable, 1};
                this_variable_source = variance_variables{i_variable, 2};
                this_variable_split = strsplit(this_variable_source, ':');
                this_variable_source_type = this_variable_split{1};
                this_variable_source_label = this_variable_split{2};
                eval(['source_data = ' this_variable_source_type '_data_to_analyze;']);
                eval(['source_labels = ' this_variable_source_type '_labels;']);
                this_variable_data = source_data(:, strcmp(source_labels, this_variable_source_label)); %#ok<IDISVAR,NODEF>
                
                % calculate variance
                this_variable_stretch_data = var(this_variable_data);
                
                % store
                variance_data_session{strcmp(variance_names_session, this_variable_label)} = [variance_data_session{strcmp(variance_names_session, this_variable_label)}, this_variable_stretch_data];
            end            
            
            % store supplementary information
            subject_list_session = [subject_list_session; subject_id]; %#ok<AGROW>
            time_point_list_session = [time_point_list_session; 'quiet stance across time']; %#ok<AGROW>
            condition_list_session = [condition_list_session; this_condition_label]; %#ok<AGROW>
            group_list_session = [group_list_session; group_label]; %#ok<AGROW>
            block_list_session = [block_list_session; 'N/A']; %#ok<AGROW>
            origin_trial_list_session = [origin_trial_list_session; this_condition_trial_numbers(1)]; %#ok<AGROW>                            
        end        
    end
    
    % determine meta information
    number_of_stretches = length(subject_list_session);
    origin_start_time_list_session = zeros(number_of_stretches, 1); % doesn't apply, but needs to be here for now
    origin_end_time_list_session = zeros(number_of_stretches, 1); % doesn't apply, but needs to be here for now
    time_list_session = zeros(number_of_stretches, 1); % doesn't apply, but needs to be here for now
    
    %% save data
    bands_per_stretch = 1;
    conditions_session = struct;
    conditions_session.subject_list = subject_list_session;
    conditions_session.condition_list = condition_list_session;
    conditions_session.time_point_list = time_point_list_session;
    conditions_session.group_list = group_list_session;
    conditions_session.block_list = block_list_session;
    if ~directoryExists('results')
        mkdir('results')
    end
    results_file_name = ['results' filesep makeFileName(collection_date, subject_id, 'resultsVarianceAcrossTime')];
    save ...
      ( ...
        results_file_name, ...
        'conditions_session', ...
        'variance_data_session', ...
        'variance_names_session', ...
        'variance_directions_session', ...
        'bands_per_stretch', ...
        'origin_trial_list_session', ...
        'origin_start_time_list_session', ...
        'origin_end_time_list_session', ...
        'time_list_session' ...
      )
end

function jacobian = determineJacobian(kinematic_tree, variable_label)
    % calculate Jacobian
    if strcmp(variable_label, 'com_ap')
        J_com = kinematic_tree.calculateCenterOfMassJacobian;
        jacobian = J_com(1, :);
    end
    if strcmp(variable_label, 'com_vert')
        J_com = kinematic_tree.calculateCenterOfMassJacobian;
        jacobian = J_com(3, :);
    end
    if strcmp(variable_label, 'com_2d')
        J_com = kinematic_tree.calculateCenterOfMassJacobian;
        jacobian = J_com([1 3], :);
    end

end

