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

function analyzeDisplacementVariables(varargin)
    [trial_type_list, trial_number_list] = parseTrialArguments(varargin{:});
    load('subjectInfo.mat', 'date', 'subject_id');
    
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
    
    % extract info from settings
    displacement_variables_header = study_settings.get('displacement_variables_header');
    displacement_variables = study_settings.get('displacement_variables');
    displacement_conditions = study_settings.get('displacement_conditions');
    group_label = subject_settings.get('group');
    block_assignment = subject_settings.get('block_assignment');
    
    % make containers to hold the data
    subject_list_session = {};
    condition_list_session = {};
    time_point_list_session = {};
    group_list_session = {};
    block_list_session = {};
    origin_trial_list_session = [];
    
    % define stretch variables from list of ucm and variance variables
    number_of_displacement_variables = size(displacement_variables, 1);
    
    displacement_data_session = cell(number_of_displacement_variables, 1);
    displacement_directions_session = cell(number_of_displacement_variables, 2);
    [displacement_directions_session{:, :}] = deal('~');
    displacement_names_session = displacement_variables(:, strcmp(displacement_variables_header, 'variable name'));
    
    
    for i_type = 1 : length(displacement_conditions)
        this_type_label = displacement_conditions{i_type};
        
        % get block information
        block_assignment_this_type = block_assignment(strcmp(block_assignment(:, 1), this_type_label), :);
        
        if any(strcmp(trial_type_list, this_type_label))
            this_type_trial_numbers = trial_number_list{strcmp(trial_type_list, this_type_label)};
            number_of_trials_this_type = length(this_type_trial_numbers);
            expected_event_labels = study_settings.get('event_labels_ramp');
            number_of_events = length(expected_event_labels);
            
            % collect data from different trials
            for i_trial = 1 : number_of_trials_this_type
                % prepare containers
                joint_data_to_analyze = cell(number_of_events, 1);
                com_data_to_analyze = cell(number_of_events, 1);
                eef_data_to_analyze = cell(number_of_events, 1);
                displacement_data_trial = cell(number_of_displacement_variables, 1);
                
                % load data
                loaded_data = load(['processed' filesep makeFileName(date, subject_id, this_type_label, this_type_trial_numbers(i_trial), 'kinematicTrajectories.mat')]);
                joint_angle_trajectories = loaded_data.joint_angle_trajectories;
                joint_angle_labels = loaded_data.joint_labels; %#ok<NASGU>
                time_mocap = loaded_data.time_mocap;
                loaded_data = load(['processed' filesep makeFileName(date, subject_id, this_type_label, this_type_trial_numbers(i_trial), 'comTrajectories.mat')]);
                com_trajectories = loaded_data.com_trajectories;
                com_labels = loaded_data.com_labels; %#ok<NASGU>
                loaded_data = load(['processed' filesep makeFileName(date, subject_id, this_type_label, this_type_trial_numbers(i_trial), 'eefTrajectories.mat')]);
                eef_trajectories = loaded_data.eef_trajectories;
                eef_labels = loaded_data.eef_labels; %#ok<NASGU>
                
                % load and check events
                loaded_data = load(['analysis' filesep makeFileName(date, subject_id, this_type_label, this_type_trial_numbers(i_trial), 'events.mat')]);
                if ~(length(loaded_data.event_data)==number_of_events)
                    error(['Trial ' num2str(this_type_trial_numbers(i_trial)) ' - expected ' num2str(number_of_events) ' events, but found ' num2str(length(event_indices_mocap))]);
                end

                % extract events
                event_times = zeros(1, number_of_events);
                for i_event = 1 : number_of_events
                    if ~(any(strcmp(expected_event_labels{i_event}, loaded_data.event_labels)))
                        error(['Trial ' num2str(this_type_trial_numbers(i_trial)) ' - event label "' expected_event_labels{i_event} '" not found']);
                    end
                    event_times(i_event) = loaded_data.event_data{strcmp(loaded_data.event_labels, expected_event_labels{i_event})};
                end
                event_indices_mocap = findClosestIndex(event_times, time_mocap);

                % extract all data for all events
                for i_event = 1 : number_of_events
                    joint_data_to_analyze{i_event} = joint_angle_trajectories(event_indices_mocap(i_event), :);
                    com_data_to_analyze{i_event} = com_trajectories(event_indices_mocap(i_event), :);
                    eef_data_to_analyze{i_event} = eef_trajectories(event_indices_mocap(i_event), :);
                end
                
                % go through variables to extract reference data and calculate displacement
                for i_variable = 1 : number_of_displacement_variables
                    % extract information for this variable
                    this_variable_name = displacement_variables{i_variable, strcmp(displacement_variables_header, 'variable name')};
                    this_source_variable_name = displacement_variables{i_variable, strcmp(displacement_variables_header, 'source variable name')};
                    this_reference_event_label = displacement_variables{i_variable, strcmp(displacement_variables_header, 'reference event label')};
                    this_reference_event_index = str2num(displacement_variables{i_variable, strcmp(displacement_variables_header, 'reference event index')});

                    % find reference event in current list
                    reference_event_indices_in_list = find(strcmp(loaded_data.event_labels, this_reference_event_label));
                    reference_event_index_in_list = reference_event_indices_in_list(this_reference_event_index);                

                    % extract variable to analyze
                    this_variable_label = displacement_variables{i_variable, 1};
                    this_variable_source = displacement_variables{i_variable, 2};
                    this_variable_split = strsplit(this_variable_source, ':');
                    this_variable_source_type = this_variable_split{1};
                    this_variable_source_label = this_variable_split{2};
                    eval(['source_data = ' this_variable_source_type '_data_to_analyze;']);
                    eval(['source_labels = ' this_variable_source_type '_labels;']);
                    
                    % extract reference
                    this_variable_reference_data = source_data{reference_event_index_in_list}(strcmp(source_labels, this_variable_source_label));

                    % calculate displacement from state at reference event for each event
                    for i_event = 1 : number_of_events
                        this_variable_raw_data = source_data{i_event}(:, strcmp(source_labels, this_variable_source_label));
                        this_variable_displacement_data = this_variable_raw_data - this_variable_reference_data;

                        % store
                        displacement_data_trial{i_variable} = [displacement_data_trial{i_variable} this_variable_displacement_data];
                    end

                    % add trial data to session data for this variable
                    displacement_data_session{i_variable} = [displacement_data_session{i_variable} displacement_data_trial{i_variable}];
                end            
                
                % create condition data for each new event
                subject_list_trial = cell(number_of_events, 1);
                [subject_list_trial{:, :}] = deal(subject_id);
                condition_list_trial = cell(number_of_events, 1);
                [condition_list_trial{:, :}] = deal(this_type_label);
                time_point_list_trial = expected_event_labels';
                group_list_trial = cell(number_of_events, 1);
                [group_list_trial{:, :}] = deal(group_label);
                this_block_label = block_assignment_this_type{strcmp(block_assignment_this_type(:, 3), num2str(this_type_trial_numbers(i_trial))), 2};
                block_list_trial = cell(number_of_events, 1);
                [block_list_trial{:, :}] = deal(this_block_label);
                origin_trial_list_trial = ones(number_of_events, 1) * this_type_trial_numbers(i_trial);
                
                % add condition data
                subject_list_session = [subject_list_session; subject_list_trial]; %#ok<AGROW>
                condition_list_session = [condition_list_session; condition_list_trial]; %#ok<AGROW>
                time_point_list_session = [time_point_list_session; time_point_list_trial]; %#ok<AGROW>
                group_list_session = [group_list_session; group_list_trial]; %#ok<AGROW>
                block_list_session = [block_list_session; block_list_trial]; %#ok<AGROW>
                origin_trial_list_session = [origin_trial_list_session; origin_trial_list_trial]; %#ok<AGROW>
            end            
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
    results_file_name = makeFileName(date, subject_id, 'resultsDisplacement');
    save ...
      ( ...
        results_file_name, ...
        'conditions_session', ...
        'displacement_data_session', ...
        'displacement_names_session', ...
        'displacement_directions_session', ...
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

