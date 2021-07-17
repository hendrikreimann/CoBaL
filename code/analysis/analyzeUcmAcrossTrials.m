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

function analyzeUcmAcrossTrials(varargin)
    [trial_type_list, trial_number_list] = parseTrialArguments(varargin{:});

    % parse arguments
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'block_assignment_suffix', '')
    parse(parser, varargin{:})
    block_assignment_suffix = parser.Results.block_assignment_suffix;
    
    % load settings
    study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
    model_data = load('subjectModel.mat');
    kinematic_tree = model_data.kinematic_tree;
    
    across_trials_conditions = study_settings.get('across_trials_conditions', 1);
    task_variables = study_settings.get('task_variables');
    
    if isempty(block_assignment_suffix)
        block_assignment = subject_settings.getTable('block_assignment');
    else
        block_assignment = subject_settings.getTable(['block_assignment_' block_assignment_suffix]);
    end
        
        
    number_of_task_variables = length(task_variables);
    group_label = subject_settings.get('group');
    
    % make containers to hold the data
    subject_list_session = {};
    condition_list_session = {};
    time_point_list_session = {};
    group_list_session = {};
    event_type_list_session = {};
    block_list_session = {};
    origin_trial_list_session = [];
    
    % define stretch variables from list of ucm and variance variables
    bands_per_stretch = 2;
    number_of_ucm_variables = number_of_task_variables * 6;
    ucm_data_session = cell(number_of_ucm_variables, 1);
    ucm_directions_session = cell(number_of_ucm_variables, 2);
    [ucm_directions_session{:, :}] = deal('~');
    ucm_names_session = {};
    for i_variable = 1 : number_of_task_variables
        this_task_variable_name = task_variables{i_variable};
        this_variable_ucm_names = ...
          { ...
            this_task_variable_name; ...
            [this_task_variable_name '_rand_ankle']; ...
            [this_task_variable_name '_rand_knee']; ...
            [this_task_variable_name '_rand_hip']; ...
            [this_task_variable_name '_rand_neck']; ...
            [this_task_variable_name '_rand_all']; ...
          };
        ucm_names_session = [ucm_names_session; this_variable_ucm_names]; %#ok<AGROW>
    end
    
    % parse event labels
    expected_event_labels = study_settings.get('event_labels_sway');
    number_of_events = length(expected_event_labels);
    event_type_list = cell(number_of_events, 1);
    time_number_list = cell(number_of_events, 1);
    for i_event = 1 : number_of_events
        this_event_label_split = strsplit(expected_event_labels{i_event}, '_');
        event_type_list{i_event} = this_event_label_split{1};
        time_number_list{i_event} = this_event_label_split{2};
    end
            
    % calculate variance across trials for specific events
    for i_condition = 1 : length(across_trials_conditions)
        % get list of trials for this condition
        this_condition_label = across_trials_conditions{i_condition};
        
        % get block information
%         block_assignment_this_condition = block_assignment(strcmp(block_assignment(:, 1), this_condition_label), :);
%         block_labels = unique(block_assignment_this_condition(:, 2));

        block_assignment_this_condition = block_assignment(strcmp(block_assignment.trial_type, this_condition_label), :);
        block_labels = unique(block_assignment_this_condition.block_label);
        
        
        if any(strcmp(trial_type_list, this_condition_label))
%             this_condition_trial_numbers = trial_number_list{strcmp(trial_type_list, this_condition_label)};
            this_condition_trial_numbers = str2double(block_assignment_this_condition.trial_number);

            % prepare containers
            number_of_trials_this_condition = length(this_condition_trial_numbers);
            joint_data_to_analyze = cell(number_of_events, 1);
            com_data_to_analyze = cell(number_of_events, 1);
            eef_data_to_analyze = cell(number_of_events, 1);
            origin_trial_numbers = [];
            
            % collect data from different trials
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
                if length(loaded_data.event_data) < number_of_events
                    error(['Trial ' num2str(this_condition_trial_numbers(i_trial)) ' - expected ' num2str(number_of_events) ' events, but found ' num2str(length(loaded_data.event_data))]);
                end

                % extract event data
                event_times = zeros(1, number_of_events);
                for i_event = 1 : number_of_events
                    if ~(any(strcmp(expected_event_labels{i_event}, loaded_data.event_labels)))
                        error(['Trial ' num2str(this_condition_trial_numbers(i_trial)) ' - event label "' expected_event_labels{i_event} '" not found']);
                    end
                    event_times(i_event) = loaded_data.event_data{strcmp(loaded_data.event_labels, expected_event_labels{i_event})};
                end
                event_indices_mocap = findClosestIndex(event_times, time_mocap);

                % extract data at time points where events happened
                for i_event = 1 : number_of_events
                    joint_angle_data_this_event = joint_angle_trajectories(event_indices_mocap(i_event), :);
                    joint_data_to_analyze{i_event} = [joint_data_to_analyze{i_event}; joint_angle_data_this_event];
                    com_data_this_event = com_trajectories(event_indices_mocap(i_event), :);
                    com_data_to_analyze{i_event} = [com_data_to_analyze{i_event}; com_data_this_event];
                    eef_data_this_event = eef_trajectories(event_indices_mocap(i_event), :);
                    eef_data_to_analyze{i_event} = [eef_data_to_analyze{i_event}; eef_data_this_event];
                    
                end
                origin_trial_numbers = [origin_trial_numbers; this_condition_trial_numbers(i_trial)]; %#ok<AGROW>
            end

            % calculate and store variance measures
            for i_block = 1 : length(block_labels)
                this_block_label = block_labels(i_block);
                this_block_trials_numbers = str2double(block_assignment_this_condition.trial_number(strcmp(block_assignment_this_condition.block_label, this_block_label)));
                number_of_trials_this_block = length(this_block_trials_numbers);
                
                trial_indices_this_block = ismember(origin_trial_numbers, this_block_trials_numbers);
                
                if number_of_trials_this_block > 0
                    for i_event = 1 : number_of_events
                        % get collected data for this event
                        joint_angle_data_to_analyze_this_event = joint_data_to_analyze{i_event}(trial_indices_this_block, :);
                        joint_angle_data_to_analyze_this_event_rand = randmat(joint_angle_data_to_analyze_this_event, 1);
                        com_data_to_analyze_this_event = com_data_to_analyze{i_event}; %#ok<NASGU>
                        eef_data_to_analyze_this_event = eef_data_to_analyze{i_event}; %#ok<NASGU>

                        % analyze joint angle variance in UCM space
                        theta_mean = mean(joint_angle_data_to_analyze_this_event, 1)';
                        kinematic_tree.jointAngles = theta_mean;
                        kinematic_tree.updateConfiguration;
                        for i_variable = 1 : number_of_task_variables
                            this_variable_label = task_variables{i_variable};
                            jacobian = determineJacobian(kinematic_tree, this_variable_label);

                            % calculate variance measures
                            [V_para_this_event_this_variable, V_perp_this_event_this_variable] = calculateUcmVariance(joint_angle_data_to_analyze_this_event', jacobian);
                            [V_para_rand_1_this_event_this_variable, V_perp_rand_1_this_event_this_variable] = calculateUcmVarianceWithRandomization(joint_angle_data_to_analyze_this_event_rand', jacobian, 1);
                            [V_para_rand_2_this_event_this_variable, V_perp_rand_2_this_event_this_variable] = calculateUcmVarianceWithRandomization(joint_angle_data_to_analyze_this_event_rand', jacobian, 2);
                            [V_para_rand_3_this_event_this_variable, V_perp_rand_3_this_event_this_variable] = calculateUcmVarianceWithRandomization(joint_angle_data_to_analyze_this_event_rand', jacobian, 3);
                            [V_para_rand_4_this_event_this_variable, V_perp_rand_4_this_event_this_variable] = calculateUcmVarianceWithRandomization(joint_angle_data_to_analyze_this_event_rand', jacobian, 4);
                            [V_para_rand_all_this_event_this_variable, V_perp_rand_all_this_event_this_variable] = calculateUcmVarianceWithRandomization(joint_angle_data_to_analyze_this_event_rand', jacobian, 'all');

                            % store calculated measures
                            ucm_data_session{strcmp(ucm_names_session, this_variable_label)} = [ucm_data_session{strcmp(ucm_names_session, this_variable_label)} [V_para_this_event_this_variable; V_perp_this_event_this_variable]];
                            ucm_data_session{strcmp(ucm_names_session, [this_variable_label '_rand_ankle'])} = [ucm_data_session{strcmp(ucm_names_session, [this_variable_label '_rand_ankle'])} [V_para_rand_1_this_event_this_variable; V_perp_rand_1_this_event_this_variable]];
                            ucm_data_session{strcmp(ucm_names_session, [this_variable_label '_rand_knee'])} = [ucm_data_session{strcmp(ucm_names_session, [this_variable_label '_rand_knee'])} [V_para_rand_2_this_event_this_variable; V_perp_rand_2_this_event_this_variable]];
                            ucm_data_session{strcmp(ucm_names_session, [this_variable_label '_rand_hip'])} = [ucm_data_session{strcmp(ucm_names_session, [this_variable_label '_rand_hip'])} [V_para_rand_3_this_event_this_variable; V_perp_rand_3_this_event_this_variable]];
                            ucm_data_session{strcmp(ucm_names_session, [this_variable_label '_rand_neck'])} = [ucm_data_session{strcmp(ucm_names_session, [this_variable_label '_rand_neck'])} [V_para_rand_4_this_event_this_variable; V_perp_rand_4_this_event_this_variable]];
                            ucm_data_session{strcmp(ucm_names_session, [this_variable_label '_rand_all'])} = [ucm_data_session{strcmp(ucm_names_session, [this_variable_label '_rand_all'])} [V_para_rand_all_this_event_this_variable; V_perp_rand_all_this_event_this_variable]];
                        end

                        % store supplementary information
                        subject_list_session = [subject_list_session; subject_id]; %#ok<AGROW>
                        time_point_list_session = [time_point_list_session; time_number_list{i_event}]; %#ok<AGROW>
                        event_type_list_session = [event_type_list_session; event_type_list{i_event}]; %#ok<AGROW>
                        condition_list_session = [condition_list_session; this_condition_label]; %#ok<AGROW>
                        group_list_session = [group_list_session; group_label]; %#ok<AGROW>
                        block_list_session = [block_list_session; this_block_label]; %#ok<AGROW>
                        origin_trial_list_session = [origin_trial_list_session; this_condition_trial_numbers(1)]; %#ok<AGROW>
                    end
                end
            end
        end
    end

    % determine meta information
    number_of_stretches = length(subject_list_session);
    origin_start_time_list_session = zeros(number_of_stretches, 1); % doesn't apply, but needs to be here for now
    origin_end_time_list_session = zeros(number_of_stretches, 1); % doesn't apply, but needs to be here for now
    time_list_session = zeros(number_of_stretches, 1); % doesn't apply, but needs to be here for now
    
    %% save data
    conditions_session = struct;
    conditions_session.subject_list = subject_list_session;
    conditions_session.condition_list = condition_list_session;
    conditions_session.time_point_list = time_point_list_session;
    conditions_session.event_type_list = event_type_list_session;
    conditions_session.group_list = group_list_session;
    conditions_session.block_list = block_list_session;
    if ~directoryExists('results')
        mkdir('results')
    end
    
    file_label = 'resultsUcmAcrossTrials';
    if ~isempty(block_assignment_suffix)
        file_label = [file_label '_' block_assignment_suffix];
    end
    results_file_name = ['results' filesep makeFileName(collection_date, subject_id, file_label)];
    save ...
      ( ...
        results_file_name, ...
        'conditions_session', ...
        'ucm_data_session', ...
        'ucm_names_session', ...
        'ucm_directions_session', ...
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

