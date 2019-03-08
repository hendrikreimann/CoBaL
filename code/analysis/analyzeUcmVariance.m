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

function analyzeUcmVariance(varargin)
    [condition_list_session, trial_number_list] = parseTrialArguments(varargin{:});
    trial_type = condition_list_session{1}; % we assume that we have a single type only for now
    load('subjectInfo.mat', 'date', 'subject_id', 'condition_list', 'trial_number_list');
    
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
    load('subjectModel.mat');
    
    across_time_conditions = study_settings.get('across_time_conditions');
    across_events_conditions = study_settings.get('across_events_conditions');
    across_trials_conditions = study_settings.get('across_trials_conditions');
    ucm_variables = study_settings.get('ucm_variables');
    number_of_ucm_variables = length(ucm_variables);
    
    % make containers to hold the data
    subject_list_session = {};
    condition_list_session = {};
    time_point_list_session = {};
    group_list_session = {};
    origin_trial_list_session = [];
    stretch_data_session = cell(number_of_ucm_variables, 1);
    stretch_directions_session = cell(number_of_ucm_variables, 2);
    [stretch_directions_session{:, :}] = deal('~');
    stretch_names_session = ucm_variables;
    group_label = subject_settings.get('group');
    
    % calculate UCM variance measures
    for i_block = 1 : number_of_blocks
        this_block_label = block_labels{i_block};
        this_block_trials = trials_by_block{i_block};
        
        if any(strcmp(this_block_label, across_time_conditions))
            % analyze across time
            number_of_trials_this_block = length(this_block_trials);
            expected_event_labels = {'trial_start_time';'trial_end_time'};
            number_of_events = length(expected_event_labels);

            for i_trial = 1 : number_of_trials_this_block
                % load data
                loaded_data = load(['processed' filesep makeFileName(date, subject_id, this_block_label, this_block_trials(i_trial), 'kinematicTrajectories.mat')]);
                joint_angle_trajectories = loaded_data.joint_angle_trajectories;
                time_mocap = loaded_data.time_mocap;
                
                % load and check events
                loaded_data = load(['analysis' filesep makeFileName(date, subject_id, this_block_label, this_block_trials(i_trial), 'events.mat')]);
                if ~(length(loaded_data.event_data)==number_of_events)
                    error(['Trial ' num2str(this_block_trials(i_trial)) ' - expected ' num2str(number_of_events) ' events, but found ' num2str(length(event_indices_mocap))]);
                end
                
                % extract events
                event_times = zeros(1, number_of_events);
                event_indices_mocap = zeros(1, number_of_events);
                for i_event = 1 : number_of_events
                    if ~(any(strcmp(expected_event_labels{i_event}, loaded_data.event_labels)))
                        error(['Trial ' num2str(this_block_trials(i_trial)) ' - event label "' expected_event_labels{i_event} '" not found']);
                    end
                    event_times(i_event) = loaded_data.event_data{strcmp(loaded_data.event_labels, expected_event_labels{i_event})};
                    event_indices_mocap(i_event) = findClosestIndex(event_times(i_event), time_mocap);
                end
                % extract data
                joint_angle_data_to_analyze = joint_angle_trajectories(event_indices_mocap(1) : event_indices_mocap(2), :)';

                % analyze data
                theta_mean = mean(joint_angle_data_to_analyze, 2);
                kinematic_tree.jointAngles = theta_mean;
                kinematic_tree.updateInternals;
                for i_variable = 1 : number_of_ucm_variables
                    % calculate Jacobian
                    if strcmp(ucm_variables{i_variable}, 'com_ap')
                        J_com = kinematic_tree.calculateCenterOfMassJacobian;
                        jacobian = J_com(1, :);
                    end
                    if strcmp(ucm_variables{i_variable}, 'com_vert')
                        J_com = kinematic_tree.calculateCenterOfMassJacobian;
                        jacobian = J_com(3, :);
                    end
                    if strcmp(ucm_variables{i_variable}, 'com_2d')
                        J_com = kinematic_tree.calculateCenterOfMassJacobian;
                        jacobian = J_com([1 3], :);
                    end

                    % calculate variance measures
                    [V_para_this_block_this_variable, V_perp_this_block_this_variable] = calculateUcmVariance(joint_angle_data_to_analyze, jacobian);
                    stretch_data_session{i_variable} = [stretch_data_session{i_variable} [V_para_this_block_this_variable; V_perp_this_block_this_variable]];
                end
                subject_list = [subject_list; subject_id]; %#ok<AGROW>
                time_point_list = [time_point_list; 'quiet stance across time']; %#ok<AGROW>
                condition_list = [condition_list; this_block_label]; %#ok<AGROW>
                origin_trial_list_session = [origin_trial_list_session; this_block_trials(1)]; %#ok<AGROW>                            
                    
                
            end
            
            
            
            
            
            
            
        end
        if any(strcmp(this_block_label, across_events_conditions))
            % analyze across events
            number_of_trials_this_block = length(this_block_trials);
            expected_event_labels = {'oscillation_peaks';'oscillation_vales'};
            number_of_events = length(expected_event_labels);
            joint_angle_data_to_analyze = cell(number_of_events, 1);
            
            for i_trial = 1 : number_of_trials_this_block
                % load data
                loaded_data = load(['processed' filesep makeFileName(date, subject_id, trial_type, this_block_trials(i_trial), 'kinematicTrajectories.mat')]);
                joint_angle_trajectories = loaded_data.joint_angle_trajectories;
                time_mocap = loaded_data.time_mocap;
                
                % load and check events
                loaded_data = load(['analysis' filesep makeFileName(date, subject_id, trial_type, this_block_trials(i_trial), 'events.mat')]);
                if ~(length(loaded_data.event_data)==number_of_events)
                    error(['Trial ' num2str(this_block_trials(i_trial)) ' - expected ' num2str(number_of_events) ' events, but found ' num2str(length(event_indices_mocap))]);
                end
                
                % extract events
                event_times = cell(1, number_of_events);
                event_indices_mocap = cell(1, number_of_events);
                for i_event = 1 : number_of_events
                    if ~(any(strcmp(expected_event_labels{i_event}, loaded_data.event_labels)))
                        error(['Trial ' num2str(this_block_trials(i_trial)) ' - event label "' expected_event_labels{i_event} '" not found']);
                    end
                    event_times{i_event} = loaded_data.event_data{strcmp(loaded_data.event_labels, expected_event_labels{i_event})};
                    event_indices_mocap{i_event} = findClosestIndex(event_times{i_event}, time_mocap);
                
                    % extract data
                    joint_angle_data_to_analyze{i_event} = joint_angle_trajectories(event_indices_mocap{i_event}, :)';
                    
                    % analyze data
                    joint_angle_data_to_analyze_this_event = joint_angle_trajectories(event_indices_mocap{i_event}, :)';
                    theta_mean = mean(joint_angle_data_to_analyze_this_event);
                    kinematic_tree.jointAngles = theta_mean;
                    kinematic_tree.updateConfiguration;
                    for i_variable = 1 : number_of_ucm_variables
                        % calculate Jacobian
                        if strcmp(ucm_variables{i_variable}, 'com_ap')
                            J_com = kinematic_tree.calculateCenterOfMassJacobian;
                            jacobian = J_com(1, :);
                        end
                        if strcmp(ucm_variables{i_variable}, 'com_vert')
                            J_com = kinematic_tree.calculateCenterOfMassJacobian;
                            jacobian = J_com(3, :);
                        end
                        if strcmp(ucm_variables{i_variable}, 'com_2d')
                            J_com = kinematic_tree.calculateCenterOfMassJacobian;
                            jacobian = J_com([1 3], :);
                        end

                        % calculate variance measures
                        [V_para_this_block_this_variable, V_perp_this_block_this_variable] = calculateUcmVariance(joint_angle_data_to_analyze_this_event, jacobian);
                        stretch_data_session{i_variable} = [stretch_data_session{i_variable} [V_para_this_block_this_variable; V_perp_this_block_this_variable]];
                    end
                    subject_list = [subject_list; subject_id]; %#ok<AGROW>
                    time_point_list = [time_point_list; expected_event_labels{i_event}]; %#ok<AGROW>
                    condition_list = [condition_list; this_block_label]; %#ok<AGROW>
                    origin_trial_list_session = [origin_trial_list_session; this_block_trials(1)]; %#ok<AGROW>                            
                    
                end
                
        
                
                
                
                
                
                
                
                
                
                
                
%                 % load indices
%                 loaded_data = load(['analysis' filesep makeFileName(date, subject_id, trial_type, this_block_trials(i_trial), 'events.mat')]);
%                 event_times = loaded_data.event_data{strcmp(loaded_data.event_labels, 'oscillation_peaks')};
%                 if ~isempty(event_times)
%                 
%                     event_indices_mocap = findClosestIndex(event_times, time_mocap);
% 
%                     % calculate variance measures
%                     joint_angle_data_to_analyze = joint_angle_trajectories(event_indices_mocap, :);
%                     theta_mean = mean(joint_angle_data_to_analyze)';
%                     kinematic_tree.jointAngles = theta_mean;
%                     kinematic_tree.updateConfiguration;
%                     for i_variable = 1 : number_of_ucm_variables
%                         % calculate Jacobian
%                         if strcmp(ucm_variables{i_variable}, 'com_ap')
%                             J_com = kinematic_tree.calculateCenterOfMassJacobian;
%                             jacobian = J_com(1, :);
%                         end
%                         if strcmp(ucm_variables{i_variable}, 'com_vert')
%                             J_com = kinematic_tree.calculateCenterOfMassJacobian;
%                             jacobian = J_com(3, :);
%                         end
%                         if strcmp(ucm_variables{i_variable}, 'com_2d')
%                             J_com = kinematic_tree.calculateCenterOfMassJacobian;
%                             jacobian = J_com([1 3], :);
%                         end
% 
%                         % calculate variance measures
%                         [V_para, V_perp] = calculateUcmVariance(joint_angle_data_to_analyze', jacobian);
%     %                     V_para_trial(i_variable, i_trial) = V_para;
%     %                     V_perp_trial(i_variable, i_trial) = V_perp;
%                         stretch_data_session{i_variable} = [stretch_data_session{i_variable} [V_para; V_perp]];
%                     end
%                     subject_list = [subject_list; subject_id]; %#ok<AGROW>
%                     time_point_list = [time_point_list; 'NA']; %#ok<AGROW>
%                     condition_list = [condition_list; this_block_label]; %#ok<AGROW>
%                     origin_trial_list_session = [origin_trial_list_session; this_block_trials(i_trial)];
%                 end
            end
        end
        if any(strcmp(this_block_label, across_trials_conditions))
            % analyze across trials
            number_of_trials_this_block = length(this_block_trials);
            expected_event_labels = {'perturbation_start';'perturbation_end';'perturbation_end_plus_one';'perturbation_end_plus_two'};
            number_of_events = length(expected_event_labels);
            joint_angle_data_to_analyze = cell(number_of_events, 1);
            for i_trial = 1 : number_of_trials_this_block
                % load data
                loaded_data = load(['processed' filesep makeFileName(date, subject_id, trial_type, this_block_trials(i_trial), 'kinematicTrajectories.mat')]);
                joint_angle_trajectories = loaded_data.joint_angle_trajectories;
                time_mocap = loaded_data.time_mocap;
                
                % load and check events
                loaded_data = load(['analysis' filesep makeFileName(date, subject_id, trial_type, this_block_trials(i_trial), 'events.mat')]);
                if ~(length(loaded_data.event_data)==number_of_events)
                    error(['Trial ' num2str(this_block_trials(i_trial)) ' - expected ' num2str(number_of_events) ' events, but found ' num2str(length(event_indices_mocap))]);
                end
                
                % extract events
                event_times = zeros(1, number_of_events);
                for i_event = 1 : number_of_events
                    if ~(any(strcmp(expected_event_labels{i_event}, loaded_data.event_labels)))
                        error(['Trial ' num2str(this_block_trials(i_trial)) ' - event label "' expected_event_labels{i_event} '" not found']);
                    end
                    event_times(i_event) = loaded_data.event_data{strcmp(loaded_data.event_labels, expected_event_labels{i_event})};
                end
                event_indices_mocap = findClosestIndex(event_times, time_mocap);

                % extract data
                for i_event = 1 : number_of_events
                    joint_angle_data_this_event = joint_angle_trajectories(event_indices_mocap(i_event), :)';
                    joint_angle_data_to_analyze{i_event} = [joint_angle_data_to_analyze{i_event} joint_angle_data_this_event];
                end
            end
            
            % calculate and store variance measures
            for i_event = 1 : number_of_events
                joint_angle_data_to_analyze_this_event = joint_angle_data_to_analyze{i_event};
                theta_mean = mean(joint_angle_data_to_analyze_this_event);
                kinematic_tree.jointAngles = theta_mean;
                kinematic_tree.updateConfiguration;
                for i_variable = 1 : number_of_ucm_variables
                    % calculate Jacobian
                    if strcmp(ucm_variables{i_variable}, 'com_ap')
                        J_com = kinematic_tree.calculateCenterOfMassJacobian;
                        jacobian = J_com(1, :);
                    end
                    if strcmp(ucm_variables{i_variable}, 'com_vert')
                        J_com = kinematic_tree.calculateCenterOfMassJacobian;
                        jacobian = J_com(3, :);
                    end
                    if strcmp(ucm_variables{i_variable}, 'com_2d')
                        J_com = kinematic_tree.calculateCenterOfMassJacobian;
                        jacobian = J_com([1 3], :);
                    end
                    
                    % calculate variance measures
                    [V_para_this_block_this_variable, V_perp_this_block_this_variable] = calculateUcmVariance(joint_angle_data_to_analyze_this_event, jacobian);
                    stretch_data_session{i_variable} = [stretch_data_session{i_variable} [V_para_this_block_this_variable; V_perp_this_block_this_variable]];
                end
                subject_list = [subject_list; subject_id]; %#ok<AGROW>
                time_point_list = [time_point_list; expected_event_labels{i_event}]; %#ok<AGROW>
                condition_list = [condition_list; this_block_label]; %#ok<AGROW>
                origin_trial_list_session = [origin_trial_list_session; this_block_trials(1)]; %#ok<AGROW>
            end
            
        end
        
        
    end
    
    number_of_stretches = length(subject_list);
    origin_start_time_list_session = zeros(number_of_stretches, 1); % doesn't apply, but needs to be here for now
    origin_end_time_list_session = zeros(number_of_stretches, 1); % doesn't apply, but needs to be here for now
    time_list_session = zeros(number_of_stretches, 1); % doesn't apply, but needs to be here for now
    

    %% save data
    bands_per_stretch = 2;
    conditions_session = struct;
    conditions_session.subject_list = subject_list;
    conditions_session.condition_list = condition_list;
    conditions_session.time_point_list = time_point_list;
    results_file_name = ['analysis' filesep makeFileName(date, subject_id, 'results')];
    save ...
      ( ...
        results_file_name, ...
        'conditions_session', ...
        'stretch_data_session', ...
        'stretch_names_session', ...
        'stretch_directions_session', ...
        'bands_per_stretch', ...
        'origin_trial_list_session', ...
        'origin_start_time_list_session', ...
        'origin_end_time_list_session', ...
        'time_list_session' ...
      )
    
    
    
    % save data
%     variable_names_session = ucm_variables; %#ok<NASGU>
%     results_file_name = ['analysis' filesep makeFileName(date, subject_id, 'results')];
%     save ...
%       ( ...
%         results_file_name, ...
%         'V_para_session', ...
%         'V_perp_session', ...
%         'variable_names_session', ...
%         'block_labels' ...
%       )
end

          

