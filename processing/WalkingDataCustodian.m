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

% input

% description

% there are different levels of variables that this object takes care of
% - basic variables can be loaded or calculated directly from loaded variables at each point in time
% - stretch variables are variables defined for each stretch. That can be a single number or a trajectory with
% normalized time. 

% order of processing is
% Construction:
% - variables are determined based on the list of variables_to_analyze, loaded from the study_settings
% For each trial
% 1. Prepare data: loads and calculates basic variables
% 2. <apply stretches and normalize in time>
% 3. calculate stretch variables

% to add a new variable, you need to do the following things:
%
% 1. Provide information about all variables this new one depends upon. Go to the determineVariables function and add a
% new conditional statement "if this.isVariableToAnalyze('<your new variable name>') ... end". Within the
% statement, add the basic variables and stretch variables that the new variable depends upon by calling addBasicVariable 
% and addStretchVariable. The order matters here, make sure for each variable, all dependencies have already been added, 
% i.e. don't call "this.addBasicVariable('lheel_y_pos')" before "this.addBasicVariable('marker_trajectories')", because 
% lheel_y_pos depends upon the marker_trajectories.
%
% 2. If you added a new basic variable, but is not already available, i.e. one of the outputs of the
% preprocess scripts, you need to specify how to calculate it. Go to prepareBasicVariables and add a new conditional
% statement "if strcmp(variable_name, '<your new variable name>')". Within the statement, calculate your new
% variable, then add it to the variable_data. Also add the corresponding time vector to the time_data.
%
% 3. If the variable you want to add is a stretch variable, you need to specify how to calculate it. Go to 
% calculateStretchVariables and add a new conditional statement "if strcmp(variable_name, '<your new variable name>')".
% Within the statement, calculate the value of your new variable for a single stretch, depending upon
% this_stretch_start_time and this_stretch_end_time. Assign that value to the variable "stretch_data".


classdef WalkingDataCustodian < handle
    properties
        date = [];
        subject_id = [];
        variables_to_analyze = {};
        basic_variable_names = {};
        stretch_variable_names = {};
        
        % these are private and should only be accessed using get functions
        number_of_time_steps_normalized;
        basic_variable_data;
        basic_variable_labels;
        stretch_variable_data;
        time_data;
        
        subject_settings;
        emg_normalization_values;
        emg_normalization_labels;
    end
    
    methods
        % constructor
%         function this = WalkingDataCustodian(date, subject_id, variables_to_analyze)
        function this = WalkingDataCustodian()
            % load this information from the subjects.mat and studySettings.txt files
            load('subjectInfo.mat', 'date', 'subject_id');
            this.subject_settings = loadSettingsFile('subjectSettings.txt');
            % load settings
            study_settings_file = '';
            if exist(['..' filesep 'studySettings.txt'], 'file')
                study_settings_file = ['..' filesep 'studySettings.txt'];
            end    
            if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
                study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
            end
            study_settings = loadSettingsFile(study_settings_file);
            emg_normalization_file_name = ['analysis' filesep makeFileName(date, subject_id, 'emgNormalization.mat')];
            if exist(emg_normalization_file_name, 'file')
                emg_normalization_data = load(emg_normalization_file_name);
                this.emg_normalization_values = emg_normalization_data.emg_normalization_values;
                this.emg_normalization_labels = emg_normalization_data.emg_variables;
            end
            
            this.date = date;
            this.subject_id = subject_id;
            this.variables_to_analyze = study_settings.variables_to_analyze;
            this.number_of_time_steps_normalized = study_settings.number_of_time_steps_normalized;
            
            this.determineVariables();
        end
        
        % initialization
        function addBasicVariable(this, variable_name)
            if ~any(strcmp(this.basic_variable_names, variable_name))
                this.basic_variable_names = [this.basic_variable_names; variable_name];
            end
        end
        function addStretchVariable(this, variable_name)
            if ~any(strcmp(this.stretch_variable_names, variable_name))
                this.stretch_variable_names = [this.stretch_variable_names; variable_name];
            end
            
            % TODO: make sure that after I added this one, all the variables are in the right order, so that
            % variables that are needed later are calculated first, e.g. c7 pos and MPSIS pos are calculated before
            % trunk angle
            % ... actually, this should not be necessary, as long as for each variable, all required variables are added
            % in the right order
            
        end
        function determineVariables(this)
            % for each possible variable to analyze, list the basic and required variables required to calculate it
            if this.isVariableToAnalyze('lheel_x_pos')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lheel_x_pos')
                this.addStretchVariable('lheel_x_pos')
            end
            if this.isVariableToAnalyze('rheel_x_pos')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('rheel_x_pos')
                this.addStretchVariable('rheel_x_pos')
            end
            if this.isVariableToAnalyze('step_length')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lheel_y_pos')
                this.addBasicVariable('rheel_y_pos')
                this.addStretchVariable('lheel_y_pos')
                this.addStretchVariable('rheel_y_pos')
                this.addStretchVariable('step_length')
            end
            if this.isVariableToAnalyze('step_width')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lheel_x_pos')
                this.addBasicVariable('rheel_x_pos')
                this.addStretchVariable('lheel_x_pos')
                this.addStretchVariable('rheel_x_pos')
                this.addStretchVariable('step_width')
            end
            if this.isVariableToAnalyze('step_time')
                this.addStretchVariable('step_time')
            end
            if this.isVariableToAnalyze('cadence')
                this.addStretchVariable('step_time')
                this.addStretchVariable('cadence')
            end
            if this.isVariableToAnalyze('velocity')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lheel_y_pos')
                this.addBasicVariable('rheel_y_pos')
                this.addStretchVariable('lheel_y_pos')
                this.addStretchVariable('rheel_y_pos')
                this.addStretchVariable('step_length')
                this.addStretchVariable('step_time')
                this.addStretchVariable('velocity')
            end
            if this.isVariableToAnalyze('trunk_angle_ap')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('trunk_angle_ap')
                this.addStretchVariable('trunk_angle_ap')
            end            
            if this.isVariableToAnalyze('trunk_angle_ml')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('trunk_angle_ml')
                this.addStretchVariable('trunk_angle_ml')
            end            
            if this.isVariableToAnalyze('copl_ap')
                this.addBasicVariable('left_foot_cop_world')
                this.addBasicVariable('copl_ap')
                this.addStretchVariable('copl_ap')
            end
            if this.isVariableToAnalyze('copl_ml')
                this.addBasicVariable('left_foot_cop_world')
                this.addBasicVariable('copl_ml')
                this.addStretchVariable('copl_ml')
            end
            if this.isVariableToAnalyze('copr_ap')
                this.addBasicVariable('right_foot_cop_world')
                this.addBasicVariable('copr_ap')
                this.addStretchVariable('copr_ap')
            end
            if this.isVariableToAnalyze('copr_ml')
                this.addBasicVariable('right_foot_cop_world')
                this.addBasicVariable('copr_ml')
                this.addStretchVariable('copr_ml')
            end
            if this.isVariableToAnalyze('cop_ap')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_ap')
                this.addStretchVariable('cop_ap')
            end
            if this.isVariableToAnalyze('body_com_to_cop_x')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_ml')
                this.addStretchVariable('cop_ml')
                this.addBasicVariable('com_trajectories')
                this.addBasicVariable('body_com_x')
                this.addStretchVariable('body_com_x')
                this.addStretchVariable('body_com_to_cop_x')
            end
            if this.isVariableToAnalyze('body_com_to_cop_y')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_ap')
                this.addStretchVariable('cop_ap')
                this.addBasicVariable('com_trajectories')
                this.addBasicVariable('body_com_y')
                this.addStretchVariable('body_com_y')
                this.addStretchVariable('body_com_to_cop_y')
            end
            if this.isVariableToAnalyze('cop_ml')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_ml')
                this.addStretchVariable('cop_ml')
            end
            if this.isVariableToAnalyze('left_arm_angle')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('left_arm_angle')
                this.addStretchVariable('left_arm_angle')
            end
            if this.isVariableToAnalyze('right_arm_angle')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('right_arm_angle')
                this.addStretchVariable('right_arm_angle')
            end
            if this.isVariableToAnalyze('left_leg_angle')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('left_leg_angle')
                this.addStretchVariable('left_leg_angle')
            end
            if this.isVariableToAnalyze('right_leg_angle')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('right_leg_angle')
                this.addStretchVariable('right_leg_angle')
            end
            if this.isVariableToAnalyze('left_arm_phase')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('left_arm_angle')
                this.addBasicVariable('left_arm_phase')
                this.addStretchVariable('left_arm_angle')
                this.addStretchVariable('left_arm_phase')
            end
            if this.isVariableToAnalyze('right_arm_phase')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('right_arm_angle')
                this.addBasicVariable('right_arm_phase')
                this.addStretchVariable('right_arm_angle')
                this.addStretchVariable('right_arm_phase')
            end
            if this.isVariableToAnalyze('left_leg_phase')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('left_leg_angle')
                this.addBasicVariable('left_leg_phase')
                this.addStretchVariable('left_leg_angle')
                this.addStretchVariable('left_leg_phase')
            end
            if this.isVariableToAnalyze('right_leg_phase')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('right_leg_angle')
                this.addBasicVariable('right_leg_phase')
                this.addStretchVariable('right_leg_angle')
                this.addStretchVariable('right_leg_phase')
            end
            if this.isVariableToAnalyze('left_arm_right_leg_relative_phase')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('left_arm_angle')
                this.addBasicVariable('left_arm_phase')
                this.addBasicVariable('right_leg_angle')
                this.addBasicVariable('right_leg_phase')
                this.addBasicVariable('left_arm_right_leg_relative_phase')
                this.addStretchVariable('left_arm_angle')
                this.addStretchVariable('left_arm_phase')
                this.addStretchVariable('right_leg_angle')
                this.addStretchVariable('right_leg_phase')
                this.addStretchVariable('left_arm_right_leg_relative_phase')
            end
            if this.isVariableToAnalyze('right_arm_left_leg_relative_phase')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('right_arm_angle')
                this.addBasicVariable('right_arm_phase')
                this.addBasicVariable('left_leg_angle')
                this.addBasicVariable('left_leg_phase')
                this.addBasicVariable('right_arm_left_leg_relative_phase')
                this.addStretchVariable('right_arm_angle')
                this.addStretchVariable('right_arm_phase')
                this.addStretchVariable('left_leg_angle')
                this.addStretchVariable('left_leg_phase')
                this.addStretchVariable('right_arm_left_leg_relative_phase')
            end
            if this.isVariableToAnalyze('left_arm_phase_at_heelstrike')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('left_arm_angle')
                this.addBasicVariable('left_arm_phase')
                this.addStretchVariable('left_arm_phase_at_heelstrike')
            end
            if this.isVariableToAnalyze('right_arm_phase_at_heelstrike')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('right_arm_angle')
                this.addBasicVariable('right_arm_phase')
                this.addStretchVariable('right_arm_phase_at_heelstrike')
            end
            if this.isVariableToAnalyze('left_glut_med')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_glut_med')
                this.addStretchVariable('left_glut_med')
            end
            if this.isVariableToAnalyze('left_delt_ant')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_delt_ant')
                this.addStretchVariable('left_delt_ant')
            end
            if this.isVariableToAnalyze('left_gastroc_med')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_gastroc_med')
                this.addStretchVariable('left_gastroc_med')
            end
            if this.isVariableToAnalyze('left_pero_lng')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_pero_lng')
                this.addStretchVariable('left_pero_lng')
            end
            if this.isVariableToAnalyze('right_glut_med')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_glut_med')
                this.addStretchVariable('right_glut_med')
            end
            if this.isVariableToAnalyze('right_delt_ant')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_delt_ant')
                this.addStretchVariable('right_delt_ant')
            end
            if this.isVariableToAnalyze('right_gastroc_med')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_gastroc_med')
                this.addStretchVariable('right_gastroc_med')
            end
            if this.isVariableToAnalyze('right_pero_lng')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_pero_lng')
                this.addStretchVariable('right_pero_lng')
            end
            if this.isVariableToAnalyze('left_glut_med_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_glut_med')
                this.addStretchVariable('left_glut_med')
                this.addStretchVariable('left_glut_med_rescaled')
            end
            if this.isVariableToAnalyze('left_delt_ant_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_delt_ant')
                this.addStretchVariable('left_delt_ant')
                this.addStretchVariable('left_delt_ant_rescaled')
            end
            if this.isVariableToAnalyze('left_gastroc_med_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_gastroc_med')
                this.addStretchVariable('left_gastroc_med')
                this.addStretchVariable('left_gastroc_med_rescaled')
            end
            if this.isVariableToAnalyze('left_pero_lng_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_pero_lng')
                this.addStretchVariable('left_pero_lng')
                this.addStretchVariable('left_pero_lng_rescaled')
            end
            if this.isVariableToAnalyze('right_glut_med_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_glut_med')
                this.addStretchVariable('right_glut_med')
                this.addStretchVariable('right_glut_med_rescaled')
            end
            if this.isVariableToAnalyze('right_delt_ant_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_delt_ant')
                this.addStretchVariable('right_delt_ant')
                this.addStretchVariable('right_delt_ant_rescaled')
            end
            if this.isVariableToAnalyze('right_gastroc_med_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_gastroc_med')
                this.addStretchVariable('right_gastroc_med')
                this.addStretchVariable('right_gastroc_med_rescaled')
            end
            if this.isVariableToAnalyze('right_pero_lng_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_pero_lng')
                this.addStretchVariable('right_pero_lng')
                this.addStretchVariable('right_pero_lng_rescaled')
            end
            if this.isVariableToAnalyze('body_com_x')
                this.addBasicVariable('com_trajectories')
                this.addBasicVariable('body_com_x')
                this.addStretchVariable('body_com_x')
            end
            if this.isVariableToAnalyze('body_com_y')
                this.addBasicVariable('com_trajectories')
                this.addBasicVariable('body_com_y')
                this.addStretchVariable('body_com_y')
            end
            if this.isVariableToAnalyze('body_com_z')
                this.addBasicVariable('com_trajectories')
                this.addBasicVariable('body_com_z')
                this.addStretchVariable('body_com_z')
            end
            
  
        end
        
        % interface
        function result = isVariableToAnalyze(this, variable_name)
            result = any(strcmp(this.variables_to_analyze, variable_name));
        end
        function result = isBasicVariable(this, variable_name)
            result = any(strcmp(this.basic_variable_names, variable_name));
        end
        function time_data = getTimeData(this, variable_name) %#ok<STOUT,INUSL>
            eval(['time_data = this.time_data.' variable_name ';']);
        end
        function variable_data = getBasicVariableData(this, variable_name) %#ok<STOUT,INUSL>
            eval(['variable_data = this.basic_variable_data.' variable_name ';']);
        end
        
        function prepareBasicVariables(this, condition, trial, variables_to_prepare)
            if nargin < 4
                variables_to_prepare = this.basic_variable_names;
            end
            
            % clear out old data
            this.basic_variable_data = struct;
            this.basic_variable_labels = struct;
            this.stretch_variable_data = struct;
            this.time_data = struct;
            
            % prepare the data by loading all the basic variables from disk and calculating the required variables
            load(['analysis' filesep makeFileName(this.date, this.subject_id, condition, trial, 'availableVariables')]);
            
            % load basic variables
            for i_variable = 1 : length(variables_to_prepare)
                variable_name = variables_to_prepare{i_variable};
                
                % try loading
                [data, time, sampling_rate, labels, success] = loadData(this.date, this.subject_id, condition, trial, variable_name, 'optional'); %#ok<ASGLU>
                
                % store
                if success
                    eval(['this.basic_variable_data.' variable_name ' = data;']);
                    eval(['this.time_data.' variable_name ' = time;']);
                    eval(['this.basic_variable_labels.' variable_name ' = labels;']);
                end
                
                % calculate variables that can't be loaded
                if strcmp(variable_name, 'lheel_x_pos')
                    LHEE_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LHEE');
                    this.basic_variable_data.lheel_x_pos = LHEE_trajectory(:, 1);
                    this.time_data.lheel_x_pos = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'rheel_x_pos')
                    RHEE_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE');
                    this.basic_variable_data.rheel_x_pos = RHEE_trajectory(:, 1);
                    this.time_data.rheel_x_pos = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'lheel_y_pos')
                    LHEE_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LHEE');
                    this.basic_variable_data.lheel_y_pos = LHEE_trajectory(:, 2);
                    this.time_data.lheel_y_pos = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'rheel_y_pos')
                    RHEE_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE');
                    this.basic_variable_data.rheel_y_pos = RHEE_trajectory(:, 2);
                    this.time_data.rheel_y_pos = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'mpsis_x_pos')
                    LPSI_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI');
                    RPSI_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI');
                    MPSI_trajectory = (LPSI_trajectory + RPSI_trajectory) * 0.5;
                    this.basic_variable_data.mpsis_x_pos = MPSI_trajectory(:, 2);
                    this.time_data.mpsis_x_pos = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'c7_x_pos')
                    c7_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7');
                    this.basic_variable_data.c7_x_pos = c7_trajectory(:, 2);
                    this.time_data.c7_x_pos = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'trunk_angle_ap')
                    c7_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7');
                    LPSI_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI');
                    RPSI_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI');
                    MPSI_trajectory = (LPSI_trajectory + RPSI_trajectory) * 0.5;
                    trunk_vector_y = c7_trajectory(:, 2) - MPSI_trajectory(:, 2);
                    trunk_vector_z = c7_trajectory(:, 3) - MPSI_trajectory(:, 3);
                    this.basic_variable_data.trunk_angle_ap = rad2deg(atan2(trunk_vector_y, trunk_vector_z));
                    this.time_data.trunk_angle_ap = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'trunk_angle_ml')
                    c7_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7');
                    LPSI_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI');
                    RPSI_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI');
                    MPSI_trajectory = (LPSI_trajectory + RPSI_trajectory) * 0.5;
                    trunk_vector_x = c7_trajectory(:, 1) - MPSI_trajectory(:, 1);
                    trunk_vector_z = c7_trajectory(:, 3) - MPSI_trajectory(:, 3);
                    this.basic_variable_data.trunk_angle_ml = rad2deg(atan2(trunk_vector_x, trunk_vector_z));
                    this.time_data.trunk_angle_ml = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'copl_ap')
                    left_foot_cop_world = this.getBasicVariableData('left_foot_cop_world');
                    this.basic_variable_data.copl_ap = left_foot_cop_world(:, 2);
                    this.time_data.copl_ap = this.time_data.left_foot_cop_world;
                end
                if strcmp(variable_name, 'copl_ml')
                    left_foot_cop_world = this.getBasicVariableData('left_foot_cop_world');
                    this.basic_variable_data.copl_ml = left_foot_cop_world(:, 1);
                    this.time_data.copl_ml = this.time_data.left_foot_cop_world;
                end
                if strcmp(variable_name, 'copr_ap')
                    right_foot_cop_world = this.getBasicVariableData('right_foot_cop_world');
                    this.basic_variable_data.copr_ap = right_foot_cop_world(:, 2);
                    this.time_data.copr_ap = this.time_data.right_foot_cop_world;
                end
                if strcmp(variable_name, 'copr_ml')
                    right_foot_cop_world = this.getBasicVariableData('right_foot_cop_world');
                    this.basic_variable_data.copr_ml = right_foot_cop_world(:, 1);
                    this.time_data.copr_ml = this.time_data.right_foot_cop_world;
                end
                if strcmp(variable_name, 'cop_ap')
                    total_forceplate_cop_world = this.getBasicVariableData('total_forceplate_cop_world');
                    this.basic_variable_data.cop_ap = total_forceplate_cop_world(:, 2);
                    this.time_data.cop_ap = this.time_data.total_forceplate_cop_world;
                end
                if strcmp(variable_name, 'cop_ml')
                    total_forceplate_cop_world = this.getBasicVariableData('total_forceplate_cop_world');
                    this.basic_variable_data.cop_ml = total_forceplate_cop_world(:, 1);
                    this.time_data.cop_ml = this.time_data.total_forceplate_cop_world;
                end
                if strcmp(variable_name, 'left_arm_right_leg_relative_phase')
                    left_arm_phase = this.getBasicVariableData('left_arm_phase');
                    right_leg_phase = this.getBasicVariableData('right_leg_phase');
                    
                    left_arm_right_leg_relative_phase = normalizeAngle(left_arm_phase - right_leg_phase);
                    
                    % store
                    this.basic_variable_data.left_arm_right_leg_relative_phase = left_arm_right_leg_relative_phase;
                    this.time_data.left_arm_right_leg_relative_phase = this.time_data.marker_trajectories;
                end             
                if strcmp(variable_name, 'right_arm_left_leg_relative_phase')
                    right_arm_phase = this.getBasicVariableData('right_arm_phase');
                    left_leg_phase = this.getBasicVariableData('left_leg_phase');
                    
                    right_arm_left_leg_relative_phase = normalizeAngle(right_arm_phase - left_leg_phase);
                    
                    % store
                    this.basic_variable_data.right_arm_left_leg_relative_phase = right_arm_left_leg_relative_phase;
                    this.time_data.right_arm_left_leg_relative_phase = this.time_data.marker_trajectories;
                end             
                if strcmp(variable_name, 'left_glut_med')
                    this.basic_variable_data.left_glut_med = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'left_glut_med'));
                    this.time_data.left_glut_med = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'left_delt_ant')
                    this.basic_variable_data.left_delt_ant = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'left_delt_ant'));
                    this.time_data.left_delt_ant = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'left_gastroc_med')
                    this.basic_variable_data.left_gastroc_med = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'left_gastroc_med'));
                    this.time_data.left_gastroc_med = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'left_pero_lng')
                    this.basic_variable_data.left_pero_lng = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'left_pero_lng'));
                    this.time_data.left_pero_lng = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'right_glut_med')
                    this.basic_variable_data.right_glut_med = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'right_glut_med'));
                    this.time_data.right_glut_med = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'right_delt_ant')
                    this.basic_variable_data.right_delt_ant = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'right_delt_ant'));
                    this.time_data.right_delt_ant = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'right_gastroc_med')
                    this.basic_variable_data.right_gastroc_med = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'right_gastroc_med'));
                    this.time_data.right_gastroc_med = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'right_pero_lng')
                    this.basic_variable_data.right_pero_lng = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'right_pero_lng'));
                    this.time_data.right_pero_lng = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'body_com_x')
                    body_com_trajectory = extractMarkerTrajectories(this.basic_variable_data.com_trajectories, this.basic_variable_labels.com_trajectories, 'BODYCOM');
                    this.basic_variable_data.body_com_x = body_com_trajectory(:, 1);
                    this.time_data.body_com_x = this.time_data.com_trajectories;
                end
                if strcmp(variable_name, 'body_com_y')
                    body_com_trajectory = extractMarkerTrajectories(this.basic_variable_data.com_trajectories, this.basic_variable_labels.com_trajectories, 'BODYCOM');
                    this.basic_variable_data.body_com_y = body_com_trajectory(:, 2);
                    this.time_data.body_com_y = this.time_data.com_trajectories;
                end
                if strcmp(variable_name, 'body_com_z')
                    body_com_trajectory = extractMarkerTrajectories(this.basic_variable_data.com_trajectories, this.basic_variable_labels.com_trajectories, 'BODYCOM');
                    this.basic_variable_data.body_com_z = body_com_trajectory(:, 3);
                    this.time_data.body_com_z = this.time_data.com_trajectories;
                end
            end
        end
        function stretch_variables = calculateStretchVariables(this, stretch_start_times, stretch_end_times, condition_stance_foot_list, variables_to_calculate)
            if nargin < 5
                variables_to_calculate = this.stretch_variable_names;
            end
            
            number_of_stretch_variables = length(variables_to_calculate);
            number_of_stretches = length(stretch_start_times);
            stretch_variables = cell(number_of_stretch_variables, 1);
            
            for i_variable = 1 : number_of_stretch_variables
                variable_name = variables_to_calculate{i_variable};
                
                % extract and normalize data from stretches
                for i_stretch = 1 : number_of_stretches
                    stretch_data = [];
                    
                    % time
                    this_stretch_start_time = stretch_start_times(i_stretch);
                    this_stretch_end_time = stretch_end_times(i_stretch);
                    
                    % calculate normalized stretch data for the basic variables
                    if this.isBasicVariable(variable_name)
                        stretch_data = this.getTimeNormalizedData(variable_name, this_stretch_start_time, this_stretch_end_time);
                    end
                
                    % calculate stretch variables that are not basic variables
                    if strcmp(variable_name, 'step_length')
                        lheel_y_pos = this.getTimeNormalizedData('lheel_y_pos', this_stretch_start_time, this_stretch_end_time);
                        rheel_y_pos = this.getTimeNormalizedData('rheel_y_pos', this_stretch_start_time, this_stretch_end_time);
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')
                            stretch_data = lheel_y_pos(end) - rheel_y_pos(end);
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')
                            stretch_data = rheel_y_pos(end) - lheel_y_pos(end);
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_BOTH')
                            stretch_data = NaN;
                        end
                    end
                    if strcmp(variable_name, 'step_width')
                        lheel_x_pos = this.getTimeNormalizedData('lheel_x_pos', this_stretch_start_time, this_stretch_end_time);
                        rheel_x_pos = this.getTimeNormalizedData('rheel_x_pos', this_stretch_start_time, this_stretch_end_time);
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')
                            stretch_data = abs(lheel_x_pos(end) - rheel_x_pos(1));
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')
                            stretch_data = abs(rheel_x_pos(end) - lheel_x_pos(1));
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_BOTH')
                            stretch_data = NaN;
                        end
                    end
                    if strcmp(variable_name, 'step_time')
                        stretch_data = this_stretch_end_time - this_stretch_start_time;
                    end
                    if strcmp(variable_name, 'cadence')
                        second_to_minute = 1/60;
                        step_time = stretch_variables{strcmp(this.stretch_variable_names, 'step_time')}(i_stretch);
                        stretch_data = (step_time * second_to_minute)^(-1);
                    end
                    if strcmp(variable_name, 'velocity')
                        step_time = stretch_variables{strcmp(this.stretch_variable_names, 'step_time')}(i_stretch);
                        step_length = stretch_variables{strcmp(this.stretch_variable_names, 'step_length')}(i_stretch);
                        stretch_data = step_length / step_time;
                    end
                    if strcmp(variable_name, 'left_arm_phase') || strcmp(variable_name, 'right_arm_phase') || strcmp(variable_name, 'left_leg_phase') || strcmp(variable_name, 'right_leg_phase')
                        % extract data
                        variable_time = this.getTimeData(variable_name);
                        variable_data = this.getBasicVariableData(variable_name);
                        [~, start_index] = min(abs(variable_time - this_stretch_start_time));
                        [~, end_index] = min(abs(variable_time - this_stretch_end_time));
                        time_extracted = variable_time(start_index : end_index);
                        data_extracted = variable_data(start_index : end_index);
                        
                        % make sure there's no leaf change within this data stretch
                        data_extracted_groomed = data_extracted;
                        for i_time = 2 : length(time_extracted)
                            if abs(data_extracted_groomed(i_time) - data_extracted_groomed(i_time-1)) > 5
                                % normalize
                                while data_extracted_groomed(i_time) - data_extracted_groomed(i_time-1) <= -pi
                                    data_extracted_groomed(i_time) = data_extracted_groomed(i_time) + 2*pi;
                                end
                                while data_extracted_groomed(i_time) - data_extracted_groomed(i_time-1) > pi
                                    data_extracted_groomed(i_time) = data_extracted_groomed(i_time) - 2*pi;
                                end
                                
                                
                            end
                        end
                        
                        % normalize time
                        time_normalized = linspace(time_extracted(1), time_extracted(end), this.number_of_time_steps_normalized)';
                        data_normalized = spline(time_extracted, data_extracted_groomed, time_normalized);
                        
                        % re-normalize angle
                        stretch_data = normalizeAngle(data_normalized);
                    end
                    if strcmp(variable_name, 'left_arm_phase_at_heelstrike')
                        left_arm_phase = stretch_variables{strcmp(variables_to_calculate, 'left_arm_phase')}(:, i_stretch);
                        stretch_data = left_arm_phase(1);
                        
                        % normalize angle, but take care that the leaf change is not too close to this value
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')
                            stretch_data = normalizeAngle(stretch_data, 0);
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')
                            stretch_data = normalizeAngle(stretch_data, pi);
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_BOTH')
                            stretch_data = NaN;
                        end
                    end
                    if strcmp(variable_name, 'right_arm_phase_at_heelstrike')
                        right_arm_phase = stretch_variables{strcmp(variables_to_calculate, 'right_arm_phase')}(:, i_stretch);
                        stretch_data = right_arm_phase(1);
                        
                        % normalize angle, but take care that the leaf change is not too close to this value
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')
                            stretch_data = normalizeAngle(stretch_data, pi);
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')
                            stretch_data = normalizeAngle(stretch_data, 0);
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_BOTH')
                            stretch_data = NaN;
                        end
                    end
                    
                    
                    if strcmp(variable_name, 'body_com_to_cop_x')
                        body_com_x = stretch_variables{strcmp(this.stretch_variable_names, 'body_com_x')}(:, i_stretch);
                        cop_ml = stretch_variables{strcmp(this.stretch_variable_names, 'cop_ml')}(:, i_stretch);
                        stretch_data = cop_ml - body_com_x;
                    end
                    if strcmp(variable_name, 'body_com_to_cop_y')
                        body_com_y = stretch_variables{strcmp(this.stretch_variable_names, 'body_com_y')}(:, i_stretch);
                        cop_ap = stretch_variables{strcmp(this.stretch_variable_names, 'cop_ap')}(:, i_stretch);
                        stretch_data = cop_ap - body_com_y;
                    end
                    
                    if strcmp(variable_name, 'left_glut_med_rescaled')
                        left_glut_med = this.getTimeNormalizedData('left_glut_med', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'left_glut_med'));
                        stretch_data = left_glut_med * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'left_delt_ant_rescaled')
                        left_delt_ant = this.getTimeNormalizedData('left_delt_ant', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'left_delt_ant'));
                        stretch_data = left_delt_ant * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'left_gastroc_med_rescaled')
                        left_gastroc_med = this.getTimeNormalizedData('left_gastroc_med', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'left_gastroc_med'));
                        stretch_data = left_gastroc_med * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'left_pero_lng_rescaled')
                        left_pero_lng = this.getTimeNormalizedData('left_pero_lng', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'left_pero_lng'));
                        stretch_data = left_pero_lng * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'right_glut_med_rescaled')
                        right_glut_med = this.getTimeNormalizedData('right_glut_med', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'right_glut_med'));
                        stretch_data = right_glut_med * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'right_delt_ant_rescaled')
                        right_delt_ant = this.getTimeNormalizedData('right_delt_ant', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'right_delt_ant'));
                        stretch_data = right_delt_ant * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'right_gastroc_med_rescaled')
                        right_gastroc_med = this.getTimeNormalizedData('right_gastroc_med', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'right_gastroc_med'));
                        stretch_data = right_gastroc_med * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'right_pero_lng_rescaled')
                        right_pero_lng = this.getTimeNormalizedData('right_pero_lng', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'right_pero_lng'));
                        stretch_data = right_pero_lng * 1 / normalization_value;
                    end
                    
                    % store in cell
                    stretch_variables{i_variable} = [stretch_variables{i_variable} stretch_data];
                end
                
            end
            
            
        end
        function data_normalized = getTimeNormalizedData(this, variable_name, start_time, end_time)
            % extract data
            variable_time = this.getTimeData(variable_name);
            variable_data = this.getBasicVariableData(variable_name);
            [~, start_index] = min(abs(variable_time - start_time));
            [~, end_index] = min(abs(variable_time - end_time));
            time_extracted = variable_time(start_index : end_index);
            data_extracted = variable_data(start_index : end_index);
                
            % normalize data in time
            if ~isempty(time_extracted)
                time_normalized = linspace(time_extracted(1), time_extracted(end), this.number_of_time_steps_normalized)';
                data_normalized = spline(time_extracted, data_extracted, time_normalized);
            else
                data_normalized = zeros(this.number_of_time_steps_normalized, 1) * NaN;
            end
        end
    end
end







