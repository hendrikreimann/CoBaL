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

% to add a new variable, you need to do two things:
%
% 1. Provide information about all variables this new one depends upon. Go to the determineVariables function and add a
% new conditional statement "if this.isVariableToAnalyze('<your new variable name>') ... end". Within the
% statement, add the basic variables and stretch variables that the new variable depends upon by calling addBasicVariable 
% and addStretchVariable. The order matters here, make sure for each variable, all dependencies have already been added, 
% i.e. don't call "this.addBasicVariable('lheel_y_pos')" before "this.addBasicVariable('marker_trajectories')", because 
% lheel_y_pos depends upon the marker_trajectories.
%
% 2A. If the variable you want to add is a basic variable, but is not already available, i.e. one of the outputs of the
% preprocess scripts, you need to specify how to calculate it. Go to prepareBasicVariables and add a new conditional
% statement "if strcmp(variable_name, '<your new variable name>')". Within the statement, calculate your new
% variable, then add it to the variable_data. Also add the corresponding time vector to the time_data.
%
% 2B. If the variable you want to add is a stretch variable, you need to specify how to calculate it. Go to 
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
    end
    
    methods
        % constructor
%         function this = WalkingDataCustodian(date, subject_id, variables_to_analyze)
        function this = WalkingDataCustodian()
            % load this information from the subjects.mat and studySettings.txt files
            load('subjectInfo.mat', 'date', 'subject_id');
            study_settings = loadSettingsFile(['..' filesep 'studySettings.txt']);
            
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
        
        function prepareBasicVariables(this, condition, trial)
            % clear out old data
            this.basic_variable_data = struct;
            this.basic_variable_labels = struct;
            this.stretch_variable_data = struct;
            this.time_data = struct;
            
            % prepare the data by loading all the basic variables from disk and calculating the required variables
            load(['analysis' filesep makeFileName(this.date, this.subject_id, condition, trial, 'availableVariables')]);
            
            % load basic variables
            for i_variable = 1 : length(this.basic_variable_names)
                variable_name = this.basic_variable_names{i_variable};
                
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
                    variable_data = LHEE_trajectory(:, 1);
                    variable_time = this.time_data.marker_trajectories;
                    this.basic_variable_data.lheel_x_pos = variable_data;
                    this.time_data.lheel_x_pos = variable_time;
                end
                if strcmp(variable_name, 'rheel_x_pos')
                    RHEE_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE');
                    variable_data = RHEE_trajectory(:, 1);
                    variable_time = this.time_data.marker_trajectories;
                    this.basic_variable_data.rheel_x_pos = variable_data;
                    this.time_data.rheel_x_pos = variable_time;
                end
                if strcmp(variable_name, 'lheel_y_pos')
                    LHEE_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LHEE');
                    variable_data = LHEE_trajectory(:, 2);
                    variable_time = this.time_data.marker_trajectories;
                    this.basic_variable_data.lheel_y_pos = variable_data;
                    this.time_data.lheel_y_pos = variable_time;
                end
                if strcmp(variable_name, 'rheel_y_pos')
                    RHEE_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE');
                    variable_data = RHEE_trajectory(:, 2);
                    variable_time = this.time_data.marker_trajectories;
                    this.basic_variable_data.rheel_y_pos = variable_data;
                    this.time_data.rheel_y_pos = variable_time;
                end
                
            end
        end
        function stretch_variables = calculateStretchVariables(this, stretch_start_times, stretch_end_times, condition_stance_foot_list)
            number_of_stretch_variables = length(this.stretch_variable_names);
            number_of_stretches = length(stretch_start_times);
            stretch_variables = cell(number_of_stretch_variables, 1);
            
            for i_variable = 1 : number_of_stretch_variables
                variable_name = this.stretch_variable_names{i_variable};
                
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
                
%                     if strcmp(variable_name, 'lheel_x_pos')
%                         stretch_data = this.getTimeNormalizedData(this.stretch_variable_names{i_variable}, this_stretch_start_time, this_stretch_end_time);
%                     end
%                     if strcmp(variable_name, 'rheel_x_pos')
%                         stretch_data = this.getTimeNormalizedData(this.stretch_variable_names{i_variable}, this_stretch_start_time, this_stretch_end_time);
%                     end
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
                            stretch_data = 0;
                        end
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
            time_normalized = linspace(time_extracted(1), time_extracted(end), this.number_of_time_steps_normalized)';
            data_normalized = spline(time_extracted, data_extracted, time_normalized);
        end
    end
end







