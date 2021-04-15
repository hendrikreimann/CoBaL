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

% this function finds the heelstrike and pushoff events

% input: 
% subjectInfo.mat
% subjectModel.mat
% markerTrajectories 
% left_touchdown_method / right_touchdown_method in subjectSettings
%
% output:
% file stepEvents.mat, containing
% - left_pushoff_times
% - left_touchdown_times
% - right_pushoff_times
% - right_touchdown_times


function calculateTimeDerivatives(varargin)

    % parse arguments
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;
    
    % load settings
    subject_settings = loadSettingsFromFile('subject');
    study_settings = loadSettingsFromFile('study');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
    variable_table = study_settings.getTable('time_derivatives');
    number_of_variables = size(variable_table, 1);
    
    for i_condition = 1 : length(condition_list)
        trial_type = condition_list{i_condition};
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            data_to_save = struct;
            save_folder = 'processed';
            save_file_name = makeFileName(collection_date, subject_id, trial_type, i_trial, 'timeDerivatives.mat');
            for i_variable = 1 : number_of_variables
                % load data
                source_variable_name = variable_table.source_variable_name{i_variable};
                filter_order = str2double(variable_table.filter_order{i_variable});
                cutoff_frequency = str2double(variable_table.cutoff_frequency{i_variable});
                [data, time, sampling_rate, labels, directions, success] = loadData(collection_date, subject_id, trial_type, i_trial, source_variable_name); %#ok<ASGLU>
                
                % filter
                [b_filter, a_filter] = butter(filter_order, cutoff_frequency/(sampling_rate/2));
                if any(~isnan(data))
                    data_filtered = deriveByTime(nanfiltfilt(b_filter, a_filter, data), 1/sampling_rate);
                    
                    
                    % save
                    variable_name = variable_table.variable_name{i_variable};
                    data_to_save.(variable_name) = data_filtered;
                    data_to_save.(['time_' variable_name]) = time;
                    data_to_save.(['sampling_rate_' variable_name]) = sampling_rate;
                    data_to_save.(['directions_' variable_name]) = directions;
                    data_to_save.(['label_' variable_name]) = ['time derivative of ' source_variable_name];
                    
                    addAvailableData ...
                      ( ...
                        variable_name, ...
                        ['time_' variable_name], ...
                        ['sampling_rate_' variable_name], ...
                        ['label_' variable_name], ...
                        ['_directions_' variable_name], ...
                        save_folder, ...
                        save_file_name ...
                      );
                else
                    error('cannot filter all-NaN timeseries')
                
                end
                
            end
            save([save_folder filesep save_file_name], '-struct', 'data_to_save');
            
            disp(['Calculating time derivatives: trial type' trial_type ', trial ' num2str(i_trial) ' completed, saved as ' save_file_name]);
        end
    end
    
end














