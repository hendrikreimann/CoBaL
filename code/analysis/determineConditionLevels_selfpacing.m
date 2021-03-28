%     This file is part of the CoBaL code base
%     Copyright (C) 2019 Hendrik Reimann <hendrikreimann@gmail.com>
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

function [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_selfpacing(subject_settings, trial_data)
    event_variables_to_save = struct;
    number_of_stretches = length(trial_data.trigger_times);
    removal_flags = false(number_of_stretches, 1);

    % extract header information
    selfpacing_settings_header = subject_settings.get('selfpacing_settings_header');
    selfpacing_settings_table = subject_settings.get('selfpacing_settings_table');
    trial_number_column_index = find(strcmp(selfpacing_settings_header, 'trial_number'));
    controller_column_index = find(strcmp(selfpacing_settings_header, 'controller'));
    parameter_set_column_index = find(strcmp(selfpacing_settings_header, 'parameter_set'));
    block_column_index = find(strcmp(selfpacing_settings_header, 'block'));

    % find row in stimulus level table for current trial
    trial_number_column = selfpacing_settings_table(:, trial_number_column_index); %#ok<FNDSB>
    trial_number_indicator = strcmp(trial_number_column, num2str(trial_data.trial_number));
    this_trial_row_index = find(trial_number_indicator);

    % extract stimulus strength for this trial
    controller = selfpacing_settings_table{this_trial_row_index, controller_column_index}; %#ok<FNDSB>
    parameter_set = selfpacing_settings_table{this_trial_row_index, parameter_set_column_index}; %#ok<FNDSB>
    block = selfpacing_settings_table{this_trial_row_index, block_column_index}; %#ok<FNDSB>
    
    % determine interval
    this_trial_type = char(trial_data.trial_type);
    interval = this_trial_type(end);
    
    % add levels
    conditions_trial.controller_list = repmat({controller}, size(removal_flags, 1), 1);
    conditions_trial.parameter_set_list = repmat({parameter_set}, size(removal_flags, 1), 1);
    conditions_trial.block_list = repmat({block}, size(removal_flags, 1), 1);
    conditions_trial.interval_list = repmat({interval}, size(removal_flags, 1), 1);





end

