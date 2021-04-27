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

function [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_visionStochastic(study_settings, subject_settings, trial_data)
    % determine conditions file
    conditions_file_name = [];
    if exist('conditions.csv', 'file')
        conditions_file_name = 'conditions.csv';
    end
    if exist(makeFileName(subject_settings.get('collection_date'), subject_settings.get('subject_id'), 'conditions.csv'), 'file')
        conditions_file_name = makeFileName(collection_date, subject_id, 'conditions.csv');
    end

    stim_frequency = loadConditionFromFile(conditions_file_name, 'frequency', trial_data.trial_number);
    stim_amplitude = loadConditionFromFile(conditions_file_name, 'SD', trial_data.trial_number);
    if trial_data.trial_number < 12
        block = 'FIRST';
    else
        block = 'SECOND';
    end

    stance_foot_data_stretch = {'STANCE_BOTH', 'STANCE_LEFT', 'STANCE_BOTH', 'STANCE_RIGHT'};
    bands_per_stretch = length(stance_foot_data_stretch);

    left_touchdown_times_relevant = ...
        trial_data.left_touchdown_times ...
          ( ...
            trial_data.left_touchdown_times > study_settings.get('analysis_start_time') ...
            & trial_data.left_touchdown_times < study_settings.get('analysis_end_time') ...
          );
    stretch_start_times = left_touchdown_times_relevant(1:end-1);
    number_of_stretches = length(stretch_start_times);
    stretch_times = zeros(number_of_stretches, bands_per_stretch+1);
    removal_flags = false(number_of_stretches, 1);
    for i_stretch = 1 : number_of_stretches
        this_stretch_start = stretch_start_times(i_stretch);
        this_right_pushoff = min(trial_data.right_pushoff_times(trial_data.right_pushoff_times > this_stretch_start));
        this_right_touchdown = min(trial_data.right_touchdown_times(trial_data.right_touchdown_times > this_stretch_start));
        this_left_pushoff = min(trial_data.left_pushoff_times(trial_data.left_pushoff_times > this_stretch_start));
        this_left_touchdown = min(trial_data.left_touchdown_times(trial_data.left_touchdown_times > this_stretch_start));
        this_stretch_times = [this_stretch_start this_right_pushoff this_right_touchdown this_left_pushoff this_left_touchdown];
        if ~issorted(this_stretch_times)
            removal_flags(i_stretch) = 1;
        end
        stretch_times(i_stretch, :) = this_stretch_times;
    end

    stance_foot_data = repmat(stance_foot_data_stretch, size(stretch_times, 1), 1);
    event_variables_to_save = struct;
    event_variables_to_save.stretch_times = stretch_times;
    event_variables_to_save.stance_foot_data = stance_foot_data;

    % conditions
    stim_frequency_list = repmat({['FRQ_' stim_frequency]}, size(stretch_times, 1), 1);
    stim_amplitude_list = repmat({['AMPL_' stim_amplitude]}, size(stretch_times, 1), 1);
    block_list = repmat({block}, size(stretch_times, 1), 1);
    conditions_trial = struct;
    conditions_trial.stim_frequency_list = stim_frequency_list;
    conditions_trial.stim_amplitude_list = stim_amplitude_list;
    conditions_trial.block_list = block_list;


    
    
end

