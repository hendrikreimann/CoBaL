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

function [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels(study_settings, subject_settings, trial_data)
    experimental_paradigm = study_settings.get('experimental_paradigm');
    
    % allocate
    conditions_trial = struct;
    event_variables_to_save = struct;
    removal_flags = zeros(size(trial_data.trigger_times));
    
    intermittent_perturbation_paradigms = {'Vision', 'CadenceVision', 'GVS', 'CadenceGVS', 'FatigueGVS', 'OculusLaneRestriction', 'CognitiveLoadVision', 'CognitiveLoadGvs'};
    if any(strcmp(experimental_paradigm, intermittent_perturbation_paradigms))
        [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_intermittentPerturbations_flex(study_settings, trial_data);
    end


    if strcmp(experimental_paradigm, 'Stochastic Resonance')

        %stim_frequency = loadConditionFromFile(conditions_file_name, 'frequency', i_trial);
        stochastic_stimulus_level_header = subject_settings.get('stochastic_stimulus_level_header');
        stochastic_stimulus_level_table = subject_settings.get('stochastic_stimulus_level_table');

        % get current trial type and number
        this_trial_type = trial_data.trial_type;
        this_trial_number = trial_data.trial_number;

        % extract header information
        trial_type_column_index = find(strcmp(stochastic_stimulus_level_header, 'trial_type'));
        trial_number_column_index = find(strcmp(stochastic_stimulus_level_header, 'trial_number'));
        stimulus_strength_column_index = find(strcmp(stochastic_stimulus_level_header, 'stimulus_strength'));

        % find row in stimulus level table for current trial
        trial_type_column = stochastic_stimulus_level_table(:, trial_type_column_index); %#ok<FNDSB>
        trial_number_column = stochastic_stimulus_level_table(:, trial_number_column_index); %#ok<FNDSB>
        trial_type_indicator = strcmp(trial_type_column, this_trial_type);
        trial_number_indicator = strcmp(trial_number_column, num2str(this_trial_number));
        this_trial_row_index = find(trial_type_indicator & trial_number_indicator);

        % extract stimulus strength for this trial
        this_trial_stimulus_strength = stochastic_stimulus_level_table{this_trial_row_index, stimulus_strength_column_index}; %#ok<FNDSB>

        stance_foot_data_stretch = {'STANCE_LEFT', 'STANCE_RIGHT'};                

        bands_per_stretch = length(stance_foot_data_stretch);

        stretch_start_times = trial_data.left_touchdown_times(1:end-1);
        number_of_stretches = length(stretch_start_times);
        stretch_times = zeros(number_of_stretches, bands_per_stretch+1);
        removal_flags = false(number_of_stretches, 1);
        for i_stretch = 1 : number_of_stretches
            this_stretch_start = stretch_start_times(i_stretch);
            this_right_touchdown = min(trial_data.right_touchdown_times(trial_data.right_touchdown_times > this_stretch_start));
            this_left_touchdown = min(trial_data.left_touchdown_times(trial_data.left_touchdown_times > this_stretch_start));
            this_stretch_times = [this_stretch_start this_right_touchdown this_left_touchdown];
            if ~issorted(this_stretch_times)
                removal_flags(i_stretch) = 1;
            end
            stretch_times(i_stretch, :) = this_stretch_times;
        end
        stretch_times(removal_flags, :) = [];

        stance_foot_data = repmat(stance_foot_data_stretch, size(stretch_times, 1), 1);
        event_variables_to_save.stretch_times = stretch_times;
        event_variables_to_save.stance_foot_data = stance_foot_data;

        % conditions
        %stim_frequency_list = repmat({['FRQ_' stim_frequency]}, size(stretch_times, 1), 1);
        stim_amplitude_list = repmat({this_trial_stimulus_strength}, size(stretch_times, 1), 1);
        conditions_trial = struct;
        conditions_trial.stim_amplitude_list = stim_amplitude_list;

        group = subject_settings.get('group');
        % add group:
        condition_group_list = cell(size(event_variables_to_save.stance_foot_data, 1), 1);
        for i_stretch = 1 : length(condition_group_list)
            condition_group_list{i_stretch} = group;
        end
        conditions_trial.group_list = condition_group_list;
    end






end