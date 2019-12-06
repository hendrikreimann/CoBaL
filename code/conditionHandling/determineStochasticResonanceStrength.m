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


function stimulus_strength = determineStochasticResonanceStrength(subject_settings, trial_type, trial_number)
    % extract header information
    stochastic_stimulus_level_header = subject_settings.get('stochastic_stimulus_level_header');
    stochastic_stimulus_level_table = subject_settings.get('stochastic_stimulus_level_table');
    trial_type_column_index = find(strcmp(stochastic_stimulus_level_header, 'trial_type'));
    trial_number_column_index = find(strcmp(stochastic_stimulus_level_header, 'trial_number'));
    stimulus_strength_column_index = find(strcmp(stochastic_stimulus_level_header, 'stimulus_strength'));

    % find row in stimulus level table for current trial
    trial_type_column = stochastic_stimulus_level_table(:, trial_type_column_index); %#ok<FNDSB>
    trial_number_column = stochastic_stimulus_level_table(:, trial_number_column_index); %#ok<FNDSB>
    trial_type_indicator = strcmp(trial_type_column, trial_type);
    trial_number_indicator = strcmp(trial_number_column, num2str(trial_number));
    this_trial_row_index = find(trial_type_indicator & trial_number_indicator);

    % extract stimulus strength for this trial
    stimulus_strength = stochastic_stimulus_level_table{this_trial_row_index, stimulus_strength_column_index}; %#ok<FNDSB>
    
end