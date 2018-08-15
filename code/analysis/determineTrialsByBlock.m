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

% This loads the condition for a single type and trial from a .csv file.

function [block_labels, trials_by_block] = determineTrialsByBlock(condition_table, header)
    trial_column = strcmp(header, 'trial');
    condition_labels = condition_table(:, ~trial_column);
    trials = cellfun(@str2num, condition_table(:, trial_column));
    block_labels = table2cell(unique(cell2table(condition_labels), 'rows'));
    number_of_trials = size(condition_table, 1);
    number_of_blocks = size(block_labels, 1);
    trials_by_block = cell(number_of_blocks, 1);
    for i_block = 1 : size(block_labels, 1)
        this_block_labels = block_labels(i_block, :);
        this_block_row_indicator = true(number_of_trials, 1);
        for i_column = 1 : length(this_block_labels)
            this_column_label = this_block_labels{i_column};
            this_block_row_indicator = this_block_row_indicator & strcmp(this_column_label, condition_labels(:, i_column));
        end
        trials_by_block{i_block} = trials(this_block_row_indicator);
    end





end