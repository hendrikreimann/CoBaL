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

function condition_indicator = getConditionIndicator(condition_combination, condition_combination_labels, condition_data, condition_labels, levels_to_remove)
%     hits = false(size(condition_data, 1), length(condition_combination_labels));

    % create indicator
    condition_indicator = ones(size(condition_data, 1), 1);
    
    % remove labels that don't fit the requested combination
    for i_label = 1 : length(condition_combination_labels)
        this_condition_label = condition_combination_labels{i_label};
        this_condition_data = condition_data(:, strcmp(condition_labels, this_condition_label));
        this_label_hits = strcmp(this_condition_data, condition_combination{i_label});
%         hits(:, i_label) = this_label_hits;
        
        condition_indicator = condition_indicator.*this_label_hits;
    end
    
    % remove specific labels as requested
    for i_row = 1 : size(levels_to_remove, 1)
        this_condition_label = levels_to_remove{i_row, 1};
        this_condition_data = condition_data(:, strcmp(condition_labels, this_condition_label));
        level_to_remove = levels_to_remove{i_row, 2};
        this_label_hits = strcmp(this_condition_data, level_to_remove);
        
        % leave only those rows in that do not fit this level
        condition_indicator = condition_indicator.*(~this_label_hits);
    end
    
    condition_indicator = logical(condition_indicator);
end