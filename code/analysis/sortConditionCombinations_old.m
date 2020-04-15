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

function [condition_combinations_stimulus_sorted, condition_combinations_control_sorted] = sortConditionCombinations(condition_combinations_stimulus, condition_combinations_control, condition_combination_labels, condition_to_compare, preferred_level_order)
    % if nothing was specified, do nothing
    condition_combinations_stimulus_sorted = condition_combinations_stimulus;
    condition_combinations_control_sorted = condition_combinations_control;
    if isempty(preferred_level_order)
        return
    end
    
    % check whether the preferred level order catches all we have
    relevant_label_column = condition_combinations_stimulus(:, strcmp(condition_combination_labels, condition_to_compare));
    unique_labels = unique(relevant_label_column);
    
    % add labels in preferred order if they're present
    level_order = {};
    for i_level = 1 : length(preferred_level_order)
        this_level = preferred_level_order(i_level);
        if any(strcmp(this_level, unique_labels))
            level_order = [level_order; this_level]; %#ok<AGROW>
        end
    end
    
    % add any levels that are in the data but not in the preferred order
    for i_level = 1 : length(unique_labels)
        this_level = unique_labels(i_level);
        if ~any(strcmp(this_level, level_order))
            level_order = [level_order; this_level]; %#ok<AGROW>
        end
    end
    
    % sort
    condition_combinations_stimulus_sorted = {};
    sortmap = [];
    for i_level = 1 : length(level_order)
        this_level = level_order(i_level);
        relevant_label_column = condition_combinations_stimulus(:, strcmp(condition_combination_labels, condition_to_compare));
        this_level_rows = strcmp(relevant_label_column, this_level);
        condition_combinations_stimulus_sorted = [condition_combinations_stimulus_sorted; condition_combinations_stimulus(this_level_rows, :)]; %#ok<AGROW>
        sortmap = [sortmap; find(this_level_rows)]; %#ok<AGROW>
    end
    
    if ~isempty(condition_combinations_control)
        condition_combinations_control_sorted = condition_combinations_control(sortmap, :);
    else
        condition_combinations_control_sorted = [];
    end
end









