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

function [comparison_indices, conditions_per_comparison_max] = determineComparisons_new(conditions_to_plot, condition_labels, plot_settings)
    % initialize
    number_of_conditions_to_plot = size(conditions_to_plot, 1);
    comparison_indices = {};
    conditions_already_compared = [];
    conditions_per_comparison_max = 0;
    condition_to_compare = plot_settings.get('condition_to_compare');
    
    % here, we go through all conditions_to_plot and group up those that go into one comparison, i.e. one figure
    while length(conditions_already_compared) < number_of_conditions_to_plot
        % start with the first available condition
        i_condition = 1;
        while ismember(i_condition, conditions_already_compared)
            i_condition = i_condition + 1;
        end

        this_comparison = i_condition; % this is the first condition in this comparison, more will be added
        % search for conditions that differ from this one only in the one we're comparing
        for j_condition = 1 : number_of_conditions_to_plot
            if i_condition ~= j_condition
                % check which conditions labels agree between these two conditions
                comparison_table = zeros(1, length(condition_labels)); % this is a table indicating equality between the two conditions in questions
                for i_label = 1 : length(condition_labels)
                    comparison_table(i_label) = strcmp(conditions_to_plot{i_condition, i_label}, conditions_to_plot{j_condition, i_label});
                end

                % look at the relevant entries of the comparison table
                comparison_table_relevant = comparison_table;
                comparison_table_relevant(strcmp(condition_labels, condition_to_compare)) = [];
                if all(comparison_table_relevant)
                    this_comparison = [this_comparison, j_condition]; %#ok<AGROW>
                end
            end
        end
        comparison_indices = [comparison_indices; this_comparison]; %#ok<AGROW>
        conditions_already_compared = [conditions_already_compared this_comparison]; %#ok<AGROW>
        
        if conditions_per_comparison_max < length(this_comparison)
            conditions_per_comparison_max = length(this_comparison);
        end
    end    
%     if ~isempty(study_settings.get('conditions_control'))
    if plot_settings.get('plot_control')
        conditions_per_comparison_max = conditions_per_comparison_max + 1;
    end
end