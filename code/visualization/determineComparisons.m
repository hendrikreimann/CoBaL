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

% this function gets information about all combinations of levels and
% groups them into "comparisons". One comparison is one list of conditions
% that will go into the same figure. 

function [comparison_indices, conditions_per_comparison_max] = determineComparisons(conditions_to_plot, condition_labels, settings)
    % initialize
    number_of_conditions_to_plot = size(conditions_to_plot, 1);
    conditions_already_compared = [];
    condition_to_compare = settings.plot_settings.get('condition_to_compare');
    condition_table = settings.conditions_settings;
    
    % separate the control conditions out
    control_row_indicator = false(number_of_conditions_to_plot, 1);
    for i_row = 1 : number_of_conditions_to_plot
        for i_factor = 1 : length(condition_labels)
            % what is the control level for this factor?
            this_factor = condition_labels{i_factor};
            control_level = condition_table(strcmp(condition_table(:, 1), this_factor), 3);
            
            % does this row match the control level?
            this_level = conditions_to_plot{i_row, i_factor};
            control_row_indicator(i_row) = control_row_indicator(i_row) | strcmp(this_level, control_level);
            
        end
    end
    conditions_to_plot_stimulus = conditions_to_plot(~control_row_indicator, :);
        
    % here, we go through all stimulus conditions and group up those that go into one comparison, i.e. one figure
    comparison_indices_stimulus = {};
    number_of_conditions_to_plot_stimulus = size(conditions_to_plot_stimulus, 1);
    while length(conditions_already_compared) < number_of_conditions_to_plot_stimulus
        % start with the first available condition
        i_condition = 1;
        while ismember(i_condition, conditions_already_compared)
            i_condition = i_condition + 1;
        end

        this_comparison = i_condition; % this is the first condition in this comparison, more will be added
        % search for conditions that differ from this one only in the one we're comparing
        for j_condition = 1 : number_of_conditions_to_plot_stimulus
            if i_condition ~= j_condition
                % check which conditions labels agree between these two conditions
                comparison_table = zeros(1, length(condition_labels)); % this is a table indicating equality between the two conditions in questions
                for i_label = 1 : length(condition_labels)
                    comparison_table(i_label) = strcmp(conditions_to_plot_stimulus{i_condition, i_label}, conditions_to_plot_stimulus{j_condition, i_label});
                end

                % look at the relevant entries of the comparison table
                comparison_table_relevant = comparison_table;
                comparison_table_relevant(strcmp(condition_labels, condition_to_compare)) = [];
                if all(comparison_table_relevant)
                    this_comparison = [this_comparison, j_condition]; %#ok<AGROW>
                end
            end
        end
        comparison_indices_stimulus = [comparison_indices_stimulus; this_comparison]; %#ok<AGROW>
        conditions_already_compared = [conditions_already_compared this_comparison]; %#ok<AGROW>
    end
    
    % now we transform the condition indices from relating to rows in conditions_to_plot_stimulus
    % to relate to rows in conditions_to_plot
    comparison_indices = cell(size(comparison_indices_stimulus));
    for i_comparison = 1 : length(comparison_indices_stimulus)
        this_comparison_stimulus = comparison_indices_stimulus{i_comparison};
        this_comparison = zeros(size(this_comparison_stimulus));
        for i_entry = 1 : length(this_comparison)
            this_entry = conditions_to_plot_stimulus(this_comparison_stimulus(i_entry), :);
            this_entry_index_in_full_list = findMatchingRow(conditions_to_plot, this_entry);
            this_comparison(i_entry) = this_entry_index_in_full_list;
        end
        comparison_indices{i_comparison} = this_comparison;
    end    
    
    % now we need to add the relevant control conditions back to each comparison
    if any(~strcmp(condition_table(:, 3), '~'))
        % there are factors with a specified control level, so add that
        for i_comparison = 1 : length(comparison_indices)
            % extract information about this comparison
            this_comparison = comparison_indices{i_comparison};
            this_comparison_conditions_stimulus = conditions_to_plot(this_comparison, :);

            % go through each factor, find relevant control condition and assemble
            control_combination_for_this_level = cell(1, size(this_comparison_conditions_stimulus, 2));
            for i_factor = 1 : length(condition_labels)
                % what is the control level for this factor?
                this_factor = condition_labels{i_factor};
                control_level = condition_table{strcmp(condition_table(:, 1), this_factor), 3};

                if strcmp(control_level, '~')
                    % no control specified for this factor, use the level from the current comparison
                    control_combination_for_this_level{i_factor} = this_comparison_conditions_stimulus{1, i_factor};
                else
                    % use specified control level
                    control_combination_for_this_level{i_factor} = control_level;
                end
            end        

            % find the assembled control condition in the list
            matching_row = findMatchingRow(conditions_to_plot, control_combination_for_this_level);

            % add the control condition to the comparison
            if ~isempty(matching_row)
                this_comparison = [this_comparison matching_row]; %#ok<AGROW>
                comparison_indices{i_comparison} = this_comparison;
            end
        end
    end
    
    % determine maximum of the number of conditions in each comparison
    comparison_sizes = zeros(length(comparison_indices), 1);
    for i_comparison = 1 : length(comparison_indices)
        this_comparison = comparison_indices{i_comparison};
        comparison_sizes(i_comparison) = length(this_comparison);
    end
    conditions_per_comparison_max = max(comparison_sizes);
        
end









