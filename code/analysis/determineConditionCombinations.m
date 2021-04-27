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

function ...
  [ ...
    condition_combination_labels, ...
    condition_combinations, ...
    condition_combinations_filtered ...
  ] ...
  = determineConditionCombinations ...
  ( ...
    condition_data_all, ...
    conditions_settings, ...
    labels_to_ignore, ...
    levels_to_remove, ...
    output_filter ...
  )
    
    % extract condition labels
    condition_labels = conditions_settings(:, 1)';
    number_of_condition_labels = length(condition_labels);
    
    % remove columns that we don't care about
    labels_to_use_indicator = true(1, number_of_condition_labels);
    for i_label = 1 : number_of_condition_labels
        this_label = condition_labels{i_label};
        if any(strcmp(this_label, labels_to_ignore))
            labels_to_use_indicator(i_label) = false;
        end
    end
    condition_combination_labels = condition_labels(labels_to_use_indicator);
    condition_data = condition_data_all(:, labels_to_use_indicator);

    % get unique entries
    condition_combinations = table2cell(unique(cell2table(condition_data), 'rows'));

    % remove levels specified in settings
    for i_row = 1 : size(levels_to_remove, 1)
        % extract info about what to remove
        this_condition_label = levels_to_remove(i_row, 1);
        this_condition_column = strcmp(condition_combination_labels, this_condition_label);
        this_level = levels_to_remove(i_row, 2);
        
        % identify and remove rows that match
        matching_rows = strcmp(condition_combinations(:, this_condition_column), this_level);
        condition_combinations(matching_rows, :) = [];
    end
    
    
    % apply output filter if required
    if nargin > 4
        % decide which filter to apply here
        if strcmp(output_filter, 'control') || strcmp(output_filter, 'stimulus')
            relevant_column_in_condition_settings = 3;
        end
        if strcmp(output_filter, 'emg')
            relevant_column_in_condition_settings = 4;
        end
        
        % make container for filter output
        condition_combinations_filtered = condition_combinations;

        % for each factor, check if a condition was specified in the relevant column
        for i_factor = 1 : number_of_condition_labels
            this_factor_condition_label = conditions_settings{i_factor, 1};
            this_factor_filter_level = conditions_settings{i_factor, relevant_column_in_condition_settings};
            if ~strcmp(this_factor_filter_level, '~')
                % there's a condition specified, so remove all others
                this_factor_column = strcmp(condition_combination_labels, this_factor_condition_label);
                this_factor_emg_match = strcmp(condition_combinations_filtered(:, this_factor_column), this_factor_filter_level);
                condition_combinations_filtered(~this_factor_emg_match, :) = [];
            end
        end

        if strcmp(output_filter, 'control') || strcmp(output_filter, 'emg')
            % assign this as output
            condition_combinations = condition_combinations_filtered;
        end
        if strcmp(output_filter, 'stimulus')
            % remove the identified conditions from the whole to get stimulus
            for i_row = 1 : size(condition_combinations_filtered)
                this_row = condition_combinations_filtered(i_row, :);
                this_row_index_in_whole = findMatchingRow(condition_combinations, this_row);
                condition_combinations(this_row_index_in_whole, :) = [];
            end
            
            
            
        end
        
    end
    

    
    
    
    
    
    
end























