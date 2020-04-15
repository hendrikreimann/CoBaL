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
    condition_combinations_emg ...
  ] ...
  = determineConditionCombinations_new_new ...
  ( ...
    condition_data_all, ...
    conditions_settings, ...
    labels_to_ignore, ...
    levels_to_remove ...
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
    
    % make EMG conditions cell
    if nargout > 2
        condition_combinations_emg = condition_combinations;
        
        % for each factor, check if an EMG condition was specified
        for i_factor = 1 : number_of_condition_labels
            this_factor_condition_label = conditions_settings{i_factor, 1};
            this_factor_emg_level = conditions_settings{i_factor, 4};
            if ~strcmp(this_factor_emg_level, '~')
                % there's a special EMG condition, so remove all others
                this_factor_column = strcmp(condition_combination_labels, this_factor_condition_label);
                this_factor_emg_match = strcmp(condition_combinations_emg(:, this_factor_column), this_factor_emg_level);
                condition_combinations_emg(~this_factor_emg_match, :) = [];
                
            end
        end
        
    end
    
    
    
    
    
    
    
% 
%         
%     % make emg conditions cell
%     if nargout > 3
%         condition_combinations_emg = {};
%         if any(~strcmp(conditions_settings(:, 4), '~'))
%             for i_combination = 1 : size(condition_combinations,1)
%                 this_combination_stimulus = condition_combinations(i_combination, :);
%                 this_combination_emg = cell(size(this_combination_stimulus));
%                 for i_label = 1 : length(condition_combination_labels)
%                     this_label = condition_combination_labels{i_label};
%                     emg_level_this_label = conditions_settings{strcmp(conditions_settings(:, 1), this_label), 4};
%                     if strcmp(emg_level_this_label, '~')
%                         this_combination_emg{i_label} = this_combination_stimulus{i_label};
%                     else
%                         this_combination_emg{i_label} = emg_level_this_label;
%                     end
%                 end
% 
%                 condition_combinations_emg = [condition_combinations_emg; this_combination_emg]; %#ok<AGROW>
%             end
%             condition_combinations_emg_unique = table2cell(unique(cell2table(condition_combinations_emg), 'rows'));
%         end
%     end
end























