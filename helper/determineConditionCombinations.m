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

function [condition_combination_labels, condition_combinations_stimulus, condition_combinations_control, condition_combinations_emg_unique] = determineConditionCombinations(condition_data_all, conditions_settings, labels_to_ignore, levels_to_remove)
    condition_labels = conditions_settings(:, 1)';
    number_of_condition_labels = length(condition_labels);
    condition_combination_labels = {};
    condition_combinations_stimulus = {};
    for i_label = 1 : number_of_condition_labels
        if ~any(strcmp(condition_labels{i_label}, labels_to_ignore))
            this_label = condition_labels{i_label};
            condition_combination_labels = [condition_combination_labels this_label]; %#ok<AGROW>
            levels_this_label = unique(condition_data_all(:, i_label));
            
            % exclude control level if there is one
            control_level_this_label = conditions_settings(strcmp(conditions_settings(:, 1), this_label), 3);
            control_level_index_in_levels = find(strcmp(levels_this_label, control_level_this_label));
            if ~isempty(control_level_index_in_levels)
                levels_this_label(control_level_index_in_levels, :) = [];
            end
            
            % remove levels specified in settings
            if numel(levels_to_remove) > 0
                levels_to_remove_this_label = levels_to_remove(strcmp(levels_to_remove(:, 1), this_label), 2);
                for i_level = 1 : length(levels_to_remove_this_label)
                    index_to_remove = strcmp(levels_to_remove_this_label(i_level), levels_this_label);
                    levels_this_label(index_to_remove) = [];
                end
            end
            
            % repeat stuff to get combinations
            if isempty(condition_combinations_stimulus)
                condition_combinations_stimulus = levels_this_label;
            else
                condition_combinations_pre_this_level = condition_combinations_stimulus;
                condition_combinations_post_this_level = {};
                for i_level = 1 : length(levels_this_label)
                    this_condition_this_level = cell(size(condition_combinations_pre_this_level, 1), 1);
                    this_condition_this_level(:) = levels_this_label(i_level);
                    condition_combinations_with_this_level = [condition_combinations_pre_this_level, this_condition_this_level];
                    condition_combinations_post_this_level = [condition_combinations_post_this_level; condition_combinations_with_this_level]; %#ok<AGROW>
                end
                condition_combinations_stimulus = condition_combinations_post_this_level;
            end
        end
    end
    
    % make control conditions cell
    condition_combinations_control = {};
    if any(~strcmp(conditions_settings(:, 3), '~'))
        % we have some conditions pointing to a control condition, so process this
        for i_combination = 1 : size(condition_combinations_stimulus,1)
            this_combination_stimulus = condition_combinations_stimulus(i_combination, :);
            this_combination_control = cell(size(this_combination_stimulus));
            for i_label = 1 : length(condition_combination_labels)
                this_label = condition_combination_labels{i_label};
                control_level_this_label = conditions_settings{strcmp(conditions_settings(:, 1), this_label), 3};
                if strcmp(control_level_this_label, '~')
                    this_combination_control{i_label} = this_combination_stimulus{i_label};
                else
                    this_combination_control{i_label} = control_level_this_label;
                end
            end
            condition_combinations_control = [condition_combinations_control; this_combination_control]; %#ok<AGROW>
        end
    end
    
    % make emg conditions cell
    if nargout > 3
        condition_combinations_emg = {};
        if any(~strcmp(conditions_settings(:, 4), '~'))
            for i_combination = 1 : size(condition_combinations_stimulus,1)
                this_combination_stimulus = condition_combinations_stimulus(i_combination, :);
                this_combination_emg = cell(size(this_combination_stimulus));
                for i_label = 1 : length(condition_combination_labels)
                    this_label = condition_combination_labels{i_label};
                    emg_level_this_label = conditions_settings{strcmp(conditions_settings(:, 1), this_label), 4};
                    if strcmp(emg_level_this_label, '~')
                        this_combination_emg{i_label} = this_combination_stimulus{i_label};
                    else
                        this_combination_emg{i_label} = emg_level_this_label;
                    end
                end

                condition_combinations_emg = [condition_combinations_emg; this_combination_emg]; %#ok<AGROW>
            end
            condition_combinations_emg_unique = table2cell(unique(cell2table(condition_combinations_emg), 'rows'));
        end
    end
end