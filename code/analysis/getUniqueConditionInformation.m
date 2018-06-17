%     This file is part of the CoBaL code base
%     Copyright (C) 2018 Hendrik Reimann <hendrikreimann@gmail.com>
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


function [condition_combinations_unique, conditions_indicators] = getUniqueConditionInformation(condition_data, condition_labels)
    number_of_stretches = size(condition_data, 1);
    number_of_conditions = size(condition_data, 2);

    % extract indicators for control
    condition_combinations_unique = table2cell(unique(cell2table(condition_data), 'rows'));
    number_of_conditions = size(condition_combinations_unique, 1);
    conditions_indicators = true(number_of_stretches, number_of_conditions);
    for i_condition = 1 : number_of_conditions
        for i_label = 1 : length(condition_labels)
            this_label_list = condition_data(:, i_label);
            this_label_indicator = strcmp(this_label_list, condition_combinations_unique(i_condition, i_label));
            conditions_indicators(:, i_condition) = conditions_indicators(:, i_condition) .* this_label_indicator;
        end
    end
end