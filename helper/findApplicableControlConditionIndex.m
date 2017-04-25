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


function applicable_control_condition_index = findApplicableControlConditionIndex(this_condition, conditions_control)
    control_conditions_with_same_stance_foot_indicator = strcmp(conditions_control(:, 1), this_condition{1});
    control_conditions_with_same_trial_type_indicator = strcmp(conditions_control(:, 5), this_condition{5});
    control_conditions_with_same_stimulus_indicator = strcmp(conditions_control(:, 6), this_condition{6});
    control_conditions_with_same_day_indicator = strcmp(conditions_control(:, 7), this_condition{7});

    applicable_control_condition_indicator = ...
      ( ...
        control_conditions_with_same_stance_foot_indicator ...
        & control_conditions_with_same_trial_type_indicator ...
        & control_conditions_with_same_stimulus_indicator ...
        & control_conditions_with_same_day_indicator ...
      );
    applicable_control_condition_index = find(applicable_control_condition_indicator, '1', 'first');
%     applicable_control_condition_index = find(strcmp(conditions_control, this_condition{1}), 1, 'first');

end