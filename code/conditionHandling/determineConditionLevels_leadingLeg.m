%     This file is part of the CoBaL code base
%     Copyright (C) 2019 Hendrik Reimann <hendrikreimann@gmail.com>
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

function conditions_trial = determineConditionLevels_leadingLeg(subject_settings, trial_data, conditions_trial)

    number_of_triggers = length(trial_data.trigger_times);
    leading_leg = subject_settings.get('leading_leg');
    if isempty(leading_leg)
        warning('"leading_leg" not specified in subject settings')
    end
    condition_leading_leg_list = cell(number_of_triggers, 1);
    for i_stretch = 1 : number_of_triggers
        condition_leading_leg_list{i_stretch} = leading_leg;
    end
    conditions_trial.leading_leg_list = condition_leading_leg_list;
end

