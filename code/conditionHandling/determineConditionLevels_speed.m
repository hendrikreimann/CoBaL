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

function conditions_trial = determineConditionLevels_speed(trial_data, conditions_trial)
    % find speed label for current trial
    protocol_info = load('protocolInfo.mat');
    this_trial_type = trial_data.trial_type;
    this_trial_number = trial_data.trial_number;
    trial_type_indicator = strcmp(protocol_info.trial_type, this_trial_type);
    trial_number_indicator = (protocol_info.trial_number == this_trial_number);
    this_trial_indicator = trial_type_indicator & trial_number_indicator;
    this_trial_speed = protocol_info.speed(this_trial_indicator);
    speed_label = num2str(this_trial_speed);

    % make list
    number_of_triggers = length(trial_data.trigger_times);
    condition_speed_list = cell(number_of_triggers, 1);
    for i_stretch = 1 : number_of_triggers
        condition_speed_list{i_stretch} = speed_label;
    end
    conditions_trial.speed_list = condition_speed_list;
end

