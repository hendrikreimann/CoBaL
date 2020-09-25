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

function [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_Ngvs(subject_settings, trial_data)
    event_variables_to_save = struct;
    number_of_stretches = length(trial_data.trigger_times);
    removal_flags = false(number_of_stretches, 1);

    % add nGVS settins
    [rms, cutoff_low, cutoff_high, repetition] = determineNgvsSettings(subject_settings, trial_data.trial_type, trial_data.trial_number);
    this_setting_label = [rms ' mA; ' cutoff_low '-' cutoff_high ' Hz'];
    if strcmp(rms, '0')
        this_setting_label = 'control';
    end
    conditions_trial.ngvs_settings_list = repmat({this_setting_label}, size(removal_flags, 1), 1);
    conditions_trial.repetition_list = repmat({repetition}, size(removal_flags, 1), 1);
end

