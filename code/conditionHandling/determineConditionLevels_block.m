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

function conditions_trial = determineConditionLevels_block(subject_settings, trial_data, conditions_trial)
    number_of_triggers = length(trial_data.trigger_times);
    block_table = subject_settings.getTable('blocks');
    condition_block_list = cell(number_of_triggers, 1);
    
    % get block
    this_trial_row = strcmp(block_table.trial_type, char(trial_data.trial_type)) & str2double(block_table.trial_number)==trial_data.trial_number;
    this_block = block_table.block{this_trial_row};
    
    for i_stretch = 1 : number_of_triggers
        condition_block_list{i_stretch} = this_block;
    end
    conditions_trial.block_list = condition_block_list;
end

