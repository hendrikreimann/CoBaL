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

function trials_to_process_table = makeTrialsToProcessTable(trial_type_list, trial_number_list)
    trials_to_process_table = table;
    for i_type = 1 : length(trial_type_list)
        this_trial_type = string(trial_type_list{i_type});
        this_trial_numbers = trial_number_list{i_type};
        for i_trial = 1 : length(this_trial_numbers)
            new_row = table(this_trial_type, this_trial_numbers(i_trial));
            trials_to_process_table = [trials_to_process_table; new_row]; %#ok<AGROW>
        end
    end
    trials_to_process_table.Properties.VariableNames = {'trial type', 'trial number'};
end