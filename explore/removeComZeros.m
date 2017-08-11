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



% corrects the already calculated angles by inverting the abduction and flexion angles for the right hip

[condition_list, trial_number_list] = parseTrialArguments();
load('subjectInfo.mat', 'date', 'subject_id');

for i_condition = 1 : length(condition_list)
    trials_to_process = trial_number_list{i_condition};
    for i_trial = trials_to_process
        %% load data
        condition = condition_list{i_condition};
        
        % load
        load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'kinematicTrajectories')]);

        % set zeros to NaN
        com_trajectories(com_trajectories==0) = NaN;
        com_trajectories_optimized(com_trajectories_optimized==0) = NaN;
        
        % save
        variables_to_save = struct;
        variables_to_save.com_trajectories = com_trajectories;
        variables_to_save.com_trajectories_optimized = com_trajectories_optimized;

        save_folder = 'processed';
        save_file_name = makeFileName(date, subject_id, condition, i_trial, 'kinematicTrajectories.mat');
        saveDataToFile([save_folder filesep save_file_name], variables_to_save);
        disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' save_folder filesep save_file_name]);        
        
    end
end








