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

function saveProblemInformation(date, subject_id, trial_type, trial_number, reason, time, bad_data_points_indicator)
    if nargin < 7
        % time is a point indicating the problem
        bad_data_points_indicator = 1;
    else
        % time is a vector, bad_data_points_indicator is a logical that indicates the problem
    end

    if ~any(bad_data_points_indicator)
        % nothing to do
        return
    end

    % create empty problem table
    bad_data_points_indices = find(bad_data_points_indicator);
    problems = table ...
      ( ...
        'Size', [0, 3], ...
        'VariableNames', {'start_time', 'end_time', 'reason'}, ...
        'VariableTypes', {'double', 'double', 'string'} ...
      );
    
    
    % add separate problems to struct one by one
    while ~isempty(bad_data_points_indices)
        addNextProblemToTable();
    end
    
    % save problem struct to file
   saveProblemsToFile();
   
function addNextProblemToTable
    % identify the end of the first problem 
    first_jump = find(diff(bad_data_points_indices) - 1);
    if isempty(first_jump)
        % no jumps in the indices, so only one problem left
        this_problem_end = length(bad_data_points_indices);
    else
        % found a jump between indices, so that's where the problem ends
        this_problem_end = first_jump;
    end
    % add this problem to the struct
    this_problem_indices = bad_data_points_indices(1 : this_problem_end);
    this_problem = {time(this_problem_indices(1)), time(this_problem_indices(end)), reason};
    problems = [problems; this_problem];
    
    % remove this problem from bad_data_points_indices
    bad_data_points_indices(1 : this_problem_end) = [];
end

function saveProblemsToFile()
    save_folder = 'analysis';
    save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'problems.mat');
    
    % load list of already available variables if existing
    if exist([save_folder filesep save_file_name], 'file')
        problems_loaded = load([save_folder filesep save_file_name], 'problems');
        
        problems = [problems; problems_loaded.problems];
    end

    % remove duplicated entries
    problems = unique(problems, 'stable');
    
    % save
    save([save_folder filesep save_file_name], 'problems');
end
end