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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.% compare the kinematic tree against the kinematic chain

function [conditions_to_analyze, trials_to_analyze] = parseTrialArguments(varargin)
    load('subjectInfo.mat', 'condition_list', 'trial_number_list');

    parser = inputParser;
    parser.KeepUnmatched = true;
    
    default_condition = 'all';
    addParameter(parser, 'condition', default_condition)
    default_trials = 0;
    addParameter(parser, 'trials', default_trials)
    
    
    parse(parser, varargin{:})
    condition = parser.Results.condition;
    trials = parser.Results.trials;
    
    % list of conditions
    if strcmp(condition, 'all')
        conditions_to_analyze = condition_list;
    else
        if iscell(condition)
            conditions_to_analyze = condition;
        else 
            conditions_to_analyze = {condition};
        end
    end
    
%     % check if a list was provided, or just a single entry
%     if iscell(condition)
%         if iscell(trials)
%            
%         else
%             
%         end
%     end
    
    
    
    % list of trials for each condition
    if iscell(trials)
        trials_to_analyze = trials;
    else
        if trials == 0;
            % find list of available trials for this condition
            trials_to_analyze = cell(size(conditions_to_analyze));
            for i_condition = 1 : length(conditions_to_analyze)
                condition_label = conditions_to_analyze{i_condition};
                condition_index_in_complete_list = find(strcmp(condition_list, condition_label));
                trials_to_analyze_in_this_condition = trial_number_list{condition_index_in_complete_list};
                trials_to_analyze{i_condition} = trials_to_analyze_in_this_condition;
            end
        else
            trials_to_analyze = {trials};
        end
    end
        
%         trials_to_analyze = trial_number_list(strcmp(condition_list, condition));
    % was there one list of trials given for several conditions?
    if length(conditions_to_analyze) > 1 && length(trials_to_analyze) == 1
        requested_trials = trials_to_analyze{1};
        trials_to_analyze = cell(size(conditions_to_analyze));
        for i_condition = 1 : length(conditions_to_analyze)
            condition_label = conditions_to_analyze{i_condition};
            condition_index_in_complete_list = find(strcmp(condition_list, condition_label));
            trials_to_analyze_in_this_condition = [];
            available_trials_in_this_condition = trial_number_list{condition_index_in_complete_list};


            % only choose trials that are not present in this condition
            for i_trial = 1 : length(requested_trials)
                if ismember(requested_trials(i_trial), available_trials_in_this_condition)
                    trials_to_analyze_in_this_condition = [trials_to_analyze_in_this_condition requested_trials(i_trial)];
                end
            end
            trials_to_analyze{i_condition} = trials_to_analyze_in_this_condition;
        end
    end
        
    
    
    % exclude calibration
    trials_to_analyze(strcmp(conditions_to_analyze, 'calibration')) = [];
    conditions_to_analyze(strcmp(conditions_to_analyze, 'calibration')) = [];

    % exclude emg
    trials_to_analyze(strcmp(conditions_to_analyze, 'emg')) = [];
    conditions_to_analyze(strcmp(conditions_to_analyze, 'emg')) = [];
    
    
    
    



%     if isempty(argin)
%         conditions_to_analyze = condition_list;
%         trials_to_analyze = trial_number_list;
%         
%         % exclude calibration
%         if any(strcmp(conditions_to_analyze, 'calibration'))
%             trials_to_analyze(strcmp(conditions_to_analyze, 'calibration')) = [];
%             conditions_to_analyze(strcmp(conditions_to_analyze, 'calibration')) = [];
%         end
%         
%         % exclude emg
%         if any(strcmp(conditions_to_analyze, 'emg'))
%             trials_to_analyze(strcmp(conditions_to_analyze, 'emg')) = [];
%             conditions_to_analyze(strcmp(conditions_to_analyze, 'emg')) = [];
%         end
%         
%         
%     end
%     if length(argin) == 1
%         if iscell(argin{1})
%             conditions_to_analyze = argin{1};
%         else
%             conditions_to_analyze = argin(1);
%         end
%         trials_to_analyze = cell(size(conditions_to_analyze));
%         for i_condition = 1 : length(conditions_to_analyze)
%             trials_to_analyze(i_condition) = trial_number_list(i_condition);
%         end
%     end
%     if length(argin) == 2
%         % if the first argument is not a cell, make it one
%         if iscell(argin{1})
%             conditions_to_analyze = argin{1};
%         else
%             conditions_to_analyze = argin(1);
%         end
%         if length(conditions_to_analyze) ~= length(argin{2})
%             error('list of conditions and list of list of trials to analyze must have the same length');
%         end
%         if iscell(argin{2})
%             trials_to_analyze = argin{2};
%         else
%             trials_to_analyze = argin(2);
%         end
%     end
%     
%     for i_condition = 1 : length(conditions_to_analyze)
%         if iscolumn(trials_to_analyze(i_condition))
%             trials_to_analyze{i_condition} = trials_to_analyze{i_condition}';
%         end
%     end
    
end