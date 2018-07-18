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

function [types_to_analyze, trials_to_analyze, types_to_exclude, trials_to_exclude] = parseTrialArguments(varargin)
    load('subjectInfo.mat', 'condition_list', 'trial_number_list', 'trial_types_to_ignore');

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
        types_to_analyze = condition_list;
    else
        if iscell(condition)
            types_to_analyze = condition;
        else 
            types_to_analyze = {condition};
        end
    end
    
    % list of trials for each condition
    if iscell(trials)
        trials_to_analyze = trials;
    else
        if trials == 0
            % find list of available trials for this condition
            trials_to_analyze = cell(size(types_to_analyze));
            for i_condition = 1 : length(types_to_analyze)
                condition_label = types_to_analyze{i_condition};
                condition_index_in_complete_list = find(strcmp(condition_list, condition_label));
                trials_to_analyze_in_this_condition = trial_number_list{condition_index_in_complete_list};
                trials_to_analyze{i_condition} = trials_to_analyze_in_this_condition;
            end
        else
            trials_to_analyze = {trials};
        end
    end
        
    % was there one list of trials given for several conditions?
    if length(types_to_analyze) > 1 && length(trials_to_analyze) == 1
        requested_trials = trials_to_analyze{1};
        trials_to_analyze = cell(size(types_to_analyze));
        for i_condition = 1 : length(types_to_analyze)
            condition_label = types_to_analyze{i_condition};
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
        
    trials_to_exclude = {};
    types_to_exclude = {};
    for i_type = 1 : length(trial_types_to_ignore)
        this_type_label = trial_types_to_ignore{i_type};
        this_type_trials = trials_to_analyze(strcmp(types_to_analyze, this_type_label));
        if ~isempty(this_type_trials)
            trials_to_exclude = [trials_to_exclude; this_type_trials];
            types_to_exclude = [types_to_exclude; this_type_label];
            trials_to_analyze(strcmp(types_to_analyze, this_type_label)) = [];
            types_to_analyze(strcmp(types_to_analyze, this_type_label)) = [];
        end
    end
    
    
    
% %     % exclude calibration
%     calibration_trials = trials_to_analyze(strcmp(types_to_analyze, 'calibration'));
%     trials_to_analyze(strcmp(types_to_analyze, 'calibration')) = [];
%     types_to_analyze(strcmp(types_to_analyze, 'calibration')) = [];
% 
%     % exclude dance
%     dance_trials = trials_to_analyze(strcmp(types_to_analyze, 'dance'));
%     trials_to_analyze(strcmp(types_to_analyze, 'dance')) = [];
%     types_to_analyze(strcmp(types_to_analyze, 'dance')) = [];
% %     
%     % exclude emg
%     emg_trials = trials_to_analyze(strcmp(types_to_analyze, 'emg'));
%     trials_to_analyze(strcmp(types_to_analyze, 'emg')) = [];
%     types_to_analyze(strcmp(types_to_analyze, 'emg')) = [];
% %     
%     % exclude metronome
%     emg_trials = trials_to_analyze(strcmp(types_to_analyze, 'metronome'));
%     trials_to_analyze(strcmp(types_to_analyze, 'metronome')) = [];
%     types_to_analyze(strcmp(types_to_analyze, 'metronome')) = [];
%     
%     % exclude washout
%     emg_trials = trials_to_analyze(strcmp(types_to_analyze, 'washout'));
%     trials_to_analyze(strcmp(types_to_analyze, 'washout')) = [];
%     types_to_analyze(strcmp(types_to_analyze, 'washout')) = [];
%     
%     % exclude adaptation
%     emg_trials = trials_to_analyze(strcmp(types_to_analyze, 'adaptation'));
%     trials_to_analyze(strcmp(types_to_analyze, 'adaptation')) = [];
%     types_to_analyze(strcmp(types_to_analyze, 'adaptation')) = [];
%     types_to_analyze(strcmp(types_to_analyze, 'adaptation')) = [];
    
end