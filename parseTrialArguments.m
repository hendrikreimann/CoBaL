function [conditions_to_analyze, trials_to_analyze] = parseTrialArguments(varargin)
    load('subjectInfo.mat', 'condition_list', 'trial_number_list');

    parser = inputParser;
    parser.KeepUnmatched = true;
    
    default_condition = 'all';
    addParameter(parser, 'condition', default_condition, @ischar)
    default_trials = 0;
    addParameter(parser, 'trials', default_trials, @isnumeric)
    
    
    parse(parser, varargin{:})
    condition = parser.Results.condition;
    trials = parser.Results.trials;
    
    % list of conditions
    if strcmp(condition, 'all')
        conditions_to_analyze = condition_list;
    else
        conditions_to_analyze = {condition};
    end
    
    % list of trials for each condition
    if trials == 0;
        % find list of available trials for this condition
        trials_to_analyze = cell(size(conditions_to_analyze));
        for i_condition = 1 : length(conditions_to_analyze)
            condition_label = conditions_to_analyze{i_condition};
            condition_index_in_complete_list = find(strcmp(condition_list, condition_label));
            trials_to_analyze_in_this_condition = trial_number_list{condition_index_in_complete_list};
            trials_to_analyze{i_condition} = trials_to_analyze_in_this_condition;
        end        
        
%         trials_to_analyze = trial_number_list(strcmp(condition_list, condition));
    else
        trials_to_analyze = {trials};
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