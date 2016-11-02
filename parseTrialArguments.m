function [conditions_to_analyze, trials_to_analyze] = parseTrialArguments(argin)

    load('subjectInfo.mat', 'condition_list', 'trial_number_list');
    if isempty(argin)
        conditions_to_analyze = condition_list;
        trials_to_analyze = trial_number_list;
        
        % exclude calibration
        if any(strcmp(conditions_to_analyze, 'calibration'))
            trials_to_analyze(strcmp(conditions_to_analyze, 'calibration')) = [];
            conditions_to_analyze(strcmp(conditions_to_analyze, 'calibration')) = [];
        end
        
        % exclude emg
        if any(strcmp(conditions_to_analyze, 'emg'))
            trials_to_analyze(strcmp(conditions_to_analyze, 'emg')) = [];
            conditions_to_analyze(strcmp(conditions_to_analyze, 'emg')) = [];
        end
        
        
    end
    if length(argin) == 1
        if iscell(argin{1})
            conditions_to_analyze = argin{1};
        else
            conditions_to_analyze = argin(1);
        end
        trials_to_analyze = cell(size(conditions_to_analyze));
        for i_condition = 1 : length(conditions_to_analyze)
            trials_to_analyze(i_condition) = trial_number_list(i_condition);
        end
    end
    if length(argin) == 2
        % if the first argument is not a cell, make it one
        if iscell(argin{1})
            conditions_to_analyze = argin{1};
        else
            conditions_to_analyze = argin(1);
        end
        if length(conditions_to_analyze) ~= length(argin{2})
            error('list of conditions and list of list of trials to analyze must have the same length');
        end
        if iscell(argin{2})
            trials_to_analyze = argin{2};
        else
            trials_to_analyze = argin(2);
        end
    end
    
    for i_condition = 1 : length(conditions_to_analyze)
        if iscolumn(trials_to_analyze(i_condition))
            trials_to_analyze{i_condition} = trials_to_analyze{i_condition}';
        end
    end
    
end