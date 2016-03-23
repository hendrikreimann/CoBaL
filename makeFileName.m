function fileName = makeFileName(date, subject_id, trial_type, trial_number, file_type)
    fileName = [date '_' subject_id '_' trial_type];
    if nargin > 3
        if isnumeric(trial_number)
            trial_number = zeroPrefixedIntegerString(trial_number, 3);
        end
        fileName = [fileName '_' trial_number];
    end
    if nargin > 4
        fileName = [fileName '_' file_type];
    end
    
%     fileName = [date '_' subject_id '_' trial_type '_' trial_number '_' file_type];
    
return
