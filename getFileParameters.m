function [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(fileName)
    elements = strsplit(fileName, '_');
    date = elements{1};
    subject_id = elements{2};
    trial_type = elements{3};
    trial_number = str2double(elements{4}(1:3));
    if length(elements) > 4
        file_type = elements{5};
    else
        file_type = '';
    end

%     underscores = find(fileName == '_' | fileName == '-'); % Location of underscores
%     date = fileName(1:underscores(1)-1); % date
%     subject_id = fileName(underscores(1)+1:underscores(2)-1);
%     trial_type = fileName(underscores(2)+1:underscores(3)-1);
%     trial_number = str2double(fileName(underscores(3)+1:underscores(4)-1));
%     file_type = fileName(underscores(4)+1:length(fileName)-4);
return
