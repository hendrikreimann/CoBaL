% prune ArmSense for duplicate files

folder = pwd;
contents = dir;
contents = contents(3:end);

file_names = {contents.name}';
bytes = [contents.bytes]';

% make list of file names and sizes for tsv and tsva
file_names_tsv = {};
file_names_tsva = {};
file_sizes_tsv = [];
file_sizes_tsva = [];
for i_file = 1 : length(file_names)
    this_file_name = file_names{i_file};
    this_file_size = bytes(i_file);
    if strcmp(this_file_name(end-1:end), 'sv')
        if strcmp(this_file_name(end-4), 'a')
            file_names_tsva = [file_names_tsva; this_file_name];
            file_sizes_tsva = [file_sizes_tsva; this_file_size];
        else
            file_names_tsv = [file_names_tsv; this_file_name];
            file_sizes_tsv = [file_sizes_tsv; this_file_size];
        end
    end
end

% go through lists and look for duplicates
if length(file_names_tsv) ~= length(file_names_tsva)
    error('Number of files does not add up, something weird is going on')
end
file_names_tsv_pruned = {};
file_names_tsva_pruned = {};
file_sizes_tsv_pruned = [];
file_sizes_tsva_pruned = [];
for i_file = 1 : length(file_names_tsv)
    this_file_name_tsv = file_names_tsv{i_file};
    this_file_size_tsv = file_sizes_tsv(i_file);
    this_file_name_tsva = file_names_tsva{i_file};
    this_file_size_tsva = file_sizes_tsva(i_file);
    
    % check if we already have files with this size
    if ~isempty(file_sizes_tsv_pruned) && any(this_file_size_tsv==file_sizes_tsv_pruned & this_file_size_tsva==file_sizes_tsva_pruned)
        disp(['Duplicate found: Trial index ' num2str(i_file)])
    else
        file_names_tsv_pruned = [file_names_tsv_pruned; this_file_name_tsv];
        file_sizes_tsv_pruned = [file_sizes_tsv_pruned; this_file_size_tsv];
        file_names_tsva_pruned = [file_names_tsva_pruned; this_file_name_tsva];
        file_sizes_tsva_pruned = [file_sizes_tsva_pruned; this_file_size_tsva];
    end
end    
    
% go through pruned lists and extract condition and trial number information
condition_list = {};
trial_number_list = [];
for i_trial = 1 : length(file_names_tsv_pruned)
    this_file_name = file_names_tsv_pruned{i_trial};
    underscores = find(this_file_name == '_'); % Location of underscores
    subject_id = this_file_name(underscores(1)+1:underscores(2)-1); % study
    condition_label = this_file_name(underscores(2)+1:underscores(3)-1); % initials
    trial_number_string = this_file_name(underscores(3)+1:end-4); % number
    condition_list = [condition_list; condition_label];
    trial_number_list = [trial_number_list; str2num(trial_number_string)];
end

conditions = unique(condition_list);
trials = {};
for i_condition = 1 : length(conditions)
    this_condition = conditions{i_condition};
    this_condition_trial_list = trial_number_list(strcmp(condition_list, this_condition));
    trials = [trials; this_condition_trial_list'];
end

% replace this information in subjectInfo.mat
subject_info = load(['..' filesep 'subjectInfo.mat']);
subject_info.condition_list = conditions;
subject_info.trial_number_list = trials;
if ~exist(['..' filesep 'subjectInfo_bak.mat'], 'file')
    copyfile(['..' filesep 'subjectInfo.mat'], ['..' filesep 'subjectInfo_bak.mat']);
end
save(['..' filesep 'subjectInfo.mat'], '-struct', 'subject_info');



















