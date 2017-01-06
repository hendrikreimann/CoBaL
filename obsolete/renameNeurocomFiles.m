
% rename all the neurocom files in the current directory to the naming scheme used by CoBaL

% hard-code information here
date = '20160311';
subject_id = 'IC';
trial_type = 'walking';

% prepare folder
if ~exist('misnamed', 'dir')
    mkdir('misnamed')
end

% find files
clear file_name_list;
data_dir = dir('*.txt');
[file_name_list{1:length(data_dir)}] = deal(data_dir.name);

number_of_header_lines = 26;
trial_number_line = 25;
for i_file = 1 : length(file_name_list)
    data_file_name = file_name_list{i_file};
    
    % read string with trial number information
    fid = fopen(data_file_name);
    for i_line = 1 : 25
        line_string = fgetl(fid);
    end
    fclose(fid);
    
    % extract trial number
    line_split = strsplit(line_string, ' ');
    trial_number = str2num(line_split{end});

    % save file
    new_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'neurocomData.txt');
    while exist([filesep new_file_name], 'file')
        disp([new_file_name ' already exists, adding ".new"']);
        new_file_name = [new_file_name '.new'];
    end
    movefile(data_file_name, ['misnamed' filesep data_file_name])
    copyfile(['misnamed' filesep data_file_name], new_file_name)
end


