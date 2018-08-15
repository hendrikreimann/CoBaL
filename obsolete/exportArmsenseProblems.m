mode = 'marker';

load_file = ['matches_' mode '.mat'];
load(load_file);
eval(['data = file_table_' mode ';']);
save_file = [load_file(1:end-3), 'txt'];

fid = fopen(save_file,'wt');
for i_entry = 1 : length(data)
    this_entry_file_names = data{i_entry, 2};
    this_entry_match_percentages = data{i_entry, 3};
    
    for i_file = 1 : length(this_entry_file_names)
        % extract file information
        this_file_name_with_path = this_entry_file_names{i_file};
        this_file_name_split = strsplit(this_file_name_with_path, filesep);
        
        this_file_subject = this_file_name_split{6};
        this_file_day = this_file_name_split{7};
        this_file_name = this_file_name_split{end};
        
        % extract match information
        this_file_matches = this_entry_match_percentages{i_file};
        
        % generate line
        line_to_save = [this_file_subject, ',' this_file_day ',' this_file_name];
        for i_match = 1 : length(this_file_matches)
            line_to_save = [line_to_save ',' num2str(this_file_matches(i_match))];
        end
        line_to_save = [line_to_save '\n'];
        fprintf(fid, line_to_save);
    end
    fprintf(fid, '\n');
end
    


fclose(fid);

