mode = 'analog';

load_file = ['timestamps_' mode '.mat'];
load(load_file);
eval(['data = file_table_' mode ';']);
save_file = [load_file(1:end-3), 'txt'];

fid = fopen(save_file,'wt');
header_line = 'date from file name, date from timestamp, time from timestamp, subject, day, file name, path\n';
fprintf(fid, header_line);

for i_entry = 1 : length(data)
    row_to_save = data(i_entry, :);
    line_to_save = row_to_save{1};
    for i_cell = 2 : length(row_to_save)
        line_to_save = [line_to_save ',' row_to_save{i_cell}];
    end
    line_to_save = [line_to_save '\n'];
    fprintf(fid, line_to_save);
end
fclose(fid);

