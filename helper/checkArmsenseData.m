% check ArmSense data for duplicates

%% get lists
% go through all available folders with data and extract the file sizes
folder_list = ...
  { ...
    '/Users/reimajbi/Downloads/Armsense_New_Upload_Check/102/day1/empty_ascii' ...
    '/Users/reimajbi/Downloads/Armsense_New_Upload_Check/102/day2/empty_ascii' ...
    '/Users/reimajbi/Downloads/Armsense_New_Upload_Check/103/day1/empty_ascii' ...
    '/Users/reimajbi/Downloads/Armsense_New_Upload_Check/103/day2/empty_ascii' ...
    '/Users/reimajbi/Downloads/Armsense_New_Upload_Check/201/day1/empty_ascii' ...
    '/Users/reimajbi/Downloads/Armsense_New_Upload_Check/201/day2/empty_ascii' ...
    '/Users/reimajbi/Downloads/Armsense_New_Upload_Check/203/day1/empty_ascii' ...
    '/Users/reimajbi/Downloads/Armsense_New_Upload_Check/203/day2/empty_ascii' ...
    '/Users/reimajbi/Downloads/Armsense_New_Upload_Check/204/day1/empty_ascii' ...
    '/Users/reimajbi/Downloads/Armsense_New_Upload_Check/204/day2/empty_ascii' ...
    '/Users/reimajbi/Downloads/Armsense_New_Upload_Check/207/day1/empty_ascii' ...
    '/Users/reimajbi/Downloads/Armsense_New_Upload_Check/207/day2/empty_ascii' ...
  };

file_names_marker = {};
file_names_analog = {};
file_sizes_marker = [];
file_sizes_analog = [];
file_srces_marker = {};
file_srces_analog = {};
for i_folder = 1 : length(folder_list)
    folder = folder_list{i_folder};
    contents = dir(folder);
    contents = contents(3:end);

    file_names = {contents.name}';
    bytes = [contents.bytes]';

    
    % make list of file names and sizes
    for i_file = 1 : length(file_names)
        this_file_name = file_names{i_file};
        this_file_size = bytes(i_file);
        if strcmp(this_file_name(end-1:end), 'sv')
            if strcmp(this_file_name(end-4), 'a')
                file_names_analog = [file_names_analog; this_file_name]; %#ok<AGROW>
                file_sizes_analog = [file_sizes_analog; this_file_size]; %#ok<AGROW>
                file_srces_analog = [file_srces_analog; folder]; %#ok<AGROW>
            else
                file_names_marker = [file_names_marker; this_file_name]; %#ok<AGROW>
                file_sizes_marker = [file_sizes_marker; this_file_size]; %#ok<AGROW>
                file_srces_marker = [file_srces_marker; folder]; %#ok<AGROW>
            end
        end
    end
end

%% gather candidates for duplicates depending on file size
% go through lists and look for duplicates - marker
% unique_file_sizes = unique(file_sizes_marker);
% file_table_marker = cell(length(unique_file_sizes), 2);
% file_table_marker(:, 1) = num2cell(unique_file_sizes);
% for i_file = 1 : length(file_sizes_marker)
%     this_file_name = file_names_marker{i_file};
%     this_file_size = file_sizes_marker(i_file);
%     this_file_srce = file_srces_marker{i_file};
%     file_name_with_path = [this_file_srce filesep this_file_name];
%     
%     this_file_index = find(this_file_size==unique_file_sizes);
%     file_table_marker{this_file_index, 2} = [file_table_marker{this_file_index, 2}; {file_name_with_path}];
% end
% % remove entries without duplicates
% duplicate_counts_marker = zeros(size(file_table_marker, 1), 1);
% for i_file = 1 : size(file_table_marker, 1)
%     duplicate_counts_marker(i_file) = length(file_table_marker{i_file, 2});
% end
% file_table_marker(duplicate_counts_marker==1, :) = [];
% duplicate_counts_marker(duplicate_counts_marker==1, :) = [];

% go through lists and look for duplicates - analog
unique_file_sizes = unique(file_sizes_analog);
file_table_analog = cell(length(unique_file_sizes), 2);
file_table_analog(:, 1) = num2cell(unique_file_sizes);
for i_file = 1 : length(file_sizes_analog)
    this_file_name = file_names_analog{i_file};
    this_file_size = file_sizes_analog(i_file);
    this_file_srce = file_srces_analog{i_file};
    file_name_with_path = [this_file_srce filesep this_file_name];
    
    this_file_index = find(this_file_size==unique_file_sizes);
    file_table_analog{this_file_index, 2} = [file_table_analog{this_file_index, 2}; {file_name_with_path}];
end    
% remove entries without duplicates
duplicate_counts_analog = zeros(size(file_table_analog, 1), 1);
for i_file = 1 : size(file_table_analog, 1)
    duplicate_counts_analog(i_file) = length(file_table_analog{i_file, 2});
end
file_table_analog(duplicate_counts_analog==1, :) = [];
duplicate_counts_analog(duplicate_counts_analog==1, :) = [];

%% check all candidates to see whether they're actually duplicates
% file_table_marker = [file_table_marker cell(size(file_table_marker, 1), 1)];
% for i_entry = 1 : size(file_table_marker, 1)
%     % load all files and extract data from this entry
%     this_file_name_list = file_table_marker{i_entry, 2};
%     this_data = cell(length(this_file_name_list), 1);
%     for i_file = 1 : length(this_file_name_list)
%         this_file = this_file_name_list{i_file};
%         if strcmp(this_file(end-2:end), 'tsv')
%             delimiter = '\t';
%             number_of_header_lines = 10;
%         end
%         if strcmp(this_file(end-2:end), 'csv')
%             delimiter = ',';
%             number_of_header_lines = 5;
%         end
%         [imported_data, delimiter, number_of_header_lines_returned] = importdata(this_file, delimiter, number_of_header_lines);
%         this_data{i_file} = imported_data.data;
%     end
%     
%     % compare each files against all others
%     this_entry_match_percentages = cell(length(this_file_name_list), 1);
%     for i_file = 1 : length(this_file_name_list)
%         this_entry_match_percentages_row = zeros(1, length(this_file_name_list));
%         
%         this_file_data = this_data{i_file};
%         this_file_data(isnan(this_file_data)) = 0;
%         for j_file = 1 : length(this_file_name_list)
%             that_file_data = this_data{j_file};
%             that_file_data(isnan(that_file_data)) = 0;
%             if numel(this_file_data) == numel(that_file_data)
%                 number_of_matches = sum(sum(this_file_data == that_file_data));
%                 match_ratio = number_of_matches / numel(that_file_data);
%                 match_percentage = match_ratio * 100;
%             else
%                 match_percentage = 0;
%             end
%             this_entry_match_percentages_row(1, j_file) = match_percentage;
%         end
%         this_entry_match_percentages{i_file} = this_entry_match_percentages_row;
%     end
%     
%     % store
%     file_table_marker{i_entry, 3} = this_entry_match_percentages;
% end
% % save
% file_name = '/Users/reimajbi/Downloads/Armsense_New_Upload_Check/matches_marker.mat'
% save(file_name, 'file_table_marker');

%%
file_table_analog = [file_table_analog cell(size(file_table_analog, 1), 1)];
for i_entry = 1 : size(file_table_analog, 1)
    % load all files and extract data from this entry
    this_file_name_list = file_table_analog{i_entry, 2};
    this_data = cell(length(this_file_name_list), 1);
    for i_file = 1 : length(this_file_name_list)
        this_file = this_file_name_list{i_file};
        if strcmp(this_file(end-2:end), 'tsv')
            delimiter = '\t';
            number_of_header_lines = 13;
        end
        if strcmp(this_file(end-2:end), 'csv')
            delimiter = ',';
            number_of_header_lines = 5;
        end
        [imported_data, delimiter, number_of_header_lines_returned] = importdata(this_file, delimiter, number_of_header_lines);
        this_data{i_file} = imported_data.data;
    end
    
    % compare each files against all others
    this_entry_match_percentages = cell(length(this_file_name_list), 1);
    for i_file = 1 : length(this_file_name_list)
        this_entry_match_percentages_row = zeros(1, length(this_file_name_list));
        
        this_file_data = this_data{i_file};
        this_file_data(isnan(this_file_data)) = 0;
        for j_file = 1 : length(this_file_name_list)
            that_file_data = this_data{j_file};
            that_file_data(isnan(that_file_data)) = 0;
            if numel(this_file_data) == numel(that_file_data)
                number_of_matches = sum(sum(this_file_data == that_file_data));
                match_ratio = number_of_matches / numel(that_file_data);
                match_percentage = match_ratio * 100;
            else
                match_percentage = 0;
            end
            this_entry_match_percentages_row(1, j_file) = match_percentage;
        end
        this_entry_match_percentages{i_file} = this_entry_match_percentages_row;
    end
    
    % store
    file_table_analog{i_entry, 3} = this_entry_match_percentages;
end

% save
file_name = '/Users/reimajbi/Downloads/Armsense_New_Upload_Check/matches_analog.mat';
save(file_name, 'file_table_analog');
















