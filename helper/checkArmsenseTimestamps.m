% check ArmSense data for duplicates

%% get lists
% go through all available folders with data and extract the file sizes
folder_list = ...
  { ...
    '/Users/hendrikreimann/Dropbox/ArmSense_intervention/101/day1/ascii', ...
    '/Users/hendrikreimann/Dropbox/ArmSense_intervention/101/day2/ascii', ...
    '/Users/hendrikreimann/Dropbox/ArmSense_intervention/102/day1/ascii', ...
    '/Users/hendrikreimann/Dropbox/ArmSense_intervention/102/day2/ascii', ...
    '/Users/hendrikreimann/Dropbox/ArmSense_intervention/103/day1/ascii', ...
    '/Users/hendrikreimann/Dropbox/ArmSense_intervention/103/day2/ascii', ...
    '/Users/hendrikreimann/Dropbox/ArmSense_intervention/107/day1/ascii', ...
    '/Users/hendrikreimann/Dropbox/ArmSense_intervention/201/day1/ascii', ...
    '/Users/hendrikreimann/Dropbox/ArmSense_intervention/201/day2/ascii', ...
    '/Users/hendrikreimann/Dropbox/ArmSense_intervention/202/day1/ascii', ...
    '/Users/hendrikreimann/Dropbox/ArmSense_intervention/202/day2/ascii', ...
    '/Users/hendrikreimann/Dropbox/ArmSense_intervention/204/day1/ascii', ... % no analog data
    '/Users/hendrikreimann/Dropbox/ArmSense_intervention/204/day2/ascii', ...
    '/Users/hendrikreimann/Dropbox/ArmSense_intervention/207/day1/ascii', ...
    '/Users/hendrikreimann/Dropbox/ArmSense_intervention/207/day2/ascii', ...
    '/Users/hendrikreimann/Dropbox/ArmSense_intervention/208/day1/ascii' ...
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
                file_srces_analog = [file_srces_analog; folder]; %#ok<AGROW>
            else
                file_names_marker = [file_names_marker; this_file_name]; %#ok<AGROW>
                file_srces_marker = [file_srces_marker; folder]; %#ok<AGROW>
            end
        end
    end
end

%% gather time stamp info for marker
file_table_marker = cell(length(file_names_marker), 7);
for i_file = 1 : length(file_names_marker)
    this_file_name = file_names_marker{i_file};
    this_file_srce = file_srces_marker{i_file};
    
    % get date info from file name
    file_name_split = strsplit(this_file_name, '_');
    date_string_from_filename = file_name_split{1};
    file_table_marker{i_file, 1} = date_string_from_filename;
    
    % get date info from time stamp
    this_file = [this_file_srce filesep this_file_name];
    if strcmp(this_file(end-2:end), 'tsv')
        delimiter = '\t';
        number_of_header_lines = 10;
        [imported_data, delimiter, number_of_header_lines_returned] = importdata(this_file, delimiter, number_of_header_lines);
        header = imported_data.textdata;
        time_stamp_line = header{8};
        time_stamp_split = strsplit(time_stamp_line, ',');
        date_string_from_ascii = time_stamp_split{1}([12 13 14 15 17 18 20 21]);
        time_string = time_stamp_split{2}(2:9);
        file_table_marker{i_file, 2} = date_string_from_ascii;
        file_table_marker{i_file, 3} = time_string;

        % get info for subject and day
        this_file_name_split = strsplit(this_file, filesep);
        this_file_subject = this_file_name_split{6};
        this_file_day = this_file_name_split{7};
        file_table_marker{i_file, 4} = this_file_subject;
        file_table_marker{i_file, 5} = this_file_day;

        file_table_marker{i_file, 6} = this_file_name;
        file_table_marker{i_file, 7} = this_file_srce;
    end    
end

% save
file_name = '/Users/hendrikreimann/Dropbox/ArmSense_intervention/timestamps_marker.mat';
save(file_name, 'file_table_marker');


%% gather time stamp info for analog
file_table_analog = cell(length(file_names_analog), 7);
for i_file = 1 : length(file_names_analog)
    this_file_name = file_names_analog{i_file};
    this_file_srce = file_srces_analog{i_file};
    
    % get date info from file name
    file_name_split = strsplit(this_file_name, '_');
    date_string_from_filename = file_name_split{1};
    file_table_analog{i_file, 1} = date_string_from_filename;
    
    % get date info from time stamp
    this_file = [this_file_srce filesep this_file_name];
    if strcmp(this_file(end-2:end), 'tsv')
        delimiter = '\t';
        number_of_header_lines = 13;
        [imported_data, delimiter, number_of_header_lines_returned] = importdata(this_file, delimiter, number_of_header_lines);
        header = imported_data.textdata;
        time_stamp_line = header{5};
        time_stamp_split = strsplit(time_stamp_line, ',');
        date_string_from_ascii = time_stamp_split{1}([12 13 14 15 17 18 20 21]);
        time_string = time_stamp_split{2}(2:9);
        file_table_analog{i_file, 2} = date_string_from_ascii;
        file_table_analog{i_file, 3} = time_string;

        % get info for subject and day
        this_file_name_split = strsplit(this_file, filesep);
        this_file_subject = this_file_name_split{6};
        this_file_day = this_file_name_split{7};
        file_table_analog{i_file, 4} = this_file_subject;
        file_table_analog{i_file, 5} = this_file_day;

        file_table_analog{i_file, 6} = this_file_name;
        file_table_analog{i_file, 7} = this_file_srce;
    end
end

% save
file_name = '/Users/hendrikreimann/Dropbox/ArmSense_intervention/timestamps_analog.mat'
save(file_name, 'file_table_analog');










