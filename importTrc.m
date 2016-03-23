% transform raw data from Cortex into matlab data
current_directory = pwd;
data_dir = dir('*.trc');
clear file_name_list;
[file_name_list{1:length(data_dir)}] = deal(data_dir.name);

number_of_files = length(file_name_list);

unit_change = 1e-3; % millimeter to meter

for i_trial = 1 : number_of_files
    % file name stuff
    trc_data_file_name = file_name_list{i_trial};
    
    % import data
    [A, delimiter, nheaderlines] = importdata(trc_data_file_name, '\t', 5);
    marker_trajectories_raw = A.data(:, 3:end) * unit_change;
    
    % save
    matlab_data_file_name = [trc_data_file_name(1 : end-4) '_rawTrajectories.mat'];
    save(matlab_data_file_name, 'marker_trajectories_raw');
    
    disp(['imported ' trc_data_file_name])
end
disp(['imported ' num2str(number_of_files) ' files'])

% save([data_root directorySeparator 'markerData_raw.mat'], 'marker_trajectories_raw', 'file_name_list');


