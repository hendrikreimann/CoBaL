
% specify data source
folder = '/Users/reimajbi/TempleDrive/Obstacle/SS/';
conditions_source_file = 'conditions.mat';

% load
load([folder conditions_source_file]);

% transform
header = {'trial', 'experimental'};
condition_map = {};
for i_trial = 1 : length(no_obstacle)
    new_entry = {no_obstacle(i_trial), 'OBS_NO'};
    condition_map = [condition_map; new_entry];
end
for i_trial = 1 : length(near_obstacle)
    new_entry = {near_obstacle(i_trial), 'OBS_NEAR'};
    condition_map = [condition_map; new_entry];
end
for i_trial = 1 : length(far_obstacle)
    new_entry = {far_obstacle(i_trial), 'OBS_FAR'};
    condition_map = [condition_map; new_entry];
end
condition_map_sorted = sortrows(condition_map,1);

% write to csv
fileID = fopen([folder 'conditions.csv'], 'w');
fprintf(fileID, '%s\n', 'trial, experimental');
for i_line = 1 : size(condition_map_sorted, 1)
    fprintf(fileID, '%i', condition_map_sorted{i_line, 1});
    fprintf(fileID, '%s', ',');
    fprintf(fileID, '%s\n', condition_map_sorted{i_line, 2});
end
fclose(fileID);