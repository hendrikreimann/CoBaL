


function deleteAvailableVariables

        data_dir = dir(['analysis' filesep '*_availableVariables.mat']);
        clear file_name_list;
        [file_name_list{1:length(data_dir)}] = deal(data_dir.name);
        number_of_files = length(file_name_list);
        
        for i_file = 1 : number_of_files
            this_file = ['analysis' filesep file_name_list{i_file}];
            delete(this_file);
        end
        
end