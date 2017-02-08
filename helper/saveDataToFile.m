%     This file is part of the CoBaL code base
%     Copyright (C) 2017 Hendrik Reimann <hendrikreimann@gmail.com>
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


function saveDataToFile(data_file_name_with_path, variables_to_save)

    if exist(data_file_name_with_path, 'file')
        variables_already_saved = load(data_file_name_with_path);
        
        variables_already_saved_names = fieldnames(variables_already_saved);
        variables_to_save_names = fieldnames(variables_to_save);
        for i_variable = 1 : length(variables_already_saved_names);
            % check if the i-th variable loaded from disk will be overwritten or not
            this_variable_name = variables_already_saved_names{i_variable};
            if ~any(strcmp(this_variable_name, variables_to_save_names))
                evalstring = ['variables_to_save.' this_variable_name ' = variables_already_saved.' this_variable_name ';'];
                eval(evalstring);
            end
        end
        
    end


    save(data_file_name_with_path, '-struct', 'variables_to_save');


end