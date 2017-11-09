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

% consolidate data from all subjects

function collectPopulationResults(varargin)

    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'subjects', [])
    parse(parser, varargin{:})
    subjects = parser.Results.subjects;

    % load settings
    if ~exist('studySettings.txt', 'file')
        error('No studySettings.txt file found. This function should be run from a study folder')
    end    
    study_settings = SettingsCustodian('studySettings_TF.txt');
    
    data_folder_list = determineDataStructure(subjects);

    %% collect data from all data folders
    variables_to_collect = study_settings.get('variables_to_collect');
    
    number_of_variables_to_collect = size(variables_to_collect, 1);
    condition_stance_foot_list = {};
    condition_perturbation_list = {};
    condition_delay_list = {};
    condition_index_list = {};
    condition_experimental_list = {};
    condition_stimulus_list = {};
    condition_day_list = {};
    subject_list = {};
    origin_trial_list = [];
    origin_start_time_list = [];
    origin_end_time_list = [];
    time_list = [];
    variable_data = cell(number_of_variables_to_collect, 1);
    response_data = cell(number_of_variables_to_collect, 1);
   
    for i_folder = 1 : length(data_folder_list)
        % load data
        data_path = data_folder_list{i_folder};
        load([data_path filesep 'subjectInfo.mat'], 'date', 'subject_id');
        disp(['Collecting from ' subject_id]);
        load([data_path filesep 'analysis' filesep date '_' subject_id '_results.mat']);

        % append data from this subject to containers for all subjects
        condition_stance_foot_list = [condition_stance_foot_list; condition_stance_foot_list_session]; %#ok<AGROW>
        condition_perturbation_list = [condition_perturbation_list; condition_perturbation_list_session]; %#ok<AGROW>
        condition_delay_list = [condition_delay_list; condition_delay_list_session]; %#ok<AGROW>
        condition_index_list = [condition_index_list; condition_index_list_session]; %#ok<AGROW>
        condition_experimental_list = [condition_experimental_list; condition_experimental_list_session]; %#ok<AGROW>
        condition_stimulus_list = [condition_stimulus_list; condition_stimulus_list_session]; %#ok<AGROW>
        condition_day_list = [condition_day_list; condition_day_list_session]; %#ok<AGROW>
        [subject_list{length(subject_list)+(1 : length(condition_stance_foot_list_session))}] = deal(subject_id); %#ok<AGROW>
%         origin_trial_list = [origin_trial_list; origin_trial_list_session]; %#ok<AGROW>
%         origin_start_time_list = [origin_start_time_list; origin_start_time_list_session]; %#ok<AGROW>
%         origin_end_time_list = [origin_end_time_list; origin_end_time_list_session]; %#ok<AGROW>
%         time_list = [time_list; time_list_session]; %#ok<AGROW>
        
        % what is this section being used for?
        for i_variable = 1 : number_of_variables_to_collect
            % load and extract data
            this_variable_name = variables_to_collect{i_variable, 1};
            index_in_saved_data = find(strcmp(variable_names_session, this_variable_name), 1, 'first');
            this_variable_data = variable_data_session{index_in_saved_data}; %#ok<USENS>
            this_response_data = response_data_session{index_in_saved_data}; %#ok<USENS>
            
            % store
            variable_data{i_variable} = [variable_data{i_variable} this_variable_data];
            response_data{i_variable} = [response_data{i_variable} this_response_data];
        end
        
    end
    population_data_absolute = variable_data;
    population_data = response_data;
    subject_list = subject_list';
    
    % add a condition_group_list
    condition_group_assignment = study_settings.get('condition_group_assignment');
    condition_group_list = cell(size(condition_stimulus_list));
    
    %% make relative illusion condition list
    condition_stimulus_list = condition_perturbation_list;
    for i_stretch = 1 : length(condition_stimulus_list)
        if ...
          (strcmp(condition_index_list{i_stretch}, 'ONE') && strcmp(condition_stimulus_list{i_stretch}, 'ILLUSION_LEFT') && strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')) ...
          || (strcmp(condition_index_list{i_stretch}, 'TWO') && strcmp(condition_stimulus_list{i_stretch}, 'ILLUSION_LEFT') && strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')) ...
          || (strcmp(condition_index_list{i_stretch}, 'THREE') && strcmp(condition_stimulus_list{i_stretch}, 'ILLUSION_LEFT') && strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')) ...
          || (strcmp(condition_index_list{i_stretch}, 'FOUR') && strcmp(condition_stimulus_list{i_stretch}, 'ILLUSION_LEFT') && strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')) ...
          || (strcmp(condition_index_list{i_stretch}, 'ONE') && strcmp(condition_stimulus_list{i_stretch}, 'ILLUSION_RIGHT') && strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')) ...
          || (strcmp(condition_index_list{i_stretch}, 'TWO') && strcmp(condition_stimulus_list{i_stretch}, 'ILLUSION_RIGHT') && strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')) ...
          || (strcmp(condition_index_list{i_stretch}, 'THREE') && strcmp(condition_stimulus_list{i_stretch}, 'ILLUSION_RIGHT') && strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')) ...
          || (strcmp(condition_index_list{i_stretch}, 'FOUR') && strcmp(condition_stimulus_list{i_stretch}, 'ILLUSION_RIGHT') && strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT'))
            % these are all the cases where the illusion is TOWARDS the first stance leg, i.e. the triggering leg
            condition_stimulus_list{i_stretch} = 'TOWARDS';
        elseif ...
          (strcmp(condition_index_list{i_stretch}, 'ONE') && strcmp(condition_stimulus_list{i_stretch}, 'ILLUSION_LEFT') && strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')) ...
          || (strcmp(condition_index_list{i_stretch}, 'TWO') && strcmp(condition_stimulus_list{i_stretch}, 'ILLUSION_LEFT') && strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')) ...
          || (strcmp(condition_index_list{i_stretch}, 'THREE') && strcmp(condition_stimulus_list{i_stretch}, 'ILLUSION_LEFT') && strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')) ...
          || (strcmp(condition_index_list{i_stretch}, 'FOUR') && strcmp(condition_stimulus_list{i_stretch}, 'ILLUSION_LEFT') && strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')) ...
          || (strcmp(condition_index_list{i_stretch}, 'ONE') && strcmp(condition_stimulus_list{i_stretch}, 'ILLUSION_RIGHT') && strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')) ...
          || (strcmp(condition_index_list{i_stretch}, 'TWO') && strcmp(condition_stimulus_list{i_stretch}, 'ILLUSION_RIGHT') && strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')) ...
          || (strcmp(condition_index_list{i_stretch}, 'THREE') && strcmp(condition_stimulus_list{i_stretch}, 'ILLUSION_RIGHT') && strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')) ...
          || (strcmp(condition_index_list{i_stretch}, 'FOUR') && strcmp(condition_stimulus_list{i_stretch}, 'ILLUSION_RIGHT') && strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT'))
            % these are all the cases where the illusion is AWAY from the first stance leg, i.e. the triggering leg
            condition_stimulus_list{i_stretch} = 'AWAY';
        elseif strcmp(condition_stimulus_list{i_stretch}, 'CONTROL')
            % do nothing
        else
            error('Something wrong with the condition: No match found')
        end
        
        % find the row that the subject of this stretch matches in the foot
        group_assignment_subject_row = find(strcmp(subject_list{i_stretch}, condition_group_assignment(:,1)));
        
           
        group_assignment_foot_col = find(strcmp(condition_stance_foot_list{i_stretch}, condition_group_assignment(group_assignment_subject_row,:)));
        if isempty(group_assignment_foot_col)
            group_assignment_foot_col = find(strcmp(condition_group_assignment(group_assignment_subject_row,:), 'BOTH'));
        end
        
        if group_assignment_foot_col == 2
            condition_group_list{i_stretch} = 'EARLY';
        elseif group_assignment_foot_col == 3
            condition_group_list{i_stretch} = 'LATE';           
        elseif group_assignment_foot_col == 4
            condition_group_list{i_stretch} = 'NO';
        else
            error('Something wrong with the group matching: No match found')
        end

    end
    


    %% make time categorical variable
%     time_category = cell(size(time_list));
%     time_category_borders = study_settings.get('time_category_borders');
%     for i_point = 1 : length(time_list)
%         if time_category_borders(1) < time_list(i_point) && time_list(i_point) < time_category_borders(2)
%             time_category{i_point} = 'TIME_ONE';
%         end
%         if time_category_borders(2) < time_list(i_point) && time_list(i_point) < time_category_borders(3)
%             time_category{i_point} = 'TIME_TWO';
%         end
%         if time_category_borders(3) < time_list(i_point) && time_list(i_point) < time_category_borders(4)
%             time_category{i_point} = 'TIME_THREE';
%         end
%         if time_category_borders(4) < time_list(i_point) && time_list(i_point) < time_category_borders(5)
%             time_category{i_point} = 'TIME_FOUR';
%         end
%     end
    

    
    
    %% save
    variables_to_save = struct;
    variables_to_save.population_data = population_data;
    variables_to_save.population_data_absolute = population_data_absolute;
    variables_to_save.variable_names = variables_to_collect(:, 1);
    variables_to_save.condition_stance_foot_list = condition_stance_foot_list;
    variables_to_save.condition_perturbation_list = condition_perturbation_list;
    variables_to_save.condition_stimulus_list = condition_perturbation_list;
    variables_to_save.condition_delay_list = condition_delay_list;
    variables_to_save.condition_index_list = condition_index_list;
    variables_to_save.condition_experimental_list = condition_experimental_list;
    variables_to_save.condition_stimulus_list = condition_stimulus_list;
    variables_to_save.condition_day_list = condition_day_list;
    variables_to_save.condition_group_list = condition_group_list;
%     variables_to_save.time_list = time_list;
%     variables_to_save.time_category = time_category;
    variables_to_save.subject_list = subject_list;
    save('results', '-struct', 'variables_to_save');





end