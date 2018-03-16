% load data and settings
function onset_classification(varargin)

    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'subjects', [])
    parse(parser, varargin{:})
    subjects = parser.Results.subjects;
    
    % set an error if run from anywhere except study folder
    
    % remove the txt settings file. The variables are set in stone...

    test_settings = SettingsCustodian('onset_classification.txt');
    % initialize
    conditions_to_test = test_settings.get('conditions_to_test');
    conditions_control = test_settings.get('conditions_control');
    number_of_conditions_to_test = size(conditions_to_test, 1);
    factor_to_analyze = test_settings.get('factor_to_analyze');
    variable_to_test = test_settings.get('variables_to_test');

    data_folder_list = determineDataStructure(subjects);
    for i_folder = 1 : length(data_folder_list)
        % load data
        study_path = cd;
        data_path = data_folder_list{i_folder};
        cd(data_path)
        load([data_path filesep 'subjectInfo.mat'], 'date', 'subject_id');
        results_file_name = ['analysis' filesep makeFileName(date, subject_id, 'results')];
        loaded_data = load(results_file_name);
        disp(['Collecting from ' subject_id]);

        % preallocate new conditions
        number_of_stretch_variables = size(loaded_data.conditions_session.condition_trigger_foot_list,1);
        condition_response_time_list_session = cell(number_of_stretch_variables, 1);
        condition_response_type_list_session = cell(number_of_stretch_variables, 1);
        
        for i_condition = 1 : number_of_conditions_to_test % in this particular analysis, there should be two conditions (one for each system/trigger foot)            
            this_determine_onset_variable = variable_to_test; % in this particular analysis, only variable_to_test should be cop_from_com_x_inverted (maybe make this a bit more clear)
            stimulus_condition = conditions_to_test(i_condition, :);
            control_condition = conditions_control(strcmp(conditions_control(:, 1), stimulus_condition{1}), :);
            
            trigger_foot_indicator_stimulus = strcmp(loaded_data.conditions_session.condition_trigger_foot_list, stimulus_condition{1});
            index_indicator_stim = strcmp(loaded_data.conditions_session.condition_index_list, stimulus_condition{4}); % the step index (i.e. step one, step two, etc)
            indicator_this_stimulus = trigger_foot_indicator_stimulus & index_indicator_stim;

            trigger_onset_variable_data = loaded_data.analysis_data_session{strcmp(loaded_data.analysis_names_session, this_determine_onset_variable)}(:, indicator_this_stimulus);

            trigger_onset_variable_mean = mean(trigger_onset_variable_data, 2);
            trigger_onset_variable_cinv = cinv(trigger_onset_variable_data, 2);

            step_time_data = loaded_data.stretch_data_session{strcmp(loaded_data.stretch_names_session, {'step_time'})}(:, indicator_this_stimulus);
            step_time_mean = mean(step_time_data, 2);
            step_time = linspace(0,step_time_mean,100)';
                
            % find the inverted heel trajectories
            if strcmp(stimulus_condition{1}, 'TRIGGER_RIGHT')
                this_swing_foot_data = loaded_data.analysis_data_session{strcmp(loaded_data.analysis_names_session, {'lheel_x_inverted'})}(:, indicator_this_stimulus);
                this_swing_foot_mean = mean(this_swing_foot_data, 2);
                this_swing_foot_cinv = cinv(this_swing_foot_data, 2);
            elseif strcmp(stimulus_condition{1}, 'TRIGGER_LEFT')
                this_swing_foot_data = loaded_data.analysis_data_session{strcmp(loaded_data.analysis_names_session, {'rheel_x_inverted'})}(:, indicator_this_stimulus);
                this_swing_foot_mean = mean(this_swing_foot_data, 2);
                this_swing_foot_cinv = cinv(this_swing_foot_data, 2);                     
            end

            trigger_foot_response_onset = estimateOnsetTime(trigger_onset_variable_mean, trigger_onset_variable_cinv, step_time);
            swing_foot_response_onset = estimateOnsetTime(this_swing_foot_mean, this_swing_foot_cinv, step_time);
            
            % hack insert known issues with onset id
            if subject_id == 'EFU' & i_condition == 2
                trigger_foot_response_onset = 0; 
                swing_foot_response_onset = 0; 
            end
            if subject_id == 'GHJ' & i_condition == 1
                trigger_foot_response_onset = 0; % seems to be problem with readings this day
                swing_foot_response_onset = 1; % definetly foot place response
            end
            if subject_id == 'IDA' & i_condition == 1
                trigger_foot_response_onset = .195; % CI leaves zero (innapropriate onset calculation)
            end
            if subject_id == 'IDA' & i_condition == 2
                
                swing_foot_response_onset = 1;
            end
            % figure out what to do with NGY i_condition == 2 cop response
            if subject_id == 'ONT' & i_condition == 2
                trigger_foot_response_onset = .355; % CI leaves zero (innapropriate onset calculation)
            end
            if subject_id == 'RON' & i_condition == 2
                trigger_foot_response_onset = 0; % CI leaves zero (innapropriate onset calculation)
            end
            if subject_id == 'RRB' & i_condition == 1
                trigger_foot_response_onset = .330; % CI leaves zero (innapropriate onset calculation)
            end
            if subject_id == 'VQN' & i_condition == 1
                trigger_foot_response_onset = .415; % CI leaves zero (innapropriate onset calculation)
            end
            if subject_id == 'WHO' & i_condition == 1
                trigger_foot_response_onset = 0; % CI leaves zero (innapropriate onset calculation)
            end
            
            if trigger_foot_response_onset == 0 || trigger_foot_response_onset < .040 || trigger_foot_response_onset > .550
                if trigger_foot_response_onset > 0 && trigger_foot_response_onset < .040 || trigger_foot_response_onset > .550 % just checking what is going on
                    disp(['check subject' subject_id 'trigger foot' stimulus_condition{1}])
                    keyboard;
                end

                if any(swing_foot_response_onset) % may need to confirm it is in appropriate direction
                    this_system_cop_response_timing = 'Late';
                    this_system_response_type = 'Step';
                else
                    this_system_cop_response_timing = 'None';
                    this_system_response_type = 'None';
                end
            elseif trigger_foot_response_onset >= .350
                this_system_cop_response_timing = 'Late';

                if any(swing_foot_response_onset)
                    this_system_response_type = 'Combo';
                else
                    this_system_response_type = 'CoP';
                end

            elseif trigger_foot_response_onset < .350
                this_system_cop_response_timing = 'Early';

                if any(swing_foot_response_onset)
                    this_system_response_type = 'Combo';
                else
                    this_system_response_type = 'CoP';
                end
            end

            % save condition_group_list
            % assign groups for all indices (steps 1-4)
            indicator_this_stimulus = trigger_foot_indicator_stimulus;
            condition_response_time_list_session(indicator_this_stimulus) = {this_system_cop_response_timing};
            condition_response_type_list_session(indicator_this_stimulus) = {this_system_response_type};
         
            perturbation_indicator_control = strcmp(loaded_data.conditions_session.condition_perturbation_list, control_condition{2});
            condition_response_time_list_session(perturbation_indicator_control) = {'Control'};
        end
        
        % save
        variables_to_save = loaded_data;
        variables_to_save.conditions_session.condition_group_list = condition_response_time_list_session;
        variables_to_save.conditions_session.condition_response_type_list = condition_response_type_list_session;
        save(results_file_name, '-struct', 'variables_to_save');
        
        cd(study_path)
    end
end