% load data and settings
function onset_classification(varargin)

    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'subjects', [])
    parse(parser, varargin{:})
    subjects = parser.Results.subjects;
    

    test_settings = SettingsCustodian('onset_classification.txt');
    % initialize
    conditions_to_test = test_settings.get('conditions_to_test');
    conditions_control = test_settings.get('conditions_control');
    number_of_conditions_to_test = size(conditions_to_test, 1);
    factor_to_analyze = test_settings.get('factor_to_analyze');
    variables_to_test = test_settings.get('variables_to_test');
    number_of_variables_to_test = size(variables_to_test, 1);

    data_folder_list = determineDataStructure(subjects);
    for i_folder = 1 : length(data_folder_list)
        % load data
        study_path = cd;
        data_path = data_folder_list{i_folder};
        cd(data_path)
        load([data_path filesep 'subjectInfo.mat'], 'date', 'subject_id');
        load([data_path filesep 'analysis' filesep date '_' subject_id '_results.mat']);
        disp(['Collecting from ' subject_id]);
        
        for i_condition = 1 : number_of_conditions_to_test
            % Separate (system) by trigger foot? consider change from stance to
            % trigger.. this will make it easier for the pushoff variable
            for i_variable = 1 : number_of_variables_to_test
                this_variable = variables_to_test{i_variable};
                stimulus_condition = conditions_to_test(i_condition, :);
                trigger_foot_indicator_stimulus = strcmp(conditions_session.condition_trigger_foot_list, stimulus_condition{1});
                index_indicator_stim = strcmp(conditions_session.condition_index_list, stimulus_condition{4}); % the step index (i.e. step one, step two, etc)
                indicator_this_stimulus = trigger_foot_indicator_stimulus & index_indicator_stim;

                trigger_foot_data = analysis_data_session{strcmp(analysis_names_session, this_variable)}(:, indicator_this_stimulus);

                trigger_foot_mean = -mean(trigger_foot_data, 2);
                trigger_foot_cinv = cinv(trigger_foot_data, 2);

                % determine average step time for this subject to place
                % into onsetCalculator
                step_time_data = stretch_data_session{strcmp(stretch_names_session, 'step_time')}(:, indicator_this_stimulus);
                step_time_mean = mean(step_time_data, 2);
                step_time = linspace(0,step_time_mean,100);
                % estimate onset time (find the time variable for the loaded variable)
                trigger_foot_response_onset = estimateOnsetTime(trigger_foot_mean, trigger_foot_cinv, step_time);
               
                if any(trigger_foot_response_onset)
                    % check for odd onsets
                    if trigger_foot_response_onset < 100 || trigger_foot_response_onset > 500

                    end
                    % if Early or Late... determine if CoP, Foot Place, or Combo
                    % run onsetEstimate for
                    if trigger_foot_response_onset >= 350

                    else trigger_foot_response_onset < 350
                        % if Early or Late... determine if CoP, Foot Place, or Combo
                        % run onsetEstimate for 
                    end
                end
                
            % classify
            % take heel positon for respective system and determine whether
            % Early, Late, No AND whether the Early/Late was CoP, Foot Place,
            % or Combo



        % assign condition_group_list

        % save condition_group_list

            end

        end
        cd(study_path)
    end
end