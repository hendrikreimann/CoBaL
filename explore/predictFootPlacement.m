function predictFootPlacement(varargin)
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;
    
     %% load and extract data
    study_settings_file = '';
    if exist('studySettings.txt', 'file')
        study_settings_file = 'studySettings.txt';
    end    
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    
    load('subjectInfo.mat', 'date', 'subject_id');
    load(['analysis' filesep date '_' subject_id '_results.mat']);
    
    conditions_control = study_settings.get('conditions_control');
    conditions_to_analyze = study_settings.get('conditions_to_analyze');
%     step_placement_x_data = variable_data_session{strcmp(variable_names_session, 'step_placement_x')};
    mpsis_x_data = variable_data_session{strcmp(variable_names_session, 'mpsis_x')};
    com_x_data = variable_data_session{strcmp(variable_names_session, 'com_x')};
    mpsis_x_vel_data = deriveByTime(mpsis_x_data, 0.01); % XXX ACHTUNG: this is a hack because I don't have the data yet
    com_x_vel_data = variable_data_session{strcmp(variable_names_session, 'com_x_vel')};
    lheelx_data = variable_data_session{strcmp(variable_names_session, 'lheel_x')};
    rheelx_data = variable_data_session{strcmp(variable_names_session, 'rheel_x')};
    lanklex_data = variable_data_session{strcmp(variable_names_session, 'lankle_x')};
    ranklex_data = variable_data_session{strcmp(variable_names_session, 'rankle_x')};
    cop_x_data = variable_data_session{strcmp(variable_names_session, 'cop_x')};
    midstance_index_data = ones(1, size(mpsis_x_data, 2)) * 65; % XXX ACHTUNG: this is a hack because I don't have the midstance index data yet
    midstance_index_data = variable_data_session{strcmp(variable_names_session, 'midstance_index')};
    
end