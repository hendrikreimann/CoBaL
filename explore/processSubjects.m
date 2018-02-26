% subjects = {'DJB';'DXT';'EFU';'FNA';'GHJ';'IDA';'MTB';'NGY';'ONT';'PAG';'RON';'RRB';'SLL';'SPA';'UJD';'VQN';'WHO';'XDY';'YMU';'ZKY'}
function processSubjects(varargin)
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'subjects', [])
    parse(parser, varargin{:})
    subjects = parser.Results.subjects;
    study_settings = SettingsCustodian('studySettings.txt');
    data_folder_list = determineDataStructure(subjects);
    for i_folder = 1 : length(data_folder_list)
        % load data
        study_path = cd;
        data_path = data_folder_list{i_folder};
        cd(data_path)
        load([data_path filesep 'subjectInfo.mat'], 'date', 'subject_id');
        disp(['Collecting from ' subject_id]);
%         load([data_path filesep 'analysis' filesep date '_' subject_id '_results.mat']);

        deleteAvailableVariables
        importAscii
        saveSubjectInfoToFile
        preprocessRawData
        createModel

        findRelevantDataStretches
        calculateKinematicTrajectories('register_without_calculating', 1)

        findEmgNormalization
        processStretchVariables

        cd(study_path)
%         processStimulusResponse
%         processGroupInformation
%         processAnalysisVariables
    end
end
