% subjects = {'DJB';'DXT';'EFU';'FNA';'GHJ';'IDA';'MTB';'NGY';'ONT';'PAG';'RON';'RRB';'SLL';'SPA';'UJD';'VQN';'WHO';'XDY';'YMU';'ZKY'}
function processSubjectsArmsense(varargin)
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

%         deleteAvailableVariables
%         importAscii
%         saveSubjectInfoToFile
%         preprocessRawData
%         createModel
%  
        findStepEvents
        findRelevantDataStretches
        calculatePhaseTrajectories
        
%         calculateKinematicTrajectories('register_without_calculating', 1)
% 
%         findEmgNormalization
%         processStretchVariables
%         processAnalysisVariables
        cd(study_path)
%         processStimulusResponse
%         processGroupInformation
%         processAnalysisVariables
    end
end
