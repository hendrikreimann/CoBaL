%     This file is part of the CoBaL code base
%     Copyright (C) 2017-2018 Hendrik Reimann <hendrikreimann@gmail.com>
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

classdef ObstacleDataProcessingTest < matlab.unittest.TestCase
    % Unit Test to process obstacle data
    %   to test, instantiate and call run
    
    properties
        tmp;
        data_root;
        subjects = {'A', 'B', 'C'};
    end

    methods
        function this = ObstacleDataProcessingTest
            % clear up temporary data
            this.tmp = [tempdir 'Obstacle'];
            if exist(this.tmp, 'dir')
                rmdir(this.tmp, 's');
            end
            mkdir(this.tmp);
            
            % locate unit test data
            this.data_root = this.locateTestData();
            this.copyDataToTmp();
        end
        function data_root = locateTestData(this)
            % locate test data
            user_dir = getUserDir();
            default_locations = ...
              { ...
                [user_dir filesep 'Dropbox' filesep 'unitTestData' filesep 'Obstacle'], ...
              };
            
            location_found = false;
            for i_location = 1 : length(default_locations)
                if exist(default_locations{i_location}, 'dir')
                    data_root = default_locations{i_location};
                    location_found = true;
                    disp(['Obstacle data found in "' data_root '"'])
                end
            end
            if ~location_found
                % TODO ask user for location of unit test data
            end
        end
        function copyDataToTmp(this)
            copyfile([this.data_root filesep 'eventGuiSettings.txt'], this.tmp)
            copyfile([this.data_root filesep 'plotSettings.txt'], this.tmp)
            copyfile([this.data_root filesep 'studySettings.txt'], this.tmp)
            copyfile([this.data_root filesep 'subjects.csv'], this.tmp)
            for i_subject = 1 : length(this.subjects)
                subject_dir = [this.tmp filesep this.subjects{i_subject} filesep];
                mkdir(subject_dir);
                copyfile([this.data_root filesep this.subjects{i_subject} filesep 'conditions.csv'], subject_dir);
                copyfile([this.data_root filesep this.subjects{i_subject} filesep 'subjectSettings.txt'], subject_dir);
                copyfile([this.data_root filesep this.subjects{i_subject} filesep 'ascii' filesep '*'], [subject_dir filesep 'ascii']);
                copyfile([this.data_root filesep this.subjects{i_subject} filesep 'neurocom' filesep '*'], [subject_dir filesep 'neurocom']);
            end
        end
    end
    methods (Test)
        function testImportAndPreprocess(this)
            current_dir = pwd;
            for i_subject = 1 : length(this.subjects)
                cd([this.tmp filesep this.subjects{i_subject}]);
                importAscii;
                saveSubjectInfoToFile;
                preprocessRawData;
            end
            cd(current_dir);
        end
        function testModelAndKinematics(this)
            current_dir = pwd;
            test_subject = 'A';
            cd([this.tmp filesep test_subject]);
            createModel;
            cd(current_dir);
        end
    end
    
end

