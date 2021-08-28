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

% this script transform raw data from ascii into matlab format

% input:
% Experimental data files generated by e.g. Nexus, QTM, labview etc.
% Files are expected to be in subfolders named <subject code>_<source type>, e.g. XYZ_labview
% any source type is viable, but if it is not in the default list, it must be specified as a name-value pair,
% e.g. importAscii('sources', 'someSource') would look for data in the subfolder XYZ_someSource
%
% output:
% Files containing the same data in .mat format, with some additional information about where they came from.
% Output files will be saved to folders "raw" and "processed".

function importProtocolData_VU(varargin)
    % parse arguments
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})

    % import protocol
    protocol_file_name = 'filenumbers.xlsx';
    if ~exist(protocol_file_name, 'file')
        error(['File ' protocol_file_name ' not found.']);
    else
        % we have a protocol file, so import that
        warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
        imported_table = readtable(protocol_file_name, 'Sheet', 'Trials');
        warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')
        
        % create containers
        trial_number = [];
        trial_type = {};
        speed = [];

        for i_line = 1 : size(imported_table, 1)
            % read line
            this_speed = imported_table.Speed(i_line);
            this_speed_first_attempt_trial_number = imported_table.firstAttempt(i_line);
            this_speed_second_attempt_trial_number = imported_table.secondAttempt(i_line);
            
            % store first attempt
            trial_number = [trial_number; this_speed_first_attempt_trial_number];
            trial_type = [trial_type; 'walking'];
            speed = [speed; this_speed];
            
            % store second attempt
            trial_number = [trial_number; this_speed_second_attempt_trial_number];
            trial_type = [trial_type; 'walking'];
            speed = [speed; this_speed];
        end

        save_file_name = 'protocolInfo.mat';
        save(save_file_name, 'trial_number', 'trial_type', 'speed');
        
        disp(['imported ' protocol_file_name])
    end

end


















  
