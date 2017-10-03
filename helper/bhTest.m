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

% t-tests with FDR limitation after Benyamini and Hochberg (1995)

function [h, p_values] = bhTest(x, varargin)

    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'tail', 'both')
    addParameter(parser, 'fdr', 0.05)
    parse(parser, varargin{:})
    tail = parser.Results.tail;
    fdr = parser.Results.fdr;

    % perform t-tests
    if iscell(x)
        number_of_data_points = length(x);
        p_values = zeros(number_of_data_points, 1);
        for i_data = 1 : number_of_data_points
            [~, p_values(i_data)] = ttest(x{i_data}, 0, 'tail', tail);
        end
    else
        number_of_data_points = size(x, 1);
        p_values = zeros(number_of_data_points, 1);
        for i_data = 1 : number_of_data_points
            [~, p_values(i_data)] = ttest(x(i_data, :), 0, 'tail', tail);
        end
    end
    h_indices = 1 : number_of_data_points;
    hypothesis_matrix = [p_values, h_indices'];
    hypothesis_matrix_ascending = sortrows(hypothesis_matrix, 1);
    decision_threshold = (1:number_of_data_points)' * 1/number_of_data_points * fdr;
    hypothesis_rejection_index = find(hypothesis_matrix_ascending(:, 1) < decision_threshold, 1, 'last');
    if isempty(hypothesis_rejection_index)
        hypothesis_rejection_index = 0;
    end
    h_results_cop = [ones(1, hypothesis_rejection_index) zeros(1, number_of_data_points - hypothesis_rejection_index)];
    result_matrix_ascending = [hypothesis_matrix_ascending h_results_cop'];
    result_matrix = sortrows(result_matrix_ascending, 2);

    h = result_matrix(:, 3);
    
%     false_discovery_rate = 0.05;
%     p_ttest_cop_pos = zeros(number_of_time_steps_twrpd, 1);
%     p_ttest_cop_neg = zeros(number_of_time_steps_twrpd, 1);
%     for i_time = 1 : number_of_time_steps_twrpd
%         [~, p_ttest_cop_pos(i_time)] = ttest(cop_response_ml_pos(i_time, :), 0, 'tail', 'right');
%         [~, p_ttest_cop_neg(i_time)] = ttest(cop_response_ml_neg(i_time, :), 0, 'tail', 'left');
%     end
%     first_time_step_to_test = 1;
%     time_steps_to_test = first_time_step_to_test : number_of_time_steps_twrpd;
%     p_values_cop = [p_ttest_cop_pos(time_steps_to_test)' p_ttest_cop_neg(time_steps_to_test)'];
%     number_of_hypotheses_cop = length(p_values_cop);
%     h_indices_cop = 1 : number_of_hypotheses_cop;
%     hypothesis_matrix_cop = [p_values_cop', h_indices_cop'];
%     hypothesis_matrix_ascending_cop = sortrows(hypothesis_matrix_cop, 1);
%     decision_threshold_cop = (1:number_of_hypotheses_cop)' * 1/number_of_hypotheses_cop * false_discovery_rate;
%     hypothesis_rejection_index_cop = find(hypothesis_matrix_ascending_cop(:, 1) < decision_threshold_cop, 1, 'last');
%     h_results_cop = [ones(1, hypothesis_rejection_index_cop) zeros(1, number_of_hypotheses_cop-hypothesis_rejection_index_cop)];
%     result_matrix_ascending_cop = [hypothesis_matrix_ascending_cop h_results_cop'];
%     result_matrix_cop = sortrows(result_matrix_ascending_cop, 2);
%     time_point_cop_ttest_result_pos = [-1*ones(1, first_time_step_to_test-1) result_matrix_cop(1:number_of_hypotheses_cop/2, 3)']';
%     time_point_cop_ttest_result_neg = [-1*ones(1, first_time_step_to_test-1) result_matrix_cop(number_of_hypotheses_cop/2+1:end, 3)']';            
%     
%     % t-tests with FDR limitation after Benyamini and Hochberg (1995) - difference of magnitude
%     false_discovery_rate = 0.05;
%     p_ttest_cop_mag = zeros(number_of_time_steps_twrpd, 1);
%     for i_time = 1 : number_of_time_steps_twrpd
%         [~, p_ttest_cop_mag(i_time)] = ttest2(cop_response_ml_pos(i_time, :), -cop_response_ml_neg(i_time, :), 'tail', 'right');
%     end
%     first_time_step_to_test = 1;
%     time_steps_to_test = first_time_step_to_test : number_of_time_steps_twrpd;
%     p_values_cop_mag = p_ttest_cop_mag(time_steps_to_test)';
%     number_of_hypotheses_cop_mag = length(p_values_cop_mag);
%     h_indices_cop_mag = 1 : number_of_hypotheses_cop_mag;
%     hypothesis_matrix_cop_mag = [p_values_cop_mag', h_indices_cop_mag'];
%     hypothesis_matrix_ascending_cop_mag = sortrows(hypothesis_matrix_cop_mag, 1);
%     decision_threshold_cop_mag = (1:number_of_hypotheses_cop_mag)' * 1/number_of_hypotheses_cop_mag * false_discovery_rate;
%     hypothesis_rejection_index_cop_mag = find(hypothesis_matrix_ascending_cop_mag(:, 1) < decision_threshold_cop_mag, 1, 'last');
%     h_results_cop_mag = [ones(1, hypothesis_rejection_index_cop_mag) zeros(1, number_of_hypotheses_cop_mag-hypothesis_rejection_index_cop_mag)];
%     result_matrix_ascending_cop_mag = [hypothesis_matrix_ascending_cop_mag h_results_cop_mag'];
%     result_matrix_cop_mag = sortrows(result_matrix_ascending_cop_mag, 2);
    
    
end
