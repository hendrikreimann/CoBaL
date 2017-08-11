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

% explore the spm1d toolbox

% load('/Users/reimajbi/Drive_UD/Vision_HY/DJB/analysis/20170227_DJB_results.mat')
load('/Users/reimajbi/Drive_UD/Vision_HY/results.mat')

% define and extract data
condition_to_analyze_stance = 'STANCE_RIGHT';
condition_to_analyze_perturbation = 'ILLUSION_RIGHT';
condition_to_analyze_index = 'ONE';
variable_to_analyze = 'cop_from_com_x';
% variable_to_analyze = 'right_ankle_eversion_angle';
% variable_to_analyze = 'right_ankle_dorsiflexion_angle';

stance_foot_indicator = strcmp(condition_stance_foot_list_all, condition_to_analyze_stance);
perturbation_indicator = strcmp(condition_perturbation_list_all, condition_to_analyze_perturbation);
index_indicator = strcmp(condition_index_list_all, condition_to_analyze_index);
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;


data_to_analyze = response_data_all{strcmp(variable_names, variable_to_analyze)};

data_to_analyze_this_condition = data_to_analyze(:, this_condition_indicator);

% figure; hold on;
% plot(data_to_analyze_this_condition);
% plot(mean(data_to_analyze_this_condition, 2), 'linewidth', 5);


% SPM1D analysis and results
spm = spm1d.stats.ttest(data_to_analyze_this_condition');
spmi = spm.inference(0.05, 'two_tailed', false);
figure;
spmi.plot();
spmi.plot_threshold_label();
spmi.plot_p_values();

% FDR-controlled t-test analysis and results
time = 1 : size(data_to_analyze_this_condition, 1);
[h, p] = bhTest(data_to_analyze_this_condition, 'tail', 'both', 'fdr', 0.5);
figure; hold on;
plot_handles = shadedErrorBar ...
  ( ...
    1 : 100, ...
    mean(data_to_analyze_this_condition, 2), ...
    cinv(data_to_analyze_this_condition, 2), ...
    { ...
      'color', [1 0 0], ...
      'linewidth', 6 ...
    }, ...
    1 ...
);
plot(time(h==1), 0, 'x', 'color', [1 0 0], 'markersize', 8, 'linewidth', 2);
















