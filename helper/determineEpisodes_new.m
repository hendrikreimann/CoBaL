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

% One episode will correspond to one figure.
% An episode is a 1-d array of comparison indices

% currently true, but about to be changed: (instead of having 2-d episodes with four columns for the indices, where each
%   column has a condition, I want to have 1-d episodes, with four entries for the indices, where each entry is a 
%   comparison index.
% One episode is a 2d-array of condition indices, where each index refers
%   to a row of condition_combinations_stimulus. The dimensions of the
%   episode array correspond to
%   - dim 1: ranges over the different levels of the comparison we're making here
%   - dim 2: ranges over the step indices (ONE to FOUR)
% episodes is a column cell array, each entry corresponds to one episode

function episodes_comparison_map = determineEpisodes_new(condition_combinations_stimulus, condition_combination_labels, comparison_indices, plot_settings)
    condition_to_compare = plot_settings.get('condition_to_compare');
    index_column = strcmp(condition_combination_labels, 'index');
    comparison_column = strcmp(condition_combination_labels, condition_to_compare);
    relevant_column_indices = find(~(index_column | comparison_column));
    index_labels = {'ONE', 'TWO', 'THREE', 'FOUR'};
    number_of_indices = length(index_labels);
    
    % use step one to figure out number of episodes and declare containers
    step_one_condition_indices = strcmp(condition_combinations_stimulus(:, index_column), 'ONE');
    step_one_conditions = condition_combinations_stimulus(step_one_condition_indices, :);
    % each unique combination of conditions is an episode, ignoring index and comparison column
    step_one_conditions_relevant = step_one_conditions;
    [step_one_conditions_relevant{:, index_column}] = deal('~');
    [step_one_conditions_relevant{:, comparison_column}] = deal('~');
    [~, unique_row_indices] = unique(cell2table(step_one_conditions_relevant), 'rows');
    unique_rows = step_one_conditions(unique_row_indices, :);
    number_of_episodes = length(unique_row_indices);
    comparison_levels = unique(condition_combinations_stimulus(:, comparison_column));
    number_of_comparison_levels = length(comparison_levels);
    episodes_condition_map = cell(number_of_episodes, 1);
    
    % populate episodes_condition_map
    for i_episode = 1 : number_of_episodes
        episodes_condition_map{i_episode} = zeros(number_of_comparison_levels, number_of_indices);
        % defining condition
        defining_condition = unique_rows(i_episode, :);
        for i_comparison = 1 : number_of_comparison_levels
            for i_index = 1 : number_of_indices
                % for each step and comparison level, find the row in the condition_combinations_stimulus that belongs to the defining condition
                index_hit = strcmp(condition_combinations_stimulus(:, index_column), index_labels{i_index});
                comparison_hit = strcmp(condition_combinations_stimulus(:, comparison_column), comparison_levels{i_comparison});
                defining_condition_hit = true(size(index_hit));
                for i_label = relevant_column_indices
                    this_label_hit = strcmp(condition_combinations_stimulus(:, i_label), defining_condition{i_label});
                    defining_condition_hit = defining_condition_hit & this_label_hit;
                end
                this_condition_index = find(index_hit & comparison_hit & defining_condition_hit);
                episodes_condition_map{i_episode}(i_comparison, i_index) = this_condition_index;
%                 this_condition = condition_combinations_stimulus(this_condition_index, :);
            end
        end
    end
    
    % transform to episodes_comparison_map, because the current version of plotResults requests that
    episodes_comparison_map = cell(number_of_episodes, 1);
    for i_episode = 1 : number_of_episodes
        this_episode_condition_map = episodes_condition_map{i_episode};
        this_episode_comparison_map = zeros(1, number_of_indices);
        for i_index = 1 : number_of_indices
            conditions_this_step = sort(this_episode_condition_map(:, i_index))';
            % go through comparisons and compare the conditions
            for i_comparison = 1 : length(comparison_indices)
                conditions_this_comparison = sort(comparison_indices{i_comparison});
                if isequal(conditions_this_step, conditions_this_comparison)
                    this_episode_comparison_map(i_index) = i_comparison;
                end
            end
        end
        episodes_comparison_map{i_episode} = this_episode_comparison_map;
        
    end
    
    
    
end
















