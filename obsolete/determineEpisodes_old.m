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

function episode_indices = determineEpisodes(study_settings, plot_settings, comparison_indices)
    conditions_to_plot = plot_settings.get('conditions_to_plot');
    condition_column_index = find(strcmp(study_settings.get('condition_labels'), 'index'));
    condition_column_stancefoot = find(strcmp(study_settings.get('condition_labels'), 'stance foot'));
%     episode_first_stretch_indices = find(strcmp(study_settings.get('conditions_to_plot(:, 4), 'ONE'));
    episode_indices = {};
    comparisons_already_used = [];
    number_of_comparisons = length(comparison_indices);
    while length(comparisons_already_used) < number_of_comparisons
        % start with the first available comparison
        i_comparison = 1;
        while ismember(i_comparison, comparisons_already_used)
            i_comparison = i_comparison + 1;
        end
        comparison_indices_in_this_episode = i_comparison; % this is the first comparison in this episode, more will be added

        % search for comparisons that differ from this one in only the step number
        base_comparison = comparison_indices{i_comparison};
        example_condition_in_base_comparison = base_comparison(1);
        example_condition_in_base_comparison_labels = conditions_to_plot(example_condition_in_base_comparison, :);
        for j_comparison = 1 : number_of_comparisons
            if i_comparison ~= j_comparison
                this_comparison = comparison_indices{j_comparison};
                example_condition_in_this_comparison = this_comparison(1);
                example_condition_in_this_comparison_labels = conditions_to_plot(example_condition_in_this_comparison, :);
                % check which conditions labels agree between these two conditions
                comparison_table = zeros(1, length(study_settings.get('condition_labels'))); % this is a table indicating equality between the two conditions in questions
                for i_label = 1 : length(study_settings.get('condition_labels'))
                    comparison_table(i_label) = strcmp(example_condition_in_base_comparison_labels{i_label}, example_condition_in_this_comparison_labels{i_label});
                end

                % look at the relevant entries of the comparison table
                comparison_table_relevant = comparison_table;
                comparison_table_relevant([condition_column_stancefoot condition_column_index plot_settings.get('comparison_to_make')]) = [];
                if all(comparison_table_relevant)
                    % check if the stance foot is alternating
                    if strcmp(example_condition_in_base_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                        if strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'TWO') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison]; %#ok<AGROW>
                        elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'THREE') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison]; %#ok<AGROW>
                        elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'FOUR') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison]; %#ok<AGROW>
                        end
                    end
                    if strcmp(example_condition_in_base_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                        if strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'TWO') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison]; %#ok<AGROW>
                        elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'THREE') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison]; %#ok<AGROW>
                        elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'FOUR') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison]; %#ok<AGROW>
                        end
                    end
                    if strcmp(example_condition_in_base_comparison_labels(condition_column_stancefoot), 'STANCE_BOTH')
                        if strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'TWO') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
        episode_indices = [episode_indices; comparison_indices_in_this_episode]; %#ok<AGROW>
        comparisons_already_used = [comparisons_already_used comparison_indices_in_this_episode]; %#ok<AGROW>

    end   
end
