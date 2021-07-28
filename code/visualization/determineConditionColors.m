function condition_colors = determineConditionColors(colors_table, comparisons, condition_to_compare)
    % find unique levels of condition to compare
    levels = unique(comparisons.condition_combinations(:, strcmp(comparisons.condition_combination_labels, condition_to_compare)));
    
    % make default color map
    default_colors = lines(length(levels));
    
    if isempty(colors_table)
        colors_table_from_settings = table('Size', [0, 3], 'VariableNames', {'condition', 'level', 'color'}, 'VariableTypes', {'string', 'string', 'string'});
    else
        colors_table_from_settings = colors_table(strcmp(colors_table.condition, condition_to_compare), :);
    end
    
    % go through levels and store default color or the one provided in the settings
    condition_colors = [levels cell(size(levels))];
    for i_level = 1 : length(levels)
        this_level = levels(i_level);
        if any(strcmp(colors_table_from_settings.level, this_level))
            % use color provided in settings
            condition_colors{i_level, 2} = hex2rgb(colors_table_from_settings.color{strcmp(colors_table_from_settings.level, this_level)});
        else
            % use default color
            condition_colors{i_level, 2} = default_colors(i_level, :);
        end
        
    end
end