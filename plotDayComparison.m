
check_out_origin                    = 0;

do_cop_std_plots                    = 1;
do_step_parameter_plots             = 1;

save_figures                        = 1;

day_directories = ...
  { ...
    '/Users/reimajbi/TempleDrive/SubconcussionW/SHW002/Day1', ...
    '/Users/reimajbi/TempleDrive/SubconcussionW/SHW002/Day2', ...
  };
figure_directory = '/Users/reimajbi/TempleDrive/SubconcussionW/SHW002';

number_of_days = length(day_directories);

colors_day = parula(number_of_days);

%% load data
lcopx_std = [];
rcopx_std = [];
stepwidth = [];
steplength = [];
conditions_day = [];
conditions_stance = {};
for i_day = 1 : number_of_days
    load([day_directories{i_day} filesep 'subjectInfo.mat']);
    load([day_directories{i_day} filesep makeFileName(date, subject_id, 'gaitParametersConditions')]);
    load([day_directories{i_day} filesep makeFileName(date, subject_id, 'gaitParametersForceplate')]);
    load([day_directories{i_day} filesep makeFileName(date, subject_id, 'gaitParametersMarker')]);
    
    lcopx_std = [lcopx_std lcop_x_std_stanceL];
    rcopx_std = [rcopx_std rcop_x_std_stanceR];
    
    stepwidth = [stepwidth; step_width_total'];
    steplength = [steplength; step_length_total'];
    conditions_day = [conditions_day; ones(size(step_width_total')) * i_day];
    conditions_stance_day = cell(size(step_width_total'));
    conditions_stance_day(conditions_stanceL) = {'stance left'};
    conditions_stance_day(conditions_stanceR) = {'stance right'};
    conditions_stance = [conditions_stance; conditions_stance_day];
end
     



if do_cop_std_plots
    figure; axes; hold on; title('left foot medial-lateral CoP STD'); set(gca, 'Fontsize', 12)
    for i_day = 1 : number_of_days
        plot(time_normalized, lcopx_std(:, i_day), 'linewidth', 3, 'color', colors_day(i_day, :), 'displayname', ['Day ' num2str(i_day)]);
    end
    legend('show')
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, [figure_directory filesep 'cop_std_left.pdf'], 'pdf')
    end
    
    figure; axes; hold on; title('right foot medial-lateral CoP STD'); set(gca, 'Fontsize', 12)
    for i_day = 1 : number_of_days
        plot(time_normalized, rcopx_std(:, i_day), 'linewidth', 3, 'color', colors_day(i_day, :), 'displayname', ['Day ' num2str(i_day)]);
    end
    legend('show')
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, [figure_directory filesep 'cop_std_right.pdf'], 'pdf')
    end
end

%% step parameter plots
if do_step_parameter_plots
    step_width_figure = figure; axes('fontsize', 12); hold on; title('step width', 'fontsize', 16);
    boxplot(stepwidth, {conditions_day, conditions_stance});
    if save_figures
        saveas(gcf, [figure_directory filesep 'step_width.pdf'], 'pdf')
    end
    
    step_length_figure = figure; axes('fontsize', 12); hold on; title('step length', 'fontsize', 16);
    boxplot(steplength, {conditions_day, conditions_stance});
    if save_figures
        saveas(gcf, [figure_directory filesep 'step_length.pdf'], 'pdf')
    end
    
    
end


