%% this script is currently choose as you go.. does not run as a functional script.. hacked together for DW and NIH grant

% load data and settings
loaded_data_OA1 = load('G:\My Drive\Vision_OA\results_grant1.mat');
loaded_data_OA2 = load('G:\My Drive\CoBaL_Test_Data\Vision_OA\results_grant1.mat');


%% invert OA and PD data
% OA1
index_indicator_step_one_OA1 =  strcmp(loaded_data_OA1.conditions.condition_index_list, 'ONE');
illusionleft_indicator_OA1 = strcmp(loaded_data_OA1.conditions.condition_perturbation_list,'ILLUSION_LEFT');
illusionright_indicator_OA1 = strcmp(loaded_data_OA1.conditions.condition_perturbation_list,'ILLUSION_RIGHT');
this_cop_data_OA1 = loaded_data_OA1.variable_data{find(strcmp(loaded_data_OA1.variable_names, 'cop_from_mpsis_x'))}(end,:)
this_cop_data_step_one_inverted_OA1 = [this_cop_data_OA1(:, illusionright_indicator_OA1 & index_indicator_step_one_OA1), this_cop_data_OA1(:, illusionleft_indicator_OA1 & index_indicator_step_one_OA1) * -1];

this_step_data_OA1 = loaded_data_OA1.variable_data{find(strcmp(loaded_data_OA1.variable_names, 'step_placement_x'))};
this_step_step_one_inverted_OA1 = [this_step_data_OA1(:, illusionright_indicator_OA1 & index_indicator_step_one_OA1), this_step_data_OA1(:, illusionleft_indicator_OA1 & index_indicator_step_one_OA1) * -1];

% OA2
% index_indicator_step_one_OA2 =  strcmp(loaded_data_OA2.conditions.condition_index_list, 'ONE');
illusionleft_indicator_OA2 = strcmp(loaded_data_OA2.conditions.stimulus_list,'STIM_LEFT');
illusionright_indicator_OA2 = strcmp(loaded_data_OA2.conditions.stimulus_list,'STIM_RIGHT');
this_cop_data_OA2 = loaded_data_OA2.variable_data{find(strcmp(loaded_data_OA2.variable_names, 'cop_from_mpsis_x'))}(100,:); % this index is for step one end
this_cop_data_step_one_inverted_OA2 = [this_cop_data_OA2(:, illusionright_indicator_OA2), this_cop_data_OA2(:, illusionleft_indicator_OA2) * -1];

this_step_data_OA2 = loaded_data_OA2.variable_data{find(strcmp(loaded_data_OA2.variable_names, 'step_placement_x'))}(1,:); % this index is for step one index
this_step_step_one_inverted_OA2 = [this_step_data_OA2(:, illusionright_indicator_OA2), this_step_data_OA2(:, illusionleft_indicator_OA2) * -1];

% concat OA1 and OA2
this_cop_data_OA = [this_cop_data_step_one_inverted_OA1 this_cop_data_step_one_inverted_OA2];
this_step_data_OA = [this_step_step_one_inverted_OA1 this_step_step_one_inverted_OA2];

% plot OA
[r , p] = corr(this_cop_data_OA', this_step_data_OA')
r2 = r^2

figure;
plot(this_cop_data_OA, this_step_data_OA,'o', 'color','k');


%% HY Data

plot_cop_step = 0;
plot_icop_step = 0;
plot_i2cop_step = 0;
trigger_leg_ankle_dorsiflexion_max = 0;
% TO DO: check for NaNs in Cop_step_end
% TO DO: integrated is less predictive than step end??
% TO DO: check cop_integrated
% TO DO: properly label the xticks on the bar graphs
% TO DO: label the x and y axes
% TO DO: transform heel trajectories into contra heel
% TO DO: check the conditions with which plotting cop and others..
% may be combining conditions that should not be combined
% TO DO: plot shadedError without occlusion..
% TO DO: add step trigger pushoff correlation with step length on step
% two
% TO DO; HAVE to uninvert triggerleg dorsiflexion for comparison to
% step length
load('D:\DataStorage\Vision_HY\results.mat')
index_indicator_step_one = strcmp(conditions.condition_index_list, 'ONE');
index_indicator_step_two = strcmp(conditions.condition_index_list, 'TWO');

stanceleft_indicator = strcmp(conditions.condition_stance_foot_list, 'STANCE_LEFT');
stanceright_indicator = strcmp(conditions.condition_stance_foot_list, 'STANCE_RIGHT');

illusionleft_indicator = strcmp(conditions.condition_perturbation_list,'ILLUSION_LEFT');
illusionright_indicator = strcmp(conditions.condition_perturbation_list,'ILLUSION_RIGHT');

index_stanceleft_one = stanceleft_indicator & index_indicator_step_one;
index_stanceright_one = stanceright_indicator & index_indicator_step_one;

group_indicator_early = strcmp(conditions.condition_group_list,'Early');
group_indicator_late = strcmp(conditions.condition_group_list, 'Late');
group_indicator_no = strcmp(conditions.condition_group_list, 'None');

index_indicator_early = index_indicator_step_one & group_indicator_early;
index_indicator_late = index_indicator_step_one & group_indicator_late;
index_indicator_no = index_indicator_step_one & group_indicator_no;

if plot_cop_step
    this_x_data = variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}(:,index_indicator_step_one);
    this_y_data = variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_step_one);

    this_x_data = variable_data{find(strcmp(variable_names, 'cop_from_com_x_integrated_inverted'))}(:,index_indicator_step_one);
    this_y_data = variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_step_one);

    this_x_data = variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_step_end_inverted'))}(:,index_indicator_step_one);
    this_y_data = variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_step_one);

    this_x_data = variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_step_end_inverted'))}(:,index_indicator_step_one);
    this_y_data = variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}(:,index_indicator_step_one);

    this_x_data = variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_inverted_max'))}(:,index_indicator_step_two);
    this_y_data = variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}(:,index_indicator_step_one);

    this_x_data = variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_step_end_inverted'))}(:,index_indicator_step_one);
    this_y_data = variable_data{find(strcmp(variable_names, 'step_length'))}(:,index_indicator_step_two);

    %         this_x_data = variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_step_end'))}(:,index_indicator_step_two);
    %         this_y_data = variable_data{find(strcmp(variable_names, 'step_length'))}(:,index_indicator_step_two);

    if any(isnan(this_x_data)) | any(isnan(this_y_data))
        indices_to_remove = find(isnan(this_x_data));
        indicies_to_remove = [indices_to_remove find(isnan(this_y_data))];
        indicies_to_remove = unique(indices_to_remove);
        this_x_data(indices_to_remove)=[]
        this_y_data(indices_to_remove)=[]
        [r , p] = corr(this_x_data', this_y_data')
        r2 = r^2
    end
    figure;
    plot(this_x_data, this_y_data,'o', 'color','k');
end

if plot_icop_step
    plot(variable_data{find(strcmp(variable_names, 'cop_from_com_x_integrated_inverted'))}(:,index_indicator_step_one), variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_step_one),'*');
end

if plot_i2cop_step
    plot(variable_data{find(strcmp(variable_names, 'cop_from_com_x_integrated_twice_inverted'))}(:,index_indicator_step_one), variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_step_one),'*');
end

if plot_pushoff_step
    plot(variable_data{find(strcmp(variable_names, 'left_ankle_dorsiflexion_step_end_inverted'))}(:,index_stanceleft_one), variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_stanceleft_one),'*');
    hold on;
    plot(variable_data{find(strcmp(variable_names, 'right_ankle_dorsiflexion_step_end_inverted'))}(:,index_stanceright_one), variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_stanceright_one),'*');
end

if plot_pushoff_step
    plot(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_step_end_inverted'))}(:,index_indicator_step_one), variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_step_one),'*');
end

time_vector = linspace(0,100,100);
bar_group_xaxis = ['Early','Late'];

if plot_group_cop
    figure; hold on;
    if strcmp(data_type, 'vision')
        mean_cop_early = mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_early),2);
        cinv_cop_early = cinv(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_early),2);
        mean_cop_late = mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_late),2);
        cinv_cop_late = cinv(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_late),2);
        mean_cop_no = mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_no),2);
        cinv_cop_no = cinv(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_no),2);

        shadedErrorBar(time_vector, mean_cop_early, cinv_cop_early, 'm');
        shadedErrorBar(time_vector, mean_cop_late, cinv_cop_late, 'c');
        shadedErrorBar(time_vector, mean_cop_no, cinv_cop_no, 'b');
    elseif strcmp(data_type, 'gvs')
        mean_cop_0 = mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_0),2);
        cinv_cop_0 = cinv(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_0),2);
        mean_cop_150 = mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_150),2);
        cinv_cop_150 = cinv(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_150),2);
        mean_cop_450 = mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_450),2);
        cinv_cop_450 = cinv(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_450),2);

        shadedErrorBar(time_vector, mean_cop_0, cinv_cop_0, 'm');
        shadedErrorBar(time_vector, mean_cop_150, cinv_cop_150, 'c');
        shadedErrorBar(time_vector, mean_cop_450, cinv_cop_450, 'b');
    end

end

if plot_group_step
    figure; hold on;
    mean_step_early = mean(variable_data{find(strcmp(variable_names, 'contra_leg_heel_x_inverted'))}(:,index_indicator_early),2);
    cinv_step_early = cinv(variable_data{find(strcmp(variable_names, 'contra_leg_heel_x_inverted'))}(:,index_indicator_early),2);
    mean_step_late = mean(variable_data{find(strcmp(variable_names, 'contra_leg_heel_x_inverted'))}(:,index_indicator_late),2);
    cinv_step_late = cinv(variable_data{find(strcmp(variable_names, 'contra_leg_heel_x_inverted'))}(:,index_indicator_late),2);
    mean_step_no = mean(variable_data{find(strcmp(variable_names, 'contra_leg_heel_x_inverted'))}(:,index_indicator_no),2);
    cinv_step_no = cinv(variable_data{find(strcmp(variable_names, 'contra_leg_heel_x_inverted'))}(:,index_indicator_no),2);

    shadedErrorBar(time_vector, mean_step_early, cinv_step_early, 'm');
    shadedErrorBar(time_vector, mean_step_late, cinv_step_late, 'c');
    shadedErrorBar(time_vector, mean_step_no, cinv_step_no, 'b');

end

if plot_group_pushoff
    figure; hold on;
    mean_pushoff_early = mean(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_inverted'))}(:,index_indicator_early),2);
    cinv_pushoff_early = cinv(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_inverted'))}(:,index_indicator_early),2);
    mean_pushoff_late = mean(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_inverted'))}(:,index_indicator_late),2);
    cinv_pushoff_late = cinv(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_inverted'))}(:,index_indicator_late),2);
    mean_pushoff_no = mean(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_inverted'))}(:,index_indicator_no),2);
    cinv_pushoff_no = cinv(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_inverted'))}(:,index_indicator_no),2);

    shadedErrorBar(time_vector, mean_pushoff_early, cinv_pushoff_early, 'm');
    shadedErrorBar(time_vector, mean_pushoff_late, cinv_pushoff_late, 'c');
    shadedErrorBar(time_vector, mean_pushoff_no, cinv_pushoff_no, 'b');
end

if plot_group_com
    figure; hold on;
    mean_com_early = mean(variable_data{find(strcmp(variable_names, 'com_x_inverted'))}(:,index_indicator_early),2);
    cinv_com_early = cinv(variable_data{find(strcmp(variable_names, 'com_x_inverted'))}(:,index_indicator_early),2);
    mean_com_late = mean(variable_data{find(strcmp(variable_names, 'com_x_inverted'))}(:,index_indicator_late),2);
    cinv_com_late = cinv(variable_data{find(strcmp(variable_names, 'com_x_inverted'))}(:,index_indicator_late),2);
    mean_com_no = mean(variable_data{find(strcmp(variable_names, 'com_x_inverted'))}(:,index_indicator_no),2);
    cinv_com_no = cinv(variable_data{find(strcmp(variable_names, 'com_x_inverted'))}(:,index_indicator_no),2);

    shadedErrorBar(time_vector, mean_com_early, cinv_com_early, 'm');
    shadedErrorBar(time_vector, mean_com_late, cinv_com_late, 'c');
    shadedErrorBar(time_vector, mean_com_no, cinv_com_no, 'b');
end

if plot_group_com_init
    figure; hold on;
    mean_com_init_early = mean(variable_data{find(strcmp(variable_names, 'com_from_com_init_x_inverted'))}(:,index_indicator_early),2);
    cinv_com_init_early = cinv(variable_data{find(strcmp(variable_names, 'com_from_com_init_x_inverted'))}(:,index_indicator_early),2);
    mean_com_init_late = mean(variable_data{find(strcmp(variable_names, 'com_from_com_init_x_inverted'))}(:,index_indicator_late),2);
    cinv_com_init_late = cinv(variable_data{find(strcmp(variable_names, 'com_from_com_init_x_inverted'))}(:,index_indicator_late),2);
    mean_com_init_no = mean(variable_data{find(strcmp(variable_names, 'com_from_com_init_x_inverted'))}(:,index_indicator_no),2);
    cinv_com_init_no = cinv(variable_data{find(strcmp(variable_names, 'com_from_com_init_x_inverted'))}(:,index_indicator_no),2);

    shadedErrorBar(time_vector, mean_com_init_early, cinv_com_init_early, 'm');
    shadedErrorBar(time_vector, mean_com_init_late, cinv_com_init_late, 'c');
    shadedErrorBar(time_vector, mean_com_init_no, cinv_com_init_no, 'b');
end

if plot_group_com_vel
    figure; hold on;
    mean_com_vel_early = mean(variable_data{find(strcmp(variable_names, 'com_x_vel_inverted'))}(:,index_indicator_early),2);
    cinv_com_vel_early = cinv(variable_data{find(strcmp(variable_names, 'com_x_vel_inverted'))}(:,index_indicator_early),2);
    mean_com_vel_late = mean(variable_data{find(strcmp(variable_names, 'com_x_vel_inverted'))}(:,index_indicator_late),2);
    cinv_com_vel_late = cinv(variable_data{find(strcmp(variable_names, 'com_x_vel_inverted'))}(:,index_indicator_late),2);
    mean_com_vel_no = mean(variable_data{find(strcmp(variable_names, 'com_x_vel_inverted'))}(:,index_indicator_no),2);
    cinv_com_vel_no = cinv(variable_data{find(strcmp(variable_names, 'com_x_vel_inverted'))}(:,index_indicator_no),2);

    shadedErrorBar(time_vector, mean_com_vel_early, cinv_com_vel_early, 'm');
    shadedErrorBar(time_vector, mean_com_vel_late, cinv_com_vel_late, 'c');
    shadedErrorBar(time_vector, mean_com_vel_no, cinv_com_vel_no, 'b');
end

if bar_group_cop
    figure; hold on;
    cop_means = [mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}(:,index_indicator_early),2), ...
        mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}(:,index_indicator_late),2)];
    cop_cinv = [cinv(variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}(:,index_indicator_early),2), ...
        cinv(variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}(:,index_indicator_late),2)];

    this_bar_early = bar([1],cop_means(1));
    this_bar_late = bar([2], cop_means(2));
    errorbar(cop_means, cop_cinv, 'LineStyle', 'none','color', 'k');
    
    this_bar_early.FaceColor = [0, 0, 1];
    this_bar_late.FaceColor = [0, 0.4470, 0.7410];
end

if bar_group_step
    figure; hold on;
    step_means = [mean(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_early),2), ...
        mean(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_late),2)];
    step_cinv = [cinv(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_early),2), ...
        cinv(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_late),2)];
    bar_step_early = bar([1],step_means(1));
    bar_step_late = bar([2], step_means(2));
    
    errorbar(step_means, step_cinv, 'LineStyle', 'none','color', 'k');
    bar_step_early.FaceColor = [0, 0, 1];
    bar_step_late.FaceColor = [0, 0.4470, 0.7410];
end

if bar_group_pushoff
    figure; hold on;
    pushoff_means = [mean(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_step_end_inverted'))}(:,index_indicator_early),2), ...
        mean(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_step_end_inverted'))}(:,index_indicator_late),2), ...
        mean(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_step_end_inverted'))}(:,index_indicator_no),2)];
    pushoff_cinv =  [cinv(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_step_end_inverted'))}(:,index_indicator_early),2), ...
        cinv(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_step_end_inverted'))}(:,index_indicator_late),2), ...
        cinv(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_step_end_inverted'))}(:,index_indicator_no),2)];
    bar(pushoff_means);
    errorbar(pushoff_means,pushoff_cinv, 'LineStyle', 'none','color', 'k');
    set(gca, 'XTick', [1, 2, 3])
    %             set(gca,'xticklabel',bar_group_xaxis)
end







