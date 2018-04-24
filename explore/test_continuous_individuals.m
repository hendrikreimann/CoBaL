plot_cop_step = 0;
plot_icop_step = 0;
plot_i2cop_step = 0;
trigger_leg_ankle_dorsiflexion_max = 0;

% load data and settings
load('D:\DataStorage\Vision_HY\results.mat')

% subjects = {'DJB';'DXT';'EFU';'FNA';'GHJ';'IDA';'MTB';'NGY';'ONT';'PAG';'RON';'RRB';'SLL';'SPA';'UJD';'VQN';'WHO';'XDY';'YMU';'ZKY'}

for i_subject = 1:length(subjects)
%     subject_indicator = ismember(conditions.subject_list,subjects(i_subject));
      
    index_indicator = strcmp(conditions.condition_index_list, 'ONE');
    stance_foot_indicator_left = strcmp(conditions.condition_stance_foot_list, 'STANCE_LEFT');
    stance_foot_indicator_right = strcmp(conditions.condition_stance_foot_list, 'STANCE_RIGHT');
    index_indicator_left = stance_foot_indicator_left & index_indicator;
    index_indicator_right = stance_foot_indicator_right & index_indicator;
    
    group_indicator_early = strcmp(conditions.condition_group_list,'Early');
    group_indicator_late = strcmp(conditions.condition_group_list, 'Late');
    group_indicator_no = strcmp(conditions.condition_group_list, 'None');
    
    index_indicator_early = index_indicator & group_indicator_early;
    index_indicator_late = index_indicator & group_indicator_late;
    index_indicator_no = index_indicator & group_indicator_no;
    
    for i_variable = 1 : number_of_variables_to_test
        
        % TO DO: check cop_inverted
        % TO DO: check cop_integrated
        % TO DO: properly label the xticks on the bar graphs
        % TO DO: label the x and y axes
        % TO DO: transform heel trajectories into contra heel
        % TO DO: check the conditions with which plotting cop and others..
        % may be combining conditions that should not be combined
        % TO DO: plot shadedError without occlusion..
        
        if plot_cop_step
            plot(variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}(:,index_indicator), variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator),'*');
        end
        
        if plot_icop_step
            plot(variable_data{find(strcmp(variable_names, 'cop_from_com_x_integrated_inverted'))}(:,index_indicator), variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator),'*');
        end
        
        if plot_i2cop_step
            plot(variable_data{find(strcmp(variable_names, 'cop_from_com_x_integrated_twice_inverted'))}(:,index_indicator), variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator),'*');
        end
        
        if plot_pushoff_step
            plot(variable_data{find(strcmp(variable_names, 'left_ankle_dorsiflexion_step_end_inverted'))}(:,index_indicator_left), variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_left),'*');
            hold on;
            plot(variable_data{find(strcmp(variable_names, 'right_ankle_dorsiflexion_step_end_inverted'))}(:,index_indicator_right), variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_right),'*');
        end
        
        if plot_pushoff_step
            plot(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_step_end_inverted'))}(:,index_indicator), variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator),'*');
        end
        
        
        
        time_vector = linspace(0,100,100);
        bar_group_xaxis = ['Early','Late', 'None'];
        
        
        
        if plot_group_cop
            figure; hold on;
            mean_cop_early = mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_early),2);
            cinv_cop_early = cinv(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_early),2);
            mean_cop_late = mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_late),2);
            cinv_cop_late = cinv(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_late),2);
            mean_cop_no = mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_no),2);
            cinv_cop_no = cinv(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_no),2);
            
            shadedErrorBar(time_vector, mean_cop_early, cinv_cop_early, 'm');
            shadedErrorBar(time_vector, mean_cop_late, cinv_cop_late, 'c');
            shadedErrorBar(time_vector, mean_cop_no, cinv_cop_no, 'b');
            
        end
        
        if plot_group_step
            figure; hold on;
%             mean_step_early = mean(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_early),2);
%             cinv_step_early = cinv(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_early),2);
%             mean_step_late = mean(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_late),2);
%             cinv_step_late = cinv(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_late),2);
%             mean_step_no = mean(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_no),2);
%             cinv_step_no = cinv(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_no),2);
            
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
                mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}(:,index_indicator_late),2), ...
            mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}(:,index_indicator_no),2)];
            cop_cinv = [cinv(variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}(:,index_indicator_early),2), ...
                cinv(variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}(:,index_indicator_late),2), ...
            cinv(variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}(:,index_indicator_no),2)];
        
            bar([1:3],cop_means);   
            errorbar(cop_means, cop_cinv, 'LineStyle', 'none','color', 'k');
            set(gca,'xticklabel',bar_group_xaxis)
        end
        
        if bar_group_step
            figure; hold on;
            step_means = [mean(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_early),2), ...
                mean(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_late),2), ...
            mean(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_no),2)];
            step_cinv = [cinv(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_early),2), ...
                cinv(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_late),2), ...
            cinv(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_no),2)];
            bar([1:3],step_means);   
            errorbar(step_means, step_cinv, 'LineStyle', 'none','color', 'k');
            set(gca,'xticklabel',bar_group_xaxis)
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
        
        
    end
end






