plot_cop_step = 0;
plot_icop_step = 0;
plot_i2cop_step = 0;
trigger_leg_ankle_dorsiflexion_max = 0;

% load data and settings
load('D:\DataStorage\Vision_HY\results.mat')

subjects = {'DJB';'DXT';'EFU';'FNA';'GHJ';'IDA';'MTB';'NGY';'ONT';'PAG';'RON';'RRB';'SLL';'SPA';'UJD';'VQN';'WHO';'XDY';'YMU';'ZKY'}

for i_subject = 1:length(subjects)
subject_indicator = ismember(conditions.subject_list,subjects(i_subject));
    
    % get condition descriptors
    stimulus_condition = conditions_to_test(i_condition, :);
    control_condition = conditions_control(strcmp(conditions_control(:, 1), stimulus_condition{1}), :);
      
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
    
   
    % get data and perform t-tests
    for i_variable = 1 : number_of_variables_to_test
        
        % TO DO: check cop_inverted
        
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
        
        
        % TO DO: add cinvs
        
        if plot_group_cop
            figure; hold on;
            plot(mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_early),2));
            plot(mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_late),2));
            plot(mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_inverted'))}(:,index_indicator_no),2));
        end
        
        if plot_group_step
            figure; hold on;
            plot(mean(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_early),2));
            plot(mean(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_late),2));
            plot(mean(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_no),2));
        end
        
        if plot_group_pushoff
            plot(mean(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_inverted'))}(:,index_indicator_early),2));
            plot(mean(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_inverted'))}(:,index_indicator_late),2));
            plot(mean(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_inverted'))}(:,index_indicator_no),2));
        end
        
        if plot_group_com
            plot(mean(variable_data{find(strcmp(variable_names, 'com_x_inverted'))}(:,index_indicator_early),2));
            plot(mean(variable_data{find(strcmp(variable_names, 'com_x_inverted'))}(:,index_indicator_late),2));
            plot(mean(variable_data{find(strcmp(variable_names, 'com_x_inverted'))}(:,index_indicator_no),2));
        end

        if plot_group_com_init
            plot(mean(variable_data{find(strcmp(variable_names, 'com_from_com_init_x_inverted'))}(:,index_indicator_early),2));
            plot(mean(variable_data{find(strcmp(variable_names, 'com_from_com_init_x_inverted'))}(:,index_indicator_late),2));
            plot(mean(variable_data{find(strcmp(variable_names, 'com_from_com_init_x_inverted'))}(:,index_indicator_no),2));
        end
        
        if plot_group_com_vel
            plot(mean(variable_data{find(strcmp(variable_names, 'com_x_vel_inverted'))}(:,index_indicator_early),2));
            plot(mean(variable_data{find(strcmp(variable_names, 'com_x_vel_inverted'))}(:,index_indicator_late),2));
            plot(mean(variable_data{find(strcmp(variable_names, 'com_x_vel_inverted'))}(:,index_indicator_no),2));
        end
        
        if bar_group_cop
            figure; hold on;
            cop_means = [mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}(:,index_indicator_early),2), ...
                mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}(:,index_indicator_late),2), ...
            mean(variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}(:,index_indicator_no),2)];
            bar(cop_means);   
        end
        
        if bar_group_step
            figure; hold on;
            step_means = [mean(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_early),2), ...
                mean(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_late),2), ...
            mean(variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,index_indicator_no),2)];
            bar(step_means);   
        end
        
        if bar_group_pushoff
            figure; hold on;
            pushoff_means = [mean(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_step_end_inverted'))}(:,index_indicator_early),2), ...
                mean(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_step_end_inverted'))}(:,index_indicator_late),2), ...
            mean(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_step_end_inverted'))}(:,index_indicator_no),2)];
            bar(pushoff_means);   
        end

    end
end






