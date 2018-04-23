plot_cop_step = 0;
plot_icop_step = 0;
plot_i2cop_step = 0;
trigger_leg_ankle_dorsiflexion_max = 0;

% load data and settings
load('D:\DataStorage\Vision_HY\results.mat')

subjects = {'DJB';'DXT';'EFU';'FNA';'GHJ';'IDA';'MTB';'NGY';'ONT';'PAG';'RON';'RRB';'SLL';'SPA';'UJD';'VQN';'WHO';'XDY';'YMU';'ZKY'}

for i_subject = 1:length(subjects)
subject_indicator = ismember(conditions.subject_list,subjects(i_subject));

% What does this script do?
    % meant to grab data depending on what you want to look at
    % 1) look at group data
    % 2) look at all mechanism data
    % 3) identify index 2
    % 4) index 3?
    % 5) all data should be inverted, so not necessary to decipher between
    % right and left, just plot all
    
    % get condition descriptors
    stimulus_condition = conditions_to_test(i_condition, :);
    control_condition = conditions_control(strcmp(conditions_control(:, 1), stimulus_condition{1}), :);
    
    % get condition indicators
    
    stance_foot_indicator_stimulus = strcmp(condition_stance_foot_list, stimulus_condition{1});
    stance_foot_indicator_control = strcmp(condition_stance_foot_list, control_condition{1});
    perturbation_indicator_stim_Right = strcmp(condition_perturbation_list, 'ILLUSION_RIGHT');
    perturbation_indicator_stim_Left = strcmp(condition_perturbation_list, 'ILLUSION_LEFT');
    perturbation_indicator_control = strcmp(condition_perturbation_list, control_condition{2});
    index_indicator_stim = strcmp(condition_index_list, stimulus_condition{4});
    index_indicator_control = strcmp(condition_index_list, control_condition{4});
    indicator_stimulus_right = stance_foot_indicator_stimulus & perturbation_indicator_stim_Right & index_indicator_stim & subject_indicator;
    indicator_stimulus_left = stance_foot_indicator_stimulus & perturbation_indicator_stim_Left & index_indicator_stim & subject_indicator;
    indicator_control = stance_foot_indicator_control & perturbation_indicator_control & index_indicator_control & subject_indicator;
    
    % get data and perform t-tests
    for i_variable = 1 : number_of_variables_to_test
        this_variable = variables_to_test{i_variable};
        this_condition = conditions_to_test(i_condition, :);
        
        if plot_cop_step
            plot(variable_data{find(strcmp(variable_names, 'cop_from_com_x_step_end_inverted'))}, variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))},'*');
        end
        
        if plot_icop_step
            plot(variable_data{find(strcmp(variable_names, 'cop_from_com_x_integrated_inverted'))}, variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))},'*');
        end
        
        if plot_i2cop_step
            plot(variable_data{find(strcmp(variable_names, 'cop_from_com_x_integrated_twice_inverted'))}, variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))},'*');
        end
        
%         if plot_pushoff_step
%             plot(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_max'))}, variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))},'*');
%         end
        if plot_pushoff_step
            plot(variable_data{find(strcmp(variable_names, 'left_ankle_dorsiflexion_step_end_inverted'))}(:,find(strcmp(conditions.condition_stance_foot_list, 'STANCE_LEFT'))), variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,find(strcmp(conditions.condition_stance_foot_list, 'STANCE_LEFT'))),'*');
            hold on;
            plot(variable_data{find(strcmp(variable_names, 'right_ankle_dorsiflexion_step_end_inverted'))}(:,find(strcmp(conditions.condition_stance_foot_list, 'STANCE_RIGHT'))), variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))}(:,find(strcmp(conditions.condition_stance_foot_list, 'STANCE_RIGHT'))),'*');
        end
        
         if plot_pushoff_step
            plot(variable_data{find(strcmp(variable_names, 'trigger_leg_ankle_dorsiflexion_step_end_inverted'))}, variable_data{find(strcmp(variable_names, 'step_placement_x_inverted'))},'*');
        end
        
        
    end
end






