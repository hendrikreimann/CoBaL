% quick hack to get some p-values


%% load and extract older adult data
load '/Users/reimajbi/Temple Drive/Vision OA/results.mat'
stepx_response_OA = response_data_all{strcmp(variable_names{:, 1}, 'step_placement_x')};

stance_foot_indicator = strcmp(condition_stance_foot_list_all, 'STANCE_RIGHT');
perturbation_indicator = strcmp(condition_perturbation_list_all, 'ILLUSION_RIGHT');
index_indicator = strcmp(condition_index_list_all, 'ONE');
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;
stepx_response_stanceR_stimR_OA = stepx_response_OA(:, this_condition_indicator);

stance_foot_indicator = strcmp(condition_stance_foot_list_all, 'STANCE_LEFT');
perturbation_indicator = strcmp(condition_perturbation_list_all, 'ILLUSION_RIGHT');
index_indicator = strcmp(condition_index_list_all, 'ONE');
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;
stepx_response_stanceL_stimR_OA = stepx_response_OA(:, this_condition_indicator);

stance_foot_indicator = strcmp(condition_stance_foot_list_all, 'STANCE_RIGHT');
perturbation_indicator = strcmp(condition_perturbation_list_all, 'ILLUSION_LEFT');
index_indicator = strcmp(condition_index_list_all, 'ONE');
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;
stepx_response_stanceR_stimL_OA = stepx_response_OA(:, this_condition_indicator);

stance_foot_indicator = strcmp(condition_stance_foot_list_all, 'STANCE_LEFT');
perturbation_indicator = strcmp(condition_perturbation_list_all, 'ILLUSION_LEFT');
index_indicator = strcmp(condition_index_list_all, 'ONE');
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;
stepx_response_stanceL_stimL_OA = stepx_response_OA(:, this_condition_indicator);

% median_stepx_response_stanceR_stimR_OA = median(stepx_response_stanceR_stimR_OA)
% median_stepx_response_stanceL_stimR_OA = median(stepx_response_stanceL_stimR_OA)
% median_stepx_response_stanceR_stimL_OA = median(stepx_response_stanceR_stimL_OA)
% median_stepx_response_stanceL_stimL_OA = median(stepx_response_stanceL_stimL_OA)

%% load and extract healthy young data
load '/Users/reimajbi/Temple Drive/Vision HY/results.mat'
stepx_response_HY = response_data_all{strcmp(variable_names{:, 1}, 'step_placement_x')};

stance_foot_indicator = strcmp(condition_stance_foot_list_all, 'STANCE_RIGHT');
perturbation_indicator = strcmp(condition_perturbation_list_all, 'ILLUSION_RIGHT');
index_indicator = strcmp(condition_index_list_all, 'ONE');
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;
stepx_response_stanceR_stimR_HY = stepx_response_HY(:, this_condition_indicator);

stance_foot_indicator = strcmp(condition_stance_foot_list_all, 'STANCE_LEFT');
perturbation_indicator = strcmp(condition_perturbation_list_all, 'ILLUSION_RIGHT');
index_indicator = strcmp(condition_index_list_all, 'ONE');
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;
stepx_response_stanceL_stimR_HY = stepx_response_HY(:, this_condition_indicator);

stance_foot_indicator = strcmp(condition_stance_foot_list_all, 'STANCE_RIGHT');
perturbation_indicator = strcmp(condition_perturbation_list_all, 'ILLUSION_LEFT');
index_indicator = strcmp(condition_index_list_all, 'ONE');
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;
stepx_response_stanceR_stimL_HY = stepx_response_HY(:, this_condition_indicator);

stance_foot_indicator = strcmp(condition_stance_foot_list_all, 'STANCE_LEFT');
perturbation_indicator = strcmp(condition_perturbation_list_all, 'ILLUSION_LEFT');
index_indicator = strcmp(condition_index_list_all, 'ONE');
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;
stepx_response_stanceL_stimL_HY = stepx_response_HY(:, this_condition_indicator);

% median_stepx_response_stanceR_stimR_HY = median(stepx_response_stanceR_stimR_HY)
% median_stepx_response_stanceL_stimR_HY = median(stepx_response_stanceL_stimR_HY)
% median_stepx_response_stanceR_stimL_HY = median(stepx_response_stanceR_stimL_HY)
% median_stepx_response_stanceL_stimL_HY = median(stepx_response_stanceL_stimL_HY)

%% collect responses

% pool different feet
stepx_response_stimR_OA = [stepx_response_stanceR_stimR_OA stepx_response_stanceL_stimR_OA];
stepx_response_stimL_OA = [stepx_response_stanceR_stimL_OA stepx_response_stanceL_stimL_OA];
stepx_response_stimR_HY = [stepx_response_stanceR_stimR_HY stepx_response_stanceL_stimR_HY];
stepx_response_stimL_HY = [stepx_response_stanceR_stimL_HY stepx_response_stanceL_stimL_HY];

% figure; hold on
% histogram(stepx_response_stimR_OA);
% histogram(stepx_response_stimL_OA);
% 
% figure; hold on
% histogram(stepx_response_stimR_HY);
% histogram(stepx_response_stimL_HY);

% pool different stimulus sides
stepx_response_OA = [stepx_response_stimR_OA -stepx_response_stimL_OA];
stepx_response_HY = [stepx_response_stimR_HY -stepx_response_stimL_HY];

figure; hold on
histogram(stepx_response_OA);
histogram(stepx_response_HY);

%% run stats
mean_OA = mean(stepx_response_OA)
mean_HY = mean(stepx_response_HY)
std_OA = std(stepx_response_OA)
std_HY = std(stepx_response_HY)
[h, p] = ttest2(stepx_response_OA, stepx_response_HY)





return






