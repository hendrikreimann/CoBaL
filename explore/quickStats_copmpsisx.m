% quick hack to get some p-values


%% load and extract older adult data
load '/Users/reimajbi/Temple Drive/Vision OA/results.mat'
copmpsisx_response_OA = response_data_all{strcmp(variable_names{:, 1}, 'cop_from_mpsis_x')};

stance_foot_indicator = strcmp(condition_stance_foot_list_all, 'STANCE_RIGHT');
perturbation_indicator = strcmp(condition_perturbation_list_all, 'ILLUSION_RIGHT');
index_indicator = strcmp(condition_index_list_all, 'ONE');
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;
copmpsisx_response_stanceR_stimR_OA = copmpsisx_response_OA(:, this_condition_indicator);

stance_foot_indicator = strcmp(condition_stance_foot_list_all, 'STANCE_LEFT');
perturbation_indicator = strcmp(condition_perturbation_list_all, 'ILLUSION_RIGHT');
index_indicator = strcmp(condition_index_list_all, 'ONE');
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;
copmpsisx_response_stanceL_stimR_OA = copmpsisx_response_OA(:, this_condition_indicator);

stance_foot_indicator = strcmp(condition_stance_foot_list_all, 'STANCE_RIGHT');
perturbation_indicator = strcmp(condition_perturbation_list_all, 'ILLUSION_LEFT');
index_indicator = strcmp(condition_index_list_all, 'ONE');
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;
copmpsisx_response_stanceR_stimL_OA = copmpsisx_response_OA(:, this_condition_indicator);

stance_foot_indicator = strcmp(condition_stance_foot_list_all, 'STANCE_LEFT');
perturbation_indicator = strcmp(condition_perturbation_list_all, 'ILLUSION_LEFT');
index_indicator = strcmp(condition_index_list_all, 'ONE');
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;
copmpsisx_response_stanceL_stimL_OA = copmpsisx_response_OA(:, this_condition_indicator);

% figure; axes; hold on
% plot(mean(copmpsisx_response_stanceR_stimR_OA, 2), 'linewidth', 3)
% plot(mean(copmpsisx_response_stanceL_stimR_OA, 2), 'linewidth', 3)
% plot(mean(copmpsisx_response_stanceR_stimL_OA, 2), 'linewidth', 3)
% plot(mean(copmpsisx_response_stanceL_stimL_OA, 2), 'linewidth', 3)

%% load and extract healthy young data
load '/Users/reimajbi/Temple Drive/Vision HY/results.mat'
copmpsisx_response_HY = response_data_all{strcmp(variable_names{:, 1}, 'cop_from_mpsis_x')};

stance_foot_indicator = strcmp(condition_stance_foot_list_all, 'STANCE_RIGHT');
perturbation_indicator = strcmp(condition_perturbation_list_all, 'ILLUSION_RIGHT');
index_indicator = strcmp(condition_index_list_all, 'ONE');
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;
copmpsisx_response_stanceR_stimR_HY = copmpsisx_response_HY(:, this_condition_indicator);

stance_foot_indicator = strcmp(condition_stance_foot_list_all, 'STANCE_LEFT');
perturbation_indicator = strcmp(condition_perturbation_list_all, 'ILLUSION_RIGHT');
index_indicator = strcmp(condition_index_list_all, 'ONE');
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;
copmpsisx_response_stanceL_stimR_HY = copmpsisx_response_HY(:, this_condition_indicator);

stance_foot_indicator = strcmp(condition_stance_foot_list_all, 'STANCE_RIGHT');
perturbation_indicator = strcmp(condition_perturbation_list_all, 'ILLUSION_LEFT');
index_indicator = strcmp(condition_index_list_all, 'ONE');
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;
copmpsisx_response_stanceR_stimL_HY = copmpsisx_response_HY(:, this_condition_indicator);

stance_foot_indicator = strcmp(condition_stance_foot_list_all, 'STANCE_LEFT');
perturbation_indicator = strcmp(condition_perturbation_list_all, 'ILLUSION_LEFT');
index_indicator = strcmp(condition_index_list_all, 'ONE');
this_condition_indicator = stance_foot_indicator & perturbation_indicator & index_indicator;
copmpsisx_response_stanceL_stimL_HY = copmpsisx_response_HY(:, this_condition_indicator);

% figure; axes; hold on
% plot(mean(copmpsisx_response_stanceR_stimR_HY, 2), 'linewidth', 3)
% plot(mean(copmpsisx_response_stanceL_stimR_HY, 2), 'linewidth', 3)
% plot(mean(copmpsisx_response_stanceR_stimL_HY, 2), 'linewidth', 3)
% plot(mean(copmpsisx_response_stanceL_stimL_HY, 2), 'linewidth', 3)

%% collect responses
% extract last data point
copmpsisx_response_stanceR_stimR_step1_OA = copmpsisx_response_stanceR_stimR_OA(end, :);
copmpsisx_response_stanceR_stimL_step1_OA = copmpsisx_response_stanceR_stimL_OA(end, :);
copmpsisx_response_stanceL_stimR_step1_OA = copmpsisx_response_stanceL_stimR_OA(end, :);
copmpsisx_response_stanceL_stimL_step1_OA = copmpsisx_response_stanceL_stimL_OA(end, :);
copmpsisx_response_stanceR_stimR_step1_HY = copmpsisx_response_stanceR_stimR_HY(end, :);
copmpsisx_response_stanceR_stimL_step1_HY = copmpsisx_response_stanceR_stimL_HY(end, :);
copmpsisx_response_stanceL_stimR_step1_HY = copmpsisx_response_stanceL_stimR_HY(end, :);
copmpsisx_response_stanceL_stimL_step1_HY = copmpsisx_response_stanceL_stimL_HY(end, :);

% figure; hold on
% histogram(copmpsisx_response_stanceR_stimR_step1_OA);
% histogram(copmpsisx_response_stanceR_stimL_step1_OA);
% histogram(copmpsisx_response_stanceL_stimR_step1_OA);
% histogram(copmpsisx_response_stanceL_stimL_step1_OA);
% 
% figure; hold on
% histogram(copmpsisx_response_stanceR_stimR_step1_HY);
% histogram(copmpsisx_response_stanceR_stimL_step1_HY);
% histogram(copmpsisx_response_stanceL_stimR_step1_HY);
% histogram(copmpsisx_response_stanceL_stimL_step1_HY);

% pool different feet
copmpsisx_response_stimR_step1_OA = [copmpsisx_response_stanceR_stimR_step1_OA copmpsisx_response_stanceL_stimR_step1_OA];
copmpsisx_response_stimL_step1_OA = [copmpsisx_response_stanceR_stimL_step1_OA copmpsisx_response_stanceL_stimL_step1_OA];
copmpsisx_response_stimR_step1_HY = [copmpsisx_response_stanceR_stimR_step1_HY copmpsisx_response_stanceL_stimR_step1_HY];
copmpsisx_response_stimL_step1_HY = [copmpsisx_response_stanceR_stimL_step1_HY copmpsisx_response_stanceL_stimL_step1_HY];

% figure; hold on
% histogram(copmpsisx_response_stimR_step1_OA);
% histogram(copmpsisx_response_stimL_step1_OA);
% 
% figure; hold on
% histogram(copmpsisx_response_stimR_step1_HY);
% histogram(copmpsisx_response_stimL_step1_HY);

% pool different stimulus sides
copmpsisx_response_step1_OA = [copmpsisx_response_stimR_step1_OA -copmpsisx_response_stimL_step1_OA];
copmpsisx_response_step1_HY = [copmpsisx_response_stimR_step1_HY -copmpsisx_response_stimL_step1_HY];

% figure; hold on
% histogram(copmpsisx_response_step1_OA);
% histogram(copmpsisx_response_step1_HY);

%% run stats
mean_OA = mean(copmpsisx_response_step1_OA)
mean_HY = mean(copmpsisx_response_step1_HY)
std_OA = std(copmpsisx_response_step1_OA)
std_HY = std(copmpsisx_response_step1_HY)
[h, p] = ttest2(copmpsisx_response_step1_OA, copmpsisx_response_step1_HY)





return






