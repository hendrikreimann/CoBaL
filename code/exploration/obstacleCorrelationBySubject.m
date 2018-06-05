subjects_to_analyze = {'A', 'AD', 'AK', 'B', 'BO', 'BR', 'C', 'D', 'DC', 'E', 'F', 'FV', 'G', 'IC', 'LH', 'RA', 'SS'}; % all obstacle
load '/Users/reimajbi/Drive_UD/Obstacle/results.mat'

data_label_x = 'cop_to_com_ml_max';
data_label_y = 'step_time';
band_index_x = 2;
band_index_y = 3;

data_label_x = 'step_time';
data_label_y = 'com_ml_min';
band_index_x = 3;
band_index_y = 2;

data_label_x = 'com_ml_vel_band_end';
data_label_y = 'step_time';
band_index_x = 2;
band_index_y = 3;

% these three are a chain, but somehow not transitive
data_label_x = 'cop_to_com_ml_min';
data_label_y = 'com_ml_vel_band_end';
band_index_x = 2;
band_index_y = 2;

% data_label_x = 'com_ml_vel_band_end';
% data_label_y = 'step_time';
% band_index_x = 2;
% band_index_y = 3;
% % 
% data_label_x = 'cop_to_com_ml_min';
% data_label_y = 'step_time';
% band_index_x = 2;
% band_index_y = 3;


number_of_subjects = length(subjects_to_analyze);

J_no = zeros(1, number_of_subjects);
J_near = zeros(1, number_of_subjects);
J_far = zeros(1, number_of_subjects);
r_no = zeros(1, number_of_subjects);
r_near = zeros(1, number_of_subjects);
r_far = zeros(1, number_of_subjects);
for i_subject = 1 : number_of_subjects
   % get indicators
    subject_indicator = strcmp(conditions.subject_list, subjects_to_analyze{i_subject});
    obs_no_indicator = strcmp(conditions.condition_experimental_list, 'OBS_NO');
    obs_near_indicator = strcmp(conditions.condition_experimental_list, 'OBS_NEAR');
    obs_far_indicator = strcmp(conditions.condition_experimental_list, 'OBS_FAR');

    % get data
    data_x_no = variable_data{strcmp(variable_names, data_label_x)}(band_index_x, obs_no_indicator & subject_indicator);
    data_y_no = variable_data{strcmp(variable_names, data_label_y)}(band_index_y, obs_no_indicator & subject_indicator);
    data_x_near = variable_data{strcmp(variable_names, data_label_x)}(band_index_x, obs_near_indicator & subject_indicator);
    data_y_near = variable_data{strcmp(variable_names, data_label_y)}(band_index_y, obs_near_indicator & subject_indicator);
    data_x_far = variable_data{strcmp(variable_names, data_label_x)}(band_index_x, obs_far_indicator & subject_indicator);
    data_y_far = variable_data{strcmp(variable_names, data_label_y)}(band_index_y, obs_far_indicator & subject_indicator); 
   
    J_no(i_subject) = (data_x_no - mean(data_x_no))'\(data_y_no - mean(data_y_no))';
    J_near(i_subject) = (data_y_near - mean(data_y_near)) * pinv(data_x_near - mean(data_x_near));
    J_far(i_subject) = (data_y_far - mean(data_y_far)) * pinv(data_x_far - mean(data_x_far));
    [correlation_c_no, correlation_p_no] = corr(data_x_no', data_y_no');
    [correlation_c_near, correlation_p_near] = corr(data_x_near', data_y_near');
    [correlation_c_far, correlation_p_far] = corr(data_x_far', data_y_far');
    rho_no(i_subject) = correlation_c_no;
    rho_near(i_subject) = correlation_c_near;
    rho_far(i_subject) = correlation_c_far;
end


disp(['X = ' data_label_x '(' num2str(band_index_x) '), Y = ' data_label_y '(' num2str(band_index_y) ')']);
disp(['No  : J = ' num2str(mean(J_no)), ', rho = ' num2str(mean(rho_no))]);
disp(['Near: J = ' num2str(mean(J_near)), ', rho = ' num2str(mean(rho_near))]);
disp(['Far : J = ' num2str(mean(J_far)), ', rho = ' num2str(mean(rho_far))]);
