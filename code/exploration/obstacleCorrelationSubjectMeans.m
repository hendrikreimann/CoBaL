subjects_to_analyze = {'A', 'AD', 'AK', 'B', 'BO', 'BR', 'C', 'D', 'DC', 'E', 'F', 'FV', 'G', 'IC', 'LH', 'RA', 'SS'}; % all obstacle
load '/Users/reimajbi/Drive_UD/Obstacle/results.mat'

data_label_x = 'cop_to_com_ml_min';
data_label_y = 'com_ml_vel_band_end';
band_index_x = 2;
band_index_y = 2;


number_of_subjects = length(subjects_to_analyze);

data_x_no = zeros(1, number_of_subjects);
data_y_no = zeros(1, number_of_subjects);
data_x_near = zeros(1, number_of_subjects);
data_y_near = zeros(1, number_of_subjects);
data_x_far = zeros(1, number_of_subjects);
data_y_far = zeros(1, number_of_subjects);
for i_subject = 1 : number_of_subjects
   % get indicators
    subject_indicator = strcmp(conditions.subject_list, subjects_to_analyze{i_subject});
    obs_no_indicator = strcmp(conditions.condition_experimental_list, 'OBS_NO');
    obs_near_indicator = strcmp(conditions.condition_experimental_list, 'OBS_NEAR');
    obs_far_indicator = strcmp(conditions.condition_experimental_list, 'OBS_FAR');

    % get data
    data_x_no(i_subject) = mean(variable_data{strcmp(variable_names, data_label_x)}(band_index_x, obs_no_indicator & subject_indicator));
    data_y_no(i_subject) = mean(variable_data{strcmp(variable_names, data_label_y)}(band_index_y, obs_no_indicator & subject_indicator));
    data_x_near(i_subject) = mean(variable_data{strcmp(variable_names, data_label_x)}(band_index_x, obs_near_indicator & subject_indicator));
    data_y_near(i_subject) = mean(variable_data{strcmp(variable_names, data_label_y)}(band_index_y, obs_near_indicator & subject_indicator));
    data_x_far(i_subject) = mean(variable_data{strcmp(variable_names, data_label_x)}(band_index_x, obs_far_indicator & subject_indicator));
    data_y_far(i_subject) = mean(variable_data{strcmp(variable_names, data_label_y)}(band_index_y, obs_far_indicator & subject_indicator));
end


J_no = (data_x_no - mean(data_x_no))'\(data_y_no - mean(data_y_no))';
J_near = (data_y_near - mean(data_y_near)) * pinv(data_x_near - mean(data_x_near));
J_far = (data_y_far - mean(data_y_far)) * pinv(data_x_far - mean(data_x_far));
[correlation_c_no, correlation_p_no] = corr(data_x_no', data_y_no');
[correlation_c_near, correlation_p_near] = corr(data_x_near', data_y_near');
[correlation_c_far, correlation_p_far] = corr(data_x_far', data_y_far');
rho_no = correlation_c_no;
rho_near = correlation_c_near;
rho_far = correlation_c_far;

linefit_no = polyfit(data_x_no, data_y_no, 1);
linefit_near = polyfit(data_x_near, data_y_near, 1);
linefit_far = polyfit(data_x_far, data_y_far, 1);

color_far = [0.2, 0.7, 0.07];
color_near = [0.7, 0.2, 0.07];
color_no = [0.2, 0.07, 0.7];


figure; axes; hold on;
xlabel([strrep(data_label_x, '_', ' ') '(' num2str(band_index_x) ')']);
ylabel([strrep(data_label_y, '_', ' ') '(' num2str(band_index_y) ')']);

plot(data_x_no, data_y_no, 'o', 'color', color_no);
plot(data_x_near, data_y_near, 'o', 'color', color_near);
plot(data_x_far, data_y_far, 'o', 'color', color_far);
xlim = get(gca, 'xlim');
ylim = get(gca, 'ylim');
plot(mean(data_x_no)+[-1 1], (mean(data_x_no)+[-1 1])*linefit_no(1)+linefit_no(2), 'color', color_no, 'linewidth', 3);
plot(mean(data_x_near)+[-1 1], (mean(data_x_near)+[-1 1])*linefit_near(1)+linefit_near(2), 'color', color_near, 'linewidth', 3);
plot(mean(data_x_far)+[-1 1], (mean(data_x_far)+[-1 1])*linefit_far(1)+linefit_far(2), 'color', color_far, 'linewidth', 3);
plot(mean(data_x_no), mean(data_y_no), 'd', 'color', [1 1 1]*0.5, 'linewidth', 3, 'markersize', 20, 'MarkerFaceColor', color_no);
plot(mean(data_x_near), mean(data_y_near), 'd', 'color', [1 1 1]*0.5, 'linewidth', 3, 'markersize', 20, 'MarkerFaceColor', color_near);
plot(mean(data_x_far), mean(data_y_far), 'd', 'color', [1 1 1]*0.5, 'linewidth', 3, 'markersize', 20, 'MarkerFaceColor', color_far);
set(gca, 'xlim', xlim);
set(gca, 'ylim', ylim);
text(xlim(2), ylim(2) - diff(ylim)*0.01, ['J = ' num2str(J_no) ', c = ' num2str(correlation_c_no) ', p = ' num2str(correlation_p_no)], 'color', color_no, 'HorizontalAlignment', 'right')
text(xlim(2), ylim(2) - diff(ylim)*0.05, ['J = ' num2str(J_near) ', c = ' num2str(correlation_c_near) ', p = ' num2str(correlation_p_near)], 'color', color_near, 'HorizontalAlignment', 'right')
text(xlim(2), ylim(2) - diff(ylim)*0.09, ['J = ' num2str(J_far) ', c = ' num2str(correlation_c_far) ', p = ' num2str(correlation_p_far)], 'color', color_far, 'HorizontalAlignment', 'right')

