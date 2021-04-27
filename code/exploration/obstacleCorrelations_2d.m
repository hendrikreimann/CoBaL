% load '/Users/hendrikreimann/Drive_UD/Obstacle/results.mat'

subject_to_analyze = 'all';
% subject_to_analyze = 'F';

data_label_x = 'cop_to_com_ml_max';
data_label_z = 'step_time';
band_index_x = 2;
band_index_z = 3;

% data_label_x = 'cop_to_com_ml_min';
% data_label_z = 'step_time';
% band_index_x = 2;
% band_index_z = 3;

% these three are a chain, but somehow not transitive
data_label_x = 'cop_to_com_ml_min';
data_label_z = 'com_ml_vel_band_end';
band_index_x = 2;
band_index_z = 2;

data_label_x = 'com_ml_vel_band_end';
data_label_z = 'step_time';
band_index_x = 2;
band_index_z = 3;

data_label_x = 'com_ml_vel_band_end';
data_label_z = 'step_time';
band_index_x = 2;
band_index_z = 3;

data_label_x = 'cop_to_com_ml_band_end';
data_label_y = 'com_ml_vel_band_end';
data_label_z = 'step_time';
band_index_x = 2;
band_index_y = 2;
band_index_z = 3;

% get indicators
subject_indicator = strcmp(conditions.subject_list, subject_to_analyze);
if strcmp(subject_to_analyze, 'all')
    subject_indicator = true(size(subject_indicator));
end
obs_no_indicator = strcmp(conditions.condition_experimental_list, 'OBS_NO');
obs_near_indicator = strcmp(conditions.condition_experimental_list, 'OBS_NEAR');
obs_far_indicator = strcmp(conditions.condition_experimental_list, 'OBS_FAR');

% get data
data_x_no = variable_data{strcmp(variable_names, data_label_x)}(band_index_x, obs_no_indicator & subject_indicator);
data_y_no = variable_data{strcmp(variable_names, data_label_y)}(band_index_y, obs_no_indicator & subject_indicator);
data_z_no = variable_data{strcmp(variable_names, data_label_z)}(band_index_z, obs_no_indicator & subject_indicator);
data_x_near = variable_data{strcmp(variable_names, data_label_x)}(band_index_x, obs_near_indicator & subject_indicator);
data_y_near = variable_data{strcmp(variable_names, data_label_y)}(band_index_y, obs_near_indicator & subject_indicator);
data_z_near = variable_data{strcmp(variable_names, data_label_z)}(band_index_z, obs_near_indicator & subject_indicator);
data_x_far = variable_data{strcmp(variable_names, data_label_x)}(band_index_x, obs_far_indicator & subject_indicator);
data_y_far = variable_data{strcmp(variable_names, data_label_y)}(band_index_y, obs_far_indicator & subject_indicator);
data_z_far = variable_data{strcmp(variable_names, data_label_z)}(band_index_z, obs_far_indicator & subject_indicator);

percentile_to_remove = 5;
[~, data_x_no_outlier_indices] = removeOutliers(data_x_no, percentile_to_remove);
[~, data_y_no_outlier_indices] = removeOutliers(data_y_no, percentile_to_remove);
[~, data_z_no_outlier_indices] = removeOutliers(data_z_no, percentile_to_remove);
data_no_outlier_indices = data_x_no_outlier_indices | data_y_no_outlier_indices | data_z_no_outlier_indices;
data_x_no(data_no_outlier_indices) = [];
data_y_no(data_no_outlier_indices) = [];
data_z_no(data_no_outlier_indices) = [];

[~, data_x_near_outlier_indices] = removeOutliers(data_x_near, percentile_to_remove);
[~, data_y_near_outlier_indices] = removeOutliers(data_y_near, percentile_to_remove);
[~, data_z_near_outlier_indices] = removeOutliers(data_z_near, percentile_to_remove);
data_near_outlier_indices = data_x_near_outlier_indices | data_y_near_outlier_indices | data_z_near_outlier_indices;
data_x_near(data_near_outlier_indices) = [];
data_y_near(data_near_outlier_indices) = [];
data_z_near(data_near_outlier_indices) = [];

[~, data_x_far_outlier_indices] = removeOutliers(data_x_far, percentile_to_remove);
[~, data_y_far_outlier_indices] = removeOutliers(data_y_far, percentile_to_remove);
[~, data_z_far_outlier_indices] = removeOutliers(data_z_far, percentile_to_remove);
data_far_outlier_indices = data_x_far_outlier_indices | data_y_far_outlier_indices | data_z_far_outlier_indices;
data_x_far(data_far_outlier_indices) = [];
data_y_far(data_far_outlier_indices) = [];
data_z_far(data_far_outlier_indices) = [];


% analyze
% J_no = (data_x_no - mean(data_x_no))'\(data_z_no - mean(data_z_no))';
[correlation_c_no, correlation_p_no] = corr([data_x_no', data_y_no'], data_z_no');


linefit_no = polyfit([data_x_no', data_y_no'], data_z_no', 1);
linefit_near = polyfit(data_x_near, data_z_near, 1);
linefit_far = polyfit(data_x_far, data_z_far, 1);

% plot
color_far = [0.2, 0.7, 0.07];
color_near = [0.7, 0.2, 0.07];
color_no = [0.2, 0.07, 0.7];

figure; axes; hold on;
xlabel([strrep(data_label_x, '_', ' ') '(' num2str(band_index_x) ')']);
ylabel([strrep(data_label_z, '_', ' ') '(' num2str(band_index_z) ')']);

plot(data_x_no, data_z_no, 'o', 'color', color_no);
plot(data_x_near, data_z_near, 'o', 'color', color_near);
plot(data_x_far, data_z_far, 'o', 'color', color_far);
xlim = get(gca, 'xlim');
ylim = get(gca, 'ylim');
plot(mean(data_x_no)+[-1 1], (mean(data_x_no)+[-1 1])*linefit_no(1)+linefit_no(2), 'color', color_no, 'linewidth', 3);
plot(mean(data_x_near)+[-1 1], (mean(data_x_near)+[-1 1])*linefit_near(1)+linefit_near(2), 'color', color_near, 'linewidth', 3);
plot(mean(data_x_far)+[-1 1], (mean(data_x_far)+[-1 1])*linefit_far(1)+linefit_far(2), 'color', color_far, 'linewidth', 3);
plot(mean(data_x_no), mean(data_z_no), 'd', 'color', [1 1 1]*0.5, 'linewidth', 3, 'markersize', 20, 'MarkerFaceColor', color_no);
plot(mean(data_x_near), mean(data_z_near), 'd', 'color', [1 1 1]*0.5, 'linewidth', 3, 'markersize', 20, 'MarkerFaceColor', color_near);
plot(mean(data_x_far), mean(data_z_far), 'd', 'color', [1 1 1]*0.5, 'linewidth', 3, 'markersize', 20, 'MarkerFaceColor', color_far);
set(gca, 'xlim', xlim);
set(gca, 'ylim', ylim);
text(xlim(2), ylim(2) - diff(ylim)*0.01, ['J = ' num2str(J_no) ', c = ' num2str(correlation_c_no) ', p = ' num2str(correlation_p_no)], 'color', color_no, 'HorizontalAlignment', 'right')
text(xlim(2), ylim(2) - diff(ylim)*0.05, ['J = ' num2str(J_near) ', c = ' num2str(correlation_c_near) ', p = ' num2str(correlation_p_near)], 'color', color_near, 'HorizontalAlignment', 'right')
text(xlim(2), ylim(2) - diff(ylim)*0.09, ['J = ' num2str(J_far) ', c = ' num2str(correlation_c_far) ', p = ' num2str(correlation_p_far)], 'color', color_far, 'HorizontalAlignment', 'right')

