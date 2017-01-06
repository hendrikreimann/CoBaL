% compare force plate data

trial_to_process = 6;


load('subjectInfo.mat');

% load and extract labview data
load(makeFileName(date, subject_id, 'walking', trial_to_process, 'labviewTrajectories'));
fxl_trajectory_labview = fxl_trajectory;
fyl_trajectory_labview = fyl_trajectory;
fzl_trajectory_labview = fzl_trajectory;
mxl_trajectory_labview = mxl_trajectory;
myl_trajectory_labview = myl_trajectory;
mzl_trajectory_labview = mzl_trajectory;
copxl_trajectory_labview = copxl_trajectory;
copyl_trajectory_labview = copyl_trajectory;
fxr_trajectory_labview = fxr_trajectory;
fyr_trajectory_labview = fyr_trajectory;
fzr_trajectory_labview = fzr_trajectory;
mxr_trajectory_labview = mxr_trajectory;
myr_trajectory_labview = myr_trajectory;
mzr_trajectory_labview = mzr_trajectory;
copxr_trajectory_labview = copxr_trajectory;
copyr_trajectory_labview = copyr_trajectory;

load(makeFileName(date, subject_id, 'walking', trial_to_process, 'forceplateTrajectories'));
fxl_trajectory_nexus = fxl_trajectory;
fyl_trajectory_nexus = fyl_trajectory;
fzl_trajectory_nexus = fzl_trajectory;
mxl_trajectory_nexus = mxl_trajectory;
myl_trajectory_nexus = myl_trajectory;
mzl_trajectory_nexus = mzl_trajectory;
copxl_trajectory_nexus = copxl_trajectory;
copyl_trajectory_nexus = copyl_trajectory;
fxr_trajectory_nexus = fxr_trajectory;
fyr_trajectory_nexus = fyr_trajectory;
fzr_trajectory_nexus = fzr_trajectory;
mxr_trajectory_nexus = mxr_trajectory;
myr_trajectory_nexus = myr_trajectory;
mzr_trajectory_nexus = mzr_trajectory;
copxr_trajectory_nexus = copxr_trajectory;
copyr_trajectory_nexus = copyr_trajectory;


figure; axes; hold on; title('fxl')
plot(time_labview, fxl_trajectory_labview, 'linewidth', 2)
plot(time_forceplate, fxl_trajectory_nexus, 'linewidth', 2)

figure; axes; hold on; title('fyl')
plot(time_labview, fyl_trajectory_labview, 'linewidth', 2)
plot(time_forceplate, fyl_trajectory_nexus, 'linewidth', 2)

figure; axes; hold on; title('fzl')
plot(time_labview, fzl_trajectory_labview, 'linewidth', 2)
plot(time_forceplate, fzl_trajectory_nexus, 'linewidth', 2)

figure; axes; hold on; title('copxl')
plot(time_labview, copxl_trajectory_labview, 'linewidth', 2)
plot(time_forceplate, copxl_trajectory_nexus, 'linewidth', 2)

figure; axes; hold on; title('copyl')
plot(time_labview, copyl_trajectory_labview, 'linewidth', 2)
plot(time_forceplate, copyl_trajectory_nexus, 'linewidth', 2)

% figure; axes; hold on; title('copxr')
% plot(time_labview, copxr_trajectory_labview, 'linewidth', 2)
% plot(time_forceplate, copxr_trajectory_nexus, 'linewidth', 2)
% 
% figure; axes; hold on; title('copyr')
% plot(time_labview, copyr_trajectory_labview, 'linewidth', 2)
% plot(time_forceplate, copyr_trajectory_nexus, 'linewidth', 2)

linkaxes(getAllAxes, 'x')
distFig
