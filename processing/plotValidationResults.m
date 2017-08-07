% plotValidationResults

function plotValidationResults(varargin)
    %% parse input
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'subjects', [])
    addParameter(parser, 'save', false)
    addParameter(parser, 'settings', 'plotSettings.txt')
    parse(parser, varargin{:})
    subjects = parser.Results.subjects;
    
     study_settings_file = '';
    if exist('studySettings.txt', 'file')
        study_settings_file = 'studySettings.txt';
    end    
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
   
    
    mean_treadmill_speed_trajectories_all = [];
    mean_coefficient_of_multiple_correlation_left_all = [];
    mean_coefficient_of_multiple_correlation_right_all = [];
    std_coefficient_of_multiple_correlation_left_all = [];
    std_coefficient_of_multiple_correlation_right_all = [];
    mean_root_mean_square_error_left_all = [];
    mean_root_mean_square_error_right_all = [];
    std_root_mean_square_error_left_all = [];
    std_root_mean_square_error_right_all = [];
    mean_angular_velocity_mocap_left_all = [];
    mean_angular_velocity_mocap_right_all = [];
    std_angular_velocity_mocap_left_all = [];
    std_angular_velocity_mocap_right_all = [];
    mean_mocap_peak_amplitude_left_all = [];
    mean_mocap_peak_amplitude_right_all = [];
    std_mocap_peak_amplitude_left_all = [];
    std_mocap_peak_amplitude_right_all = [];
    std_treadmill_speed_trajectories_all = [];
    mean_angular_acceleration_mocap_left_all = [];
    mean_angular_acceleration_mocap_right_all = [];
    std_angular_acceleration_mocap_left_all = [];
    std_angular_acceleration_mocap_right_all = [];
    percent_error_left_all = [];
    percent_error_right_all = [];
    
    data_folder_list = determineDataStructure(subjects);
    for i_folder = 1 : length(data_folder_list)
        % load data
        data_path = data_folder_list{i_folder};
        load([data_path filesep 'subjectInfo.mat'], 'date', 'subject_id');
        load([data_path filesep 'analysis' filesep date '_' subject_id '_validation.mat']);    
        
        mean_treadmill_speed_trajectories_all = [mean_treadmill_speed_trajectories_all; mean_treadmill_speed_trajectories];
        mean_coefficient_of_multiple_correlation_left_all = [mean_coefficient_of_multiple_correlation_left_all, mean_coefficient_of_multiple_correlation_left];
        mean_coefficient_of_multiple_correlation_right_all = [mean_coefficient_of_multiple_correlation_right_all, mean_coefficient_of_multiple_correlation_right];
        std_coefficient_of_multiple_correlation_left_all = [std_coefficient_of_multiple_correlation_left_all, std_coefficient_of_multiple_correlation_left];
        std_coefficient_of_multiple_correlation_right_all = [std_coefficient_of_multiple_correlation_right_all, std_coefficient_of_multiple_correlation_right];
        mean_root_mean_square_error_left_all = [mean_root_mean_square_error_left_all, mean_root_mean_square_error_left];
        mean_root_mean_square_error_right_all = [mean_root_mean_square_error_right_all, mean_root_mean_square_error_right];
        std_root_mean_square_error_left_all = [std_root_mean_square_error_left_all, std_root_mean_square_error_left];
        std_root_mean_square_error_right_all = [std_root_mean_square_error_right_all, std_root_mean_square_error_right];
        mean_angular_velocity_mocap_left_all = [mean_angular_velocity_mocap_left_all; mean_angular_velocity_mocap_left];
        mean_angular_velocity_mocap_right_all = [mean_angular_velocity_mocap_right_all; mean_angular_velocity_mocap_right];
        std_angular_velocity_mocap_left_all = [std_angular_velocity_mocap_left_all; std_angular_velocity_mocap_left];
        std_angular_velocity_mocap_right_all = [std_angular_velocity_mocap_right_all; std_angular_velocity_mocap_right];
        mean_mocap_peak_amplitude_left_all = [mean_mocap_peak_amplitude_left_all; mean_mocap_peak_amplitude_left];
        mean_mocap_peak_amplitude_right_all = [mean_mocap_peak_amplitude_right_all; mean_mocap_peak_amplitude_right];
        std_mocap_peak_amplitude_left_all = [std_mocap_peak_amplitude_left_all; std_mocap_peak_amplitude_left];
        std_mocap_peak_amplitude_right_all = [std_mocap_peak_amplitude_right_all; std_mocap_peak_amplitude_right];
        std_treadmill_speed_trajectories_all = [std_treadmill_speed_trajectories_all; std_treadmill_speed_trajectories];
        mean_angular_acceleration_mocap_left_all = [mean_angular_acceleration_mocap_left_all; mean_angular_acceleration_mocap_left];
        mean_angular_acceleration_mocap_right_all = [mean_angular_acceleration_mocap_right_all; mean_angular_acceleration_mocap_right];
        std_angular_acceleration_mocap_left_all = [std_angular_acceleration_mocap_left_all; std_angular_acceleration_mocap_left];
        std_angular_acceleration_mocap_right_all = [std_angular_acceleration_mocap_right_all; std_angular_acceleration_mocap_right];
        percent_error_left_all = [percent_error_left_all, percent_error_left];
        percent_error_right_all = [percent_error_right_all, percent_error_right];
    end
    
%     subjects = {'1', '2', '3', '4', '5', '6'};
    mean_coefficient_of_multiple_correlation_left_all  = mean_coefficient_of_multiple_correlation_left_all';
    mean_coefficient_of_multiple_correlation_right_all = mean_coefficient_of_multiple_correlation_right_all';
    std_coefficient_of_multiple_correlation_left_all = std_coefficient_of_multiple_correlation_left_all';
    std_coefficient_of_multiple_correlation_right_all = std_coefficient_of_multiple_correlation_right_all';
    mean_root_mean_square_error_left_all = mean_root_mean_square_error_left_all';
    mean_root_mean_square_error_right_all = mean_root_mean_square_error_right_all';
    std_root_mean_square_error_left_all = std_root_mean_square_error_left_all';
    std_root_mean_square_error_right_all = std_root_mean_square_error_right_all';
    percent_error_left_all = percent_error_left_all';
    percent_error_right_all = percent_error_right_all';
    
    % Create Subject descriptive table
    T1 = table;
    
    Treadmill_Speed = [mean_treadmill_speed_trajectories_all; std_treadmill_speed_trajectories_all];
    Treadmill_Speed(1:2:end-1) = mean_treadmill_speed_trajectories_all;
    Treadmill_Speed(2:2:end) = std_treadmill_speed_trajectories_all; 
    T1.Treadmill_Speed = Treadmill_Speed;
    
    Armswing_Amplitude_Left = [mean_mocap_peak_amplitude_left_all; std_mocap_peak_amplitude_left_all];
    Armswing_Amplitude_Left(1:2:end-1) = mean_mocap_peak_amplitude_left_all;
    Armswing_Amplitude_Left(2:2:end) = std_mocap_peak_amplitude_left_all;
    T1.Armswing_Amplitude_Left = Armswing_Amplitude_Left;
    
    Armswing_Amplitude_Right = [mean_mocap_peak_amplitude_right_all; std_mocap_peak_amplitude_right_all];
    Armswing_Amplitude_Right(1:2:end-1) = mean_mocap_peak_amplitude_right_all;
    Armswing_Amplitude_Right(2:2:end) = std_mocap_peak_amplitude_right_all;
    T1.Armswing_Amplitude_Right = Armswing_Amplitude_Right;
    
    Armswing_Velocity_Left = [mean_angular_velocity_mocap_left_all; std_angular_velocity_mocap_left_all];
    Armswing_Velocity_Left(1:2:end-1) = mean_angular_velocity_mocap_left_all;
    Armswing_Velocity_Left(2:2:end) = std_angular_velocity_mocap_left_all;
    T1.Armswing_Velocity_Left = Armswing_Velocity_Left;
    
    Armswing_Velocity_Right = [mean_angular_velocity_mocap_right_all; std_angular_velocity_mocap_right_all];
    Armswing_Velocity_Right(1:2:end-1) = mean_angular_velocity_mocap_right_all;
    Armswing_Velocity_Right(2:2:end) = std_angular_velocity_mocap_right_all;
    T1.Armswing_Velocity_Right = Armswing_Velocity_Right;
    
    Armswing_Acceleration_Left = [mean_angular_acceleration_mocap_left_all; std_angular_acceleration_mocap_left_all];
    Armswing_Acceleration_Left(1:2:end-1) = mean_angular_acceleration_mocap_left_all;
    Armswing_Acceleration_Left(2:2:end) = std_angular_acceleration_mocap_left_all;
    T1.Armswing_Acceleration_Left = Armswing_Acceleration_Left;
    
    Armswing_Acceleration_Right = [mean_angular_acceleration_mocap_right_all; std_angular_acceleration_mocap_right_all];
    Armswing_Acceleration_Right(1:2:end-1) = mean_angular_acceleration_mocap_right_all; 
    Armswing_Acceleration_Right(2:2:end) = std_angular_acceleration_mocap_right_all;
    T1.Armswing_Acceleration_Right = Armswing_Acceleration_Right;
    
%     T1.Properties.RowNames = subjects;
    
%     (mean_treadmill_speed_trajectories_all, mean_mocap_peak_amplitude_left_all, mean_mocap_peak_amplitude_right_all, mean_angular_velocity_mocap_left_all, ...
%         mean_angular_velocity_mocap_right_all, mean_angular_acceleration_mocap_left_all, mean_angular_acceleration_mocap_right_all, ...
%         'RowNames', subjects); %'VariableNames',{'Treadmill Speed' 'Armswing Amplitude Left' 'Armswing Amplitude Right' 'Armswing Velocity Left', ...
        %'Armswing Velocity Right', 'Armswing Acceleration Left', 'Armswing Acceleration Right'});
    
        
    % Create RMS table
    T2 = table
    
    Beta1 = [mean_root_mean_square_error_left_all(:,1); std_root_mean_square_error_left_all(:,1); percent_error_left_all(:,1)]; 
    Beta1(1:3:end-2) = mean_root_mean_square_error_left_all(:,1);
    Beta1(2:3:end-1) = std_root_mean_square_error_left_all(:,1);
    Beta1(3:3:end) = percent_error_left_all(:,1);
    T2.Beta1 = Beta1;
    
    Beta2 = [mean_root_mean_square_error_left_all(:,2); std_root_mean_square_error_left_all(:,2); percent_error_left_all(:,2)];
    Beta2(1:3:end-2) = mean_root_mean_square_error_left_all(:,2);
    Beta2(2:3:end-1) = std_root_mean_square_error_left_all(:,2);
    Beta2(3:3:end) = percent_error_left_all(:,2);
    T2.Beta2 = Beta2;
   
    Beta3 = [mean_root_mean_square_error_left_all(:,3); std_root_mean_square_error_left_all(:,3); percent_error_left_all(:,3)];
    Beta3(1:3:end-2) = mean_root_mean_square_error_left_all(:,3);
    Beta3(2:3:end-1) = std_root_mean_square_error_left_all(:,3);
    Beta3(3:3:end) = percent_error_left_all(:,3);
    T2.Beta3 = Beta3;
    
    Beta4 = [mean_root_mean_square_error_left_all(:,4); std_root_mean_square_error_left_all(:,4); percent_error_left_all(:,4)];
    Beta4(1:3:end-2) = mean_root_mean_square_error_left_all(:,4);
    Beta4(2:3:end-1) = std_root_mean_square_error_left_all(:,4);
    Beta4(3:3:end) = percent_error_left_all(:,4);
    T2.Beta4 = Beta4;
    

    Beta5 = [mean_root_mean_square_error_left_all(:,5); std_root_mean_square_error_left_all(:,5); percent_error_left_all(:,5)];
    Beta5(1:3:end-2) = mean_root_mean_square_error_left_all(:,5);
    Beta5(2:3:end-1) = std_root_mean_square_error_left_all(:,5);
    Beta5(3:3:end) = percent_error_left_all(:,5);
    T2.Beta5 = Beta5;
    
    Beta6 = [mean_root_mean_square_error_left_all(:,6); std_root_mean_square_error_left_all(:,6); percent_error_left_all(:,6)];
    Beta6(1:3:end-2) = mean_root_mean_square_error_left_all(:,6);
    Beta6(2:3:end-1) = std_root_mean_square_error_left_all(:,6);
    Beta6(3:3:end) = percent_error_left_all(:,6);
    T2.Beta6 = Beta6;
    
    Beta7 = [mean_root_mean_square_error_left_all(:,7); std_root_mean_square_error_left_all(:,7); percent_error_left_all(:,7)];
    Beta7(1:3:end-2) = mean_root_mean_square_error_left_all(:,7);
    Beta7(2:3:end-1) = std_root_mean_square_error_left_all(:,7);
    Beta7(3:3:end) = percent_error_left_all(:,7);
    T2.Beta7 = Beta7;
    
    Beta8 = [mean_root_mean_square_error_left_all(:,8); std_root_mean_square_error_left_all(:,8); percent_error_left_all(:,8)];
    Beta8(1:3:end-2) = mean_root_mean_square_error_left_all(:,8);
    Beta8(2:3:end-1) = std_root_mean_square_error_left_all(:,8);
    Beta8(3:3:end) = percent_error_left_all(:,8);
    T2.Beta8 = Beta8;
    
    % Create CMC table
    T3 = table;
    
    Beta1 = [mean_coefficient_of_multiple_correlation_left_all(:,1); std_coefficient_of_multiple_correlation_left_all(:,1)]; 
    Beta1(1:2:end-1) = mean_coefficient_of_multiple_correlation_left_all(:,1);
    Beta1(2:2:end) = std_coefficient_of_multiple_correlation_left_all(:,1);
    T3.Beta1 = Beta1;
    
    Beta2 = [mean_coefficient_of_multiple_correlation_left_all(:,2); std_coefficient_of_multiple_correlation_left_all(:,2)];
    Beta2(1:2:end-1) = mean_coefficient_of_multiple_correlation_left_all(:,2);
    Beta2(2:2:end) = std_coefficient_of_multiple_correlation_left_all(:,2);
    T3.Beta2 = Beta2;
   
    Beta3 = [mean_coefficient_of_multiple_correlation_left_all(:,3); std_coefficient_of_multiple_correlation_left_all(:,3)];
    Beta3(1:2:end-1) = mean_coefficient_of_multiple_correlation_left_all(:,3);
    Beta3(2:2:end) = std_coefficient_of_multiple_correlation_left_all(:,3);
    T3.Beta3 = Beta3;
    
    Beta4 = [mean_coefficient_of_multiple_correlation_left_all(:,4); std_coefficient_of_multiple_correlation_left_all(:,4)];
    Beta4(1:2:end-1) = mean_coefficient_of_multiple_correlation_left_all(:,4);
    Beta4(2:2:end) = std_coefficient_of_multiple_correlation_left_all(:,4);
    T3.Beta4 = Beta4;
    

    Beta5 = [mean_coefficient_of_multiple_correlation_left_all(:,5); std_coefficient_of_multiple_correlation_left_all(:,5)];
    Beta5(1:2:end-1) = mean_coefficient_of_multiple_correlation_left_all(:,5);
    Beta5(2:2:end) = std_coefficient_of_multiple_correlation_left_all(:,5);
    T3.Beta5 = Beta5;
    
    Beta6 = [mean_coefficient_of_multiple_correlation_left_all(:,6); std_coefficient_of_multiple_correlation_left_all(:,6)];
    Beta6(1:2:end-1) = mean_coefficient_of_multiple_correlation_left_all(:,6);
    Beta6(2:2:end) = std_coefficient_of_multiple_correlation_left_all(:,6);
    T3.Beta6 = Beta6;
    
    Beta7 = [mean_coefficient_of_multiple_correlation_left_all(:,7); std_coefficient_of_multiple_correlation_left_all(:,7)];
    Beta7(1:2:end-1) = mean_coefficient_of_multiple_correlation_left_all(:,7);
    Beta7(2:2:end) = std_coefficient_of_multiple_correlation_left_all(:,7);
    T3.Beta7 = Beta7;
    
    Beta8 = [mean_coefficient_of_multiple_correlation_left_all(:,8); std_coefficient_of_multiple_correlation_left_all(:,8)];
    Beta8(1:2:end-1) = mean_coefficient_of_multiple_correlation_left_all(:,8);
    Beta8(2:2:end) = std_coefficient_of_multiple_correlation_left_all(:,8);
    T3.Beta8 = Beta8;
    
    % Create velocity vs RMS
    figure
%     mean_angular_velocity_mocap_left_all = sort(mean_angular_velocity_mocap_left_all);
    shapes = {'-^','-o','-*','-x','-+','-p','-d','-v'};
    Vel_RMS_matrix = [mean_angular_velocity_mocap_left_all, mean_root_mean_square_error_left_all];
    [X,I] = sort(Vel_RMS_matrix(:,1));
    Vel_RMS_matrix = Vel_RMS_matrix(I,:);
    for i = 1:8
        plot(Vel_RMS_matrix(:,1), Vel_RMS_matrix(:,i+1),shapes{i},'MarkerSize',12,'LineWidth',2)
        hold on
    end
    xlim([20 160])
    xlabel('Average Armswing Angular Velocity (deg/sec)')
    ylabel('RMS (sensor vs mocap)')
    legend('alpha1','alpha2','alpha3','alpha4','alpha5','alpha6','alpha7','alpha8')
    
    figure
    % create amplitude vs RMS
    Amp_RMS_matrix = [mean_mocap_peak_amplitude_left_all  ,mean_root_mean_square_error_left_all];
    [X,I] = sort(Amp_RMS_matrix(:,1));
    Amp_RMS_matrix = Amp_RMS_matrix(I,:);
    for i = 1:8
       plot( Amp_RMS_matrix(:,1), Amp_RMS_matrix(:,i+1),shapes{i},'MarkerSize',12,'LineWidth',2)
       hold on
    end
    xlim([15 55]);
    xlabel('Average Armswing Amplitude (deg)')
    ylabel('RMS (sensor vs mocap)')
    legend('alpha1','alpha2','alpha3','alpha4','alpha5','alpha6','alpha7','alpha8')
    
%     save('tables','T1', 'T2');
    
    writetable(T1,'descriptives')
    writetable(T2,'RMS');
    writetable(T3,'CMC')
%     
%     if ~exist('figures', 'dir')
%             mkdir('figures')
%     end
%     filename = ['figures' filesep 'noLabels' filesep get(figure_handles(i_figure), 'UserData')];
%                 saveas(figure_handles(i_figure), filename, parser.Results.format)

end