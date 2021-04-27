function [pos_paced, pace] = selfPacer_spline(pos, time, heelstrike_indices)
    pelvis_pos_belt_splined = spline(heelstrike_indices, pos(heelstrike_indices), (1 : length(time))');

    belt_position_timeseries = pelvis_pos_belt_splined;
    
    pos_paced = pos - belt_position_timeseries;
    pace = deriveByTime(belt_position_timeseries, time);
    
    vel = deriveByTime(pos, time);
    pelvis_vel_belt_splined = deriveByTime(pelvis_pos_belt_splined, time);
    
%     colors = lines(3);
%     markersize = 12;
%     
%     pos_fig = figure; axes('fontsize', 16); hold on; xlim([14, 16])
%     plot(time, pos, 'linewidth', 2, 'color', colors(1, :))
%     plot(time(heelstrike_indices), pos(heelstrike_indices), 'v', 'MarkerFaceColor', colors(2, :), 'MarkerSize', markersize)
%     plot(time, pelvis_pos_belt_splined, 'linewidth', 2, 'color', colors(3, :))
%     xlabel('time (s)'); ylabel('pelvis position (m)')
%     l = legend('pelvis position', 'heel-strikes', 'splined position');
%     set(l, 'location', 'northwest')
%     
%     vel_fig = figure; axes('fontsize', 16); hold on; xlim([10, 20])
%     plot(time, vel, 'linewidth', 2, 'color', colors(1, :))
%     plot(time(heelstrike_indices), vel(heelstrike_indices), 'v', 'MarkerFaceColor', colors(2, :), 'MarkerSize', markersize)
%     plot(time, pelvis_vel_belt_splined, 'linewidth', 2, 'color', colors(3, :))
%     xlabel('time (s)'); ylabel('pelvis velocity (m/s)')
%     l = legend('pelvis velocity', 'heel-strikes', 'splined position time derivative');
%     set(l, 'location', 'northwest')
%     
%     print(pos_fig, 'posBeltSpace_spline', '-djpeg', '-r150')
%     print(vel_fig, 'velBeltSpace_spline', '-djpeg', '-r150')
end