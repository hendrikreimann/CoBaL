%     This file is part of the CoBaL code base
%     Copyright (C) 2017 Hendrik Reimann <hendrikreimann@gmail.com>
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.% compare the kinematic tree against the kinematic chain

% class for the display of a stick figure

classdef stickFigure < handle
    properties
        scene_figure;
        scene_axes;
        scene_bound;
        
        % data
        number_of_markers;
%         number_of_virtual_markers;
        segment_com_plots;
        marker_plots;
%         virtual_marker_plots;
        marker_labels;
%         virtual_marker_labels;
        
        line_plots;
        line_marker_indices = [];
        
%         virtual_line_plots;
%         virtual_line_marker_indices = [];
        
        marker_positions = [];
%         virtual_marker_positions = [];
        
        % colors
        marker_color = [1 1 1]*0.5;
        
        % flags
        show_segment_coms = false;
        show_connection_lines = true;
        show_markers = true;
        show_marker_labels = false;
    end 
    methods
%         function this = stickFigure(markerPositions, markerLabels, virtualMarkerPositions, virtualMarkerLabels, scene_bound, axesHandle)
        function this = stickFigure(markerPositions, markerLabels, scene_bound, axesHandle)
            if nargin < 3
                this.scene_bound = [-1 1; -1 1; -1 1];
            else
                this.scene_bound = scene_bound;
            end
            if nargin < 4
                this.scene_figure = figure( 'Position', [500, 500, 600, 600], 'Name', 'scene' );
                this.scene_axes = axes( 'Position', [0.1 0.1 0.8 0.8]);
            else
                this.scene_figure = get(axesHandle, 'parent');
                this.scene_axes = axesHandle;
            end
            
            % prepare axes
            axis equal; hold on;
            plot3([this.scene_bound(1, 1), this.scene_bound(1, 2)], [0, 0]+mean(this.scene_bound(2, :)), [0, 0]+mean(this.scene_bound(3, :)), 'color', 'k', 'Linewidth', 1, 'Linestyle', ':');
            plot3([0, 0]+mean(this.scene_bound(1, :)), [this.scene_bound(2, 1), this.scene_bound(2, 2)], [0, 0]+mean(this.scene_bound(3, :)), 'color', 'k', 'Linewidth', 1, 'Linestyle', ':');
            plot3([0, 0]+mean(this.scene_bound(1, :)), [0, 0]+mean(this.scene_bound(2, :)), [this.scene_bound(3, 1), this.scene_bound(3, 2)], 'color', 'k', 'Linewidth', 1, 'Linestyle', ':');

            % set up marker plots
            this.number_of_markers = length(markerLabels);
            this.marker_labels = markerLabels;
            this.marker_plots = cell(this.number_of_markers, 1);
            for i_marker = 1 : this.number_of_markers
                this.marker_plots{i_marker} = plot3(0, 0, 0, 'color', this.marker_color, 'markersize', 10, 'linewidth', 2, 'Marker', 'o');
            end
%             this.number_of_virtual_markers = length(virtualMarkerLabels);
%             this.virtual_marker_labels = virtualMarkerLabels;
%             this.virtual_marker_plots = cell(this.number_of_virtual_markers, 1);
%             for i_marker = 1 : this.number_of_virtual_markers
%                 this.virtual_marker_plots{i_marker} = plot3(0, 0, 0, 'color', this.marker_color, 'markersize', 10, 'linewidth', 3, 'Marker', 'x');
%             end
            
            % groom
            set(gca,'xlim',[this.scene_bound(1, 1), this.scene_bound(1, 2)],'ylim',[this.scene_bound(2, 1), this.scene_bound(2, 2)], 'zlim',[this.scene_bound(3, 1), this.scene_bound(3, 2)]);
            xlabel('x');
            ylabel('y');
            zlabel('z');
            
%             daspect([5 5 1])
%             axis tight
            view(-50,30)
%             camlight left
%             camlight('headlight')
            camlight(30, 30)
%             colormap cool
%             colormap autumn
            alpha(.4)
            
%             this.update(markerPositions, virtualMarkerPositions);
            this.update(markerPositions);
        end
        
%         function update(this, marker_positions, virtual_marker_positions)
        function update(this, marker_positions)
            if nargin > 1
                this.marker_positions = marker_positions;
            end
%             if nargin > 2
%                 this.virtual_marker_positions = virtual_marker_positions;
%             end
            
            % update marker positions
            for i_marker = 1 : this.number_of_markers
                set ...
                  ( ...
                    this.marker_plots{i_marker}, ...
                    'XData', this.marker_positions((i_marker-1)*3 + 1), ...
                    'YData', this.marker_positions((i_marker-1)*3 + 2), ...
                    'ZData', this.marker_positions((i_marker-1)*3 + 3) ...
                  )
            end
            
            % update line positions
            for i_line = 1 : length(this.line_plots)
                from_marker_index = this.line_marker_indices(i_line, 1);
                to_marker_index = this.line_marker_indices(i_line, 2);
                
                set ...
                  ( ...
                    this.line_plots(i_line), ...
                    'Xdata', [this.marker_positions((from_marker_index-1)*3 + 1), this.marker_positions((to_marker_index-1)*3 + 1)], ...
                    'Ydata', [this.marker_positions((from_marker_index-1)*3 + 2), this.marker_positions((to_marker_index-1)*3 + 2)], ...
                    'Zdata', [this.marker_positions((from_marker_index-1)*3 + 3), this.marker_positions((to_marker_index-1)*3 + 3)] ...
                  )
            end
            
            
        end
        function setMiscellaneousPlotColor(this, index, color)
%             if index > length(this.miscellaneousPlots)
%                 error('index larger than number of miscellaneous plots')
%             else
%                 this.miscellaneousPlots(index) = plot3(this.scene_axes, [0 1], [0 1], [0 0], 'color', color, 'Linewidth', 1, 'Linestyle', '-');
%                 this.update();
%             end
        end
        
        function setColors(this, marker_set)
            if strcmp(marker_set, 'default')
                warning('marker set "default" not implemented yet');
            elseif strcmp(marker_set, 'extended plug-in gait')
                red = [1 0 0];
                green = [0 1 0];
                blue = [0 0 1];
                
                marker_marker_style = '*';
                marker_marker_size = 3;
                
                joint_center_marker_style = 'd';
                joint_center_marker_size = 9;
                joint_center_color = 'c';
                
                com_marker_style = 'o';
                com_marker_size = 8;
                com_color = lightenColor([1 0.5 0], 0.5);

                body_com_marker_style = 'o';
                body_com_marker_size = 12;
                body_com_color = [1 0.5 0];
                
                % label, markercolor, markerfacecolor, style, size
                marker_candidates = ...
                  { ...
                    'LFHD', red, red, marker_marker_style, marker_marker_size; ...
                    'LBHD', red, red, marker_marker_style, marker_marker_size; ...
                    'LFHD', red, red, marker_marker_style, marker_marker_size; ...
                    'RFHD', green, green, marker_marker_style, marker_marker_size; ...
                    'LBHD', red, red, marker_marker_style, marker_marker_size; ...
                    'RBHD', green, green, marker_marker_style, marker_marker_size; ...
                    'C7', blue, blue, marker_marker_style, marker_marker_size; ...
                    'T10', blue, blue, marker_marker_style, marker_marker_size; ...
                    'RBAK', blue, blue, marker_marker_style, marker_marker_size; ...
                    'CLAV', blue, blue, marker_marker_style, marker_marker_size; ...
                    'STRN', blue, blue, marker_marker_style, marker_marker_size; ...
                    'LSHO', red, red, marker_marker_style, marker_marker_size; ...
                    'RSHO', green, green, marker_marker_style, marker_marker_size; ...
                    'LUPA', red, red, marker_marker_style, marker_marker_size; ...
                    'LELB', red, red, marker_marker_style, marker_marker_size; ...
                    'LFRM', red, red, marker_marker_style, marker_marker_size; ...
                    'LFRA', red, red, marker_marker_style, marker_marker_size; ...
                    'LWRA', red, red, marker_marker_style, marker_marker_size; ...
                    'LWRB', red, red, marker_marker_style, marker_marker_size; ...
                    'LFIN', red, red, marker_marker_style, marker_marker_size; ...
                    'RUPA', green, green, marker_marker_style, marker_marker_size; ...
                    'RELB', green, green, marker_marker_style, marker_marker_size; ...
                    'RFRM', green, green, marker_marker_style, marker_marker_size; ...
                    'RFRA', green, green, marker_marker_style, marker_marker_size; ...
                    'RWRA', green, green, marker_marker_style, marker_marker_size; ...
                    'RWRB', green, green, marker_marker_style, marker_marker_size; ...
                    'RFIN', green, green, marker_marker_style, marker_marker_size; ...
                    'LASI', red, red, marker_marker_style, marker_marker_size; ...
                    'RASI', green, green, marker_marker_style, marker_marker_size; ...
                    'LPSI', red, red, marker_marker_style, marker_marker_size; ...
                    'RPSI', green, green, marker_marker_style, marker_marker_size; ...
                    'LTHI', red, red, marker_marker_style, marker_marker_size; ...
                    'LTHIA', red, red, marker_marker_style, marker_marker_size; ...
                    'LKNE', red, red, marker_marker_style, marker_marker_size; ...
                    'LTIB', red, red, marker_marker_style, marker_marker_size; ...
                    'LTIBA', red, red, marker_marker_style, marker_marker_size; ...
                    'LANK', red, red, marker_marker_style, marker_marker_size; ...
                    'LHEE', red, red, marker_marker_style, marker_marker_size; ...
                    'LTOE', red, red, marker_marker_style, marker_marker_size; ...
                    'LTOEL', red, red, marker_marker_style, marker_marker_size; ...
                    'RTHI', green, green, marker_marker_style, marker_marker_size; ...
                    'RTHIA', green, green, marker_marker_style, marker_marker_size; ...
                    'RKNE', green, green, marker_marker_style, marker_marker_size; ...
                    'RTIB', green, green, marker_marker_style, marker_marker_size; ...
                    'RTIBA', green, green, marker_marker_style, marker_marker_size; ...
                    'RANK', green, green, marker_marker_style, marker_marker_size; ...
                    'RHEE', green, green, marker_marker_style, marker_marker_size; ...
                    'RTOE', green, green, marker_marker_style, marker_marker_size; ...
                    'RTOEL', green, green, marker_marker_style, marker_marker_size; ...
                    'CERVIXCOR', joint_center_color, joint_center_color, joint_center_marker_style, joint_center_marker_size; ...
                    'LSHOULDERCOR', joint_center_color, joint_center_color, joint_center_marker_style, joint_center_marker_size; ...
                    'RSHOULDERCOR', joint_center_color, joint_center_color, joint_center_marker_style, joint_center_marker_size; ...
                    'LELBOWCOR', joint_center_color, joint_center_color, joint_center_marker_style, joint_center_marker_size; ...
                    'RELBOWCOR', joint_center_color, joint_center_color, joint_center_marker_style, joint_center_marker_size; ...
                    'LWRISTCOR', joint_center_color, joint_center_color, joint_center_marker_style, joint_center_marker_size; ...
                    'RWRISTCOR', joint_center_color, joint_center_color, joint_center_marker_style, joint_center_marker_size; ...
                    'LUMBARCOR', joint_center_color, joint_center_color, joint_center_marker_style, joint_center_marker_size; ...
                    'LHIPCOR', joint_center_color, joint_center_color, joint_center_marker_style, joint_center_marker_size; ...
                    'RHIPCOR', joint_center_color, joint_center_color, joint_center_marker_style, joint_center_marker_size; ...
                    'LKNEECOR', joint_center_color, joint_center_color, joint_center_marker_style, joint_center_marker_size; ...
                    'RKNEECOR', joint_center_color, joint_center_color, joint_center_marker_style, joint_center_marker_size; ...
                    'LANKLECOR', joint_center_color, joint_center_color, joint_center_marker_style, joint_center_marker_size; ...
                    'RANKLECOR', joint_center_color, joint_center_color, joint_center_marker_style, joint_center_marker_size; ...
                    'HEADCOM', com_color, com_color, com_marker_style, com_marker_size; ...
                    'TORSOCOM', com_color, com_color, com_marker_style, com_marker_size; ...
                    'LUPPERARMCOM', com_color, com_color, com_marker_style, com_marker_size; ...
                    'RUPPERARMCOM', com_color, com_color, com_marker_style, com_marker_size; ...
                    'LFOREARMCOM', com_color, com_color, com_marker_style, com_marker_size; ...
                    'RFOREARMCOM', com_color, com_color, com_marker_style, com_marker_size; ...
                    'LHANDCOM', com_color, com_color, com_marker_style, com_marker_size; ...
                    'RHANDCOM', com_color, com_color, com_marker_style, com_marker_size; ...
                    'PELVISCOM', com_color, com_color, com_marker_style, com_marker_size; ...
                    'LTHIGHCOM', com_color, com_color, com_marker_style, com_marker_size; ...
                    'RTHIGHCOM', com_color, com_color, com_marker_style, com_marker_size; ...
                    'LSHANKCOM', com_color, com_color, com_marker_style, com_marker_size; ...
                    'RSHANKCOM', com_color, com_color, com_marker_style, com_marker_size; ...
                    'LFOOTCOM', com_color, com_color, com_marker_style, com_marker_size; ...
                    'RFOOTCOM', com_color, com_color, com_marker_style, com_marker_size; ...
                    'BODYCOM', body_com_color, body_com_color, body_com_marker_style, body_com_marker_size; ...
                  };
              
                for i_label = 1 : size(marker_candidates, 1)
                    try
                        set ...
                          ( ...
                            this.marker_plots{find(strcmp(this.marker_labels, marker_candidates{i_label, 1}))}, ...
                            'color', marker_candidates{i_label, 2}, ...
                            'markerfacecolor', marker_candidates{i_label, 3}, ...
                            'marker', marker_candidates{i_label, 4}, ...
                            'markersize', marker_candidates{i_label, 5} ...
                          );
                    catch exception
                        % label not found, but we don't care
                    end
                end
              
              
                

%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'CERVIXCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'LSHOULDERCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'RSHOULDERCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'LELBOWCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'RELBOWCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'LWRISTCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'RWRISTCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'LUMBARCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'LHIPCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'RHIPCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'LKNEECOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'RKNEECOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'LANKLECOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'RANKLECOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
%     
%                 % segment CoMs
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'HEADCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'TORSOCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'LUPPERARMCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'RUPPERARMCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'LFOREARMCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'RFOREARMCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'LHANDCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'RHANDCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'PELVISCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'LTHIGHCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'RTHIGHCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'LSHANKCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'RSHANKCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'LFOOTCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'RFOOTCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
%                 
%                 % body CoM
%                 set(this.marker_plots{find(strcmp(this.marker_labels, 'BODYCOM'))}, 'color', com_color, 'marker', com_marker_style, 'markersize', round(com_marker_size*1.5), 'markerfacecolor', com_color);
            else
                error('color scheme not recognized, options are "default", "extended plug-in gait".')
            end
        end
        function addLines(this, marker_set)
            if strcmp(marker_set, 'default')
                warning('marker set "default" not implemented yet');
            elseif strcmp(marker_set, 'extended plug-in gait')
                marker_line_width = 1;
                marker_line_color = [1, 1, 1] * 0.7;
                bone_line_width = 5;
                bone_line_color = lightenColor([0, 0, 1], 0.3);
                
                
                % from-label, to_label, line_width, color
                line_candidates = ...
                  { ...
                    'LFHD', 'RFHD', marker_line_width, marker_line_color; ...
                    'LFHD', 'RFHD', marker_line_width, marker_line_color; ...
                    'RFHD', 'RBHD', marker_line_width, marker_line_color; ...
                    'RBHD', 'LBHD', marker_line_width, marker_line_color; ...
                    'LBHD', 'LFHD', marker_line_width, marker_line_color; ...
                    'C7', 'T10', marker_line_width, marker_line_color; ...
                    'T10', 'STRN', marker_line_width, marker_line_color; ...
                    'STRN', 'CLAV', marker_line_width, marker_line_color; ...
                    'CLAV', 'C7', marker_line_width, marker_line_color; ...
                    'C7', 'LSHO', marker_line_width, marker_line_color; ...
                    'CLAV', 'LSHO', marker_line_width, marker_line_color; ...
                    'LSHO', 'LELB', marker_line_width, marker_line_color; ...
                    'LELB', 'LWRA', marker_line_width, marker_line_color; ...
                    'LELB', 'LWRB', marker_line_width, marker_line_color; ...
                    'LWRA', 'LWRB', marker_line_width, marker_line_color; ...
                    'LWRA', 'LFIN', marker_line_width, marker_line_color; ...
                    'LWRB', 'LFIN', marker_line_width, marker_line_color; ...
                    'C7', 'RSHO', marker_line_width, marker_line_color; ...
                    'CLAV', 'RSHO', marker_line_width, marker_line_color; ...
                    'RSHO', 'RELB', marker_line_width, marker_line_color; ...
                    'RELB', 'RWRA', marker_line_width, marker_line_color; ...
                    'RELB', 'RWRB', marker_line_width, marker_line_color; ...
                    'RWRA', 'RWRB', marker_line_width, marker_line_color; ...
                    'RWRA', 'RFIN', marker_line_width, marker_line_color; ...
                    'RWRB', 'RFIN', marker_line_width, marker_line_color; ...
                    'LASI', 'RASI', marker_line_width, marker_line_color; ...
                    'RASI', 'RPSI', marker_line_width, marker_line_color; ...
                    'RPSI', 'LPSI', marker_line_width, marker_line_color; ...
                    'LPSI', 'LASI', marker_line_width, marker_line_color; ...
                    'LASI', 'LKNE', marker_line_width, marker_line_color; ...
                    'LPSI', 'LKNE', marker_line_width, marker_line_color; ...
                    'LKNE', 'LANK', marker_line_width, marker_line_color; ...
                    'LKNE', 'LHEE', marker_line_width, marker_line_color; ...
                    'LANK', 'LHEE', marker_line_width, marker_line_color; ...
                    'LANK', 'LTOE', marker_line_width, marker_line_color; ...
                    'LTOE', 'LTOEL', marker_line_width, marker_line_color; ...
                    'LANK', 'LTOEL', marker_line_width, marker_line_color; ...
                    'LHEE', 'LTOE', marker_line_width, marker_line_color; ...
                    'RASI', 'RKNE', marker_line_width, marker_line_color; ...
                    'RPSI', 'RKNE', marker_line_width, marker_line_color; ...
                    'RKNE', 'RANK', marker_line_width, marker_line_color; ...
                    'RKNE', 'RHEE', marker_line_width, marker_line_color; ...
                    'RANK', 'RHEE', marker_line_width, marker_line_color; ...
                    'RANK', 'RTOE', marker_line_width, marker_line_color; ...
                    'RTOE', 'RTOEL', marker_line_width, marker_line_color; ...
                    'RANK', 'RTOEL', marker_line_width, marker_line_color; ...
                    'RHEE', 'RTOE', marker_line_width, marker_line_color; ...
                    'CERVIXCOR', 'LUMBARCOR', bone_line_width, bone_line_color; ... % spine
                    'CERVIXCOR', 'LSHOULDERCOR', bone_line_width, bone_line_color; ... % left shoulder
                    'CERVIXCOR', 'RSHOULDERCOR', bone_line_width, bone_line_color; ... % right shoulder
                    'LUMBARCOR', 'LSHOULDERCOR', bone_line_width, bone_line_color; ... % left torso
                    'LUMBARCOR', 'RSHOULDERCOR', bone_line_width, bone_line_color; ... % right torso
                    'LSHOULDERCOR', 'LELBOWCOR', bone_line_width, bone_line_color; ... % left upper arm
                    'RSHOULDERCOR', 'RELBOWCOR', bone_line_width, bone_line_color; ... % right upper arm
                    'LELBOWCOR', 'LWRISTCOR', bone_line_width, bone_line_color; ... % left lower arm
                    'RELBOWCOR', 'RWRISTCOR', bone_line_width, bone_line_color; ... % right lower arm
                    'LUMBARCOR', 'LHIPCOR', bone_line_width, bone_line_color; ... % left hip
                    'LUMBARCOR', 'RHIPCOR', bone_line_width, bone_line_color; ... % right hip
                    'LHIPCOR', 'LKNEECOR', bone_line_width, bone_line_color; ... % left thigh
                    'RHIPCOR', 'RKNEECOR', bone_line_width, bone_line_color; ... % right thigh
                    'LKNEECOR', 'LANKLECOR', bone_line_width, bone_line_color; ... % left shank
                    'RKNEECOR', 'RANKLECOR', bone_line_width, bone_line_color; ... % right shank
                  };

                for i_label = 1 : size(line_candidates, 1)
                    try
                        this.addLine ...
                          ( ...
                            line_candidates{i_label, 1}, ...
                            line_candidates{i_label, 2}, ...
                            line_candidates{i_label, 3}, ...
                            line_candidates{i_label, 4} ...
                          );
                    catch exception
                        % label not found, but we don't care
                    end
                end
                
                
                this.update();
            end
        end
        function addLine(this, from_label, to_label, width, color)
            if nargin < 5
                color = [0, 0, 0];
            end
            if nargin < 4
                width = 1;
            end
            label_indices = [find(strcmp(this.marker_labels, from_label)), find(strcmp(this.marker_labels, to_label))];
            if size(label_indices, 1) == 1 && size(label_indices, 2) == 2
                this.line_marker_indices(end+1, :) = [find(strcmp(this.marker_labels, from_label)), find(strcmp(this.marker_labels, to_label))];
                this.line_plots(end+1) = plot3([0 0], [0 0], [0 0], 'color', color, 'linewidth', width);
            end
        end
        function setMarkerColor(this, marker_index, color)
            set(this.marker_plots{marker_index}, 'color', color);
        end
        function setMarkerStyle(this, marker_index, style)
            set(this.marker_plots{marker_index}, 'marker', style);
        end
        function color = getMarkerColor(this, joint_index, marker_index)
            if joint_index == 0
                color = get(this.fixedMarkerPlots(marker_index), 'color');
            else
                color = get(this.marker_plots{joint_index}(marker_index), 'color');
            end
        end
        
        function setLinkPlotsColor(this, color)
            for i_joint = 1 : this.kinematicTree.numberOfJoints
                this.linkPlotsColor = color;
                set(this.linkPlots(i_joint), 'color', color);
            for i_eef = 1 : this.kinematicTree.numberOfBranches
                set(this.endEffectorLinkPlots(i_eef), 'color', color);
            end
        end
        end
        function setLinkPlotsLinewidth(this, linewidth)
            for i_joint = 1 : this.kinematicTree.numberOfJoints
                this.linkPlotsLinewidth = linewidth;
                set(this.linkPlots(i_joint), 'linewidth', linewidth);
            end
            for i_eef = 1 : this.kinematicTree.numberOfBranches
                set(this.endEffectorLinkPlots(i_eef), 'linewidth', linewidth);
            end
        
        end
    end
end