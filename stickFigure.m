% class for the display of a stick figure

classdef stickFigure < handle
    properties
        scene_figure;
        scene_axes;
        scene_bound;
        
        % data
        number_of_markers;
        number_of_virtual_markers;
        segment_com_plots;
        marker_plots;
        virtual_marker_plots;
        marker_labels;
        virtual_marker_labels;
        
        line_plots;
        line_marker_indices = [];
        
        virtual_line_plots;
        virtual_line_marker_indices = [];
        
        marker_positions = [];
        virtual_marker_positions = [];
        
        % colors
        marker_color = [1 1 1]*0.5;
        
        % flags
        show_segment_coms = false;
        show_connection_lines = true;
        show_markers = true;
        show_marker_labels = false;
    end 
    methods
        function this = stickFigure(markerPositions, markerLabels, virtualMarkerPositions, virtualMarkerLabels, scene_bound, axesHandle)
            if nargin < 5
                this.scene_bound = [-1 1; -1 1; -1 1];
            else
                this.scene_bound = scene_bound;
            end
            if nargin < 6
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
            this.number_of_virtual_markers = length(virtualMarkerLabels);
            this.virtual_marker_labels = virtualMarkerLabels;
            this.virtual_marker_plots = cell(this.number_of_virtual_markers, 1);
            for i_marker = 1 : this.number_of_virtual_markers
                this.virtual_marker_plots{i_marker} = plot3(0, 0, 0, 'color', this.marker_color, 'markersize', 10, 'linewidth', 3, 'Marker', 'x');
            end
            
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
            
            this.update(markerPositions, virtualMarkerPositions);
        end
        
        function update(this, marker_positions, virtual_marker_positions)
            if nargin > 1
                this.marker_positions = marker_positions;
            end
            if nargin > 2
                this.virtual_marker_positions = virtual_marker_positions;
            end
            
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
            
            % update virtual marker positions
            for i_marker = 1 : this.number_of_virtual_markers
                set ...
                  ( ...
                    this.virtual_marker_plots{i_marker}, ...
                    'XData', this.virtual_marker_positions((i_marker-1)*3 + 1), ...
                    'YData', this.virtual_marker_positions((i_marker-1)*3 + 2), ...
                    'ZData', this.virtual_marker_positions((i_marker-1)*3 + 3) ...
                  )
            end
            
            % update virtual line positions
            for i_line = 1 : length(this.virtual_line_plots)
                from_marker_index = this.virtual_line_marker_indices(i_line, 1);
                to_marker_index = this.virtual_line_marker_indices(i_line, 2);
                
                set ...
                  ( ...
                    this.virtual_line_plots(i_line), ...
                    'Xdata', [this.virtual_marker_positions((from_marker_index-1)*3 + 1), this.virtual_marker_positions((to_marker_index-1)*3 + 1)], ...
                    'Ydata', [this.virtual_marker_positions((from_marker_index-1)*3 + 2), this.virtual_marker_positions((to_marker_index-1)*3 + 2)], ...
                    'Zdata', [this.virtual_marker_positions((from_marker_index-1)*3 + 3), this.virtual_marker_positions((to_marker_index-1)*3 + 3)] ...
                  )
            end
            
            
%             for i_joint = 1 : this.kinematicTree.numberOfJoints
%                 if this.show_connection_lines
%                     set( ...
%                          this.jointPlots(i_joint), ...
%                          'Xdata', this.kinematicTree.jointPositions{i_joint}(1), ...
%                          'Ydata', this.kinematicTree.jointPositions{i_joint}(2), ...
%                          'Zdata', this.kinematicTree.jointPositions{i_joint}(3) ...
%                        )
%                 else
%                     set(this.jointPlots(i_joint), 'visible', 'off');
%                 end
%                 if this.kinematicTree.linkMasses(i_joint) > 0 && this.show_segment_coms
%                     set( ...
%                          this.linkCenterPlots(i_joint), ...
%                          'Xdata', this.kinematicTree.linkTransformations{i_joint}(1, 4), ...
%                          'Ydata', this.kinematicTree.linkTransformations{i_joint}(2, 4), ...
%                          'Zdata', this.kinematicTree.linkTransformations{i_joint}(3, 4) ...
%                        )
%                 else
%                     set(this.linkCenterPlots(i_joint), 'visible', 'off');
%                 end
%                 if this.show_marker_labels
%                     set(this.jointLabels(i_joint), 'position', this.kinematicTree.jointPositions{i_joint}, 'visible', 'on');
%                 else
%                     set(this.jointLabels(i_joint), 'visible', 'off');
% 
%                 end
%             end
%             for i_joint = 2 : this.kinematicTree.numberOfJoints
%                 if this.show_connection_lines
%                     set( ...
%                          this.linkPlots(i_joint), ...
%                          'Xdata', [this.kinematicTree.jointPositions{this.kinematicTree.jointParents(i_joint)}(1), this.kinematicTree.jointPositions{i_joint}(1)], ...
%                          'Ydata', [this.kinematicTree.jointPositions{this.kinematicTree.jointParents(i_joint)}(2), this.kinematicTree.jointPositions{i_joint}(2)], ...
%                          'Zdata', [this.kinematicTree.jointPositions{this.kinematicTree.jointParents(i_joint)}(3), this.kinematicTree.jointPositions{i_joint}(3)] ...
%                        )
%                 else
%                     set(this.linkPlots(i_joint), 'visible', 'off');
%                 end
%             end
            
%             for i_eef = 1 : this.kinematicTree.numberOfBranches
%                 if this.show_connection_lines
%                     set(this.endEffectorPlots(i_eef), ...
%                             'Xdata', this.kinematicTree.endEffectorPositions{i_eef}(1), ...
%                             'Ydata', this.kinematicTree.endEffectorPositions{i_eef}(2), ...
%                             'Zdata', this.kinematicTree.endEffectorPositions{i_eef}(3))
%                     set( ...
%                          this.endEffectorLinkPlots(i_eef), ...
%                          'Xdata', [this.kinematicTree.jointPositions{this.kinematicTree.endEffectorParents(i_eef)}(1), this.kinematicTree.endEffectorPositions{i_eef}(1)], ...
%                          'Ydata', [this.kinematicTree.jointPositions{this.kinematicTree.endEffectorParents(i_eef)}(2), this.kinematicTree.endEffectorPositions{i_eef}(2)], ...
%                          'Zdata', [this.kinematicTree.jointPositions{this.kinematicTree.endEffectorParents(i_eef)}(3), this.kinematicTree.endEffectorPositions{i_eef}(3)] ...
%                        )
% %                     set( ...
% %                          this.endEffectorVelocityPlots(i_eef), ...
% %                          'Xdata', [this.kinematicTree.endEffectorPositions{i_eef}(1), this.kinematicTree.endEffectorPositions{i_eef}(1) + this.kinematicTree.endEffectorVelocities{i_eef}(1)], ...
% %                          'Ydata', [this.kinematicTree.endEffectorPositions{i_eef}(2), this.kinematicTree.endEffectorPositions{i_eef}(2) + this.kinematicTree.endEffectorVelocities{i_eef}(2)], ...
% %                          'Zdata', [this.kinematicTree.endEffectorPositions{i_eef}(3), this.kinematicTree.endEffectorPositions{i_eef}(3) + this.kinematicTree.endEffectorVelocities{i_eef}(3)] ...
% %                        )
%                 else
%                     set(this.endEffectorPlots(i_eef), 'visible', 'off');
%                     set(this.endEffectorLinkPlots(i_eef), 'visible', 'off');
%                     set(this.endEffectorVelocityPlots(i_eef), 'visible', 'off');
%                 end
%             end
            
            % update segment mass shapes
%             for i_joint = 1 : this.kinematicTree.numberOfJoints
%                 if this.showLinkMassEllipsoids && (this.kinematicTree.linkMasses(i_joint) > 0) && ~isempty(this.linkMassShapeData{i_joint})
%                     % transform ellipsoid
%                     link_transformation = this.kinematicTree.linkTransformations{i_joint};
%                     inertia_transformation = [this.linkFrameToLinkInertiaRotations{i_joint} zeros(3, 1); 0 0 0 1];
%                     shape_transformed = zeros(size(this.linkMassShapeData{i_joint}));
%                     mass_shape_data = this.linkMassShapeData{i_joint};
%                     for i_point = 1 : size(mass_shape_data, 2)
%                         for j_point = 1 : size(mass_shape_data, 3)
%                             point = squeeze(mass_shape_data(:, i_point, j_point));
%                             shape_transformed(:, i_point, j_point) = link_transformation * inertia_transformation * point;
%                         end
%                     end
%                     
%                     set( ...
%                          this.linkMassShapeSurfs(i_joint), ...
%                          'Xdata', squeeze(shape_transformed(1, :, :)), ...
%                          'Ydata', squeeze(shape_transformed(2, :, :)), ...
%                          'Zdata', squeeze(shape_transformed(3, :, :)) ...
%                        );
%                 else
%                     set(this.linkMassShapeSurfs(i_joint), 'visible', 'off');
%                 end
%             end
                
            % update marker plots
%             for i_marker = 1 : this.kinematicTree.getNumberOfMarkers(0)
%                 if this.show_markers
%                     position = this.kinematicTree.getMarkerPosition(0, i_marker);
%                     set( ...
%                          this.fixedMarkerPlots(i_marker), ...
%                             'Xdata', position(1), ...
%                             'Ydata', position(2), ...
%                             'Zdata', position(3) ...
%                        )
%                 else
%                     set(this.fixedMarkerPlots(i_marker), 'visible', 'off');
%                 end
%             end
%             for i_joint = 1 : this.kinematicTree.numberOfJoints
%                 for i_marker = 1 : this.kinematicTree.getNumberOfMarkers(i_joint)
%                     if this.show_markers
%                         position = this.kinematicTree.getMarkerPosition(i_joint, i_marker);
%                         set( ...
%                              this.marker_plots{i_joint}(i_marker), ...
%                                 'Xdata', position(1), ...
%                                 'Ydata', position(2), ...
%                                 'Zdata', position(3) ...
%                            )
%                     else
%                         set(this.marker_plots{i_joint}(i_marker), 'visible', 'off');
%                     end
%                 end
%             end
%             for i_line = 1 : length(this.marker_connection_plots)
%                 if this.show_markers
%                     joint_index_marker_one = this.kinematicTree.markerExportMap(1, this.kinematicTree.markerConnectionLineIndices(i_line, 1));
%                     marker_index_marker_one = this.kinematicTree.markerExportMap(2, this.kinematicTree.markerConnectionLineIndices(i_line, 1));
%                     joint_index_marker_two = this.kinematicTree.markerExportMap(1, this.kinematicTree.markerConnectionLineIndices(i_line, 2));
%                     marker_index_marker_two = this.kinematicTree.markerExportMap(2, this.kinematicTree.markerConnectionLineIndices(i_line, 2));
%                     position_marker_one = this.kinematicTree.getMarkerPosition(joint_index_marker_one, marker_index_marker_one);
%                     position_marker_two = this.kinematicTree.getMarkerPosition(joint_index_marker_two, marker_index_marker_two);
%                     set( ...
%                          this.marker_connection_plots(i_line), ...
%                             'Xdata', [position_marker_one(1) position_marker_two(1)], ...
%                             'Ydata', [position_marker_one(2) position_marker_two(2)], ...
%                             'Zdata', [position_marker_one(3) position_marker_two(3)] ...
%                        )
%                 else
%                     set(this.marker_connection_plots(i_line), 'visible', 'off');
%                 end
%             end
            
            
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
                
                % head
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LFHD')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RFHD')), green);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LBHD')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RBHD')), green);
                
                % torso
                this.setMarkerColor(find(strcmp(this.marker_labels, 'C7')), blue);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'T10')), blue);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RBAK')), blue);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'CLAV')), blue);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'STRN')), blue);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LSHO')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RSHO')), green);
                
                % left arm
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LUPA')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LELB')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LFRM')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LWRA')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LWRB')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LFIN')), red);
                
                % right arm
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RUPA')), green);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RELB')), green);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RFRM')), green);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RWRA')), green);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RWRB')), green);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RFIN')), green);
                
                % pelvis
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LASI')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RASI')), green);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LPSI')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RPSI')), green);
                
                % left leg
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LTHI')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LTHIA')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LKNE')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LTIB')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LTIBA')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LANK')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LHEE')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LTOE')), red);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'LTOEL')), red);
                
                % right leg
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RTHI')), green);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RTHIA')), green);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RKNE')), green);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RTIB')), green);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RTIBA')), green);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RANK')), green);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RHEE')), green);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RTOE')), green);
                this.setMarkerColor(find(strcmp(this.marker_labels, 'RTOEL')), green);
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
                virtual_line_width = 5;
                virtual_line_color = lightenColor([0, 0, 1], 0.3);
                
                % head
                this.addLine('LFHD', 'RFHD', marker_line_width, marker_line_color);
                this.addLine('RFHD', 'RBHD', marker_line_width, marker_line_color);
                this.addLine('RBHD', 'LBHD', marker_line_width, marker_line_color);
                this.addLine('LBHD', 'LFHD', marker_line_width, marker_line_color);
                
                % torso
                this.addLine('C7', 'T10', marker_line_width, marker_line_color);
                this.addLine('T10', 'STRN', marker_line_width, marker_line_color);
                this.addLine('STRN', 'CLAV', marker_line_width, marker_line_color);
                this.addLine('CLAV', 'C7', marker_line_width, marker_line_color);

                % left arm
                this.addLine('C7', 'LSHO', marker_line_width, marker_line_color);
                this.addLine('CLAV', 'LSHO', marker_line_width, marker_line_color);
                this.addLine('LSHO', 'LELB', marker_line_width, marker_line_color);
                this.addLine('LELB', 'LWRA', marker_line_width, marker_line_color);
                this.addLine('LELB', 'LWRB', marker_line_width, marker_line_color);
                this.addLine('LWRA', 'LWRB', marker_line_width, marker_line_color);
                this.addLine('LWRA', 'LFIN', marker_line_width, marker_line_color);
                this.addLine('LWRB', 'LFIN', marker_line_width, marker_line_color);
                
                % right arm
                this.addLine('C7', 'RSHO', marker_line_width, marker_line_color);
                this.addLine('CLAV', 'RSHO', marker_line_width, marker_line_color);
                this.addLine('RSHO', 'RELB', marker_line_width, marker_line_color);
                this.addLine('RELB', 'RWRA', marker_line_width, marker_line_color);
                this.addLine('RELB', 'RWRB', marker_line_width, marker_line_color);
                this.addLine('RWRA', 'RWRB', marker_line_width, marker_line_color);
                this.addLine('RWRA', 'RFIN', marker_line_width, marker_line_color);
                this.addLine('RWRB', 'RFIN', marker_line_width, marker_line_color);

                % pelvis
                this.addLine('LASI', 'RASI', marker_line_width, marker_line_color);
                this.addLine('RASI', 'RPSI', marker_line_width, marker_line_color);
                this.addLine('RPSI', 'LPSI', marker_line_width, marker_line_color);
                this.addLine('LPSI', 'LASI', marker_line_width, marker_line_color);
                
                % left leg
                this.addLine('LASI', 'LKNE', marker_line_width, marker_line_color);
                this.addLine('LPSI', 'LKNE', marker_line_width, marker_line_color);
                this.addLine('LKNE', 'LANK', marker_line_width, marker_line_color);
                this.addLine('LKNE', 'LHEE', marker_line_width, marker_line_color);
                this.addLine('LANK', 'LHEE', marker_line_width, marker_line_color);
                this.addLine('LANK', 'LTOE', marker_line_width, marker_line_color);
                this.addLine('LTOE', 'LTOEL', marker_line_width, marker_line_color);
                this.addLine('LANK', 'LTOEL', marker_line_width, marker_line_color);
                this.addLine('LHEE', 'LTOE', marker_line_width, marker_line_color);
                
                % right leg
                this.addLine('RASI', 'RKNE', marker_line_width, marker_line_color);
                this.addLine('RPSI', 'RKNE', marker_line_width, marker_line_color);
                this.addLine('RKNE', 'RANK', marker_line_width, marker_line_color);
                this.addLine('RKNE', 'RHEE', marker_line_width, marker_line_color);
                this.addLine('RANK', 'RHEE', marker_line_width, marker_line_color);
                this.addLine('RANK', 'RTOE', marker_line_width, marker_line_color);
                this.addLine('RTOE', 'RTOEL', marker_line_width, marker_line_color);
                this.addLine('RANK', 'RTOEL', marker_line_width, marker_line_color);
                this.addLine('RHEE', 'RTOE', marker_line_width, marker_line_color);
                
                % skeleton
                this.addLine('CERVIXCOR', 'LUMBARCOR', virtual_line_width, virtual_line_color, 'virtual'); % spine
                this.addLine('CERVIXCOR', 'LSHOULDERCOR', virtual_line_width, virtual_line_color, 'virtual'); % left shoulder
                this.addLine('CERVIXCOR', 'RSHOULDERCOR', virtual_line_width, virtual_line_color, 'virtual'); % right shoulder
                this.addLine('LUMBARCOR', 'LSHOULDERCOR', virtual_line_width, virtual_line_color, 'virtual'); % left torso
                this.addLine('LUMBARCOR', 'RSHOULDERCOR', virtual_line_width, virtual_line_color, 'virtual'); % right torso
                this.addLine('LSHOULDERCOR', 'LELBOWCOR', virtual_line_width, virtual_line_color, 'virtual'); % left upper arm
                this.addLine('RSHOULDERCOR', 'RELBOWCOR', virtual_line_width, virtual_line_color, 'virtual'); % right upper arm
                this.addLine('LELBOWCOR', 'LWRISTCOR', virtual_line_width, virtual_line_color, 'virtual'); % left lower arm
                this.addLine('RELBOWCOR', 'RWRISTCOR', virtual_line_width, virtual_line_color, 'virtual'); % right lower arm
                this.addLine('LUMBARCOR', 'LHIPCOR', virtual_line_width, virtual_line_color, 'virtual'); % left hip
                this.addLine('LUMBARCOR', 'RHIPCOR', virtual_line_width, virtual_line_color, 'virtual'); % right hip
                this.addLine('LHIPCOR', 'LKNEECOR', virtual_line_width, virtual_line_color, 'virtual'); % left thigh
                this.addLine('RHIPCOR', 'RKNEECOR', virtual_line_width, virtual_line_color, 'virtual'); % right thigh
                this.addLine('LKNEECOR', 'LANKLECOR', virtual_line_width, virtual_line_color, 'virtual'); % left shank
                this.addLine('RKNEECOR', 'RANKLECOR', virtual_line_width, virtual_line_color, 'virtual'); % right shank
                
                
                this.update();
            end
        end
        function addLine(this, from_label, to_label, width, color, mode)
            if nargin < 6
                mode = 'default';
            end
            if nargin < 5
                color = [0, 0, 0];
            end
            if nargin < 4
                width = 1;
            end
            if strcmp(mode, 'default')
                this.line_plots(end+1) = plot3([0 0], [0 0], [0 0], 'color', color, 'linewidth', width);
                this.line_marker_indices(end+1, :) = [find(strcmp(this.marker_labels, from_label)), find(strcmp(this.marker_labels, to_label))];
            elseif strcmp(mode, 'virtual')
                this.virtual_line_plots(end+1) = plot3([0 0], [0 0], [0 0], 'color', color, 'linewidth', width);
                this.virtual_line_marker_indices(end+1, :) = [find(strcmp(this.virtual_marker_labels, from_label)), find(strcmp(this.virtual_marker_labels, to_label))];
            else
                error('mode must be "default" or "virtual"')
            end
        end
        function setMarkerColor(this, marker_index, color)
            set(this.marker_plots{marker_index}, 'color', color);
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