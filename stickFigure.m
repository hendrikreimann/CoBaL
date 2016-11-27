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
            
            if false
%             % update virtual marker positions
%             for i_marker = 1 : this.number_of_virtual_markers
%                 set ...
%                   ( ...
%                     this.virtual_marker_plots{i_marker}, ...
%                     'XData', this.virtual_marker_positions((i_marker-1)*3 + 1), ...
%                     'YData', this.virtual_marker_positions((i_marker-1)*3 + 2), ...
%                     'ZData', this.virtual_marker_positions((i_marker-1)*3 + 3) ...
%                   )
%             end
%             
%             % update virtual line positions
%             for i_line = 1 : length(this.virtual_line_plots)
%                 from_marker_index = this.virtual_line_marker_indices(i_line, 1);
%                 to_marker_index = this.virtual_line_marker_indices(i_line, 2);
%                 
%                 set ...
%                   ( ...
%                     this.virtual_line_plots(i_line), ...
%                     'Xdata', [this.virtual_marker_positions((from_marker_index-1)*3 + 1), this.virtual_marker_positions((to_marker_index-1)*3 + 1)], ...
%                     'Ydata', [this.virtual_marker_positions((from_marker_index-1)*3 + 2), this.virtual_marker_positions((to_marker_index-1)*3 + 2)], ...
%                     'Zdata', [this.virtual_marker_positions((from_marker_index-1)*3 + 3), this.virtual_marker_positions((to_marker_index-1)*3 + 3)] ...
%                   )
%             end
            
            
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
                com_color = [1 0.5 0];
                
                
                % head
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LFHD'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RFHD'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LBHD'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RBHD'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                
                % torso
                set(this.marker_plots{find(strcmp(this.marker_labels, 'C7'))}, 'color', blue, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'T10'))}, 'color', blue, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RBAK'))}, 'color', blue, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'CLAV'))}, 'color', blue, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'STRN'))}, 'color', blue, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LSHO'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RSHO'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);

                % left arm
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LUPA'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LELB'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LFRM'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LWRA'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LWRB'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LFIN'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);

                % right arm
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RUPA'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RELB'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RFRM'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RWRA'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RWRB'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RFIN'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);

                % pelvis
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LASI'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RASI'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LPSI'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RPSI'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);

                % left leg
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LTHI'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LTHIA'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LKNE'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LTIB'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LTIBA'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LANK'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LHEE'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LTOE'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LTOEL'))}, 'color', red, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                
                % right leg
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RTHI'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RTHIA'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RKNE'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RTIB'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RTIBA'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RANK'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RHEE'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RTOE'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RTOEL'))}, 'color', green, 'marker', marker_marker_style, 'markersize', marker_marker_size);
                
                % joint centers
                set(this.marker_plots{find(strcmp(this.marker_labels, 'CERVIXCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LSHOULDERCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RSHOULDERCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LELBOWCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RELBOWCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LWRISTCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RWRISTCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LUMBARCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LHIPCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RHIPCOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LKNEECOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RKNEECOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LANKLECOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RANKLECOR'))}, 'color', joint_center_color, 'marker', joint_center_marker_style, 'markersize', joint_center_marker_size);
    
                % segment CoMs
                set(this.marker_plots{find(strcmp(this.marker_labels, 'HEADCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
                set(this.marker_plots{find(strcmp(this.marker_labels, 'TORSOCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LUPPERARMCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RUPPERARMCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LFOREARMCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RFOREARMCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LHANDCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RHANDCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
                set(this.marker_plots{find(strcmp(this.marker_labels, 'PELVISCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LTHIGHCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RTHIGHCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LSHANKCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RSHANKCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
                set(this.marker_plots{find(strcmp(this.marker_labels, 'LFOOTCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
                set(this.marker_plots{find(strcmp(this.marker_labels, 'RFOOTCOM'))}, 'color', lightenColor(com_color, 0.5), 'marker', com_marker_style, 'markersize', com_marker_size, 'markerfacecolor', lightenColor(com_color, 0.5));
                
                % body CoM
                set(this.marker_plots{find(strcmp(this.marker_labels, 'BODYCOM'))}, 'color', com_color, 'marker', com_marker_style, 'markersize', round(com_marker_size*1.5), 'markerfacecolor', com_color);
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
                this.addLine('CERVIXCOR', 'LUMBARCOR', bone_line_width, bone_line_color); % spine
                this.addLine('CERVIXCOR', 'LSHOULDERCOR', bone_line_width, bone_line_color); % left shoulder
                this.addLine('CERVIXCOR', 'RSHOULDERCOR', bone_line_width, bone_line_color); % right shoulder
                this.addLine('LUMBARCOR', 'LSHOULDERCOR', bone_line_width, bone_line_color); % left torso
                this.addLine('LUMBARCOR', 'RSHOULDERCOR', bone_line_width, bone_line_color); % right torso
                this.addLine('LSHOULDERCOR', 'LELBOWCOR', bone_line_width, bone_line_color); % left upper arm
                this.addLine('RSHOULDERCOR', 'RELBOWCOR', bone_line_width, bone_line_color); % right upper arm
                this.addLine('LELBOWCOR', 'LWRISTCOR', bone_line_width, bone_line_color); % left lower arm
                this.addLine('RELBOWCOR', 'RWRISTCOR', bone_line_width, bone_line_color); % right lower arm
                this.addLine('LUMBARCOR', 'LHIPCOR', bone_line_width, bone_line_color); % left hip
                this.addLine('LUMBARCOR', 'RHIPCOR', bone_line_width, bone_line_color); % right hip
                this.addLine('LHIPCOR', 'LKNEECOR', bone_line_width, bone_line_color); % left thigh
                this.addLine('RHIPCOR', 'RKNEECOR', bone_line_width, bone_line_color); % right thigh
                this.addLine('LKNEECOR', 'LANKLECOR', bone_line_width, bone_line_color); % left shank
                this.addLine('RKNEECOR', 'RANKLECOR', bone_line_width, bone_line_color); % right shank
                
                
                this.update();
            end
        end
        function addLine(this, from_label, to_label, width, color)
%         function addLine(this, from_label, to_label, width, color, mode)
%             if nargin < 6
%                 mode = 'default';
%             end
            if nargin < 5
                color = [0, 0, 0];
            end
            if nargin < 4
                width = 1;
            end
            this.line_plots(end+1) = plot3([0 0], [0 0], [0 0], 'color', color, 'linewidth', width);
            this.line_marker_indices(end+1, :) = [find(strcmp(this.marker_labels, from_label)), find(strcmp(this.marker_labels, to_label))];
%             if strcmp(mode, 'default')
%                 this.line_plots(end+1) = plot3([0 0], [0 0], [0 0], 'color', color, 'linewidth', width);
%                 this.line_marker_indices(end+1, :) = [find(strcmp(this.marker_labels, from_label)), find(strcmp(this.marker_labels, to_label))];
%             elseif strcmp(mode, 'virtual')
%                 this.virtual_line_plots(end+1) = plot3([0 0], [0 0], [0 0], 'color', color, 'linewidth', width);
%                 this.virtual_line_marker_indices(end+1, :) = [find(strcmp(this.virtual_marker_labels, from_label)), find(strcmp(this.virtual_marker_labels, to_label))];
%             else
%                 error('mode must be "default" or "virtual"')
%             end
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