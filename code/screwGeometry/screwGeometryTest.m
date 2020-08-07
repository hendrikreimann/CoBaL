%     This file is part of the ScrewGeometry library
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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


classdef screwGeometryTest < matlab.unittest.TestCase
    % Unit Test for screwGeometry
    %   to test, instantiate and call run
    
    properties
    end
    
    methods (Test)
        function testWedgeAxis(testCase)
            axis = [3; 5; 7];
            actualSolution = wedgeAxis(axis);
            expectedSolution = [0 -7 5; 7 0 -3; -5 3 0];
            testCase.verifyEqual(actualSolution, expectedSolution);
        end
        function testVeeAxis(testCase)
            matrix = [0 -7 5; 7 0 -3; -5 3 0];
            actualSolution = veeAxis(matrix);
            expectedSolution = [3; 5; 7];
            testCase.verifyEqual(actualSolution, expectedSolution);
        end
        function testWedgeTwist(testCase)
            twistCoordinates = [2; 3; 5; 7; 11; 13];
            actualSolution = wedgeTwist(twistCoordinates);
            expectedSolution = [0 -13 11 2; 13 0 -7 3; -11 7 0 5; 0 0 0 0];
            testCase.verifyEqual(actualSolution, expectedSolution);
        end
        function testVeeTwist(testCase)
            twistMatrix = [0 -13 11 2; 13 0 -7 3; -11 7 0 5; 0 0 0 0];
            actualSolution = veeTwist(twistMatrix);
            expectedSolution = [2; 3; 5; 7; 11; 13];
            testCase.verifyEqual(actualSolution, expectedSolution);
        end
        function testExpAxis(testCase)
            axis = [1; 2; 3];
            axisNormed = axis * 1 / norm(axis);
            angle = 1;
            actualSolution = expAxis(axisNormed, angle);
            expectedSolution = [0.57313785544898699 -0.60900664213739331 0.54829180960859991; 0.74034884046078198 0.67164450419152844 -0.027879282947946255; -0.35127851212351696 0.42190587791811218 0.83582225209576411];
            testCase.verifyEqual(actualSolution, expectedSolution, 'AbsTol', sqrt(eps));
        end
        function testLogAxis(testCase)
            matrix = [0.57313785544898699 -0.60900664213739331 0.54829180960859991; 0.74034884046078198 0.67164450419152844 -0.027879282947946255; -0.35127851212351696 0.42190587791811218 0.83582225209576411];
            [actualSolutionAxis, actualSolutionAngle] = logAxis(matrix);
            expectedSolutionAxis = [1; 2; 3] * 1 / norm([1; 2; 3]);
            expectedSolutionAngle = 1;
            testCase.verifyEqual(actualSolutionAxis, expectedSolutionAxis, 'AbsTol', sqrt(eps));
            testCase.verifyEqual(actualSolutionAngle, expectedSolutionAngle, 'AbsTol', sqrt(eps));
        end
        function testExpTwist(testCase)
            axis = [1; 2; 3] * 1 / norm([1; 2; 3]);
            twistCoordinates = [-axis(2)+axis(3); -axis(3)+axis(1); -axis(1)+axis(2); axis(1); axis(2); axis(3);];
            actualSolution = expTwist(twistCoordinates, pi);
            expectedSolution = ...
              [ ...
                -0.85714285714285676 0.28571428571428564 0.4285714285714286 1.1428571428571423; ...
                0.28571428571428586 -0.42857142857142816 0.8571428571428571 0.28571428571428537; ...
                0.42857142857142849 0.8571428571428571 0.28571428571428559 -0.57142857142857117; ...
                0 0 0 1 ...
              ];
            testCase.verifyEqual(actualSolution, expectedSolution, 'AbsTol', sqrt(eps));
                
              
        end
        function testLogTwist(testCase)
            % normal case
            transformation = [1 0 0 0; 0 0.540302305868140 -0.841470984807897 3.443808342687410; 0 0.841470984807897 0.540302305868140 -0.303848887220212;  0 0 0 1];
            [actualSolutionTwist, actualSolutionAngle] = logTwist(transformation);
            expectedSolutionTwist = [0;3;-2;1;0;0];
            expectedSolutionAngle = 1;
            testCase.verifyEqual(actualSolutionTwist, expectedSolutionTwist, 'AbsTol', eps^(1/2));
            testCase.verifyEqual(actualSolutionAngle, expectedSolutionAngle, 'AbsTol', eps^(1/4));
            
            % border case
            transformation = ...
              [ ...
                -0.85714285714285676 0.28571428571428564 0.4285714285714286 1.1428571428571423; ...
                0.28571428571428586 -0.42857142857142816 0.8571428571428571 0.28571428571428537; ...
                0.42857142857142849 0.8571428571428571 0.28571428571428559 -0.57142857142857117; ...
                0 0 0 1 ...
              ];
            [actualSolutionTwist, actualSolutionAngle] = logTwist(transformation);
            axis = [1; 2; 3] * 1 / norm([1; 2; 3]);
            expectedSolutionTwist = [-axis(2)+axis(3); -axis(3)+axis(1); -axis(1)+axis(2); axis(1); axis(2); axis(3);];
            expectedSolutionAngle = pi;
            testCase.verifyEqual(actualSolutionTwist, expectedSolutionTwist, 'AbsTol', eps^(1/2));
            testCase.verifyEqual(actualSolutionAngle, expectedSolutionAngle, 'AbsTol', eps^(1/4));
        end
        function testRigidToAdjointTransformation(testCase)
            transformation = ...
              [ ...
                -0.85714285714285676 0.28571428571428564 0.4285714285714286 1.1428571428571423; ...
                0.28571428571428586 -0.42857142857142816 0.8571428571428571 0.28571428571428537; ...
                0.42857142857142849 0.8571428571428571 0.28571428571428559 -0.57142857142857117; ...
                0 0 0 1 ...
              ];
            expectedSolution = ...
              [ ...
                -0.85714285714285676 0.28571428571428564 0.4285714285714286 0.28571428571428559 0 0.57142857142857106; ...
                0.28571428571428586 -0.42857142857142816 0.8571428571428571 0 -1.1428571428571423 -0.57142857142857106; ...
                0.42857142857142849 0.8571428571428571 0.28571428571428559 0.57142857142857106 -0.57142857142857062 0.85714285714285676; ...
                0 0 0 -0.85714285714285676 0.28571428571428564 0.4285714285714286; ...
                0 0 0 0.28571428571428586 -0.42857142857142816 0.8571428571428571; ...
                0 0 0 0.42857142857142849 0.8571428571428571 0.28571428571428559; ...
              ];
            actualSolution = rigidToAdjointTransformation(transformation);
            testCase.verifyEqual(actualSolution, expectedSolution, 'AbsTol', sqrt(eps));
            

        end
        function testAdjointToRigidTransformation(testCase)
             adjoint = ...
              [ ...
                -0.85714285714285676 0.28571428571428564 0.4285714285714286 0.28571428571428559 0 0.57142857142857106; ...
                0.28571428571428586 -0.42857142857142816 0.8571428571428571 0 -1.1428571428571423 -0.57142857142857106; ...
                0.42857142857142849 0.8571428571428571 0.28571428571428559 0.57142857142857106 -0.57142857142857062 0.85714285714285676; ...
                0 0 0 -0.85714285714285676 0.28571428571428564 0.4285714285714286; ...
                0 0 0 0.28571428571428586 -0.42857142857142816 0.8571428571428571; ...
                0 0 0 0.42857142857142849 0.8571428571428571 0.28571428571428559; ...
              ];
            expectedSolution = ...
              [ ...
                -0.85714285714285676 0.28571428571428564 0.4285714285714286 1.1428571428571423; ...
                0.28571428571428586 -0.42857142857142816 0.8571428571428571 0.28571428571428537; ...
                0.42857142857142849 0.8571428571428571 0.28571428571428559 -0.57142857142857117; ...
                0 0 0 1 ...
              ];
            actualSolution = adjointToRigidTransformation(adjoint);
            testCase.verifyEqual(actualSolution, expectedSolution, 'AbsTol', sqrt(eps));
        end
        function testInvertAdjointTransformation(testCase)
             adjoint = ...
              [ ...
                -0.85714285714285676 0.28571428571428564 0.4285714285714286 0.28571428571428559 0 0.57142857142857106; ...
                0.28571428571428586 -0.42857142857142816 0.8571428571428571 0 -1.1428571428571423 -0.57142857142857106; ...
                0.42857142857142849 0.8571428571428571 0.28571428571428559 0.57142857142857106 -0.57142857142857062 0.85714285714285676; ...
                0 0 0 -0.85714285714285676 0.28571428571428564 0.4285714285714286; ...
                0 0 0 0.28571428571428586 -0.42857142857142816 0.8571428571428571; ...
                0 0 0 0.42857142857142849 0.8571428571428571 0.28571428571428559; ...
              ];
            expectedSolution = inv(adjoint);
            actualSolution = invertAdjointTransformation(adjoint);
            testCase.verifyEqual(actualSolution, expectedSolution, 'AbsTol', sqrt(eps));
        end
        function testGenerateTwistCoordinates(testCase)
            supportPoint = [0; 0; 1];
            axis = [1; 0; 0];
            actualSolution = generateTwistCoordinates(supportPoint, axis);
            expectedSolution = [0; 1; 0; 1; 0; 0];
            testCase.verifyEqual(actualSolution, expectedSolution);
            
        end
    end
    
end

