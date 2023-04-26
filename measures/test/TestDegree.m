classdef TestDegree < matlab.unittest.TestCase
    methods (Test)
        function testCalculateBinary(testCase)
            % Create a Degree object.
            degreeObj = Degree();
            
            % Define a binary adjacency matrix.
            adjMatrix = [0 1 1;
                         1 0 0;
                         1 0 0];
            
            % Expected degree values for each node.
            expectedValues = [2; 1; 1];
            
            % Call the calculateBinary method.
            calculatedValues = degreeObj.calculateBinary(adjMatrix);
            
            % Verify that the expected and calculated degree values match.
            testCase.verifyEqual(calculatedValues, expectedValues);
        end
        
        function testCalculateWeighted(testCase)
            % Create a Degree object.
            degreeObj = Degree();
            
            % Define a weighted adjacency matrix.
            adjMatrix = [0 2 3;
                         2 0 0;
                         3 0 0];
            
            % Expected degree values for each node.
            expectedValues = [5; 2; 3];
            
            % Call the calculateWeighted method.
            calculatedValues = degreeObj.calculateWeighted(adjMatrix);
            
            % Verify that the expected and calculated degree values match.
            testCase.verifyEqual(calculatedValues, expectedValues);
        end
    end
end