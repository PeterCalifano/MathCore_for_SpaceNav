classdef testCInterpolator < matlab.unittest.TestCase
        properties (Access = private)
        charSavedPath   % char vector of original path
        charPathFile
    end

    methods (TestClassSetup)
        function setupPath(testCase)
            % Save current path (memory + file, if you want a record)
            testCase.charSavedPath = path;
            testCase.charPathFile  = fullfile(fileparts(mfilename("fullpath")), 'path_before_tests.mat');
            charOrigPath = testCase.charSavedPath;
            save(testCase.charPathFile, 'charOrigPath');

            restoredefaultpath;

            % Add the folder you want to test (plus subs if needed)
            charSrcPath = "/home/peterc/devDir/MathCore_for_SpaceNav";
            addpath(genpath(charSrcPath), '-begin');             % or addpath(genpath(p), '-begin');
            testCase.addTeardown(@() teardownPath(testCase));
        end
    end

    methods (Access = private)
        function teardownPath(testCase)
            % Restore original path
            if ~isempty(testCase.charSavedPath)
                path(testCase.charSavedPath);
            end

            % Clean up temp file
            if ~isempty(testCase.charPathFile) && isfile(testCase.charPathFile)
                delete(testCase.charPathFile);
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TEST
    methods (Test)
        function testConstructorValid(test)
            dInterpDomain = linspace(0,1,10);
            ui8PolyDeg = uint8(3);
            enumType = EnumInterpType.VECTOR;

            % Provide explicit bounds and output size
            objHelper = CInterpolatorTestHelper(dInterpDomain, ui8PolyDeg, enumType, true, [0,1], int32(3));
            
            test.verifyEqual(objHelper.enumInterpType, enumType);
            test.verifyEqual(objHelper.ui8PolyDeg, ui8PolyDeg);
            test.verifyEqual(objHelper.dDomainBounds, [0,1]);
            
            % Buffer should be empty initially
            test.verifyEmpty(objHelper.dInterpCoeffsBuffer);
        end
        
        function testConstructorInvalidDegree(test)

            dInterpDomain = 1:5;
            ui8PolyDeg = uint8(1);
            
            % Degree <= 1 should assert
            test.verifyError(@() CInterpolatorTestHelper(dInterpDomain, ui8PolyDeg, ...
                    EnumInterpType.VECTOR, true, [0,1], int32(4)), 'MATLAB:assert:failed');
        end
        
        function testFixQuatNoSwitch(test)

            % Create constant quaternion data (no sign changes)
            dQuat = ones(5,4);
            objHelper = CInterpolatorTestHelper(1:5, uint8(2), EnumInterpType.QUAT, true, [1,5], int32(4));
            
            [dModified, dSwitchIntervals, bIsSwitched, ui8SwitchCount] = objHelper.testFixQuatSignDiscontinuity(dQuat);
            
            % Should simply transpose without modification
            test.verifyEqual(dModified, dQuat');
            
            % No switches detected
            test.verifyEqual(ui8SwitchCount, uint8(0));
            test.verifyEmpty(dSwitchIntervals);
            test.verifyFalse(any(bIsSwitched));
        end
    end
end


