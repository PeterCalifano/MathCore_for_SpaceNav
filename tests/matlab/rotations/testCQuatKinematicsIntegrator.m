classdef testCQuatKinematicsIntegrator < matlab.unittest.TestCase
    % Unit tests for CQuatKinematicsIntegrator

    properties
        Integr  CQuatKinematicsIntegrator
        Tolerance = 1e-10
    end

    methods (TestMethodSetup)
        function createIntegrator(testCase)
            testCase.Integr = CQuatKinematicsIntegrator();
        end
    end

    methods (Test)

        function testNormalizeSeq(testCase)
            % Norm of each column should be 1
            Q = [2 0; 0 3; 0 4; 0 0];
            Qn = CQuatKinematicsIntegrator.NormalizeSeq(Q);
            norms = vecnorm(Qn);
            testCase.verifyLessThanOrEqual(abs(norms - 1), testCase.Tolerance);
        end

        function testExpMapIdentity(testCase)
            % Zero rotation should map to identity quaternion
            dq = [0;0;0];
            q = CQuatKinematicsIntegrator.ExpMap(dq);
            testCase.verifyEqual(q, [1;0;0;0], 'AbsTol', testCase.Tolerance);
        end

        function testExpMapPiRotation(testCase)
            % Rotation of pi about x-axis
            dq = [pi;0;0];
            q = CQuatKinematicsIntegrator.ExpMap(dq);
            expected = [cos(pi/2); sin(pi/2); 0; 0];
            testCase.verifyEqual(q, expected, 'AbsTol', testCase.Tolerance);
        end

        function testZeroOmegaRK4(testCase)
            % Integrating zero angular velocity yields constant quaternion
            q0   = [1;0;0;0];
            omega = [0;0;0];
            tgrid = linspace(0,1,11);

            % Integrate
            [qEnd, qSeq] = testCase.Integr.integrate(q0, omega, tgrid, 'rk4');
            % All outputs should equal the initial quaternion
            testCase.verifyEqual(qSeq, repmat(q0,1,numel(tgrid)), 'AbsTol', testCase.Tolerance);
            testCase.verifyEqual(qEnd, q0, 'AbsTol', testCase.Tolerance);
        end

        function testConstantOmegaRK4(testCase)
            % Constant rotation about Z at pi rad/s for 1 second
            q0   = [1;0;0;0];
            omega = [0;0;pi];
            dt    = 0.1;
            tgrid = 0:dt:1;

            % Integrate
            [qEnd, ~] = testCase.Integr.integrate(q0, omega, tgrid, 'rk4', dt);
            % Final rotation is by pi radians about Z:
            expected = [cos(pi/2); 0;0;sin(pi/2)];
            testCase.verifyEqual(qEnd, expected, 'AbsTol', 1e-3);
        end

        function testConstantOmegaLieGroupEuler(testCase)
            % Constant rotation about Z at pi rad/s for 1 second
            q0   = [1;0;0;0];
            omega = [0;0;pi];
            dt    = 0.1;
            tgrid = 0:dt:1;

            [qEnd, ~] = testCase.Integr.integrate(q0, omega, tgrid, 'lie_euler', dt);
            % Final rotation is by pi radians about Z:
            expected = [cos(pi/2); 0;0;sin(pi/2)];
            testCase.verifyEqual(qEnd, expected, 'AbsTol', 1e-3);
        end

        function testConstantOmegaRKMK4(testCase)
            % Constant rotation about Z at pi rad/s for 1 second
            q0   = [1;0;0;0];
            omega = [0;0;pi];
            dt    = 0.1;
            tgrid = 0:dt:1;

            % Integrate
            [qEnd, ~] = testCase.Integr.integrate(q0, omega, tgrid, 'rkmk4', dt);

            % Final rotation is by pi radians about Z:
            expected = [cos(pi/2); 0;0;sin(pi/2)];
            testCase.verifyEqual(qEnd, expected, 'AbsTol', 1e-3);
        end

        function testLieEulerVsRKMK4(testCase)
            % For small dt, Lie‐Euler and RKMK4 should agree to first order
            q0 = [1;0;0;0];
            omega = [0.1; 0.2; 0.3];
            dt = 1e-3;
            tgrid = [0 dt];

            % Integrate and compare
            [qLE, ~] = testCase.Integr.integrate(q0, omega, tgrid, 'lie_euler', dt);
            [qRK, ~] = testCase.Integr.integrate(q0, omega, tgrid, 'rkmk4',     dt);
            testCase.verifyLessThan(norm(qLE - qRK), 1e-3);
        end

        function testOmegaAngVelProfileRKMK4(testCase)
            % TODO
            error('Not implemented yet')
        end

        function testOmegaFcnHandleRKMK4(testCase)
            
            % Constant rotation about Z at pi rad/s for 1 second
            q0   = [1;0;0;0];
            omega = @(dT) [0;0;pi];
            dt    = 0.1;
            tgrid = 0:dt:1;

            % Integrate
            [qEnd, ~] = testCase.Integr.integrate(q0, omega, tgrid, 'rkmk4', dt);

            % Final rotation is by pi radians about Z:
            expected = [cos(pi/2); 0;0;sin(pi/2)];
            testCase.verifyEqual(qEnd, expected, 'AbsTol', 1e-3);

        end

    end
end
