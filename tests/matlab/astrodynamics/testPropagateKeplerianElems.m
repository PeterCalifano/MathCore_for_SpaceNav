classdef testPropagateKeplerianElems < matlab.unittest.TestCase
    % Test suite for PropagateKeplerianElems function
    methods (Test)
        function testZeroDtReturnsInitial(testCase)
            dx0 = [7000; 0.1; 0; 0; 0; 1.0];
            mu  = 398600.4418;
            dt  = 0;
            out = PropagateKeplerianElems(dx0, mu, dt);
            testCase.verifyEqual(out(:,1), dx0, 'AbsTol', 1e-12);
        end

        function testCircularQuarterPeriod(testCase)
            % Circular orbit: true anomaly = M
            a   = 7000;
            dx0 = [a; 0; 0; 0; 0; 0];
            mu  = 398600.4418;
            n   = sqrt(mu/(a^3));
            % Quarter period dt = (1/4)*(2*pi/n)
            dt  = (pi/2)/n;
            out = PropagateKeplerianElems(dx0, mu, dt);
            expectedNu = mod(n*dt, 2*pi);
            testCase.verifyEqual(out(6,2), expectedNu, 'AbsTol', 1e-9);
        end

        function testEllipticalFullPeriod(testCase)
            % After one full period, anomaly returns to initial
            a   = 10000;
            e   = 0.5;
            dx0 = [a; e; 0; 0; 0; 0];
            mu  = 398600.4418;
            T   = 2*pi*sqrt(a^3/mu);
            out = PropagateKeplerianElems(dx0, mu, T);
            % True anomaly at end should match initial (mod 2*pi)
            nuEnd = out(6,end);
            testCase.verifyEqual(mod(nuEnd,2*pi), 0, 'AbsTol', 1e-6);
        end

        function testHyperbolicZeroDt(testCase)
            dx0 = [-10000; 1.5; 0; 0; 0; 0.5];
            mu  = 398600.4418;
            dt  = 0;
            out = PropagateKeplerianElems(dx0, mu, dt);
            testCase.verifyEqual(out(:,1), dx0, 'AbsTol', 1e-12);
        end

        function testMultiStepOutputSize(testCase)
            dx0 = [7000; 0.1; 0; 0; 0; 0.2];
            mu  = 398600.4418;
            dt  = [100, 200, 300];
            out = PropagateKeplerianElems(dx0, mu, dt);
            testCase.verifySize(out, [6, numel(dt)+1]);
        end

        function testInvalidEccentricityError(testCase)
            % Negative eccentricity should error
            dx0 = [7000; -0.1; 0; 0; 0; 0];
            mu  = 398600.4418;
            testCase.verifyError(@() PropagateKeplerianElems(dx0, mu, 0), 'MATLAB:assertion:failed', '?');
        end
    end
end
