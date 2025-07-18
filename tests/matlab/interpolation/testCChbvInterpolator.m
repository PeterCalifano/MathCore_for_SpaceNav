classdef testCChbvInterpolator < matlab.unittest.TestCase
    methods (Test)
        function testConstructorScaling(test)

            dInterpDomain = linspace(0,2,5);
            dDomainBounds = [0,2];
            ui8PolyDeg = uint8(3);
            
            self = CChbvInterpolator(dInterpDomain, ...
                                    ui8PolyDeg, ...
                                    EnumInterpType.VECTOR, ...
                                    false, ...
                                    dDomainBounds, ...
                                    int32(1), ...
                                    true);

            % Expected scaled domain: (2*x - (UB+LB)) / (UB-LB)
            expectedScaled = (2.*dInterpDomain - sum(dDomainBounds)) ./ diff(dDomainBounds);
            test.verifyEqual(self.dScaledInterpDomain, expectedScaled, 'AbsTol', 1e-12);
            test.verifyTrue(self.bUSE_MEX);
        
        end
        
        function testEvalPolyChebyshev(test)
            interpDomain = linspace(0,1,5);

            % Degree 3: returns 4 terms
            self = CChbvInterpolator(interpDomain, uint8(3));
            [~, dPoly] = self.evalPoly(0.5);

            test.verifySize(dPoly, [4,1]);

            % Check known Chebyshev values: T0=0, T1=1, T2=2*x*T1-T0=1, T3=2*x*T2-T1=0
            test.verifyEqual(dPoly(1), 0);
            test.verifyEqual(dPoly(2), 1);
            test.verifyEqual(dPoly(3), 1);
            test.verifyEqual(dPoly(4), 0);
        end
        
        function testFitAndEvalInterpolantConstant(test)
            
            % Simple constant data: f(x)=1
            dInterpDomain = [0,0.5,1];
            dDataMatrix = ones(1, numel(dInterpDomain));

            self = CChbvInterpolator(dInterpDomain, uint8(2), EnumInterpType.VECTOR, false);
            
            % Fit data matrix
            [self, dCoeffsT, stats] = self.fitDataMatrix(dDataMatrix);
            
            test.verifyNotEmpty(dCoeffsT);
            test.verifyTrue(isstruct(stats));
            
            % Evaluate at an arbitrary point
            [~, dInterp] = self.evalInterpolant(0.2, true);
            
            % Should reproduce constant value
            test.verifyEqual(dInterp, 1, 'AbsTol', 1e-10);
        end


        function testFitSinusoidal(test)

            % Sinusoidal function f(x)=sin(2*pi*x)
            dInterpDomain = linspace(0,1,50);
            
            dDataMatrix = sin(2*pi*dInterpDomain);
            
            self = CChbvInterpolator(dInterpDomain, uint8(8), EnumInterpType.VECTOR, false, [0,1], int32(1));
            
            % Fit data and retrieve statistics
            [self, dCoeffsT, strFitStats] = self.fitDataMatrix(dDataMatrix);
            
            % Evaluate at test points
            dEvalPoints = linspace(0,1,20);
            
            for id = 1:numel(dEvalPoints)

                [~, dInterp] = self.evalInterpolant(dEvalPoints(id), true);
                test.verifyLessThan(abs(dInterp - sin(2*pi*dEvalPoints(id))), 1e-2);
            
            end
        end

        function testSlerpQuaternionFitWithDiscontinuity(test)
            
            dInterpDomain = linspace(0,1,1000);
            dInterpDomainNormalized = (dInterpDomain - dInterpDomain(1)) / (dInterpDomain(end) - dInterpDomain(1));
            
            dQuat0 = [1;0;0;0];
            dQuat1 = [-1;0.5;0.9;1];
            dQuat1 = dQuat1/norm(dQuat1);

            [dQuatSequence] = InterpolateSlerp(dQuat0, dQuat1, dInterpDomainNormalized);

            % Introduce a sign change
            ui32SwitchIdx = floor(size(dQuatSequence, 2)/2);
            dQuatSequence(:, ui32SwitchIdx:end) = -dQuatSequence(:, ui32SwitchIdx:end);

            dDataMatrix = dQuatSequence;

            % Instantiate quaternion interpolator
            self = CChbvInterpolator(dInterpDomain, uint8(5), EnumInterpType.QUAT, true, [0,1], int32(4));
            % Fit data matrix (should fix discontinuity internally)
            [self, dCoeffsT, strFitStats] = self.fitDataMatrix(dDataMatrix);
            
            % Evaluate at points before and after discontinuity
            testPoints = [dInterpDomain(ui32SwitchIdx-5), dInterpDomain(ui32SwitchIdx+5)];
        
            ui32GridIndex = [ui32SwitchIdx-5, ui32SwitchIdx+5];
            for id = 1:length(ui32GridIndex)
        
                [~, dInterp] = self.evalInterpolant(testPoints(id), true);
                dExpected = dQuatSequence(:, ui32GridIndex(id));
                
                test.verifyLessThan(dot(dInterp, dExpected) - 1, 1e-5);
            
            end
        end

        

    end
end
