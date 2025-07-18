classdef CInterpolatorTestHelper < CInterpolator
    %% Helper subclass to expose protected methods for testing
    methods
        function self = CInterpolatorTestHelper(dInterpDomain, ...
                                            ui8PolyDeg, ...
                                            enumInterpType, ...
                                            bENABLE_AUTO_CHECK, ...
                                            dDomainBounds, ...
                                            i32OutputVectorSize)
            % Call superclass constructor
            self@CInterpolator(dInterpDomain, ...
                                ui8PolyDeg, ...
                                enumInterpType, ...
                                bENABLE_AUTO_CHECK, ...
                                dDomainBounds, ...
                                i32OutputVectorSize);
        end
        function [self, dInterpVector] = evalInterpolant(self, ~, ~)
            % Dummy implementation for abstract method
            dInterpVector = zeros(self.i32OutputVectorSize, 1);
        end
        function [self, dPolyTermsValues] = evalPoly(self, ~, ~)
            % Dummy implementation for abstract method
            dPolyTermsValues = zeros(self.ui8PolyDeg + 1, 1);
        end
        function [self, dInterpCoeffsMatrix, strFitStats] = fitDataMatrix(self, dDataMatrix)
            % Dummy implementation for abstract method
            dInterpCoeffsMatrix = zeros(self.ui8PolyDeg, size(dDataMatrix, 2));
            strFitStats = struct();
        end
        function [dModifiedDataMatrix, dSwitchIntervals, bIsSignSwitched, ui8HowManySwitches] = testFixQuatSignDiscontinuity(self, dQuatMatrix)
            % Public wrapper for protected fixQuatSignDiscontinuity
            [~, dModifiedDataMatrix, dSwitchIntervals, ...
                bIsSignSwitched, ui8HowManySwitches] = fixQuatSignDiscontinuity(self, dQuatMatrix);
        end
    end
end
