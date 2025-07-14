classdef QuatKinematicsIntegrator < handle & matlab.mixin.Copyable
    %% DESCRIPTION
    % What the class represent
    % -------------------------------------------------------------------------------------------------------------
    %% CHANGELOG
    % 14-07-2025    Pietro Califano, GPT o4-mini-high     First prototype implementation
    % -------------------------------------------------------------------------------------------------------------
    %% DEPENDENCIES
    % qCross function
    
    properties

    end

    methods
        function self = QuatKinematicsIntegrator()
            % Default constructor
        end

        function [dTmpQuatOut, dQuatOutSeq] = integrate(self, ...
                                                        dQuat0, ...
                                                        dOmegaAngVel, ...
                                                        dTimegrid, ...
                                                        enumMethod, ...
                                                        dDeltaT)
            arguments
                self            (1,1) QuatKinematicsIntegrator
                dQuat0          (4,1) double {mustBeFinite}
                dOmegaAngVel    (3,:) double
                dTimegrid       (1,:) double {mustBeFinite, mustBeVector}
                enumMethod      (1,:) char {mustBeMember(enumMethod, {'lie', 'rkmk4'})}
                dDeltaT         (1,1) double {mustBeScalarOrEmpty} = 0.0;
            end
        
            % Initialize quaternion and timestamps
            dTmpQuatOut = self.NormalizeSeq(dQuat0);

            dT0 = dTimegrid(1);
            dTf = dTimegrid(end);

            idT = 2;            
            if size(dTimegrid,2) == 2
                assert(dDeltaT > 0.0, 'ERROR: timestep must be provided if timegrid is given as [t0,tf].')
                dTimegrid = dT0:dDeltaT:dTf;
            end

            % Allocate output sequence array
            ui32NumSteps = size(dTimegrid, 2);
            dQuatOutSeq = zeros(4, ui32NumSteps);

            % Determine angular velocity mode
            if size(dOmegaAngVel, 2) == 1
                bUseAngularVelProfile = false;
            else
                % TODO add interpolation of omega angular velocity to perform step (use splining)
                bUseAngularVelProfile = true;
                assert(size(dOmegaAngVel, 2) == ui32NumSteps, 'ERROR: angular velocity profile size must match size of timegrid.');

            end            
            dOmegaAngVel_ = dOmegaAngVel(:,1);

            % Select implementation of step
            switch enumMethod
                case 'lie'
                    dQuatIntegr = @(dQuatOut, dOmegaAngVel, dTmpStep) self.IntegrStep_LieGroup(dQuatOut, dOmegaAngVel, dTmpStep);

                case 'rkmk4'
                    dQuatIntegr = @(dQuatOut, dOmegaAngVel, dTmpStep) self.IntegrStep_LieGroup(dQuatOut, dOmegaAngVel, dTmpStep);

            end

            % Integration cycle
            while dCurrentTime < dTf - 1e-12
                % Update integration time step
                if idT < ui32NumSteps
                    dDeltaT = dTimegrid(idT) - dTimegrid(idT-1);
                end

                dTmpDeltaStep = min(dDeltaT, dTf - dCurrentTime);
            
                % Determine angular velocity
                if bUseAngularVelProfile
                    dOmegaAngVel_ = dOmegaAngVel(:,ui32CurrentIntegrStepIdx);
                end

                % Integrate over dTmpDeltaStep time
                dTmpQuatOut = dQuatIntegr(dTmpQuatOut, dOmegaAngVel_, dTmpDeltaStep);
                    
                % Store in sequence
                dTmpQuatOut = self.NormalizeSeq(dTmpQuatOut);
                dQuatOutSeq(:, ui32CurrentIntegrStepIdx) = dTmpQuatOut;

                % Update current index and timegrid
                dCurrentTime = dCurrentTime + dTmpDeltaStep;
                ui32CurrentIntegrStepIdx = ui32CurrentIntegrStepIdx + uint32(1);
            end

        end
    end

    methods (Static)

        function dQuatOut = ExpMap(dDeltaAngle)
            arguments
                dDeltaAngle (3,1) double {mustBeFinite}
            end

            dCurrentTimeheta = norm(dDeltaAngle);

            if dCurrentTimeheta < 1e-12
                dQuatOut = [1.0; 0.0; 0.0; 0.0];
                return;
            end
            % Function applying delta rotation to a quaternion using ExpMap operation.
            % TODO: extend to support sequences of quaternions.

            dvTmpAxis = dDeltaAngle / dCurrentTimeheta;
            dTmpHalfTheta = 0.5 * dCurrentTimeheta;
            dQuatOut = [cos(dTmpHalfTheta); sin(dTmpHalfTheta) * dvTmpAxis];
        
        end

        function dQuatSeqOut = QuatSeqCross(dQuatSeq1, dQuatSeq2, bIS_VSRPplus)
            arguments
                dQuatSeq1 (4,:) double {mustBeFinite}
                dQuatSeq2 (4,:) double {mustBeFinite}
                bIS_VSRPplus (1,1) logical {isscalar, islogical} = false
            end
            % Quaternion convention definition
            % (SV) Scalar first, Vector last
            % (P) Passive
            % (R) Successive coordinate transformations have the unmodified quaternion chain on the Right side of
            %     the triple product.
            % (plus) Right-Handed Rule for the imaginary numbers i, j, k. (aka Hamilton)

            dQuatSeqOut = zeros(size(dQuatSeq1));

            for idQ = 1:size(dQuatSeq1, 2)
                dQuatSeqOut(:,idQ) = qCross(dQuatSeq1(:,idQ), ...
                                        dQuatSeq2(:,idQ), ...
                                        bIS_VSRPplus);
            end

        end

        function dQuatOut = NormalizeSeq(dQuat)
            arguments
                dQuat (4,:) double {mustBeFinite, mustBeNumeric}
            end
            dQuatOut = dQuat ./ vecnorm(dQuat);
        end

        function dQuatOut = IntegrStep_LieGroup(dQuat0, ...
                                                dOmegaAngVel, ...
                                                dDeltaTime)
            arguments
                dQuat0          (4,1) double {mustBeFinite}
                dOmegaAngVel    (3,1) double {mustBeFinite}
                dDeltaTime      (1,1) double {mustBePositive}
            end
            % Function performing an integration setp assuming constant angular velocity and computing the
            % delta quaternion accordingly over the dDeltaTime timestep. ExpMap operation is then used to
            % rotate the initial quaternion.

            dTmpOmegaDelta = 0.5 * dOmegaAngVel * dDeltaTime;
            dQuatOut = self.QuatSeqCross( self.ExpMap(dTmpOmegaDelta), dQuat0 );

        end

        function dQuatOut = IntegrStep_RKMK4(self, dQuat, fOmega, dT, dDeltaT)
            arguments
                self (1,1) QuatKinematicsIntegrator
                dQuat (4,1) double {mustBeFinite}
                fOmega function_handle
                dT (1,1) double {mustBeFinite}
                dDeltaT (1,1) double {mustBePositive}
            end

            dvTmpK1 = 0.5 * fOmega(dT);
            dvTmpK2 = 0.5 * fOmega(dT + dDeltaT/2);
            dvTmpK3 = 0.5 * fOmega(dT + dDeltaT/2);
            dvTmpK4 = 0.5 * fOmega(dT + dDeltaT);

            dvTmpU = (dvTmpK1 + 2*dvTmpK2 + 2*dvTmpK3 + dvTmpK4) * (dDeltaT/6);
            dQuatOut = self.QuatSeqCross(self.ExpMap(dvTmpU), dQuat);

        end
    end
end

