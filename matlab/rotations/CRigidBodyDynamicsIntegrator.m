%% RigidBodyDynamicsIntegrator.m
classdef CRigidBodyDynamicsIntegrator < handle & matlab.mixin.Copyable
    properties
        dInertiaMatrix          (3,3) double {mustBeFinite}     = zeros(3,3)
        objQuatKinIntegrator    QuaternionIntegrator            = QuaternionIntegrator()
    end

    methods (Access = public)
         % CONSTRUCTOR
        function self = CRigidBodyDynamicsIntegrator(dInertiaMatrix, objQuatInt)
            arguments
                dInertiaMatrix (3,3) double {mustBeFinite} = eye(3)
                objQuatInt (1,1) QuaternionIntegrator = QuaternionIntegrator()
            end
            % Assign attributes
            self.dInertiaMatrix = dInertiaMatrix;
            self.objQuatKinIntegrator = objQuatInt;
        end

        % PUBLIC METHODS
        function [dQuatSeq, dOmegaSeq] = integrate(self, ...
                                                   dTimegrid, ...
                                                   dQuat0, ...
                                                   dOmega0, ...
                                                   varTorque, ...
                                                   dDeltaT, ...
                                                   dDefaultMaxDeltaT)
            arguments
                self                (1,1) CRigidBodyDynamicsIntegrator
                dTimegrid           (1,:) double {mustBeFinite, mustBeVector}
                dQuat0              (4,1) double {mustBeFinite}
                dOmega0             (3,1) double {mustBeFinite}
                varTorque           {mustBeA(varTorque, ["function_handle", "double"])}
                dDeltaT             (1,1) double {mustBeScalarOrEmpty, mustBeGreaterThanOrEqual(dDeltaT, 0.0)} = 0.0;
                dDefaultMaxDeltaT   (1,1) double {mustBeScalarOrEmpty, mustBeGreaterThanOrEqual(dDefaultMaxDeltaT, 0.0)} = 1.0;
            end

            % TODO implement integration step from t0 to tf
            dQuatSeq = zeros(4,  ui32NumTimes);
            dOmegaSeq = zeros(3, ui32NumTimes);

        end

    end

    methods (Static)
        %%% RHS components
        % Angular acceleration RHS
        function dAngAccel = EvalRHS_AngAccel_(dInertiaMatrix, ...
                                               dOmegaAngVel, ...
                                               dExtTorque) %#codegen
            arguments
                dInertiaMatrix  (3,3) double {mustBeFinite}
                dOmegaAngVel    (3,1) double {mustBeFinite}
                dExtTorque      (3,1) double {mustBeFinite} = zeros(3,1);
            end
            arguments

            end

            % Compute angular acceleration alpha = I^-1 * (T - w x Iw)
            dAngAccel =  dInertiaMatrix \ (dExtTorque - cross(dOmegaAngVel, dInertiaMatrix * dOmegaAngVel));
        end

        function [dQuatOut, dOmegaOut] = integrateStep_(dTimestamp, ...
                                                        dInertiaMatrix, ...
                                                        dQuat0, ...
                                                        dOmega0, ...
                                                        dDeltaT, ...
                                                        varTorque) %#codegen
            arguments
                dTimestamp      (1,1) double {mustBeGreaterThanOrEqual(dTimestamp, 0.0)}
                dInertiaMatrix  (3,3) double {ismatrix, isnumeric}
                dQuat0          (4,1) double {mustBeFinite}
                dOmega0         (3,1) double {mustBeFinite}
                dDeltaT         (1,1) double {mustBePositive}
                varTorque       {mustBeA(varTorque, ["function_handle", "double"])} = zeros(3,1)
            end
            % Combined integration step (decoupled or coupled)
            persistent fcnEvalOmegaRHS 

            if isa(varTorque, "function_handle")
                dvTmpTorque_ = varTorque(dQuat0, dTimestamp);
            else
                dvTmpTorque_ = @(dTstamp, dOmega, dQuat0) varTorque; % Constant, defaults to zero
            end

            if isempty(fcnEvalOmegaRHS)
                fcnEvalOmegaRHS = @(dTstamp, dOmega, dQuat0) CRigidBodyDynamicsIntegrator.EvalRHS_AngAccel_(dInertiaMatrix, ...
                                                                                                           dOmegaAngVel, ...
                                                                                                           dvTmpTorque_(dTstamp, dOmega, dQuat0));
            end


            if bDecoupledKinoDynamics
                
                % Integrate RK4 step for angular velocity
                [dOmegaOut, dOmegaStages, dStagesTimes] = CRigidBodyDynamicsIntegrator.IntegrStep_RK4(dOmega0, ...
                                                                                                      dQuat0, ...
                                                                                                      fcnEvalOmegaRHS, ...
                                                                                                      dTimestamp, ...
                                                                                                      dDeltaT);

                % Integrate attitude quaternion
                dQuatOut = CQuatKinematicsIntegrator.IntegrStep_RKMK4(dQuat0, @(dTime) [interp1(dStagesTimes, dOmegaStages(1,:), dTime), ...
                                                                                       interp1(dStagesTimes, dOmegaStages(2,:), dTime), ...
                                                                                       interp1(dStagesTimes, dOmegaStages(3,:), dTime)], dDeltaT);
            else
                error('Not yet implemented')
            end
        end

        function [dNextOmega, dStagesVals, dStagesTimes] = IntegrStep_RK4(dOmega0, ...
                                                                          dQuat0, ...
                                                                          fcnEvalOmegaRHS, ...
                                                                          dTstamp, ...
                                                                          dDeltaTime)
            arguments
                dOmega0         (3,1) double {mustBeFinite}
                dQuat0          (4,1) double {mustBeFinite}
                fcnEvalOmegaRHS  function_handle % Expected signatfcnEvalOmegaRHS(dTstamp, dOmega, dQuat0)
                dTstamp         (1,1) double {mustBeFinite}
                dDeltaTime      (1,1) double {mustBePositive}
            end
            %%% IntegrStep_RK4  Single RK4 step for attitude dynamics of a rigid body
            %
            % Butcher table:
            %    A = [0    0    0    0;
            %         1/2  0    0    0;
            %         0    1/2  0    0;
            %         0    0    1    0];
            %    b = [1/6; 1/3; 1/3; 1/6];
            %    c = [0; 1/2; 1/2; 1];

            % DEVNOTE: current implementation only supports torque-free dynamics integration
            dQuat0   = [1;0;0;0];
            dStagesTimes = dTstamp + [0.0, 0.5 * dDeltaTime, 0.5 * dDeltaTime, dDeltaTime];

            % Stage 1
            dTmpRHS_K1 = fcnEvalOmegaRHS(dStagesTimes(1), dOmega, dQuat0);
            dOmega_K1 = dInitAngVel + 0.5 * dDeltaTime * dTmpRHS_K1;
            dQuat_K1 = dQuat0;

            % Stage 2
            dTmpRHS_K2 = fcnEvalOmegaRHS(dStagesTimes(2), dOmega_K1, dQuat_K1); 
            dOmega_K2 = dInitAngVel + 0.5 * dDeltaTime * dTmpRHS_K2;
            dQuat_K2 = dQuat0;

            % Stage 3
            dTmpRHS_K3 = fcnEvalOmegaRHS(dStagesTimes(3), dOmega_K2, dQuat_K2);
            dOmega_K3 = dInitAngVel + dDeltaTime * dTmpRHS_K3;
            dQuat_K3 = dQuat0;

            % Stage 4
            dTmpRHS_K4 = fcnEvalOmegaRHS(dStagesTimes(4), dOmega_K3, dQuat_K3);
            
            % Define output
            dStagesVals = [dOmega0; dOmega_K1; dOmega_K2; dOmega_K3];

            % Update solution
            dNextOmega = dQuat0 + dDeltaTime * (dTmpRHS_K1 + 2*dTmpRHS_K2 + 2*dTmpRHS_K3 + dTmpRHS_K4)/6;
        end
    end

    
        
end

% Example usage:
% objQuatInt  = QuaternionIntegrator();
% objDynInt   = RigidBodyDynamicsIntegrator(diag([1.2,1.5,2.0]), objQuatInt);
% fTorque     = @() [0.01; 0; -0.02];
% dQuat0      = [1; 0; 0; 0];
% dvOmega0    = [0.1; 0.2; 0.3];
% [dQEnd, dvWEnd] = objDynInt.Integrate(dQuat0, dvOmega0, fTorque, 0, 5, 0.01);
% disp(dQEnd); disp(dvWEnd);
