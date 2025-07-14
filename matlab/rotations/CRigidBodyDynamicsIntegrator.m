%% RigidBodyDynamicsIntegrator.m
classdef CRigidBodyDynamicsIntegrator < handle & matlab.mixin.Copyable
    properties
        dInertiaMatrix          (3,3) double {mustBeFinite} = eye(3)
        objQuatKinIntegrator    QuaternionIntegrator = QuaternionIntegrator()
    end

    methods
        function self = CRigidBodyDynamicsIntegrator(bus_Inertia, objQuatInt)
            arguments
                bus_Inertia (3,3) double {mustBeFinite} = eye(3)
                objQuatInt (1,1) QuaternionIntegrator = QuaternionIntegrator()
            end
            self.dInertiaMatrix = bus_Inertia;
            self.objQuatKinIntegrator = objQuatInt;
        end

        function dvAlpha = AngularAcceleration(self, dvOmega, dvTorque)
            arguments
                self (1,1) CRigidBodyDynamicsIntegrator
                dvOmega (3,1) double {mustBeFinite}
                dvTorque (3,1) double {mustBeFinite}
            end
            dvTmpInertiaOmega = self.dInertiaMatrix * dvOmega;
            dvAlpha = self.bus_Inertia \
                      (dvTorque - cross(dvOmega, dvTmpInertiaOmega));
        end

        function [dQuatOut, dvOmegaOut] = JointEulerStep(self, dQuat, dvOmega, fTorque, dDeltaT)
            arguments
                self (1,1) CRigidBodyDynamicsIntegrator
                dQuat (4,1) double {mustBeFinite}
                dvOmega (3,1) double {mustBeFinite}
                fTorque function_handle
                dDeltaT (1,1) double {mustBePositive}
            end
            dvTmpTorque = fTorque();
            dvTmpAlpha  = self.AngularAcceleration(dvOmega, dvTmpTorque);
            dvOmegaOut  = dvOmega + dvTmpAlpha * dDeltaT;
            dQuatTmp    = self.objQuatKinIntegrator.LieEulerStep(dQuat, dvOmega, dDeltaT);
            dQuatOut    = self.objQuatKinIntegrator.Normalize(dQuatTmp);
        end

        function [dQuatOut, dvOmegaOut] = Integrate(self, dQuat0, dvOmega0, fTorque, dT0, dT1, dDeltaT)
            arguments
                self (1,1) CRigidBodyDynamicsIntegrator
                dQuat0 (4,1) double {mustBeFinite}
                dvOmega0 (3,1) double {mustBeFinite}
                fTorque function_handle
                dT0 (1,1) double {mustBeFinite}
                dT1 (1,1) double {mustBeGreaterThan(dT1,dT0)}
                dDeltaT (1,1) double {mustBePositive}
            end
            dQuat = self.objQuatKinIntegrator.Normalize(dQuat0);
            dvOmega = dvOmega0;
            dTmpT = dT0;
            while dTmpT < dT1 - 1e-12
                dTmpStep = min(dDeltaT, dT1 - dTmpT);
                [dQuat, dvOmega] = self.JointEulerStep(dQuat, dvOmega, fTorque, dTmpStep);
                dTmpT = dTmpT + dTmpStep;
            end
            dQuatOut = dQuat;
            dvOmegaOut = dvOmega;
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
