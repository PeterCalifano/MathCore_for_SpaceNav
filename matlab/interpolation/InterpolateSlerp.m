function dQuatOut = InterpolateSlerp(dQ1, dQ2, dTGrid)%#codegen
%SLERP Spherical linear interpolation between two quaternions
%   dQuatOut = Slerp(dQ1, dQ2, dTGrid) returns a 4xN array of interpolated
%   quaternions between dQ1 and dQ2 at the fractional times specified in
%   dTGrid (1xN vector in [0,1]).

arguments
    dQ1     (4,1) double
    dQ2     (4,1) double
    dTGrid  (1,:) double 
end

assert(dTGrid(1) >= 0, 'ERROR: timegrid must be normalized to [0,1] values.')
assert(dTGrid(end) <= 1, 'ERROR: timegrid must be normalized to [0,1] values.')

% Compute dot product for angle
dTmpDot = dot(dQ1, dQ2);

% Shortest path correction
if dTmpDot < 0
    dTmpDot = -dTmpDot;
    dQ2   = -dQ2;
end

% Preallocate output
ui32TmpNumTimes = uint32(numel(dTGrid));
dQuatOut        = zeros(4, ui32TmpNumTimes);

% Check for nearly linear case
if dTmpDot > 0.9995
    % Linear interpolation (LERP) with renormalization
    dTmpRatioA = 1 - dTGrid;
    dTmpRatioB = dTGrid;
    
    for ui32TmpIdx = 1:ui32TmpNumTimes
        dTmpQuat = dTmpRatioA(ui32TmpIdx)*dQ1 + dTmpRatioB(ui32TmpIdx)*dQ2;
        dTmpNorm = norm(dTmpQuat);
        dQuatOut(:,ui32TmpIdx) = dTmpQuat / dTmpNorm;
    end

else
    % Full SLERP
    dTmpTheta    = acos(dTmpDot);
    dTmpSinTheta = sin(dTmpTheta);
    
    for ui32TmpIdx = 1:ui32TmpNumTimes
        dTmpRatioA = sin((1 - dTGrid(ui32TmpIdx))*dTmpTheta) / dTmpSinTheta;
        dTmpRatioB = sin(dTGrid(ui32TmpIdx)*dTmpTheta)       / dTmpSinTheta;
        dQuatOut(:,ui32TmpIdx) = dTmpRatioA*dQ1 + dTmpRatioB*dQ2;
    end
end

end
