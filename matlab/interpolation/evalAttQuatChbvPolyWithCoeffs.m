function dChbvInterpVector = evalAttQuatChbvPolyWithCoeffs(ui32PolyDeg, ui32OutputSize, ...
    dEvalPoint, dChbvCoeffs, dDomainLB, dDomainUB, ui32PolyMaxDeg) %#codegen
%% SIGNATURE
% dChbvInterpVector = evalAttQuatChbvPolyWithCoeffs(ui32PolyDeg, ui32OutputSize, ...
%     dEvalPoint, dChbvCoeffs, dDomainLB, dDomainUB, ui32PolyMaxDeg)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Evaluate quaternion component series with a runtime degree within fixed capacity.
% Pack the degree+1 active coefficients of each component consecutively, followed by unused
% storage. Reuse the vector evaluator's bounds checks and allocation. The output is not normalized.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% ui32PolyDeg       Runtime active degree, at least 2 and no greater than ui32PolyMaxDeg.
% ui32OutputSize    Component count; normally 4 for quaternions. Constant for static codegen.
% dEvalPoint        Evaluation time in [dDomainLB,dDomainUB].
% dChbvCoeffs       Packed active coefficients, optionally followed by unused capacity.
% dDomainLB        Lower interpolation bound.
% dDomainUB        Upper interpolation bound, strictly greater than dDomainLB.
% ui32PolyMaxDeg   Fixed workspace degree bound; defaults to the active degree.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dChbvInterpVector Interpolated components without normalization or a sign change.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-05-2024    Pietro Califano     First version, modified from general purpose utility. Validated.
% 18-07-2025    Pietro Califano     Fix basis and fitting problem errors.
% 10-09-2026    Pietro Califano, Codex gpt-6    Reuse vector evaluation for runtime quaternion degrees.
% 11-09-2026  Pietro Califano, Codex gpt-6    Remove unused runtime sign-switch metadata.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalChbvPolyWithCoeffs.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    ui32PolyDeg      (1, 1) uint32
    ui32OutputSize   (1, 1) uint32
    dEvalPoint       (1, 1) double
    dChbvCoeffs      (:, 1) double
    dDomainLB        (1, 1) double
    dDomainUB        (1, 1) double
    ui32PolyMaxDeg   (1, 1) uint32 {coder.mustBeConst} = ui32PolyDeg
end
arguments (Output)
    dChbvInterpVector (:, 1) double
end

% The shared evaluator reads only the packed active prefix, including for padded storage.
ui32ActiveCoeffCount = ui32OutputSize * (ui32PolyDeg + 1);
dChbvInterpVector = evalChbvPolyWithCoeffs(ui32PolyDeg, ui32OutputSize, dEvalPoint, ...
    dChbvCoeffs, dDomainLB, dDomainUB, ui32ActiveCoeffCount, ui32PolyMaxDeg);
end
