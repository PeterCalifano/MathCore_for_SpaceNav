function dChbvInterpVector = evalChbvPolyWithCoeffs(ui32PolyDeg, ui32OutputSize, ...
    dEvalPoint, dChbvCoeffs, dDomainLB, dDomainUB, ui32PtrToLastCoeff, ui32PolyMaxDeg) %#codegen
%% SIGNATURE
% dChbvInterpVector = evalChbvPolyWithCoeffs(ui32PolyDeg, ui32OutputSize, ...
%     dEvalPoint, dChbvCoeffs, dDomainLB, dDomainUB, ui32PtrToLastCoeff, ui32PolyMaxDeg)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Evaluate a vector Chebyshev series with a runtime degree and a fixed maximum workspace.
% Pack active coefficients by component: all degree+1 coefficients of component 1, then component 2,
% and so on. Ignore storage after the active prefix. For static codegen, output size and maximum
% degree must be constant; active degree, coefficients, time and domain bounds may change at runtime.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% ui32PolyDeg           Active degree, at least 2 and no greater than ui32PolyMaxDeg.
% ui32OutputSize        Number of vector components.
% dEvalPoint            Evaluation time in [dDomainLB,dDomainUB].
% dChbvCoeffs           Packed active coefficients, optionally followed by unused storage.
% dDomainLB             Lower bound of the interpolation interval.
% dDomainUB             Upper bound, strictly greater than dDomainLB.
% ui32PtrToLastCoeff    Active coefficient count; defaults to (degree+1)*output size.
% ui32PolyMaxDeg        Workspace degree bound; defaults to the active degree.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dChbvInterpVector     Interpolated vector, with ui32OutputSize components.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-04-2024    Pietro Califano     First version, verified in unit test.
% 01-02-2025    Pietro Califano     Upgrade of functions for codegen with static-sized arrays.
% 18-07-2025    Pietro Califano     Fix basis and fitting problem errors.
% 09-09-2026    Pietro Califano, Codex gpt-6    Evaluate active degrees within fixed capacity.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% EvalRecursiveChbv.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    ui32PolyDeg         (1,1) uint32
    ui32OutputSize      (1,1) uint32
    dEvalPoint          (1,1) double
    dChbvCoeffs         (:,1) double
    dDomainLB           (1,1) double
    dDomainUB           (1,1) double
    ui32PtrToLastCoeff  (1,1) uint32 = ui32OutputSize * (ui32PolyDeg + 1)
    ui32PolyMaxDeg      (1,1) uint32 {coder.mustBeConst} = ui32PolyDeg
end
arguments (Output)
    dChbvInterpVector   (:,1) double
end

% Validate the active prefix before evaluating or indexing its storage.
if coder.target('MATLAB') || coder.target('MEX')

    assert(ui32PolyDeg <= ui32PolyMaxDeg, 'evalChbvPolyWithCoeffs:DegreeExceedsMaximum', ...
        'Active polynomial degree exceeds the configured maximum.');
    assert(ui32PtrToLastCoeff == (ui32PolyDeg+1)*ui32OutputSize && ...
        ui32PtrToLastCoeff <= numel(dChbvCoeffs), 'evalChbvPolyWithCoeffs:CoefficientCount', ...
        'Coefficient storage does not contain the requested active series.');
    assert(dDomainUB > dDomainLB, 'evalChbvPolyWithCoeffs:InvalidDomain', ...
        'Interpolation upper bound must exceed its lower bound.');
    assert(dEvalPoint >= dDomainLB && dEvalPoint <= dDomainUB, ...
        'ERROR: invalid evaluation point. Out of interpolation bound.');
end

% The recurrence returns the full workspace; only the active prefix contributes to the result.
dScaledPoint = (2*dEvalPoint - (dDomainLB + dDomainUB)) / (dDomainUB - dDomainLB);
dChbvPolynomial = EvalRecursiveChbv(ui32PolyDeg, dScaledPoint, ui32PolyMaxDeg);
dChbvInterpVector = coder.nullcopy(zeros(ui32OutputSize, 1));

% Evaluate each component by its active coefficients and the Chebyshev basis.
for ui32Component = uint32(1):ui32OutputSize
    
    ui32CoeffOffset = (ui32Component-1)*(ui32PolyDeg+1);
    dComponentValue = 0.0;

    for ui32Term = uint32(1):ui32PolyDeg+1
        dComponentValue = dComponentValue + ...
            dChbvCoeffs(ui32CoeffOffset+ui32Term)*dChbvPolynomial(ui32Term);
    end

    dChbvInterpVector(ui32Component) = dComponentValue;
end
end
