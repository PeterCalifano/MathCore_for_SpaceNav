function dChbvPolynomial = EvalRecursiveChbv(ui32PolyDeg, dScaledPoint, ui32PolyMaxDeg) %#codegen
%% SIGNATURE
% dChbvPolynomial = EvalRecursiveChbv(ui32PolyDeg, dScaledPoint, ui32PolyMaxDeg)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Evaluate the Chebyshev basis up to the active degree within a fixed-size workspace.
% Entries beyond the active degree remain zero. Keep the maximum degree constant for static codegen.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% ui32PolyDeg       Active degree, at least two and no greater than ui32PolyMaxDeg.
% dScaledPoint      Scalar evaluation coordinate in [-1,1].
% ui32PolyMaxDeg    Workspace degree bound; defaults to the active degree.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dChbvPolynomial  Basis values T0 through Tdegree, followed by zero-filled unused capacity.
%                  The vector always contains ui32PolyMaxDeg+1 entries.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-04-2024    Pietro Califano         First version. Validated.
% 01-02-2025    Pietro Califano         Function modified for compatibility with static size requirements
% 18-07-2025    Pietro Califano         Fix basis and fitting problem errors
% 09-09-2026    Pietro Califano, Codex gpt-6    Enforce the fixed workspace degree bound.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    ui32PolyDeg       (1,1) uint32 {mustBeNumeric}
    dScaledPoint      (1,1) double {mustBeNumeric}
    ui32PolyMaxDeg    (1,1) uint32 {coder.mustBeConst, mustBeNumeric} = ui32PolyDeg
end
arguments (Output)
    dChbvPolynomial  (:,1) double
end

%% Function code

if coder.target('MATLAB') || coder.target('MEX')
    assert(ui32PolyDeg >= 2, 'Error: selected degree is too low!')
    assert(ui32PolyDeg <= ui32PolyMaxDeg, 'EvalRecursiveChbv:DegreeExceedsMaximum', ...
        'Active polynomial degree exceeds the configured maximum.');
end

dChbvPolynomial = coder.nullcopy(zeros(ui32PolyMaxDeg + 1, 1));

% Initialize recursion
dChbvPolynomial(1) = 1.0;
dChbvPolynomial(2) = dScaledPoint;

for idN = 3:ui32PolyDeg + 1
    dChbvPolynomial(idN) = 2.0 * dScaledPoint * dChbvPolynomial(idN-1) - dChbvPolynomial(idN-2);
end

end
