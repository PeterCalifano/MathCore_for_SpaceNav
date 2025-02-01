function dChbvPolynomial = EvalRecursiveChbv(ui32PolyDeg, dScaledPoint, ui32PtrToLastCoeff, ui32PolyMaxDeg) %#codegen
arguments
    ui32PolyDeg         (1,1) uint32 {isscalar, isnumeric}
    dScaledPoint        (1,1) double {isscalar, isnumeric}
    ui32PtrToLastCoeff  (1,1) uint32 {isscalar, isnumeric}  
    ui32PolyMaxDeg      (1,1) uint32 {isscalar, isnumeric} = ui32PolyDeg
end
%% PROTOTYPE
% TODO: update doc
% dChbvPolynomial = EvalRecursiveChbv(ui32PolyDeg, dScaledPoint, ui32PolyMaxDeg)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function constructing the Chebyshev Polynomial evaluated at dScaledPoint point in [-1, 1]
% interval up to dPolyDeg degree. No coefficient is applied (assumed as ones).
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% ui32PolyDeg:      [1]   Degree of the Chebyshev polynomial
% dScaledPoint:     [1]   Chebyshev polynomial evaluation point (scalar only). Must be in [-1,1]
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dChbvPolynomial: [ui32PolyDeg+1, 1]  Vector of evaluation Chebyshev polynomials (no coefficients) 
%                                      defined over [-1,1] domain. Max size is determined by ui32PolyMaxDeg.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-04-2024    Pietro Califano         First version. Validated.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

assert(ui32PolyDeg > 2, 'Error: selected degree is too low!')
dChbvPolynomial = coder.nullcopy(zeros(ui32PolyDeg+1, 1));

% Initialize recursion
dChbvPolynomial(1) = 0.0;
dChbvPolynomial(2) = 1.0;

for idN = 3:ui32PolyDeg+1
    dChbvPolynomial(idN) = 2.0 * dScaledPoint * dChbvPolynomial(idN-1) - dChbvPolynomial(idN-2);
end

end
