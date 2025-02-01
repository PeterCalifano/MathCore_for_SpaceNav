function [outVal, convFlag] = erf_custom(z) %#codegen
% Validated against erf() function. Error in the order of 1e-16.

assert(length(z) == 1, 'Input must be a scalar.')

if z >= 1
    warning('erf function may have bad convergence property if z >= 1.')
end

% Initialize variables
convFlag = 0;
outVal = 0;
nFactorial = 1;
n = 0;
err = 10;

while err > 1e-16

    % Evaluate next term in the MacLaurin expansion
    nextInSeries = ((-1)^n * z^(2*n+1))./(nFactorial * (2*n+1));
    
    % Sum to current value of the expansion
    prevVal = outVal;
    outVal = outVal + nextInSeries;

    if n == 1e3 % Safe exit 
        break;     
    end

    % Evaluate current increase of the series expansion
    err = abs(outVal-prevVal);

    % Increase the term order
    n = n + 1;

    % Update the factorial for next term
    nFactorial = n * nFactorial;

end

if n < 1e4
    convFlag = 1;
end

% Multiply by scale factor
outVal = 2/(sqrt(pi)) * outVal;


end
