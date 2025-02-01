function [outMatrix] = mldivideMex(A, B)

% A: RHS matrix
% B: LHS vector or matrix
% outMatrix: solution of the Linear system (exact or in LS sense)

outMatrix = mldivide(A, B);

end