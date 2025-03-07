function dFiniteDiffJac = ComputeFiniteDiffJacobian(fcn_handle, dX0diff, dEps, varargin)
arguments
    fcn_handle  (1,1) {mustBeA(fcn_handle, ["function_handle", "string", "char"])}
    dX0diff      (:,1) double {isvector, isnumeric} 
    dEps        (1,1) double = 1e-9
end
arguments (Repeating)
    varargin % Varargin passed as additional args to function
end

% Compute expansion point map
dF0 = fcn_handle(dX0diff, varargin{:});
dFiniteDiffJac = zeros(numel(dF0), length(dX0diff));

d2ndOrdFiniteDiffDenom = 1/(2*dEps);

% Loop over input states to evaluate function 
for ii = 1:length(dX0diff)
   
    dEpsVec = zeros(size(dX0diff));
    dEpsVec(ii) = dEps;

    dfPlus = fcn_handle(dX0 + dEpsVec, varargin{:});
    dfMinus = fcn_handle(dX0 - dEpsVec, varargin{:});

    dFiniteDiffJac(:, ii) = d2ndOrdFiniteDiffDenom * (dfPlus - dfMinus);
end

end
