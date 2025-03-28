function [xk_tf] = RK_evalfinalstate(F, x0, h, t0, tf)

% Evaluate how many nsteps are required to cover tspan
h = abs(h);
nstep = floor( (tf-t0)/h ); 

if h*nstep > tf
    warning('tf has been exceeded')
end

% Compute final state of integration
if h*nstep == tf
    xk_tf = (F(h)^nstep)*x0;
else
    
    xk_tf = F(tf - h*(nstep))*(F(h)^nstep)*x0;
end

end