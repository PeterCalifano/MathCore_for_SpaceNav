function [state, evalcounter, starter_fk] = RKN_integrator(RHS, tspan, state0, order) %#codegen

ifnotstarter = 1;

h_fixed = tspan(2) - tspan(1);
h_vec = h_fixed*ones(length(tspan), 1);

% Modify last stepsize to solve up to the correct interval 
if tspan(end) - tspan(end-1) ~= h_fixed
    h_vec(end) = tspan(end) - tspan(end-1);
end


if nargout > 2
    ifnotstarter = 0; % If the method is used as starter, fk is evaluated at all point of tspan. Otherwise last cycle is skipper.
end

% ACHTUNG: RHS MUST BE SPECIFIED WITH VARS AS (x, t)
% NOTE: RHS must be a function with dxdt as output

switch order
    case 2

        % ADD CHECK on tspan and h

        % ADD CHECK on x0 and RHS size and adjust input to RHS appropriately to
        % handle multivariable functions

        % Static memory allocation
        xk = nan(length(tspan)+1, length(x0));
        evalcounter = 0;
        % Initial condition
        xk(1, :) = x0;

        h = tspan(2) - tspan(1);

        for idt = 1:length(tspan)

            fk = RHS(xk(idt, :), tspan(idt));
            evalcounter = evalcounter + 1;

            if nargout > 2
                starter_fk(idt, :) = fk;
            end

            % Predictor step
            x_p = xk(idt, :) + h_vec(idt)*fk;
            % Corrector step to get x at k+1 step
            xk(idt + 1, :) = xk(idt, :) + h_vec(idt)* (0.5*fk + 0.5* RHS(x_p, tspan(idt)+h_vec(idt)));
            evalcounter = evalcounter + 1;
        end
        state = xk(1:end-1, :);


    case 4

        %% RK4 Integrator scheme
%         RK4timer = tic;
        evalcounter = 0;

        time_vec = tspan; %et0:h:etf;
%         h = tspan(2) - tspan(1); % Fixed step derived from tspan

        alfa = [1/2, 1/2, 1, 1]';

        beta = [1/2, 0, 0, 0;
            0, 1/2, 0, 0;
            0, 0, 1, 0;
            1/6, 1/3, 1/3, 1/6];


        state = nan(length(time_vec), length(state0));
        % Initial condition assignment
        state(1, :) = state0; %x0;

        for timestep = 1:length(time_vec) - ifnotstarter

            x_k = state(timestep, :)';
            t_k = time_vec(timestep);
            
            % 1st Predictor stage
            fk = RHS(x_k, t_k);
            evalcounter = evalcounter + 1;
            x_P1 = x_k + alfa(1)* h_vec(timestep) * fk;  %f(x_k, t_k)

            if nargout > 2
                if timestep == 1
                    starter_fk = nan(length(tspan), length(fk));
                end

                starter_fk(timestep, :) = fk;
            end

            % 2nd Predictor stage
            fP1 = RHS(x_P1, t_k + alfa(1)*h_vec(timestep));
            evalcounter = evalcounter + 1;
            x_P2 = x_k + alfa(2)*h_vec(timestep) *fP1;  %f(x_P1, t_k + alfa(1)*h)

            % 3rd Predictor stage
            fP2 = RHS(x_P2,t_k + alfa(2)*h_vec(timestep));
            evalcounter = evalcounter + 1;
            x_P3 = x_k + alfa(3)*h_vec(timestep) *fP2; %f(x_P2, t_k + alfa(2)*h)

            % 4th Corrector stage
            fP3 = RHS(x_P3, t_k + alfa(3)*h_vec(timestep));
            evalcounter = evalcounter + 1;

            x_next = x_k + alfa(4)*h_vec(timestep) * (beta(4, 1)*fk  + beta(4, 2)* fP1 + beta(4, 3)* fP2 + beta(4, 4)* fP3);

            state(timestep + 1, :) = x_next;

        end

        if ifnotstarter == 0
            state = state(1:end-1, :);
        end

%         toc(RK4timer)

end

