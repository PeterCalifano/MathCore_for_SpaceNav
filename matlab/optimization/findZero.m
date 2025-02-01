function [x_zero, feval, info] = findZero(fun, interval, method, tol) %#codegen

if nargin < 4
    tol = 1e-8;
    disp(['Default tolerance: ', num2str(tol)])
end

% Initialization of counters and variables
evalcounter = 0;
MaxIter = 200;
safety_var = 0;

a = interval(1);
b = interval(2);

f_a = fun(a);
evalcounter = evalcounter + 1;
f_b = fun(b);
evalcounter = evalcounter + 1;



if b < a
    warning('Left value of the interval seems to be greater than the right value. Zero finding procedure may fail.')
end

switch method

    case 1 % Bi-section
        c = b;
        err = 1;
        while abs(err) > tol
            if safety_var == MaxIter
                warning('Solution not found in specified Max Iteration number')
                break;
            end

            if f_a * f_b > 0
                error('Zero not existent in the selected interval');
            end

            c_prev = c;
            % Compute middle point of the interval
            c = (a+b)/2;
            % Compute difference with respect to 
            err = abs(c_prev - c);
            % Evaluate function in middle point
            
            f_c = fun(c);
            evalcounter = evalcounter + 1;

            % Determine next search interval
            if f_a*f_c < 0
                b = c;
                f_b = f_c;
            else
                a = c;
                f_a = f_c;
            end
            safety_var = safety_var + 1;
            
        end

        fprintf('\nNumber of cycles for bi-section method: %i\n', safety_var);

        % Return output
        x_zero = c;
        feval = f_c;
        info.evalcounter = evalcounter;

    case 2 % Secants

        x_k = b;
        xprev = a;

        f_xprev = f_b;
        f_k = f_a;
        err = 1;

        while abs(err) > tol
            if safety_var == MaxIter
                warning('Solution not found in specified Max Iteration number')
                break;
            end

            x_new = x_k - (x_k - xprev)*f_k./(f_k-f_xprev);

            % Update x and f
            xprev = x_k;
            f_xprev = f_k;

            err = abs(x_new - x_k);

            if err < tol
                break;
            end

            x_k = x_new;       
            f_k = fun(x_k);         
            evalcounter = evalcounter + 1;

            safety_var = safety_var + 1;
        end

        fprintf('\nNumber of cycles for secants method: %i\n', safety_var);

        % Return output
        x_zero = x_k;
        feval = f_k;
        info.evalcounter = evalcounter;

    case 3 % Regula falsi

        x_k = b;
        xprev = a;

        % Initilize variables
        f_xk = f_b;
        f_xprev = f_a;
        err = 1;

        while err > tol

            if safety_var == MaxIter
                break;
            end

            x_new = x_k - (x_k - xprev)*f_xk./(f_xk-f_xprev);

            err = abs(x_k - x_new);

            if err < tol
                break;
            end

            f_xnew = fun(x_new);
            evalcounter = evalcounter + 1;

            if f_xprev*f_xnew < 0
                x_k = x_new;
                f_xk = f_xnew;
            else
                xprev = x_new;
                f_xprev = f_xnew;
            end

            safety_var = safety_var + 1;
        end

        fprintf('\nNumber of cycles for Regula Falsi method: %i \n', safety_var);

        % Return Output
        x_zero = x_new;
        feval = f_xnew;
        info.evalcounter = evalcounter;
end
end