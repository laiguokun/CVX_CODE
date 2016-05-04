% The funciton is f(x) = exp(x1 + 3*x2 - 0.1) + exp(x1-3*x2-0.1) +
% exp(-x1-0.1)
function [] = unconstrained()

x0 = [1,1];
x = x0;
gap = f(x);
res = 2.55927;
disp(gap - res);
cnt = 0;

while (cnt < 20)
    df = fg(x);
    ddf = Hession(x);
    dx = -(ddf^-1 * df')';
    lambda = df * ddf ^ -1 * df';
    if (lambda < 1e-6)
        break;
    end;
    x = backtrack_line_search(x, df, dx, 0.01, 0.5);
    gap = f(x) - res; 
    disp(gap);
    cnt = cnt + 1;
end
disp(x);
end

function value = f(x)
    value = exp(x(1) + 3*x(2) - 0.1) + exp(x(1)-3*x(2)-0.1) + exp(-x(1) - 0.1);
end

function g = fg(x)
    g = [exp(x(1) + 3*x(2) - 0.1) + exp(x(1)-3*x(2)-0.1) - exp(-x(1) - 0.1), 3 * exp(x(1) + 3*x(2) - 0.1) - 3 * exp(x(1)-3*x(2)-0.1)];
end

function x = backtrack_line_search(x0, df, dx, alpha, beta)
    t = 1;
    while (f(x0 + t * dx) > f(x0) + alpha * t * df * dx')
        t = t * beta;
    end;
    x = x0 + t * dx;
end

function P = Hession(x)
    P = size(2,2);
    P(1,1) = exp(x(1) + 3*x(2) - 0.1) + exp(x(1)-3*x(2)-0.1) + exp(-x(1) - 0.1);
    P(1,2) = 3 * exp(x(1) + 3*x(2) - 0.1) - 3 * exp(x(1)-3*x(2)-0.1);
    P(2,1) = P(1,2);
    P(2,2) = 9 * exp(x(1) + 3*x(2) - 0.1) + 9 * exp(x(1)-3*x(2)-0.1);
end
