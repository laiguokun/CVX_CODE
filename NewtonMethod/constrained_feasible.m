% The funciton is f(x) = exp(x1 + 3*x2 - 0.1) + exp(x1-3*x2-0.1) +
% exp(-x1-0.1)
% under equality Ax = b A and b is random
% feasible start
function [] = constrained_feasible()

x0 = [1,1];
x = x0;
gap = f(x);
res = 2.55927;
n = 2; % n is the number of variable
m = 1; % m is the number of equality
A = rand(m,n);
b = rand(m,1);
x(1) = rand();
x(2) = (b - A(1,1) * x(1))/A(1,2);
disp(gap - res);
cnt = 0;

while (cnt < 10)
    df = fg(x);
    ddf = Hession(x);
    P = [ddf, A'; A, zeros(m,m)];
    tmp = P ^ -1 * [-df'; zeros(m,1)];
    dx = tmp(1:n)';
    lambda = dx * ddf * dx';
    if (lambda < 1e-6)
        break;
    end;
    x = backtrack_line_search(x, df, dx, 0.01, 0.5);
    gap = f(x) - res; 
    disp('gap');
    disp(gap);
%    disp('equlity');
%    disp(A * x' - b);
    cnt = cnt + 1;
end

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
