% The funciton is f(x) = exp(x1 + 3*x2 - 0.1) + exp(x1-3*x2-0.1) +
% exp(-x1-0.1)
% under equality Ax = b A and b is random
% infeasible start
function [] = constrained_infeasible()

x0 = rand(1,2);;
x = x0';
gap = f(x);
res = 2.55927;
n = 2; % n is the number of variable
m = 1; % m is the number of equality
A = rand(m,n);
b = rand(m,1);
v = zeros(m,1);
disp(gap - res);
disp(A * x - b);
cnt = 0;

while (cnt < 10)
    df = fg(x);
    ddf = Hession(x);
    P = [ddf, A'; A, zeros(m,m)];
    delta = P ^ -1 * -[df' + A' * v; (A*x - b)];
    xv = [x;v];
    lambda = residual(xv, n, m, A, b);
    if (lambda < 1e-6)
        break;
    end;
    xv = backtrack_line_search(xv, delta, 0.01, 0.5, n, m, A, b);
    x = xv(1:n);
    v = xv(n+1:n+m);
    gap = f(x) - res; 
    disp('gap');
    disp(gap);
    disp('equlity');
    disp(A*x - b);
    cnt = cnt + 1;
end

end

function value = f(x)
    value = exp(x(1) + 3*x(2) - 0.1) + exp(x(1)-3*x(2)-0.1) + exp(-x(1) - 0.1);
end

function g = fg(x)
    g = [exp(x(1) + 3*x(2) - 0.1) + exp(x(1)-3*x(2)-0.1) - exp(-x(1) - 0.1), 3 * exp(x(1) + 3*x(2) - 0.1) - 3 * exp(x(1)-3*x(2)-0.1)];
end

function x = backtrack_line_search(xv, delta, alpha, beta, n, m, A, b)
    t = 1;
    while (residual(xv + t * delta, n, m, A, b) > (1 - alpha * t) * residual(xv, n, m, A, b))
        t = t * beta;
    end;
    x = xv + t * delta;
end

function P = Hession(x)
    P = size(2,2);
    P(1,1) = exp(x(1) + 3*x(2) - 0.1) + exp(x(1)-3*x(2)-0.1) + exp(-x(1) - 0.1);
    P(1,2) = 3 * exp(x(1) + 3*x(2) - 0.1) - 3 * exp(x(1)-3*x(2)-0.1);
    P(2,1) = P(1,2);
    P(2,2) = 9 * exp(x(1) + 3*x(2) - 0.1) + 9 * exp(x(1)-3*x(2)-0.1);
end

function res = residual(xv, n, m, A, b)
    x = xv(1:n);
    v = xv(n+1:n+m);
    df = fg(x)';
    rdual = df + A'*v;
    rpri = A*x - b;
    res = norm([rdual;rpri], 2);
end

