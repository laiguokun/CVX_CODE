% The funciton is f(x) = exp(x1 + 3*x2 - 0.1) + exp(x1-3*x2-0.1) +
% exp(-x1-0.1)
% under inequality Ax + b > 0 A and b is random
function [] = interiorpoint_feasible()

n = 2;
p = 1;
res = 2.55927;
global A;
%A = rand(p,n);
A = [1,1];
global b;
%b = rand(p,1);
b = [0];
x = zeros(2,1);
%%find feasible start point
x(1) = rand();
x(2) = (b - A(1,1) * x(1))/A(1,2);
if (A*(x+1) + b > 0)
    x = x + 1;
else
    x = x - 1;
end
t = abs(origin_f(x));
mu = 10;
while (p/t > 1e-10)
    x = newton(t,x);
    t = mu * t;
end
disp('result');
disp(origin_f(x));
disp(x);
end

function x = newton(t,x0)
    x = x0;
    lambda = 1;
    cnt = 0;
    while(lambda > 1e-10)
        df = fg(x,t);
        ddf = Hession(x,t);
        P = ddf;
        dx = P^-1 * -df;
        lambda = (dx' * ddf * dx) / 2;
        x = backtrack_line_search(x, df, dx, 0.01, 0.5, t);  
    end
end

function value = origin_f(x)
    value = (exp(x(1) + 3*x(2) - 0.1) + exp(x(1)-3*x(2)-0.1) + exp(-x(1) - 0.1));
end
function value = f(x, t)
    value = (exp(x(1) + 3*x(2) - 0.1) + exp(x(1)-3*x(2)-0.1) + exp(-x(1) - 0.1));
    value = value * t - sum(fie(x));
end

function value = fie(x)  %function of inequality
    global A;
    global b;
    value = log(A * x + b);
end

function g = fg(x, t)
    g = [exp(x(1) + 3*x(2) - 0.1) + exp(x(1)-3*x(2)-0.1) - exp(-x(1) - 0.1); 3 * exp(x(1) + 3*x(2) - 0.1) - 3 * exp(x(1)-3*x(2)-0.1)];
    g = t * g - fgie(x);
end

function g = fgie(x) %gradient function of inequality
    global A;
    global b;
    F = A * x + b;
    g = F.^-1;
    g = g * sum(A, 1)';
end

function x = backtrack_line_search(x0, df, dx, alpha, beta, tt)
    global A;
    global b;
    t = 1;
    while ((min(A * (x0 + t * dx) + b) < 0) || f(x0 + t * dx, tt) > f(x0, tt) + alpha * t * df' * dx)
        t = t * beta;
    end;
    x = x0 + t * dx;
end

function P = Hession(x, t)
    P = zeros(2,2);
    P(1,1) = exp(x(1) + 3*x(2) - 0.1) + exp(x(1)-3*x(2)-0.1) + exp(-x(1) - 0.1);
    P(1,2) = 3 * exp(x(1) + 3*x(2) - 0.1) - 3 * exp(x(1)-3*x(2)-0.1);
    P(2,1) = P(1,2);
    P(2,2) = 9 * exp(x(1) + 3*x(2) - 0.1) + 9 * exp(x(1)-3*x(2)-0.1);
    P = t * P - Hessionie(x);
end

function P = Hessionie(x) %Hession of inequality
    global A;
    global b;
    F = A * x + b;
    df = fgie(x);
    P = -(df/F) * sum(A,1);
end

