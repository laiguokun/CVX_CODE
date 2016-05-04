% The funciton is f(x) = exp(x1 + 3*x2 - 0.1) + exp(x1-3*x2-0.1) +
% exp(-x1-0.1)
% under inequality Ax + b < 0 A and b is random
function [] = interiorpoint()

n = 2;
p = 1;
res = 2.55927;
global A;
A = rand(p,n);
global b;
b = rand(p,1);
%%find feasible start point
x = A^-1*(-b);
if (A*(x+1) + b < 0)
    x = x + 1;
else
    x = x - 1;
end
t = f(x,0);
disp(f(x,0) - res);
mu = 10;
while (p/t > 1e-10)
    x = newton(t,x);
    t = mu * t;
    break;
end
end

function newton(t,x)
    lambda = 1;
    disp(f(x,t));
    while(lambda > 1e-10)
        df = fg(x,t);
        ddf = Hession(x,t);
        P = ddf;
        dx = P^-1 * -df;
        lambda = dx' * ddf * dx;
        x = backtrack_line_search(x, df, dx, 0.01, 0.5, t);   
        disp(f(x));
        break;
    end
end

function value = f(x, t)
    value = (exp(x(1) + 3*x(2) - 0.1) + exp(x(1)-3*x(2)-0.1) + exp(-x(1) - 0.1));
    value = value * t + sum(fie(x));
end

function value = fie(x) %function of inequality
    global A;
    global b;
    value = A * x + b;
end

function g = fg(x, t)
    g = [exp(x(1) + 3*x(2) - 0.1) + exp(x(1)-3*x(2)-0.1) - exp(-x(1) - 0.1); 3 * exp(x(1) + 3*x(2) - 0.1) - 3 * exp(x(1)-3*x(2)-0.1)];
    g = t * g + fgie(x);
end

function g = fgie(x) %gradient function of inequality
    global A;
    global b;
    g = f(x).^-1;
    g = g .* sum(A)';
end

function x = backtrack_line_search(x0, df, dx, alpha, beta, tt)
    t = 1;
    while (f(x0 + t * dx, tt) > f(x0, tt) + alpha * t * df' * dx)
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
    P = t * P + Hessionie(x);
end

function P = Hessionie(x) %Hession of inequality
    P = zeros(2,2);
    f = fie(x);
    df = fgie(x);
    P(1,1) = -df(1)/f(1) * df(1);
    P(1,2) = -df(1)/f(1) * df(2);
    P(2,1) = -df(2)/f(2) * df(1);
    P(2,2) = -df(2)/f(2) * df(2);
end

