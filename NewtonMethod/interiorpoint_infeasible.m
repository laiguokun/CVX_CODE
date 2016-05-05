% The funciton is f(x) = exp(x1 + 3*x2 - 0.1) + exp(x1-3*x2-0.1) +
% exp(-x1-0.1)
% under inequality Ax + b > 0 A and b is random
% under equality A1x = b1
function [] = interiorpoint_infeasible()

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
x = [-1;-1];

%stage 1 find feasible start point
%maximize s == minimize -s

s = min(A * x + b) - 1;
disp(s);
x = [s;x];
t = abs(s);
mu = 10;
while (x(1) < 1)
    x = newton_stage1(t,x);
    t = mu * t;
end
disp(A * x(2:3) + b);
%stage 2 

res = 2.55927;
x = x(2:3);
disp(origin_f(x) - res);
t = abs(origin_f(x));
mu = 10;
while (p/t > 1e-10)
    x = newton(t,x);
    t = mu * t;
end
disp(origin_f(x) - res);
end

function x = newton_stage1(t,x0)
    x = x0;
    lambda = 1;
    while(lambda > 1e-10 || x(1) < 1)
        df = fg_stage1(x,t);
        ddf = Hession_stage1(x,t);
        P = ddf;
        if (det(P) == 0)
           dx = -df;
        else
           dx = P \ -df;
        end
        lambda = (dx' * ddf * dx) / 2;
        x = backtrack_line_search_stage1(x, df, dx, 0.01, 0.5, t);  
    end
end

function value = f_stage1(x, t)
    value = -x(1);
    value = value * t - sum(fie_stage1(x));
end

function value = fie_stage1(x)  %function of inequality
    global A;
    global b;
    value = log(A * x(2:3) + b - x(1));
end

function g = fg_stage1(x, t)
    g = [-1;0;0];
    g = t * g - fgie_stage1(x);
end

function g = fgie_stage1(x) %gradient function of inequality
    global A;
    global b;
    F = A * x(2:3) + b - x(1);
    g = F.^-1;
    g = g * [-1;sum(A, 1)'];
end

function x = backtrack_line_search_stage1(x0, df, dx, alpha, beta, tt)
    global A;
    global b;
    t = 1;
    while (((A * (x0(2:3) + t * dx(2:3)) + b - (x0(1) + t * dx(1))) < 0) || f_stage1(x0 + t * dx, tt) > f_stage1(x0, tt) + alpha * t * df' * dx)
        t = t * beta;
    end;
    x = x0 + t * dx;
end

function P = Hession_stage1(x, t)
    P = zeros(3,3);
    P = t * P - Hessionie_stage1(x);
end

function P = Hessionie_stage1(x) %Hession of inequality
    global A;
    global b;
    F = A*x(2:3) + b - x(1);
    df = fgie_stage1(x);
    P = -(df/F) * [-1;sum(A,1)']';
end


function x = newton(t,x0)
    x = x0;
    lambda = 1;
    cnt = 0;
    while(lambda > 1e-10)
        df = fg(x,t);
        ddf = Hession(x,t);
        P = ddf;
        dx = P \ -df;
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

