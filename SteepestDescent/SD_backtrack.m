% The funciton is f(x) = 0.5 * (x(1) ^ 2 + gamma * x(2) ^ 2)
function [] = SD_backtrack()

x0 = [10,1];
gamma = 10;
x = x0;
gap = f(x, gamma);
disp(gap);
P = [1,0;0,100];
cnt = 0;
while (gap > 1e-6 && cnt < 20)
    df = fg(x, gamma);
    dx = -(P^-1 * df')'; 
    x = backtrack_line_search(x, df, dx, gamma, 0.3, 0.5);
    gap = f(x,gamma); 
    disp(gap);
    cnt = cnt + 1;
end
end

function value = f(x,gamma)
    value = 0.5 * (x(1)^2 + gamma * x(2)^2);
end

function g = fg(x,gamma)
    g = [x(1), gamma * x(2)];
end

function x = backtrack_line_search(x0, df, dx, gamma, alpha, beta)
    t = 1;
    while (f(x0 + t * dx, gamma) > f(x0, gamma) + alpha * t * df * dx')
        t = t * beta;
    end;
    x = x0 + t * dx;
end

