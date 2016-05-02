% The funciton is f(x) = 0.5 * (x(1) ^ 2 + gamma * x(2) ^ 2)
function [] = GD_exact()

x0 = [10,1];
gamma = 3;
x = x0;
gap = f(x, gamma);
disp(gap);
y = [];
y = [y,gap];
cnt = 0;
while (gap > 1e-6 && cnt < 20)
    df = -fg(x, gamma);
    x = backtrack_line_search(x, df, gamma, 0.3, 0.5);
    gap = f(x,gamma); 
    disp(gap);
    y = [y,gap];
    cnt = cnt + 1;
end
end

function value = f(x,gamma)
    value = 0.5 * (x(1)^2 + gamma * x(2)^2);
end

function g = fg(x,gamma)
    g = [x(1), gamma * x(2)];
end

function x = backtrack_line_search(x0, df, gamma, alpha, beta)
    t = 1;
    while (f(x0 + t * df, gamma) > f(x0, gamma) + alpha * t * (- df * df'))%% dx = -df
        t = t * beta;
    end;
    x = x0 + t * df;
end

