% The funciton is f(x) = 0.5 * (x(1) ^ 2 + gamma * x(2) ^ 2)
function [] = GD_exact()

x0 = [10,1];
gamma = 10;
x = x0;
gap = f(x, gamma);
disp(gap);
cnt = 0;
while (gap > 1e-6 && cnt < 20)
    df = fg(x, gamma);
    x = exact_line_search(x, -df, gamma);
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

function x = exact_line_search(x0, dx, gamma)
    t = (-gamma * dx(2) * x0(2) - dx(1) * x0(1))/ (dx(1)^2 + gamma * dx(2)^2);
    x = x0 + t * dx;
end

