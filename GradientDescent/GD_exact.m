% The funciton is f(x) = 0.5 * (x(1) ^ 2 + gamma * x(2) ^ 2)
function [] = GD_exact()

x0 = [10,1];
gamma = 10;
x = x0;
gap = f(x, gamma);
disp(gap);
y = [];
y = [y,gap];
cnt = 0;
while (gap > 1e-6 && cnt < 20)
    deltax = fg(x, gamma);
    x = exact_line_search(x, deltax, gamma);
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

function x = exact_line_search(x0, deltax, gamma)
    t = (-gamma * deltax(2) * x0(2) - deltax(1) * x0(1))/ (deltax(1)^2 + gamma * deltax(2)^2);
    x = x0 + t * deltax;
end

