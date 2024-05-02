function V = simpson2(x, y, g, al, omega, bound)
    n = length(x);
    vs = @(x, y, al, omega) sin(omega * (x * cos(al) + y * sin(al)));
    f = @(x, y, al, omega, g) vs(x, y, al, omega) * g;
    h = bound(2) - bound(1);
    s = f(x(1), y(1), al, omega, g(1)) + f(x(end), y(end), al, omega, g(end)); 
    
    for i = 3:2:(n-1)
        s = s + 4 * f(x(i), y(i), al, omega, g(i)); 
    end
    for i = 2:2:(n-1)
        s = s + 2 * f(x(i), y(i), al, omega, g(i));
    end
    
    V = h / 3 * s;
end
