function I_c = simpson_Ic(x, y, g, aa, w, bound)
    n = length(x);
    
    vc = @(x, y, aa, w) cos(w * (x * cos(aa) + y * sin(aa)));
    f = @(x, y, aa, w, g) vc(x, y, aa, w) .* g;
    
    h = bound(2) - bound(1);
    
    s = f(x(1), y(1), aa, w, g(1)) + f(x(end), y(end), aa, w, g(end)); 
    for i = 3:2:(n-1)
        s = s + 4 * f(x(i), y(i), aa, w, g(i)); 
    end
    for i = 2:2:(n-1)
        s = s + 2 * f(x(i), y(i), aa, w, g(i)); 
    end
    I_c = h / 3 * s;
end

