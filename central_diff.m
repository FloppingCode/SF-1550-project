function ic_diff = central_diff(x, y, g, al, omega, bound, h)
% Funktiom som uppskattar derivatan med hjälp av I_c med hjälp av
% centraldifferens
    ic_plus = simpson(x, y, g, al + h, omega, bound);
    ic_minus = simpson(x, y, g, al - h, omega, bound);
    ic_diff = (ic_plus - ic_minus) / (2*h);
end
