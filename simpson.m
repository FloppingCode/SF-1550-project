function integral_value = simpson(x, y, g, a, omega, bound)
% Funktion för att beräkna integralen av en funktion över ett givet intervall
% med Simpsons regel.

    n = length(x);
    vc = @(x, y, a, omega) cos(omega * (x * cos(a) + y * sin(a)));
    f = @(x, y, a, omega, g) vc(x, y, a, omega) * g;
    h = bound(2) - bound(1);
    
    % Beräkning av de första och sista termerna i Simpsons regel
    s = f(x(1), y(1), a, omega, g(1)) + f(x(end), y(end), a, omega, g(end)); 
    
    % Loopa över de ojämna termera 
    for i = 3:2:(n-1)
        s = s + 4 * f(x(i), y(i), a, omega, g(i)); 
    end
    
    % Loopa över de jämna termerna
    for i = 2:2:(n-1)
        s = s + 2 * f(x(i), y(i), a, omega, g(i));
    end
        integral_value = h / 3 * s;
end
