function integral_value = simpson2(x, y, g, al, omega, bound)
% Funktion för att beräkna integralen av en funktion över ett givet intervall
% med Simpsons regel.

    n = length(x);
    vs = @(x, y, al, omega) sin(omega * (x * cos(al) + y * sin(al)));
    f = @(x, y, al, omega, g) vs(x, y, al, omega) * g;
    h = bound(2) - bound(1);
    s = f(x(1), y(1), al, omega, g(1)) + f(x(end), y(end), al, omega, g(end)); 
    
    % Loopa över de ojämna termerna
    for i = 3:2:(n-1)
        s = s + 4 * f(x(i), y(i), al, omega, g(i)); 
    end
    
    % Loopa över de jämna termerna
    for i = 2:2:(n-1)
        s = s + 2 * f(x(i), y(i), al, omega, g(i));
    end
    integral_value = h / 3 * s;
end
