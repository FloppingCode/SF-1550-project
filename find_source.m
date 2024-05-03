function [x, y, a, it, res, MSE] = find_source(B, omega, eta)
% Funktion som uppskattar position och styrka hos en ljudkälla med hjälp av Gauss-Newton-metoden.

[x0_improved, y0_improved, a0_improved] = generate_guesses(B, omega, eta);

M = 15;
aa_vals = linspace(0, 2*pi, M)';
ic_vals = zeros(M,1);
i = 1;
for al = linspace(0, 2*pi, M)
    ic_vals(i,1) = simpson(B.x, B.y, B.un, al, omega, B.s);
    aa_vals(i,1) = al;
    i = i + 1;
end

[x, y, a, it, res] = gaussnewton(eta, ic_vals, aa_vals, omega, x0_improved, y0_improved, a0_improved);

predicted_ic_vals = a * cos(omega * (x * cos(aa_vals) + y * sin(aa_vals)));
final_residual = ic_vals - predicted_ic_vals;

MSE = calculate_MSE(final_residual);
end

function MSE = calculate_MSE(residual)
% Funktion för att beräkna medelkvadratfelet (MSE) med hjälp av residualen.

n = length(residual);
squared_residual = residual .^ 2;
mean_squared_residual = mean(squared_residual);
MSE = sqrt(mean_squared_residual);
end
