function [x_tilde, y_tilde, a_tilde, it, residual] = gaussnewton(eta, ic_vals, aa_vals, omega, x0, y0, a0)
% Funktion som anpassar position och styrka genom att ta in startgissningar
% och iterera mha Gauss-Newton-metoden.

maxiter = 1000;
tol_x_y = 1e-8; % Tolerans för felet i x och y
a_range = max(ic_vals) - min(ic_vals); 
tol_a = 1e-8 * a_range; % Tolerans för felet i a
it = 0;

for i = 1:maxiter
    old = [omega*x0; omega*y0; eta*a0];
    cos_vals = cos(old(1) .* cos(aa_vals) + old(2) .* sin(aa_vals));
    I_c = old(3) .* cos_vals - ic_vals;
    J = jac(old(1), old(2), old(3), aa_vals);
    d = J \ (-I_c);

    new = old + d;
    it = it + 1;
    if norm(old(1:2) - new(1:2)) < tol_x_y && norm(old(3) - new(3)) < tol_a % Om tol är mindre än både gränserna för x, y och a
        break
    end
end

residual = norm(I_c);

x_tilde = new(1, 1) / omega;
y_tilde = new(2, 1) / omega;
a_tilde = new(3, 1) / eta;
end

function J = jac(x0, y0, a0, aa_vals)
% Funktion för att beräkna Jacobianen för Gauss-Newton-metoden

dx = (-sin(x0 .* cos(aa_vals) + y0 .* sin(aa_vals)) .* a0 .* cos(aa_vals))';
dy = (-sin(x0 .* cos(aa_vals) + y0 .* sin(aa_vals)) .* a0 .* sin(aa_vals))';
da = (cos(x0 .* cos(aa_vals) + y0 .* sin(aa_vals)))';

J = [dx', dy', da'];
end

