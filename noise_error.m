function [noiselevel, err_total, err_x, err_y, err_a] = noise_error(Bound, g, omega, eta, xs, ys, a, x_tilde, y_tilde, a_tilde)
% Funktion som beräknar felen, både totalt och för varje parameter, när
% brus läggs till 

    % Intervall för brusnivån
    noiselevel = linspace(0.01, 1);
    
    err_total = zeros(length(noiselevel), 1);
    err_x = zeros(length(noiselevel), 1);
    err_y = zeros(length(noiselevel), 1);
    err_a = zeros(length(noiselevel), 1);

    M = 15;
    aa_vals = linspace(0, 2*pi, M)';
    ic_vals = zeros(M, 1);

    for j = 1:length(noiselevel)
        % Lägg till brus till ljudkällan
        gnoise = g + max(abs(g)) * randn(size(g)) * noiselevel(j);

        % Beräkna Ic för brusig ljudkälla
        for k = 1:M
            ic_vals(k) = simpson(Bound.x, Bound.y, gnoise, aa_vals(k), omega, Bound.s);
        end

        [x0, y0, a0] = gaussnewton(eta, ic_vals, aa_vals, omega, x_tilde, y_tilde, a_tilde);

        err_x(j) = abs(x0 - xs);
        err_y(j) = abs(y0 - ys);
        err_a(j) = abs(a0 - a);
        err_total(j) = sqrt(err_x(j)^2 + err_y(j)^2 + err_a(j)^2);
    end
end
