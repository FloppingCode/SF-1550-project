function [x, y, a, it] = find_source(B, omega, eta)
g = B.un;
M = 15;

[x0_improved, y0_improved, a0_improved] = generate_guesses(B, omega, eta);


aa_vals = zeros(M,1);
ic_vals = zeros(M,1);
i = 1;
for al = linspace (0, 2*pi, M)
    ic_vals(i,1) = simpson(B.x, B.y, g, al, omega, B.s);
    aa_vals(i,1)= al;
    i = i + 1;
end

[x, y, a, it] = gaussnewton(eta, ic_vals, aa_vals, omega, x0_improved, y0_improved, a0_improved);

end
