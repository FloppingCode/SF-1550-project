function [x, y, a, it] = find_source(B, omega, eta)
g = B.un;
M = 15;

ic1_diff = central_diff(B.x, B.y, g, pi/2, omega, B.s, 1e-8);
ic2_diff = central_diff(B.x, B.y, g, 0, omega, B.s, 1e-8);
is1 = simpson2(B.x, B.y, g, pi/2, omega, B.s);
is2 = simpson2(B.x, B.y, g, 0, omega, B.s);
ic1 = simpson(B.x, B.y, g, 0, omega, B.s);

x0_improved = ic1_diff / (omega*is1);
y0_improved =  -1 * (ic2_diff / (omega*is2));
a0_improved = sqrt(ic1^2+is2^2) * 1/eta;

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

