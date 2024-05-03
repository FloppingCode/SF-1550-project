%% Uppgift 3.1: Beräkna η
format long;
omega = 19;

S0 = @(x,y) cos(20.*sqrt(x.^2+y.^2)).*exp(-1000.*(x.^2+y.^2));
vc = @(x,y, al) cos(omega*(x.*cos(al) + y.*sin(al)));
f2 = @(x) cos(omega*x);
integrand = @(x,y,al) S0(x,y).*vc(x, y, al);

for al = linspace(0, 2*pi, 5) 
    eta = trapets2d(integrand, -0.5, 0.5, -0.5, 0.5, 500, al);
    disp(eta)
end


%% Uppgift 3.2: Beräkna g

N = 500;
xs  = 0.6;
ys  = 0.2;
a = 9;
omega = 19;

S0 = @(x,y) cos(20.*sqrt(x.^2+y.^2)).*exp(-1000.*(x.^2+y.^2));
S = @(x,y) a * S0(x-xs,y-ys);


[Bound,Sol]=hhsolver(omega,S,N); 
g = Bound.un;


%disp(g)

% 3.3a: Beräkna Ic mha numerisk integration
S0 = @(x,y) cos(20.*sqrt(x.^2+y.^2)).*exp(-1000.*(x.^2+y.^2));
vc = @(x,y,al) cos(omega*(x.*cos(al) + y.*sin(al)));
f2 = @(x) cos(omega*x);

integrand = @(x, y, al) vc(x, y, al) .* S0(x, y);

M = 13;
aa_vals = linspace(0, 2*pi, M)';
ic_vals = zeros(M, 1);
for i = 1:M
    ic_vals(i) = simpson(Bound.x, Bound.y, g, aa_vals(i), omega, Bound.s);
end

al = 2*pi;
eta = trapets2d(integrand, 0, 1, 0, 0.75, 700, al);

% Anpassa med Gauss-Newton
x0 = xs + 0.03;
y0 = ys - 0.02;
a0 = a - 0.01;
correct_values = [xs, ys, a];

eta = trapets2d(integrand, -0.5, 0.5, -0.5, 0.5, 500, a);
disp(eta)

[x_tilde, y_tilde, a_tilde, iterations] = gaussnewton(eta, ic_vals, aa_vals, omega, x0, y0, a0);
result = [x_tilde, y_tilde, a_tilde];

fprintf("x: %.4f\n", x_tilde)
fprintf("y: %.4f\n", y_tilde)
fprintf("a: %.4f\n", a_tilde)

ic_real = eta * a * cos(omega * (xs * cos(aa_vals) + ys * sin(aa_vals)));

figure;
plot(aa_vals, ic_vals, '-o', aa_vals, ic_real, '--')
title("Verifiering")
legend("Anpassade värden", "Förväntade värden")
xlabel("0 < \alpha < 2\pi")
ylabel("VL/HL")

% Uppgift 3.4: Brus

noise_error(Bound, g, omega, eta, x0, y0, a0, x_tilde, y_tilde, a_tilde);

% Uppgift 3.5: 
omega = 19;
[x0_improved, y0_improved, a0_improved] = generate_guesses(Bound, omega, eta);
fprintf("Bättre x0: %f\n", x0_improved)
fprintf("Bättre y0: %f\n", y0_improved)
fprintf("Bättre a0: %f\n", a0_improved)

noise_error(Bound, g, omega, eta, x0_improved, y0_improved, a0_improved, x_tilde, y_tilde, a_tilde);


%% 3.6 Hitta källor
eta = 0.00237;

integrand = @(x,y,al) S0(x,y).*f2(x, al);
for i = 1:5
    fprintf("Källa %d\n", i);
    
    filename = sprintf("source%d.mat", i);
    load(filename, "B", "omega");
    
    [x, y, a, counter] = find_source(B, omega, eta);
    
    fprintf("x: %f\n", x)
    fprintf("y: %f\n", y)
    fprintf("a: %f\n", a)
    fprintf("\n")
end
