%% Uppgift 3.1: Beräkna η
fprintf("Uppgift 1: Beräkna η\n")

format long;
omega = 19;

S0 = @(x,y) cos(20.*sqrt(x.^2+y.^2)).*exp(-1000.*(x.^2+y.^2));
f2 = @(x) cos(omega*x);
integrand = @(x,y,al) S0(x,y).*f2(x);

% for al = linspace(0, 2*pi, 5) 
%     eta = trapets2d(integrand, -0.5, 0.5, -0.5, 0.5, 500, al);
% end

eta = trapets2d(integrand, -0.5, 0.5, -0.5, 0.5, 700)

disp("η: " + eta)

%% Uppgift 3.2: Beräkna g
fprintf("Uppgift 2: Beräkna g\n")

N = 500;
xs  = 0.6;
ys  = 0.2;
a = 9;
omega = 19;

S0 = @(x,y) cos(20.*sqrt(x.^2+y.^2)).*exp(-1000.*(x.^2+y.^2));
S = @(x,y) a * S0(x-xs,y-ys);

[Bound,Sol]=hhsolver(omega,S,N); 
g = Bound.un;
%% 3.3: Beräkna och anpassa Ic
fprintf("Uppgift 3: Beräkna och anpassa Ic\n")

% Funktioner
S0 = @(x, y) cos(20.*sqrt(x.^2+y.^2)).*exp(-1000.*(x.^2+y.^2));
f2 = @(x) cos(omega*x);
integrand = @(x,y,al) S0(x,y).* f2(x);

M = 15;
aa_vals = linspace(0, 2*pi, M)';
ic_vals = zeros(M, 1);

% Beräkna ic_vals
for i = 1:M
    ic_vals(i) = simpson(Bound.x, Bound.y, g, aa_vals(i), omega, Bound.s);
end

eta = trapets2d(integrand, 0, 1, 0, 0.75, 700);

% 3.3b Anpassa med Gauss-Newton
x0 = xs + 0.03;
y0 = ys - 0.02;
a0 = a - 0.01;
correct_values = [xs, ys, a];

eta = trapets2d(integrand, -0.5, 0.5, -0.5, 0.5, 500);

[x_tilde, y_tilde, a_tilde, iterations, res] = gaussnewton(eta, ic_vals, aa_vals, omega, x0, y0, a0);

% Skriv ut anpassade värden
fprintf("Anpassat x: %.4f\n", x_tilde)
fprintf("Anpassat y: %.4f\n", y_tilde)
fprintf("Anpassat a: %.4f\n", a_tilde)
fprintf("\n")

% Verifiering
ic_real = eta * a * cos(omega * (xs * cos(aa_vals) + ys * sin(aa_vals)));

figure;
plot(aa_vals, ic_vals, '-o', aa_vals, ic_real, '--')
title("Verifiering")
legend("Anpassade värden", "Förväntade värden")
xlabel("\alpha")
ylabel("VL/HL")
%% Uppgift 3.4: Brus

[noiselevel, err_total, err_x, err_y, err_a] = noise_error(Bound, g, omega, eta, x0, y0, a0, x_tilde, y_tilde, a_tilde);

figure;
subplot(2, 1, 1);
plot(noiselevel, err_total, '-');
xlabel('Noiselevel');
ylabel('Totalt fel (euk.)');
title('Totalt fel');

subplot(2, 1, 2);
plot(noiselevel, err_x, '-', noiselevel, err_y, '--', noiselevel, err_a, '-.');
legend('Fel i x', 'Fel i y', 'Fel i a', 'Location', 'best');
xlabel('Noiselevel');
ylabel('Fel');
title('Felet i varje parameter');


%% Uppgift 3.5: 
omega = 19;
[x0_improved, y0_improved, a0_improved] = generate_guesses(Bound, omega, eta);
fprintf("Bättre x0: %f\n", x0_improved)
fprintf("Bättre y0: %f\n", y0_improved)
fprintf("Bättre a0: %f\n", a0_improved)

[noiselevel, err_total, err_x, err_y, err_a] = noise_error(Bound, g, omega, eta, x0_improved, y0_improved, a0_improved, x_tilde, y_tilde, a_tilde);

figure;
subplot(2, 1, 1);
plot(noiselevel, err_total, '-');
xlabel('Noiselevel');
ylabel('Totalt fel (euk.)');
title('Totalt fel');

subplot(2, 1, 2);
plot(noiselevel, err_x, '-', noiselevel, err_y, '--', noiselevel, err_a, '-.');
legend('Fel i x', 'Fel i y', 'Fel i a', 'Location', 'best');
xlabel('Noiselevel');
ylabel('Fel');
title('Felet i varje parameter');


%% 3.6 Hitta källor
fprintf("Uppgift 3.6: Hitta källors posisition och styrka\n")
format long;

eta = 0.002378611784500; % Hämtat från 3.1

for i = 1:2
    fprintf("Källa %d\n", i);
    
    filename = sprintf("source%d.mat", i);
    load(filename, "B", "omega");
    
    [x, y, a, counter, res] = find_source(B, omega, eta);
    
    fprintf("x: %.10f\n", x)
    fprintf("y: %.10f\n", y)
    fprintf("a: %.10f\n", a)
    fprintf("res: %.10f\n", res)
    fprintf("\n")
end
