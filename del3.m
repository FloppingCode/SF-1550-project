%% Uppgift 3.1: Beräkna η
format long;
omega = 19;

S0 = @(x,y) cos(20.*sqrt(x.^2+y.^2)).*exp(-1000.*(x.^2+y.^2));
vc = @(x,y,a) cos(omega*(x.*cos(a) + y.*sin(a)));
f2 = @(x, a) cos(a*x);
integrand = @(x,y,a) S0(x,y).*f2(x, a);

for a = linspace(0, 2*pi, 5) 
    eta = trapets2d(integrand, -0.5,0.5,-0.5,0.5, 500, a);
    disp(eta)
end

%% Uppgift 3.2: Beräkna g

N = 500;
xs  = 0.6;
ys  = 0.2;
a = 10;
omega = 19;

S = @(x,y) a * S0(x-xs,y-ys);

[Bound,Sol]=hhsolver(omega,S,N); 
g = Bound.un;

%disp(g)

% 3.3a: Beräkna Ic mha numerisk integration
% Definiera funktioner
S0 = @(x,y) cos(20.*sqrt(x.^2+y.^2)).*exp(-1000.*(x.^2+y.^2));
vc = @(x,y,a) cos(omega*(x.*cos(a) + y.*sin(a)));
f2 = @(x, a) cos(a*x);

integrand = @(x, y, a) f2(x, a) .* S0(x, y);

M = 13;
aa_vals = linspace(0, 2*pi, M)';
ic_vals = zeros(M, 1);
for i = 1:M
    ic_vals(i) = simpson(Bound.x, Bound.y, g, aa_vals(i), omega, Bound.s);
end

eta = trapets2d(integrand, -0.05, 0.05, -0.05, 0.05, 700, 10);

% Anpassa med Gauss-Newton
x0 = xs + 0.03;
y0 = ys - 0.02;
a0 = a - 0.04;
correct_values = [xs, ys, a];

eta = trapets2d(integrand, -0.05, 0.05, -0.05, 0.05, 500, a);

[x_tilde, y_tilde, a_tilde, iterations] = gaussnewton(eta, ic_vals, aa_vals, omega, x0, y0, a0);
result = [x_tilde, y_tilde, a_tilde];

fprintf("x: %.4f\n", x_tilde)
fprintf("y: %.4f\n", y_tilde)
fprintf("a: %.4f\n", a_tilde)

% Uppgift 3.4: Brus
M = 15;
noiselevel = linspace(0.01, 1);
eta = trapets2d(integrand, -0.05, 0.05, -0.05, 0.05, 500, a);
err_total = zeros(length(noiselevel), 1);
err_x = zeros(length(noiselevel), 1);
err_y = zeros(length(noiselevel), 1);
err_a = zeros(length(noiselevel), 1);

aa_vals = linspace(0, 2*pi, M)';
ic_vals = zeros(M, 1);

for j = 1:length(noiselevel)
    gnoise = g + max(abs(g)) * randn(size(g)) * noiselevel(j);
    
    for k = 1:M
        ic_vals(k) = simpson(Bound.x, Bound.y, gnoise, aa_vals(k), omega, Bound.s);
    end
    
    [x0, y0, a0] = gaussnewton(eta, ic_vals, aa_vals, omega, x_tilde, y_tilde, a_tilde);
    
    err_x(j) = abs(x0 - xs);
    err_y(j) = abs(y0 - ys);
    err_a(j) = abs(a0 - a);
    
    err_total(j) = sqrt(err_x(j)^2 + err_y(j)^2 + err_a(j)^2);
end

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

% Uppgift 3.5: 
omega = 19;

ic1_diff = central_diff(Bound.x, Bound.y, g, pi/2, omega, Bound.s, 1e-8);
ic2_diff = central_diff(Bound.x, Bound.y, g, 0, omega, Bound.s, 1e-5);
is1 = simpson2(Bound.x, Bound.y, g, pi/2, omega, Bound.s);
is2 = simpson2(Bound.x, Bound.y, g, 0, omega, Bound.s);
ic1 = simpson(Bound.x, Bound.y, g, 0, omega, Bound.s);

x0_improved = ic1_diff / (omega*is1);
y0_improved =  -1 * (ic2_diff / (omega*is2));
a0_improved = sqrt(ic1^2+is2^2) * 1/eta;
fprintf("Bättre x0: %f\n", x0_improved)
fprintf("Bättre y0: %f\n", y0_improved)
fprintf("Bättre a: %f\n", a0_improved)
