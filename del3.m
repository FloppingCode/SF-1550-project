S0 = @(x, y) cos(20*(x.^2 + y.^2)) .* exp(-1000*(x.^2 + y.^2).^2);

% Steg 1: Beräkna eta mha trapetsregeln i 2D
format long; 
num_points = 500;
M = 10;
alpha_list = linspace(0, 2*pi, M);
eta_values = zeros(M, 1);
x_min = -1; x_max = 1;
y_min = -1; y_max = 1;
omega = 19;

for i = 1:M
    alpha = alpha_list(i);
    integrand = @(x, y) S0(x, y) .* cos(omega * (x * cos(alpha) + y * sin(alpha)));
    eta_values(i) = trapets2d(integrand, x_min, x_max, y_min, y_max, num_points);
end

%disp(eta_values);

% Steg 2: Beräkna g mha hhsolver
omega = 19;
N = 400;
a = 10; 
x0 = 0.5; 
y0 = 0.2;

S = @(x, y) a * S0(x - x0, y - y0);
[B, Sol] = hhsolver(omega, S, N);
g = B.un;

% Steg 3.1: Beräkna I_c mha Simpson
M = 10;
alpha_list = linspace(0, 2*pi, M)';
Ic_list = zeros(M, 1);

for i = 1:M
    Ic_list(i) = simpson_ic(B.x, B.y, g, alpha_list(i), omega, B.s);
end


%disp(Ic_values);

% Steg 3.2: Hitta konstanterna a˜, x˜0, y˜0 mha Gauss-Newton
params = [x0, y0, a];
x_tilde = x0 - 0.02;
y_tilde = y0 + 0.04;
a_tilde = a + 0.03;

v_c = @(x,y, alpha) cos(w*(x.*cos(alpha) + y.*sin(alpha)));

[x, y, a] = gaussnewton(omega, eta_values, Ic_list, alpha_list, x_tilde, y_tilde, a_tilde);

fprintf('a = ', a);
fprintf('x = ', x_0);
fprintf('y = ', y_0);

%Fortsätt... 
