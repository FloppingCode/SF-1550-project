S0 = @(x, y) cos(20*(x.^2 + y.^2)) .* exp(-1000*(x.^2 + y.^2).^2);

% Steg 1: Beräkna eta mha trapetsregeln i 2D
x_min = -10; x_max = 10;
y_min = -10; y_max = 10;
num_points = 1000;
integrand = @(x, y) S0(x, y) .* cos(19 * x);
eta = trapets2d(integrand, x_min, x_max, y_min, y_max, num_points);

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

% Steg 3.1: Hitta konstanterna a˜, x˜0, y˜0 mha Gauss-Newton
params = [x0, y0, a];
x_tilde = x0 - 0.02;
x-tilde = y0 + 0.04;
a_tilde = a + 0.03;


%Fortsätt... 
