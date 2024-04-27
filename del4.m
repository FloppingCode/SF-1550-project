%% grovlokaliser 4.
format long
xs = 0.6;
% arrayfunc stuff for terseness
runs = 200
xs_arr = linspace(0.6, 1.0, runs);
A_arr = arrayfun(@(xs) A(xs,0.6,100), xs_arr)
% plotting
% verkar vara runt 0.77-0.78
figure(1)
plot(xs_arr,A_arr,LineStyle="--")
%% bestäm optimala placeringen i en variabel
N = 100; % välja N
x_opt = gyllnesnitt(@(xs) A(xs,0.6,N),0.74,0.8,10^(-4))
A(x_opt,0.6,N)

% N = 100 gav:  
% x_opt: 0.785871508588702
% A(x_opt): 0.067753585710342
% N = 1000 gav: 
% x_opt: 0.777117626563683
% A(x_opt): 0.022311382629263
% N = 2000 gav:
% x_opt: 0.777117626563683 (gav exakt samma? lite konstigt)
% A(x_opt): 0.022346745275933

%% visualisering av bra vs dålig Tv-placering
N = 200;
xs = x_opt;
ys = 0.6;
 % bra placering
S0 = @(x, y) cos(20*(x.^2 + y.^2)) .* exp(-1000*(x.^2 + y.^2).^2);
S = @(x,y) aa*S0(x-xs,y-ys);
[Bound,Sol] = hhsolver(omega,S,N);
%surf(Sol.x,Sol.y,Sol.S)
figure(2)
mesh(Sol.x,Sol.y,Sol.u)

xs = 0.6 % x dålig placering
figure(3)

S0 = @(x, y) cos(20*(x.^2 + y.^2)) .* exp(-1000*(x.^2 + y.^2).^2);
S = @(x,y) aa*S0(x-xs,y-ys);
[Bound,Sol] = hhsolver(omega,S,N);
mesh(Sol.x,Sol.y,Sol.u)

%% Nelson-Meads part 1

% vi börjar med grovlokalisering
N = 100;
% generate 2d equidis point grid
runs = 50;
x = linspace(0.6,1.0,runs);
y = linspace(0.6,1.0,runs);
[X, Y] = meshgrid(x,y);
Z = arrayfun(@(xs,ys) A(xs,ys,N),X,Y);
%plotting
figure(2)
surf(X,Y,Z)
% grovskattning, kan kanske välja en bättre
% x = 0.828571, y = 0.828571, z = 0.0950528

%% Nealson-Meads part 2
% välj N
N = 200;
X0 = [0.828571, 0.828571]; % bra
option = optimset('Display','iter');
[Xopt,fval] = fminsearch(@(x) A(x(1),x(2),N), X0, option);

%% visualisering av placering
xs = Xopt(1);
ys = Xopt(2);
S0 = @(x, y) cos(20*(x.^2 + y.^2)) .* exp(-1000*(x.^2 + y.^2).^2);
S = @(x,y) aa*S0(x-xs,y-ys);
[Bound,Sol] = hhsolver(omega,S,N);
figure(4)
mesh(Sol.x,Sol.y,Sol.u)

% lite oskäker om det här räcker, vi hittade en som uppfyller toleransen.
% Vi kan inte få så mycket lägre ljudstyrka

%Prova med olika startgissningar till du hittar ett minimum inne
%i rummet. Hur mycket lägre ljudstyrka kan du uppnå? Visualisera resultatet!

%% funktioner
function kvotLjud = A(xs,ys,N)
    % konstanter
    omega = 30;
    aa = 1.0;
    % funktionen
    S0 = @(x, y) cos(20*(x.^2 + y.^2)) .* exp(-1000*(x.^2 + y.^2).^2);
    S = @(x,y) aa*S0(x-xs,y-ys);
    [Bound,Sol] = hhsolver(omega,S,N);
    w = find(Sol.x <= 0.25 & Sol.y>=0.5);
    kvotLjud = max(abs(Sol.u(w)))/max(abs(Sol.u(:)));
end