%% Vonalak módszere beépített ODE-megoldókkal (az 1. példa feladatán)
%
%   A térbeli diszkretizáció után  du/dt = A u  merev (stiff) rendszer.
%   Az ode45 (explicit Runge–Kutta) a stabilitás miatt rengeteg apró lépést
%   tesz; az ode15s (implicit, változó rendű BDF/NDF) a pontossághoz igazítja
%   a lépésközt, és - ha megadjuk a ritka Jacobi-mátrixot - nagyon gyors.

clear; close all;

alpha = 1;  L = 1;  Nx = 100;  T = 0.2;
dx = L/Nx;  x = (0:Nx)'*dx;  xi = x(2:end-1);  m = Nx - 1;
e  = ones(m, 1);
A  = alpha * spdiags([e -2*e e], -1:1, m, m) / dx^2;

u0    = sin(pi*xi) + 0.5*sin(3*pi*xi);
upont = @(x,t) exp(-pi^2*alpha*t)*sin(pi*x) + 0.5*exp(-9*pi^2*alpha*t)*sin(3*pi*x);
rhs   = @(t, u) A*u;

opt45 = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'Refine', 1);  % Refine=1: minden lépés egy kimenet
opt15 = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'Jacobian', A);

tic; [t45, u45] = ode45(rhs,  [0 T], u0, opt45);  ido45 = toc;
tic; [t15, u15] = ode15s(rhs, [0 T], u0, opt15);  ido15 = toc;

fprintf('ode45 : %6d lépés, %.3f s, hiba = %.2e\n', numel(t45)-1, ido45, ...
        max(abs(u45(end,:)' - upont(xi, T))));
fprintf('ode15s: %6d lépés, %.3f s, hiba = %.2e\n', numel(t15)-1, ido15, ...
        max(abs(u15(end,:)' - upont(xi, T))));

% A merevség oka: A sajátértékei -4 alpha/dx^2 sin^2(k pi dx/2), k = 1..m
lam = -4*alpha/dx^2 * sin((1:m)*pi*dx/2).^2;
fprintf('merevségi hányados |lambda_max/lambda_min| = %.0f\n', ...
        max(abs(lam))/min(abs(lam)));

figure(1); clf;
semilogy(t45(2:end), diff(t45), '.', t15(2:end), diff(t15), 'o');
xlabel('t'); ylabel('lépésköz \Delta t'); grid on;
legend('ode45 (explicit)', 'ode15s (implicit)', 'Location', 'best');

save('pelda1_ode15s_eredmenyek.mat', 't45', 't15', '-v7');
