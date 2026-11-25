%% CSAPDA: a Crank–Nicolson-módszer oszcillációja nem sima kezdeti adatnál
%
%   u_t = u_xx,  u(0,t) = u(1,t) = 0,  u(x,0) = 1  (0 < x < 1)
%   A kezdeti és a peremfeltétel x = 0-ban és x = 1-ben nem illeszkedik
%   (szakadás), így a megoldásban a nagy hullámszámú módusok is jelen vannak.
%
%   Pontos megoldás:  u(x,t) = 4/pi * sum_{k páratlan} sin(k pi x)/k * exp(-k^2 pi^2 t)
%
%   Nagy időlépésnél (r >> 1) a CN-módszer ezeket a módusokat alig csillapítja,
%   és előjelet vált minden lépésben -> fűrészfog-szerű oszcilláció.
%   Gyógymód: implicit Euler-indítás (Rannacher), vagy kisebb dt.

clear; close all;

Nx = 50;  dt = 0.01;  T = 0.05;        % r = dt/dx^2 = 25
nulla = @(t) 0;
p = struct('alpha', 1, 'L', 1, 'Nx', Nx, 'dt', dt, 'T', T, ...
           'u0', @(x) ones(size(x)), 'gL', nulla, 'gR', nulla, 'tsave', T);

p.theta = 0.5;  p.nBE = 0;  [x, ~, uCN]  = hovez1d_theta(p);   % Crank–Nicolson
p.theta = 1;    p.nBE = 0;  [~, ~, uBE]  = hovez1d_theta(p);   % implicit Euler
p.theta = 0.5;  p.nBE = 2;  [~, ~, uRAN] = hovez1d_theta(p);   % CN + Rannacher

% pontos megoldás (a Fourier-sor 2000 tagjával)
xf = linspace(0, 1, 500)';
upont = zeros(size(xf));
for k = 1:2:3999
    upont = upont + 4/pi * sin(k*pi*xf)/k * exp(-k^2*pi^2*T);
end

figure(1); clf;
plot(xf, upont, 'k-', x, uCN, 'o-', x, uBE, 's-', x, uRAN, 'd-');
legend('pontos', 'Crank–Nicolson', 'implicit Euler', 'CN + Rannacher-indítás', ...
       'Location', 'south');
xlabel('x'); ylabel('u(x, 0.05)'); grid on;
title(sprintf('r = \\Delta t/\\Delta x^2 = %.0f', dt*Nx^2));

% a hibák a rácspontokban
up = zeros(size(x));
for k = 1:2:3999
    up = up + 4/pi * sin(k*pi*x)/k * exp(-k^2*pi^2*T);
end
fprintf('max. hiba:  CN = %.3e,  BE = %.3e,  CN+Rannacher = %.3e\n', ...
        max(abs(uCN-up)), max(abs(uBE-up)), max(abs(uRAN-up)));

save('cn_oszcillacio_eredmenyek.mat', 'x', 'uCN', 'uBE', 'uRAN', 'xf', 'upont', '-v7');
