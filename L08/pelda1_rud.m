%% 1. PÉLDA: hővezetés egy vékony rúdban, ismert analitikus megoldással
%
%   u_t = alpha u_xx,  0 < x < 1,  u(0,t) = u(1,t) = 0,
%   u(x,0) = sin(pi x) + 0.5 sin(3 pi x)
%
%   Pontos megoldás (Fourier-sor két tagja):
%   u(x,t) = exp(-pi^2 alpha t) sin(pi x) + 0.5 exp(-9 pi^2 alpha t) sin(3 pi x)
%
%   A szkript bemutatja:
%     (a) a Crank–Nicolson-megoldást és összevetését a pontos megoldással,
%     (b) az időbeli konvergenciarendet (implicit Euler: 1, CN: 2),
%     (c) a térbeli konvergenciarendet (2),
%     (d) az explicit módszer instabilitását r > 1/2 esetén.

clear; close all;

alpha = 1;  L = 1;
u0    = @(x) sin(pi*x) + 0.5*sin(3*pi*x);
upont = @(x,t) exp(-pi^2*alpha*t)*sin(pi*x) + 0.5*exp(-9*pi^2*alpha*t)*sin(3*pi*x);
nulla = @(t) 0;

%% (a) Crank–Nicolson-megoldás néhány időpontban
p = struct('alpha', alpha, 'L', L, 'Nx', 50, 'dt', 0.005, 'T', 0.2, ...
           'u0', u0, 'gL', nulla, 'gR', nulla, 'theta', 0.5, ...
           'tsave', [0 0.01 0.03 0.1 0.2]);
[x, ts, U] = hovez1d_theta(p);

xf = linspace(0, L, 400)';
figure(1); clf; hold on;
for k = 1:numel(ts)
    plot(xf, upont(xf, ts(k)), 'k-');
    plot(x, U(:,k), 'o', 'MarkerSize', 4);
end
xlabel('x'); ylabel('u(x,t)'); grid on;
title('1. példa: Crank–Nicolson (pontok) és pontos megoldás (vonal)');
pelda1a = struct('x', x, 'ts', ts, 'U', U);

%% (b) Időbeli konvergencia: finom térbeli rács, változó dt
Tv  = 0.1;
dts = 0.02 ./ 2.^(0:5);
hibaBE = zeros(size(dts));  hibaCN = zeros(size(dts));
p.Nx = 1000;  p.T = Tv;  p.tsave = Tv;
for i = 1:numel(dts)
    p.dt = dts(i);
    p.theta = 1;    [x, ~, U] = hovez1d_theta(p);
    hibaBE(i) = max(abs(U - upont(x, Tv)));
    p.theta = 0.5;  [x, ~, U] = hovez1d_theta(p);
    hibaCN(i) = max(abs(U - upont(x, Tv)));
end
rendBE = log2(hibaBE(1:end-1) ./ hibaBE(2:end));
rendCN = log2(hibaCN(1:end-1) ./ hibaCN(2:end));
fprintf('\nIdőbeli konvergencia (t = %.2f, Nx = %d)\n', Tv, p.Nx);
fprintf('   dt         hiba BE     rend    hiba CN     rend\n');
for i = 1:numel(dts)
    if i == 1
        fprintf('%9.6f   %9.3e      -    %9.3e      -\n', dts(i), hibaBE(i), hibaCN(i));
    else
        fprintf('%9.6f   %9.3e   %5.2f   %9.3e   %5.2f\n', dts(i), ...
                hibaBE(i), rendBE(i-1), hibaCN(i), rendCN(i-1));
    end
end

figure(2); clf;
loglog(dts, hibaBE, 's-', dts, hibaCN, 'o-', ...
       dts, 0.5*hibaBE(end)*dts/dts(end), 'k--', ...
       dts, 0.5*hibaCN(end)*(dts/dts(end)).^2, 'k:');
xlabel('\Delta t'); ylabel('max. hiba'); grid on;
legend('implicit Euler', 'Crank–Nicolson', 'O(\Delta t)', 'O(\Delta t^2)', ...
       'Location', 'southeast');

%% (c) Térbeli konvergencia: nagyon kicsi dt, változó Nx
Nxs = [10 20 40 80 160];
hibaX = zeros(size(Nxs));
p.theta = 0.5;  p.dt = 1e-5;
for i = 1:numel(Nxs)
    p.Nx = Nxs(i);
    [x, ~, U] = hovez1d_theta(p);
    hibaX(i) = max(abs(U - upont(x, Tv)));
end
dxs = L ./ Nxs;
fprintf('\nTérbeli konvergencia (CN, dt = %g)\n', p.dt);
fprintf('   dx        hiba        rend\n');
fprintf('%8.5f   %9.3e      -\n', dxs(1), hibaX(1));
for i = 2:numel(Nxs)
    fprintf('%8.5f   %9.3e   %5.2f\n', dxs(i), hibaX(i), log2(hibaX(i-1)/hibaX(i)));
end

%% (d) Az explicit módszer stabilitása: r = 0.48 és r = 0.60
Nx = 20;  dx = L/Nx;  Tst = 0.2;
rs = [0.48 0.60];
pelda1d = struct('r', {}, 'x', {}, 'U', {}, 't', {});
figure(3); clf;
for i = 1:2
    q = struct('alpha', alpha, 'L', L, 'Nx', Nx, 'dt', rs(i)*dx^2/alpha, ...
               'T', Tst, 'u0', u0, 'gL', nulla, 'gR', nulla, 'theta', 0);
    q.tsave = round(Tst/q.dt) * q.dt;      % az utolsó lépés
    [x, ts, U] = hovez1d_theta(q);
    subplot(1, 2, i);
    plot(xf, upont(xf, ts), 'k-', x, U, 'o-');
    title(sprintf('FTCS, r = %.2f, t = %.3f', rs(i), ts));
    xlabel('x'); ylabel('u'); grid on;
    pelda1d(i).r = rs(i);  pelda1d(i).x = x;  pelda1d(i).U = U;  pelda1d(i).t = ts;
end

%% eredmények mentése (a dokumentum ábráihoz)
save('pelda1_eredmenyek.mat', 'pelda1a', 'dts', 'hibaBE', 'hibaCN', ...
     'Nxs', 'hibaX', 'pelda1d', '-v7');
