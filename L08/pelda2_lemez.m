%% 2. PÉLDA: egy alumíniumlemez felmelegedése (2D, vegyes peremfeltételek)
%
%   10 cm x 10 cm-es alumíniumlemez, kezdetben 20 °C, egy forró folttal.
%   Bal perem:   100 °C (fűtött),        jobb perem: 20 °C (hűtött)
%   Alsó és felső perem: hőszigetelt (du/dy = 0)
%
%   A szkript:
%     (a) Crank–Nicolson (+ Rannacher-indítás) megoldás, pillanatképek,
%     (b) a stacionárius állapot ellenőrzése (lineáris profil x-ben),
%     (c) futásidő- és pontosságösszehasonlítás öt módszerrel.

clear; close all;

% --- fizikai paraméterek -------------------------------------------------
p.alpha = 9.7e-5;          % hődiffúziós együttható [m^2/s] (alumínium)
p.Lx = 0.10;  p.Ly = 0.10; % méretek [m]
p.TL = 100;   p.TR = 20;   % peremhőmérsékletek [°C]

% --- numerikus paraméterek -----------------------------------------------
p.Nx = 100;   p.Ny = 100;  % 1 mm-es rács -> 99 x 101 = 9999 ismeretlen
op = hovez2d_operatorok(p);
[X, Y] = ndgrid(op.xi, op.y);

% kezdeti feltétel: 20 °C + Gauss-alakú forró folt (150 °C csúcs)
U0 = 20 + 130 * exp(-((X - 0.065).^2 + (Y - 0.035).^2) / (2 * 0.008^2));

tau = p.Lx^2 / p.alpha;    % jellemző diffúziós idő
fprintf('Jellemző idő: L^2/alpha = %.0f s,  explicit dt_max = %.2e s\n', ...
        tau, op.dtmax);

% a teljes lemez (a Dirichlet-peremekkel együtt) a rajzoláshoz
teljes = @(U) [op.TL * ones(1, op.my); U; op.TR * ones(1, op.my)];

%% (a) Crank–Nicolson + Rannacher, dt = 0.5 s, 300 s-ig
dt = 0.5;  T = 300;  tkep = [0 1 4 15 60];
tsave = union([0:dt:20, 20:2:T], tkep);   % az elején sűrűbben mentünk
tic;
[ts, Us] = hovez2d_megold(op, U0, 'CN_LU', dt, T, tsave, 2);
fprintf('CN-megoldás (%d lépés): %.2f s\n', round(T/dt), toc);


figure(1); clf;
for i = 1:numel(tkep)
    k = find(abs(ts - tkep(i)) < 1e-9);
    subplot(2, 3, i);
    contourf(100*op.x, 100*op.y, teljes(Us(:,:,k)).', 20, 'LineColor', 'none');
    axis equal tight; colorbar; caxis([20 150]);
    title(sprintf('t = %g s', ts(k)));  xlabel('x [cm]'); ylabel('y [cm]');
end

% hőmérséklet két szondapontban az idő függvényében
[~, i1] = min(abs(op.xi - 0.065));  [~, j1] = min(abs(op.y - 0.035));
[~, i2] = min(abs(op.xi - 0.050));  [~, j2] = min(abs(op.y - 0.050));
szonda1 = squeeze(Us(i1, j1, :));   szonda2 = squeeze(Us(i2, j2, :));
subplot(2, 3, 6);
plot(ts, szonda1, ts, szonda2);  grid on;
xlabel('t [s]'); ylabel('T [°C]'); legend('folt közepe', 'lemez közepe');

%% (b) stacionárius állapot: u(x) = TL + (TR - TL) x / Lx
Ust = teljes(Us(:,:,end));
ustac = p.TL + (p.TR - p.TL) * op.x / p.Lx;
fprintf('t = %g s: max |u - u_stac| = %.2e °C\n', T, ...
        max(max(abs(Ust - ustac))));

%% (c) Módszerek összehasonlítása t = 60 s-ig
T = 60;
fprintf('\nReferenciamegoldás (CN, dt = 0.02 s) számítása...\n');
[~, Uref] = hovez2d_megold(op, U0, 'CN_LU', 0.02, T, T, 2);

esetek = { 'explicit',  0.95*op.dtmax;
           'BE_perjel', 0.5;
           'BE_LU',     0.5;
           'CN_LU',     0.5;
           'ADI',       0.5;
           'BE_LU',     2.0;
           'CN_LU',     2.0;
           'ADI',       2.0 };
ne = size(esetek, 1);
ido = zeros(ne, 1);  hiba = zeros(ne, 1);  lepes = zeros(ne, 1);
fprintf('\n%-10s %10s %8s %12s %14s\n', 'módszer', 'dt [s]', 'lépés', ...
        'idő [s]', 'max hiba [°C]');
for i = 1:ne
    m = esetek{i, 1};  dti = esetek{i, 2};
    tic;
    [~, Ui] = hovez2d_megold(op, U0, m, dti, T, T, 2);
    ido(i)   = toc;
    lepes(i) = round(T / dti);
    hiba(i)  = max(abs(Ui(:) - Uref(:)));
    fprintf('%-10s %10.4g %8d %12.3f %14.2e\n', m, dti, lepes(i), ido(i), hiba(i));
end

%% eredmények mentése (a dokumentum ábráihoz)
op_x = op.x;  op_y = op.y;  modszerek = esetek(:, 1);  dtk = cell2mat(esetek(:, 2));
Ukep = zeros(numel(op.x), numel(op.y), numel(tkep));
for i = 1:numel(tkep)
    Ukep(:, :, i) = teljes(Us(:, :, abs(ts - tkep(i)) < 1e-9));
end
save('pelda2_eredmenyek.mat', 'op_x', 'op_y', 'tkep', 'Ukep', 'ts', ...
     'szonda1', 'szonda2', 'modszerek', 'dtk', 'ido', 'hiba', 'lepes', '-v7');
