function [x, ts, U] = hovez1d_theta(p)
%HOVEZ1D_THETA  Az 1D hővezetési egyenlet megoldása theta-módszerrel.
%
%   u_t = alpha * u_xx + f(x,t),   0 < x < L,   0 < t <= T
%   u(0,t) = gL(t),  u(L,t) = gR(t)            (Dirichlet-perem)
%   u(x,0) = u0(x)
%
%   Térben: másodrendű centrális differencia, Nx egyenlő részintervallum.
%   Időben: theta-módszer
%       theta = 0   : explicit Euler (FTCS)
%       theta = 1   : implicit Euler (BTCS)
%       theta = 1/2 : Crank–Nicolson
%
%   Bemenet: p struktúra a következő mezőkkel
%       alpha, L, Nx, dt, T   - fizikai és numerikus paraméterek
%       u0(x)                 - kezdeti feltétel (függvény handle)
%       gL(t), gR(t)          - peremértékek (függvény handle)
%       tsave                 - azok az időpontok, ahol a megoldást eltároljuk
%       theta   (opcionális)  - alapértelmezés: 0.5
%       nBE     (opcionális)  - az első nBE lépést két fél lépésnyi implicit
%                               Euler-lépéssel helyettesíti (Rannacher-indítás)
%       f(x,t)  (opcionális)  - forrástag, alapértelmezés: 0
%
%   Kimenet:
%       x   - rácspontok (oszlopvektor, a peremekkel együtt)
%       ts  - a ténylegesen eltárolt időpontok
%       U   - a megoldás: U(:,k) a ts(k) időpontban

    % --- alapértelmezett értékek ---------------------------------------
    if ~isfield(p, 'theta'), p.theta = 0.5; end
    if ~isfield(p, 'nBE'),   p.nBE   = 0;   end
    if ~isfield(p, 'f'),     p.f     = @(x, t) zeros(size(x)); end

    % --- rács ------------------------------------------------------------
    dx = p.L / p.Nx;
    x  = (0:p.Nx)' * dx;
    xi = x(2:end-1);                  % belső (ismeretlen) pontok
    m  = p.Nx - 1;
    Nt = round(p.T / p.dt);

    r = p.alpha * p.dt / dx^2;        % rácsparaméter (a rács Fourier-száma)
    if p.theta < 0.5 && r > 1 / (2 * (1 - 2 * p.theta))
        warning('hovez1d_theta:instabil', ...
            'r = %.3f > a stabilitási határ, a megoldás instabil lehet!', r);
    end

    % --- ritka mátrix: a második derivált közelítése (D_xx) -------------
    e = ones(m, 1);
    A = p.alpha * spdiags([e -2*e e], -1:1, m, m) / dx^2;
    I = speye(m);

    % a peremfeltételekből és a forrásból származó vektor: b(t)
    b = @(t) [p.alpha * p.gL(t) / dx^2; zeros(m-2, 1); ...
              p.alpha * p.gR(t) / dx^2] + p.f(xi, t);

    % --- a lépésenként felhasznált mátrixok (egyszer építjük fel) ---------
    th  = p.theta;
    ML  = I - th * p.dt * A;          % bal oldal (implicit rész)
    MR  = I + (1 - th) * p.dt * A;    % jobb oldal (explicit rész)
    MBE = I - 0.5 * p.dt * A;         % fél lépésnyi implicit Euler (indítás)
    % Megjegyzés: ML tridiagonális, ezt a "\" operátor felismeri és sávos
    % (Thomas-típusú) megoldót használ, így O(m) műveletigényű.

    % --- kezdeti feltétel és tárolás ---------------------------------------
    u = p.u0(xi);
    nsave = round(p.tsave / p.dt);
    ts = nsave * p.dt;
    U  = zeros(p.Nx + 1, numel(nsave));
    k  = 1;
    while k <= numel(nsave) && nsave(k) == 0
        U(:, k) = [p.gL(0); u; p.gR(0)];
        k = k + 1;
    end

    % --- időléptetés -----------------------------------------------------
    t = 0;
    for n = 1:Nt
        if n <= p.nBE
            % Rannacher-indítás: két fél lépés implicit Euler-rel
            u = MBE \ (u + 0.5 * p.dt * b(t + 0.5 * p.dt));
            u = MBE \ (u + 0.5 * p.dt * b(t + p.dt));
        elseif th == 0
            % explicit lépés: csak mátrix–vektor szorzás
            u = u + p.dt * (A * u + b(t));
        else
            rhs = MR * u + p.dt * (th * b(t + p.dt) + (1 - th) * b(t));
            u = ML \ rhs;
        end
        t = n * p.dt;

        while k <= numel(nsave) && nsave(k) == n
            U(:, k) = [p.gL(t); u; p.gR(t)];
            k = k + 1;
        end
    end
end
