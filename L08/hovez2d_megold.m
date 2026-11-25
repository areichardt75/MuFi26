function [ts, Us] = hovez2d_megold(op, U0, modszer, dt, T, tsave, nBE)
%HOVEZ2D_MEGOLD  A 2D hővezetési egyenlet időléptetése különböző módszerekkel.
%
%   [ts, Us] = hovez2d_megold(op, U0, modszer, dt, T, tsave, nBE)
%
%   op      - a hovez2d_operatorok függvény kimenete
%   U0      - kezdeti hőmérséklet (mx x my mátrix)
%   modszer - 'explicit'   : explicit Euler (FTCS), mátrixos alakban
%             'BE_perjel'  : implicit Euler, minden lépésben  M \ b
%             'BE_LU'      : implicit Euler, egyszeri LU-felbontással
%             'CN_LU'      : Crank–Nicolson, egyszeri LU-felbontással
%             'ADI'        : Peaceman–Rachford ADI (csak tridiagonális rendszerek)
%   dt, T   - időlépés és végidő
%   tsave   - a tárolandó időpontok
%   nBE     - (CN_LU, ADI) ennyi kezdő lépést két fél implicit lépéssel
%             helyettesítünk (Rannacher-indítás), alapértelmezés: 2
%
%   ts      - a tárolt időpontok,  Us(:,:,k) - a megoldás ts(k)-ban

    if nargin < 7, nBE = 2; end
    Nt = round(T / dt);
    nsave = round(tsave / dt);
    ts = nsave * dt;
    Us = zeros(op.mx, op.my, numel(nsave));
    U  = U0;   u = U0(:);   k = 1;
    [Us, k] = tarol(Us, k, nsave, 0, U);
    N  = op.mx * op.my;   I = speye(N);

    switch modszer
        % -----------------------------------------------------------------
        case 'explicit'
            if dt > op.dtmax
                warning('dt = %.3g > dt_max = %.3g: instabil!', dt, op.dtmax);
            end
            for n = 1:Nt
                % U*Dy.' + Dx*U ugyanaz, mint A*u, de nem kell a nagy mátrix
                U = U + dt * (op.Dx * U + U * op.Dy.' + op.Bx);
                [Us, k] = tarol(Us, k, nsave, n, U);
            end

        % -----------------------------------------------------------------
        case 'BE_perjel'
            M = I - dt * op.A;
            for n = 1:Nt
                u = M \ (u + dt * op.b);      % minden lépésben új felbontás!
                [Us, k] = tarol(Us, k, nsave, n, reshape(u, op.mx, op.my));
            end

        % -----------------------------------------------------------------
        case 'BE_LU'
            % P*M*Q = L*R  (ritka LU sor- és oszlopcserével, kevés kitöltéssel)
            % MATLAB-ban egyszerűbben:  dM = decomposition(M);  u = dM \ rhs;
            [L, R, P, Q] = lu(I - dt * op.A);
            for n = 1:Nt
                u = Q * (R \ (L \ (P * (u + dt * op.b))));
                [Us, k] = tarol(Us, k, nsave, n, reshape(u, op.mx, op.my));
            end

        % -----------------------------------------------------------------
        case 'CN_LU'
            % I - dt/2*A: a CN-lépés és a fél BE-lépés közös bal oldali mátrixa
            [L, R, P, Q] = lu(I - 0.5 * dt * op.A);
            Mjobb = I + 0.5 * dt * op.A;                 % a CN jobb oldali mátrixa
            for n = 1:Nt
                if n <= nBE      % Rannacher-indítás: 2 x fél implicit Euler-lépés
                    u = Q * (R \ (L \ (P * (u + 0.5 * dt * op.b))));
                    u = Q * (R \ (L \ (P * (u + 0.5 * dt * op.b))));
                else
                    u = Q * (R \ (L \ (P * (Mjobb * u + dt * op.b))));
                end
                [Us, k] = tarol(Us, k, nsave, n, reshape(u, op.mx, op.my));
            end

        % -----------------------------------------------------------------
        case 'ADI'
            % Peaceman–Rachford: fél lépés x-ben implicit, y-ban explicit,
            % majd fordítva. Csak tridiagonális rendszereket kell megoldani,
            % ráadásul egyszerre az összes rácsvonalra (több jobb oldal).
            c  = 0.5 * dt;
            Mx = op.Ix - c * op.Dx;      % tridiagonális (mx x mx)
            My = op.Iy - c * op.Dy;      % tridiagonális (my x my)
            for n = 1:Nt
                if n <= nBE
                    % indítás: két fél lépés "lokálisan egydimenziós" (LOD)
                    % implicit Euler-rel (erősen csillapít, mint a BE)
                    for s = 1:2
                        U = Mx \ (U + c * op.Bx);
                        U = (My \ U.').';
                    end
                else
                    Ucs = Mx \ (U + c * (U * op.Dy.') + c * op.Bx);   % x-implicit
                    U   = (My \ (Ucs + c * (op.Dx * Ucs) + c * op.Bx).').';  % y-implicit
                end
                [Us, k] = tarol(Us, k, nsave, n, U);
            end

        otherwise
            error('Ismeretlen módszer: %s', modszer);
    end
end

function [Us, k] = tarol(Us, k, nsave, n, U)
% eltárolja U-t, ha az n-edik lépés a mentendők között van
    while k <= numel(nsave) && nsave(k) == n
        Us(:, :, k) = U;
        k = k + 1;
    end
end
