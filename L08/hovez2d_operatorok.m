function op = hovez2d_operatorok(p)
%HOVEZ2D_OPERATOROK  Ritka véges differenciás operátorok a 2D hővezetéshez.
%
%   u_t = alpha (u_xx + u_yy)   a [0,Lx] x [0,Ly] téglalapon
%   Peremfeltételek:
%       x = 0  : u = TL   (Dirichlet, adott hőmérséklet)
%       x = Lx : u = TR   (Dirichlet)
%       y = 0 és y = Ly : du/dy = 0  (Neumann, hőszigetelt perem)
%
%   Ismeretlenek: U(i,j) ~ u(x_i, y_j), ahol i = 1..Nx-1 (belső x-pontok),
%   j = 0..Ny (az y-irányú peremeken is ismeretlen a hőmérséklet).
%   Mátrixos alakban a szemidiszkrét rendszer:
%
%       dU/dt = Dx*U + U*Dy.' + Bx          (Dx, Dy már alpha-val szorozva)
%
%   vektoros alakban (u = U(:)):
%
%       du/dt = A*u + b,   A = kron(Iy, Dx) + kron(Dy, Ix),   b = Bx(:)

    dx = p.Lx / p.Nx;   dy = p.Ly / p.Ny;
    mx = p.Nx - 1;      my = p.Ny + 1;

    % x-irány: Dirichlet mindkét végen -> csak belső pontok
    ex = ones(mx, 1);
    Dx = p.alpha * spdiags([ex -2*ex ex], -1:1, mx, mx) / dx^2;

    % y-irány: Neumann mindkét végen, fiktív pont (ghost point) módszerrel:
    % u_{-1} = u_{1}  =>  (u_xx)_0 ~ (2 u_1 - 2 u_0)/dy^2   (másodrendű!)
    ey = ones(my, 1);
    Dy = spdiags([ey -2*ey ey], -1:1, my, my);
    Dy(1, 2)        = 2;
    Dy(my, my - 1)  = 2;
    Dy = p.alpha * Dy / dy^2;

    % a Dirichlet-peremértékek hozzájárulása
    Bx = zeros(mx, my);
    Bx(1, :)   = p.alpha * p.TL / dx^2;
    Bx(end, :) = p.alpha * p.TR / dx^2;

    op.Dx = Dx;  op.Dy = Dy;  op.Bx = Bx;
    op.Ix = speye(mx);  op.Iy = speye(my);
    op.A  = kron(op.Iy, Dx) + kron(Dy, op.Ix);
    op.b  = Bx(:);
    op.x  = (0:p.Nx)' * dx;   op.xi = op.x(2:end-1);
    op.y  = (0:p.Ny)' * dy;
    op.dx = dx;  op.dy = dy;  op.mx = mx;  op.my = my;
    op.TL = p.TL;  op.TR = p.TR;
    % az explicit módszer stabilitási határa
    op.dtmax = 1 / (2 * p.alpha * (1/dx^2 + 1/dy^2));
end
