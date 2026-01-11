% xfem_1d_fixed_jump_verbose.m
% 1D Poisson: -(k u')' = f on (0,L), Dirichlet at x=0 and x=L.
% CG P_p (p=1..3) + XFEM Heaviside on the single cut element to enforce a
% prescribed value jump [u](xJ) = J. No static condensation: the two enriched
% DOFs are appended to the global system.
%
% Verbose printing of each assembly step and local matrices.

clear; clc; close all; format short g;

%% ===== Controls =====
L       = 1.0;     % domain [0, L]
Ne      = 3;       % *** keep only 3 elements ***
pDeg    = 2;       % polynomial degree: 1,2,3
kappa   = 1.0;     % diffusion coefficient
fVal    = 1.0;     % constant source
u0      = 0.0;     % Dirichlet at x=0
uL      = 0.0;     % Dirichlet at x=L

xJ      = 0.43;    % jump location (strictly inside (0,L))
Jjump   = 0.75;    % prescribed jump [u](xJ) = u(+)-u(-)
gammaJ  = 50;      % interface penalty strength

%% ===== Mesh, space, quadrature =====
xe   = linspace(0, L, Ne+1).';              % uniform mesh
p    = pDeg;  nLoc = p+1;
ndof_std = Ne*p + 1;                        % standard CG DOFs
xiNodes = lagrange_nodes_1D(p);
[~, gp, gw] = gauss_rule_1D(max(5,p+3));
evalShapes = @(xL,xR,xi) shapes_Pk_on_interval(xiNodes, xi, xL, xR);
mapSub     = @(gp,gw,aa,bb) deal(((bb-aa)/2).*gp + (aa+bb)/2, ((bb-aa)/2).*gw);
e2g        = @(e) ((e-1)*p + (1:nLoc));

% find cut element
assert(xJ>0 && xJ<L, 'xJ must be strictly inside (0,L).');
eJ = find(xe(1:end-1) < xJ & xJ < xe(2:end), 1);
assert(~isempty(eJ), 'xJ must not coincide with a mesh node.');

fprintf('\n=== SETUP ===\n');
fprintf('L=%g, Ne=%d, p=%d, ndof_std=%d\n', L, Ne, p, ndof_std);
fprintf('Jump location xJ=%.6f lies in element eJ=%d, [xL,xR]=[%.6f, %.6f]\n', ...
        xJ, eJ, xe(eJ), xe(eJ+1));

%% ===== Base CG assembly (standard DOFs only) =====
Kuu = sparse(ndof_std, ndof_std);
fu  = zeros(ndof_std,1);

fprintf('\n--- STEP 1: Base CG element matrices (Ke, Fe) and assembly ---\n');
for e = 1:Ne
    xL = xe(e); xR = xe(e+1); hE = xR-xL; Jmap = hE/2;
    [N, dNdx] = evalShapes(xL, xR, gp);
    Ke = (dNdx .* (gw.' * (Jmap*kappa))) * dNdx.';      % \int k u' v'
    Fe = N * (gw * (Jmap * fVal));                      % \int f v
    idx = e2g(e);

    fprintf('\nElement e=%d, [xL,xR]=[%.6f, %.6f], global dofs idx=%s\n', e, xL, xR, mat2str(idx));
    fprintf('Ke (nLoc x nLoc):\n'); disp(full(Ke));
    fprintf('Fe (nLoc x 1):\n');  disp(full(Fe));

    Kuu(idx,idx) = Kuu(idx,idx) + Ke;
    fu(idx)      = fu(idx)      + Fe;
end

%% ===== XFEM on the cut element: build local blocks Kua, Kau, Kaa, Fa (no condensation) =====
fprintf('\n--- STEP 2: XFEM on the cut element eJ=%d ---\n', eJ);
xL = xe(eJ); xR = xe(eJ+1); hE = xR-xL; Jmap = hE/2; idx = e2g(eJ);
xiJ = 2*(xJ - xL)/hE - 1;
hL  = max(xJ - xL, eps);
hR  = max(xR - xJ, eps);

% split the cut element into left/right subpieces in reference coords
xiBreaks = [-1, xiJ, 1];

% local enriched blocks: sizes (nLoc x 2), (2 x nLoc), (2 x 2)
Kua = zeros(nLoc,2); Kau = zeros(2,nLoc); Kaa = zeros(2,2); Fa = zeros(2,1);

for j = 1:2
    aR = xiBreaks(j); bR = xiBreaks(j+1);
    [gpS, gwS] = mapSub(gp, gw, aR, bR);
    [N, dNdx]  = evalShapes(xL, xR, gpS);

    onLeft  = (bR <= xiJ + 1e-14);
    onRight = (aR >= xiJ - 1e-14);

    % enriched gradients G = [G1; G2]
    % left subpiece:  a1' = 0,          a2' = -N_right'
    % right subpiece: a1' = +N_left',   a2' = 0
    G1 = zeros(1,numel(gpS)); G2 = zeros(1,numel(gpS));
    if onLeft
        G1(:) = 0;
        G2(:) = -dNdx(end,:);           % derivative of right end-node shape
        tag = 'LEFT';
    elseif onRight
        G1(:) =  dNdx(1,:);             % derivative of left end-node shape
        G2(:) =  0;
        tag = 'RIGHT';
    else
        error('Unexpected subpiece classification.');
    end
    G = [G1; G2];

    % subpiece contributions
    Kua_sub = (dNdx .* (gwS.' * (Jmap*kappa))) * G.';   % (nLoc x 2)
    Kaa_sub = (G     .* (gwS.' * (Jmap*kappa))) * G.';  % (2 x 2)
    Kau_sub = Kua_sub.';                                % symmetry

    fprintf('\n  Subpiece %d (%s): [xi_a, xi_b] = [%.6f, %.6f]\n', j, tag, aR, bR);
    fprintf('  Kua_sub (nLoc x 2):\n'); disp(full(Kua_sub));
    fprintf('  Kaa_sub (2 x 2):\n');    disp(full(Kaa_sub));

    Kua = Kua + Kua_sub;
    Kaa = Kaa + Kaa_sub;
    Kau = Kau + Kau_sub;
end

% interface penalty to enforce [u](xJ) = Jjump
[NJ, ~] = evalShapes(xL, xR, xiJ);
N1J = NJ(1); NnJ = NJ(end);
JJ  = [N1J; NnJ];
etaJ = gammaJ * kappa * (1/hL + 1/hR);

Kaa_pen = etaJ * (JJ*JJ.');   % only touches enriched DOFs (jump is purely in a)
Fa_pen  = etaJ * Jjump * JJ;

fprintf('\n--- STEP 3: Interface penalty to enforce [u](xJ)=J ---\n');
fprintf('etaJ = gammaJ*kappa*(1/hL+1/hR) = %g\n', etaJ);
fprintf('JJ = [N_left(xJ); N_right(xJ)] = [%g; %g]\n', N1J, NnJ);
fprintf('Kaa_pen (2x2):\n'); disp(full(Kaa_pen));
fprintf('Fa_pen  (2x1):\n');  disp(full(Fa_pen));

% add penalty to local enriched block
Kaa = Kaa + Kaa_pen;
Fa  = Fa  + Fa_pen;

fprintf('\n--- STEP 4: Local XFEM blocks (AFTER penalty) ---\n');
fprintf('Kua (nLoc x 2):\n'); disp(full(Kua));
fprintf('Kau (2 x nLoc):\n'); disp(full(Kau));
fprintf('Kaa (2 x 2):\n');    disp(full(Kaa));
fprintf('Fa  (2 x 1):\n');    disp(full(Fa));

%% ===== Build the full global system with enriched DOFs appended =====
ndof_tot = ndof_std + 2;
K = sparse(ndof_tot, ndof_tot);
F = zeros(ndof_tot,1);

% insert base CG blocks
K(1:ndof_std, 1:ndof_std) = Kuu;
F(1:ndof_std)             = fu;

% global indices for enriched DOFs
idA = ndof_std + (1:2);

% assemble the XFEM local blocks into global K/F
K(idx, idA) = K(idx, idA) + Kua;      % Kua
K(idA, idx) = K(idA, idx) + Kau;      % Kau = Kua'
K(idA, idA) = K(idA, idA) + Kaa;      % Kaa
F(idA)      = F(idA)      + Fa;       % Fa

fprintf('\n--- STEP 5: Global system BEFORE Dirichlet ---\n');
fprintf('Global DOFs: standard=%d, enriched=%d (total=%d)\n', ndof_std, 2, ndof_tot);
fprintf('Global K (showing full matrix):\n'); disp(full(K));
fprintf('Global F:\n'); disp(full(F));

%% ===== Strong Dirichlet at x=0 and x=L on standard DOFs =====
fprintf('\n--- STEP 6: Apply strong Dirichlet at x=0 (u=%g) and x=L (u=%g) ---\n', u0, uL);
[K,F] = apply_dirichlet_rc(K,F, 1,        u0);   % node at x=0
[K,F] = apply_dirichlet_rc(K,F, ndof_std, uL);   % node at x=L

fprintf('Global K AFTER Dirichlet:\n'); disp(full(K));
fprintf('Global F AFTER Dirichlet:\n'); disp(full(F));

%% ===== Solve =====
fprintf('\n--- STEP 7: Solve the global linear system ---\n');
sol = K \ F;

u_std = sol(1:ndof_std);
a_enr = sol(idA);
fprintf('Solution: enriched DOFs a = [a1; a2] = [%g; %g]\n', a_enr(1), a_enr(2));

%% ===== Evaluate jump and plot =====
% measure the jump just left/right of xJ
dx = 1e-8 * (xe(eJ+1)-xe(eJ));
uLm = eval_piecewise_with_enrich(xJ - dx, xe, u_std, eJ, a_enr, evalShapes, xJ);
uRp = eval_piecewise_with_enrich(xJ + dx, xe, u_std, eJ, a_enr, evalShapes, xJ);
measJ = uRp - uLm;
fprintf('Measured jump near xJ: u(+)-u(-) = %.8f (target %.8f)\n', measJ, Jjump);

% plot
figure('Color','w'); hold on; grid on;
xlabel('x'); ylabel('u(x)');
title(sprintf('1D XFEM jump (no condensation): Ne=%d, p=%d, [u](x_J)=%.3g at x_J=%.3f',Ne,pDeg,Jjump,xJ));
% draw mesh lines
for i=1:numel(xe), xline(xe(i), ':', 'Color',[.8 .8 .8]); end
% element-wise curve with enrichment on the cut element
colors = lines(Ne);
for e = 1:Ne
    xL = xe(e); xR = xe(e+1);
    xs = linspace(xL, xR, max(50, 18*p));
    ts = 2*(xs - xL)/(xR - xL) - 1;
    [N, ~] = evalShapes(xL, xR, ts);
    ue = u_std(e2g(e));
    us = (N.'*ue).';
    if e == eJ
        tJ = 2*(xJ - xL)/(xR - xL) - 1;
        leftMask  = ts <= tJ;
        rightMask = ts >= tJ;
        Nleft  = N(1,:); Nright = N(end,:);
        us(leftMask)  = us(leftMask)  + (-1)*Nright(leftMask) .* a_enr(2);
        us(rightMask) = us(rightMask) + (+1)*Nleft(rightMask)  .* a_enr(1);
    end
    plot(xs, us, '-', 'Color', colors(e,:), 'LineWidth', 1.8);
end
yl = ylim; plot([xJ xJ], yl, 'k--', 'LineWidth', 1.2);
legend(sprintf('measured [u]=%.6f', measJ), 'Location','southoutside');

%% ===== Helper functions =====
function xi = lagrange_nodes_1D(p)
    switch p
        case 1, xi = [-1, 1];
        case 2, xi = [-1, 0, 1];
        case 3, xi = [-1, -1/3, 1/3, 1];
        otherwise, error('pDeg must be 1..3');
    end
end

function [nG, xg, wg] = gauss_rule_1D(n)
    nG = n;
    switch n
      case 2
        xg=[-1/sqrt(3);1/sqrt(3)]; wg=[1;1];
      case 3
        xg=[-sqrt(3/5);0;sqrt(3/5)]; wg=[5/9;8/9;5/9];
      case 4
        xg=[-0.8611363116;-0.3399810436;0.3399810436;0.8611363116];
        wg=[ 0.3478548451;  0.6521451549;0.6521451549;0.3478548451];
      otherwise
        xg=[-0.9061798459;-0.5384693101;0;0.5384693101;0.9061798459];
        wg=[ 0.2369268851;  0.4786286705;0.5688888889;0.4786286705;0.2369268851];
    end
end

function [N,dNdx] = shapes_Pk_on_interval(xiNodes, xi, xL, xR)
    xi = xi(:).';  nG = numel(xi);  nLoc = numel(xiNodes);
    N = zeros(nLoc,nG); dN_dxi = zeros(nLoc,nG);
    for a = 1:nLoc
        La = ones(1,nG); denom = 1;
        for m = 1:nLoc
            if m==a, continue; end
            La = La .* (xi - xiNodes(m));
            denom = denom * (xiNodes(a) - xiNodes(m));
        end
        dLa = zeros(1,nG);
        for m = 1:nLoc
            if m==a, continue; end
            prod = ones(1,nG);
            for b = 1:nLoc
                if b==a || b==m, continue; end
                prod = prod .* (xi - xiNodes(b));
            end
            dLa = dLa + prod;
        end
        N(a,:)      = La  / denom;
        dN_dxi(a,:) = dLa / denom;
    end
    J = (xR - xL)/2;
    dNdx = (1/J) * dN_dxi;
end

function [K,F] = apply_dirichlet_rc(K,F,i,uD)
    F = F - K(:,i) * uD;
    K(:,i) = 0;  K(i,:) = 0;  K(i,i) = 1;  F(i) = uD;
end

function uval = eval_piecewise_with_enrich(x, xe, u_std, eJ, a_enr, evalShapes, xJ)
    % Evaluate u(x) with enriched piecewise contribution on the cut element
    Ne = numel(xe)-1; p = (numel(u_std)-1)/Ne; nLoc = p+1;
    e = find(xe(1:end-1) - 1e-14 <= x & x < xe(2:end) + 1e-14, 1, 'first');
    e = min(max(e,1), Ne);
    xL = xe(e); xR = xe(e+1);
    t = 2*(x - xL)/(xR - xL) - 1;
    [N,~] = evalShapes(xL, xR, t);
    idx = ((e-1)*p + (1:nLoc)).';
    uval = (N.' * u_std(idx)).';
    if e == eJ
        tJ = 2*(xJ - xL)/(xR - xL) - 1;
        if t < tJ
            uval = uval + (-1)*N(end).*a_enr(2);
        else
            uval = uval + (+1)*N(1).*a_enr(1);
        end
    end
end
