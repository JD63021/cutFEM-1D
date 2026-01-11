% cg1d_unit_domain_marked.m  (self-contained)
% Classical CG FEM on [0,1]:  - (k u')' = f,  u(0)=1, u(1)=0
% k = kappa (const), f = fVal (const). No cuts, no Nitsche, no ghost.
% Plots: element edges (vertical lines), per-element colored FEM curve,
% markers at element endpoints, and the exact solution.

clear; clc; close all;

%% ===== Controls =====
L        = 1.0;          % domain [0,1]
Ne       = 3;           % number of elements (kept at 10 as requested)
pDeg     = 3;            % polynomial degree: 1, 2, or 3
kappa    = 1.0;          % diffusion
fVal     = 10.0;         % constant source
u0       = 1.0;          % Dirichlet at x=0
uL       = 0.0;          % Dirichlet at x=1
% =====================

% Mesh: exactly 10 elements on [0,1]
xe = linspace(0, L, Ne+1).';    % nodes
ne = Ne;  nd_nodes = numel(xe);

% Space
p    = pDeg;
nLoc = p+1;
ndof = ne*p + 1;

% Reference nodes & quadrature
xiNodes      = lagrange_nodes_1D(p);
[~, gp, gw]  = gauss_rule_1D(max(5,p+3));
evalShapes   = @(xL,xR,xi) shapes_Pk_on_interval(xiNodes, xi, xL, xR);
e2g          = @(e) ((e-1)*p + (1:nLoc));

%% --- Assemble K and F (standard CG) ---
K = sparse(ndof, ndof);
F = zeros(ndof,1);

for e = 1:ne
    xL = xe(e); xR = xe(e+1); hE = xR - xL; J = hE/2;
    [N, dNdx] = evalShapes(xL, xR, gp);

    % Element stiffness and load
    Ke = (dNdx .* (gw.' * (J*kappa))) * dNdx.';   % ∫ kappa * dNdx^T dNdx
    Fe = N * (gw * (J * fVal));                   % ∫ f * N

    % Scatter
    idx = e2g(e);
    K(idx,idx) = K(idx,idx) + Ke;
    F(idx)     = F(idx)     + Fe;
end

%% --- Strong Dirichlet at x=0 and x=1 ---
[K,F] = apply_dirichlet_rc(K,F, 1,    u0);
[K,F] = apply_dirichlet_rc(K,F, ndof, uL);

%% --- Solve ---
u = K \ F;

%% --- Exact solution on [0,1] ---
u_exact = @(s) -5*s.^2 + 4*s + 1;

%% --- Plot: edges, exact, per-element colored FEM, endpoint markers ---
figure('Color','w'); ax = gca; hold(ax,'on'); grid(ax,'on');
xlabel(ax,'x'); ylabel(ax,'u(x)'); set(ax,'LineWidth',1.25);
ylim(ax,[-0.1, 2.0]); xlim(ax,[0, L]);

% 1) vertical element edges
for i = 1:nd_nodes
    xline(ax, xe(i), ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
end

% 2) exact solution (black)
xx = linspace(0,L,1601);
plot(ax, xx, u_exact(xx), 'k-', 'LineWidth', 1.8, 'DisplayName','exact');

% 3) FEM curve per element, colored, mark endpoints
colors = lines(ne);
for e = 1:ne
    xL = xe(e); xR = xe(e+1);
    xs = linspace(xL, xR, max(12, 6*p));
    t  = 2*(xs - xL)/(xR - xL) - 1;
    [Ns,~] = evalShapes(xL, xR, t);
    ue = u(e2g(e));
    us = (Ns.'*ue).';
    plot(ax, xs, us, '-', 'Color', colors(e,:), 'LineWidth', 1.6);

    % endpoints on FEM curve
    [Nend,~] = evalShapes(xL, xR, [-1, 1]);
    uEnds = (Nend.' * ue).';
    plot(ax, [xL, xR], uEnds, 'o', 'Color', colors(e,:), ...
         'MarkerFaceColor', colors(e,:), 'MarkerSize', 4);
end

title(ax, sprintf('CG FEM on [0,1] | Ne=10, p=%d, kappa=%.3g, f=%.3g', pDeg, kappa, fVal));
legend(ax,'Location','southoutside','NumColumns',3);

%% ========= Local helper functions =========
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
    % Strong Dirichlet at node i by row/column zeroing.
    F = F - K(:,i) * uD;
    K(:,i) = 0;  K(i,:) = 0;  K(i,i) = 1;  F(i) = uD;
end
