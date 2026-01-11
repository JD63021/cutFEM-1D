% cutfem_1d_marked_driver_boxes.m  (self-contained)
% 3 elements on [0,1.1], P3 CG, symmetric Nitsche at two internal
% Dirichlet points c1 (u=1) and c2=c1+1 (u=0). Source f=10 only on (c1,c2].
% Strong Dirichlet at x=0 and x=L. Ghost penalty optional.
% Plots: main FEM view + two Gauss "boxes" figures on a cut element.

clear; clc; close all;

%% ===== Controls =====
L        = 1.1;          % domain [0, L]
Ne       = 3;            % exactly 3 elements
pDeg     = 3;            % polynomial degree: 1,2,3
kappa    = 1.0;          % diffusion (uniform)
fVal     = 10.0;         % source magnitude inside the active segment
segLen   = 1.0;          % active segment length

c0       = 0.05;         % base (left) pin
shift    = 0.0;          % translation
gammaN   = 6;            % Nitsche penalty
useGhost = true;         % ghost penalty near each cut
gammaG   = 0.5;          % ghost coefficient
tiny     = 1e-12;        % avoid exact-on-node issues
% =====================

% Derived pins
c1 = c0 + shift;              % u(c1)=1
c2 = c1 + segLen;             % u(c2)=0
assert(c1>0 && c2<L, 'Segment must lie strictly inside (0,L).');

% Mesh: exactly 3 elements
xe = linspace(0, L, Ne+1).';  % nodes: [0, L/3, 2L/3, L]
ne = Ne;  nd_nodes = numel(xe);

% Space
p    = pDeg;
nLoc = p+1;
ndof = ne*p + 1;

% Reference nodes & quadrature
xiNodes      = lagrange_nodes_1D(p);
[~, gp, gw]  = gauss_rule_1D(max(5,p+3));
evalShapes   = @(xL,xR,xi) shapes_Pk_on_interval(xiNodes, xi, xL, xR);
mapSub       = @(gp,gw,aa,bb) deal(((bb-aa)/2).*gp + (aa+bb)/2, ((bb-aa)/2).*gw);
e2g          = @(e) ((e-1)*p + (1:nLoc));

%% --- Assemble K,F with piecewise (subpiece) integration ---
K = sparse(ndof, ndof);
F = zeros(ndof,1);

for e = 1:ne
    xL = xe(e); xR = xe(e+1); hE = xR - xL; J = hE/2;
    idx = e2g(e);
    Ke = zeros(nLoc); Fe = zeros(nLoc,1);

    % reference breaks for cuts (if they fall inside this element)
    xiBreaks = [-1, 1];
    if (xL < c1 && c1 < xR), xiBreaks(end+1) = 2*(c1 - xL)/hE - 1; end
    if (xL < c2 && c2 < xR), xiBreaks(end+1) = 2*(c2 - xL)/hE - 1; end
    xiBreaks = sort(xiBreaks);

    for j = 1:numel(xiBreaks)-1
        aR = xiBreaks(j); bR = xiBreaks(j+1);
        if bR-aR < tiny, continue; end

        [gpS, gwS] = mapSub(gp, gw, aR, bR);
        [N, dNdx]  = evalShapes(xL, xR, gpS);

        % Is this subpiece inside (c1,c2]? (midpoint test)
        xiMid = 0.5*(aR+bR);
        xMid  = xL + (xiMid+1)*J;
        f_here = (xMid > c1 && xMid <= c2) * fVal;

        % stiffness & load on this subpiece
        Ke = Ke + (dNdx .* (gwS.' * (J*kappa))) * dNdx.';
        Fe = Fe + N * (gwS * (J * f_here));
    end

    K(idx,idx) = K(idx,idx) + Ke;
    F(idx)     = F(idx)     + Fe;
end

%% --- Internal Dirichlet pins via symmetric Nitsche ---
[K,F] = add_nitsche_point(K, F, xe, e2g, evalShapes, c1, 1.0, gammaN, kappa);
[K,F] = add_nitsche_point(K, F, xe, e2g, evalShapes, c2, 0.0, gammaN, kappa);

%% --- Ghost penalty near each cut (optional) ---
if useGhost
    eCut1 = find(xe(1:end-1) < c1 & c1 < xe(2:end), 1);
    eCut2 = find(xe(1:end-1) < c2 & c2 < xe(2:end), 1);
    if ~isempty(eCut1), K = add_ghost_faces(K, xe, eCut1, gammaG); end
    if ~isempty(eCut2), K = add_ghost_faces(K, xe, eCut2, gammaG); end
end

%% --- STRONG Dirichlet at domain ends (hard enforce) ---
[K,F] = apply_dirichlet_rc(K,F, 1,      1.0);   % u(0)=1
[K,F] = apply_dirichlet_rc(K,F, ndof,   0.0);   % u(L)=0

%% --- Solve ---
u = K \ F;

%% --- Main plot: edges, exact, per-element FEM, endpoint markers, pins ---
figure('Color','w'); ax = gca; hold(ax,'on'); grid(ax,'on');
xlabel(ax,'x'); ylabel(ax,'u(x)'); set(ax,'LineWidth',1.25);
ylim(ax,[-0.1, 2.0]); xlim(ax,[0, L]);

% element edges
for i = 1:nd_nodes
    xline(ax, xe(i), ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
end

% exact reference on the active segment (black)
u_exact = @(x) (x<=c1).*1 ...
             + (x>c1 & x<=c2).*(-5*(x-c1).^2 + 4*(x-c1) + 1) ...
             + (x>c2).*0;
xx = linspace(0,L,800);
plot(ax, xx, u_exact(xx), 'k-', 'LineWidth', 1.8, 'DisplayName','exact');

% FEM per element (colored) + endpoint markers
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

% moving boundary markers
plot(ax, c1, 1.0, 's', 'MarkerSize', 9, 'MarkerFaceColor',[0.2 0.7 0.2], ...
     'MarkerEdgeColor','k', 'DisplayName','c_1 (u=1)');
plot(ax, c2, 0.0, 'd', 'MarkerSize', 9, 'MarkerFaceColor',[0.8 0.2 0.8], ...
     'MarkerEdgeColor','k', 'DisplayName','c_2 (u=0)');

title(ax, sprintf('CutFEM 1D | Ne=3, p=%d, Nitsche=%g, ghost=%d (\\gamma_g=%g)', ...
      pDeg, gammaN, useGhost, gammaG));
legend(ax,'Location','southoutside','NumColumns',3);

%% --- Gauss boxes illustration on one cut element (left & right subpieces) ---
% Pick the cut at c1 (first cut element)
eCut1 = find(xe(1:end-1) < c1 & c1 < xe(2:end), 1);
xL1 = xe(eCut1); xR1 = xe(eCut1+1); c = c1;

% standard 2-pt Gauss on [-1,1]
xi2 = [-1/sqrt(3), 1/sqrt(3)];
w2  = [1, 1];

% helper to map [-1,1] -> [a,b]
mapRef = @(a,b,xi) a + 0.5*(b-a)*(xi+1);

% ---- Left subpiece [xL1, c] ----
aL = xL1; bL = c; lenL = bL - aL;
xL_pts = mapRef(aL,bL,xi2);         % two Gauss points (left subpiece)
boxW_L  = 0.5*lenL;                 % each box width = (subLen)/2 (since w=1)
boxH_L  = 1.0;                      % schematic height=1

figure('Color','w'); axL = gca; hold(axL,'on'); grid(axL,'on');
xlabel(axL,'x'); ylabel(axL,'schematic integrand'); title(axL,'Left subpiece [x_L, c]: 2-pt Gauss boxes');
xlim(axL,[xL1-0.05, xR1+0.05]); ylim(axL,[0, 1.3]);

% draw subpiece span
plot(axL, [aL bL], [0 0], 'k-', 'LineWidth',1.2); yline(axL,1,'k:');
xline(axL, xL1, ':', 'Color',[0.5 0.5 0.5]);
xline(axL, c,   ':', 'Color',[0.5 0.5 0.5]);

% draw the two boxes centered at Gauss points (width=boxW_L, height=1)
for k = 1:2
    x0 = xL_pts(k) - boxW_L/2;
    rectangle('Position',[x0, 0, boxW_L, boxH_L],'FaceColor',[0.6 0.8 1],'EdgeColor','b');
    plot(axL, xL_pts(k), 1, 'bo', 'MarkerFaceColor','b');
end
legend(axL, {'subpiece','h=1'}, 'Location','northwest');

% ---- Right subpiece [c, xR1] ----
aR = c; bR = xR1; lenR = bR - aR;
xR_pts = mapRef(aR,bR,xi2);
boxW_R  = 0.5*lenR;
boxH_R  = 1.0;

figure('Color','w'); axR = gca; hold(axR,'on'); grid(axR,'on');
xlabel(axR,'x'); ylabel(axR,'schematic integrand'); title(axR,'Right subpiece [c, x_R]: 2-pt Gauss boxes');
xlim(axR,[xL1-0.05, xR1+0.05]); ylim(axR,[0, 1.3]);

plot(axR, [aR bR], [0 0], 'k-', 'LineWidth',1.2); yline(axR,1,'k:');
xline(axR, c,   ':', 'Color',[0.5 0.5 0.5]);
xline(axR, xR1, ':', 'Color',[0.5 0.5 0.5]);

for k = 1:2
    x0 = xR_pts(k) - boxW_R/2;
    rectangle('Position',[x0, 0, boxW_R, boxH_R],'FaceColor',[1 0.85 0.6],'EdgeColor',[0.85 0.3 0]);
    plot(axR, xR_pts(k), 1, 'o', 'Color',[0.85 0.3 0], 'MarkerFaceColor',[0.85 0.3 0]);
end
legend(axR, {'subpiece','h=1'}, 'Location','northwest');

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

function [K,F] = add_nitsche_point(K, F, xe, e2g, evalShapes, c, g, gammaN, kappa)
    % Symmetric Nitsche at x=c: enforce u(c)=g (c strictly inside an element).
    tiny = 1e-12;
    eCut = find(xe(1:end-1) < c & c < xe(2:end), 1);
    if isempty(eCut)
        c = c + sign(0.5-rand)*tiny;  % nudge off a node if necessary
        eCut = find(xe(1:end-1) < c & c < xe(2:end), 1);
    end
    xL = xe(eCut); xR = xe(eCut+1); hE = xR - xL;
    xi = 2*(c - xL)/hE - 1;

    [N_c, dNdx_c] = evalShapes(xL, xR, xi);
    N_c    = N_c(:);  dNdx_c = dNdx_c(:);

    % 1D "sides" (nL=-1, nR=+1)
    hL = max(c - xL, tiny);  nL = -1;
    hR = max(xR - c, tiny);  nR = +1;

    % consistency + adjoint-consistency (sum both sides)
    Kc = - nL*kappa*( dNdx_c*N_c.' + N_c*dNdx_c.' ) ...
         - nR*kappa*( dNdx_c*N_c.' + N_c*dNdx_c.' );
    % penalty (sum both sides)
    Kp = (gammaN*kappa/hL)*(N_c*N_c.') + (gammaN*kappa/hR)*(N_c*N_c.');

    % RHS with u(c)=g
    Fc = ( -nL*kappa*dNdx_c + -nR*kappa*dNdx_c ) * g;
    Fp = ( (gammaN*kappa/hL) + (gammaN*kappa/hR) ) * (N_c*g);

    idx = e2g(eCut);
    K(idx,idx) = K(idx,idx) + Kc + Kp;
    F(idx)     = F(idx)     + Fc + Fp;
end

function K = add_ghost_faces(K, xe, eCut, gammaG)
    % Add ghost penalty on the two faces adjacent to element eCut
    nd = numel(xe);
    faces = [];
    if eCut>1,           faces(end+1) = eCut;     end
    if eCut+1 < nd,      faces(end+1) = eCut+1;   end
    for node = faces
        if node<=1 || node>=nd, continue; end
        hL = xe(node)   - xe(node-1);
        hR = xe(node+1) - xe(node);
        hF = 0.5*(hL+hR);
        b  = [-1/hL; 1/hL + 1/hR; -1/hR];     % slope-jump coeffs
        d  = [node-1, node, node+1];
        Jface = gammaG * hF * (b*b.');        % 3x3 block
        K(d,d) = K(d,d) + Jface;
    end
end

function [K,F] = apply_dirichlet_rc(K,F,i,uD)
    % Strong Dirichlet at node i (row/column zeroing).
    F = F - K(:,i) * uD;
    K(:,i) = 0;  K(i,:) = 0;  K(i,i) = 1;  F(i) = uD;
end
