function out = figure3(varargin)
%FIGURE3 Reproduce the three panels of Figure 3 (search for two donors) and
%   assemble the lettered figure output/figs/fig3.pdf.
%
% Usage
%   figure3
%   figure3('outdir','output/figure3')
%
% Panel map
%   (a) fig2d   schematic of two donors on the ring and the split of searchers
%   (b) fig2e   averaged splitting probability against reversal rate, two
%               donors of unequal absorptivity, arc solution and simulation
%   (c) fig2f   arrival rates and conditional times with one advertised donor
%
% The two-target solver is atwo. The Monte Carlo ensembles of (b) and (c) are
% cached in outdir under stamped keys.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'src')));

ip = inputParser;
addParameter(ip,'outdir',fullfile('output','figure3'));
parse(ip,varargin{:});
OUTDIR = ip.Results.outdir;
if ~isfolder(OUTDIR), mkdir(OUTDIR); end

A   = afpt();
prm = A.mkprm('kbar',0.0314,'sig',0.8,'beta0',1,'sm',1,'v',1);
panel_a(OUTDIR);
out.e = panel_b(prm, OUTDIR);
out.f = panel_c(prm, OUTDIR);

% Lettered figure, one row of three panel canvases at the panel size.
figdir = fullfile('output','figs');
if ~isfolder(figdir), mkdir(figdir); end
T = [22 19];
tiles = [0 0 T; T(1) 0 T; 2*T(1) 0 T];
panels = fullfile(OUTDIR, strcat({'fig2d','fig2e','fig2f'}, '.fig'));
assemble_figure(fullfile(figdir,'fig3'), panels, tiles, [3*T(1) T(2)]);
end

% =========================================================================
% Panel (a), fig2d. Two-donor schematic.
% =========================================================================
function out = panel_a(outdir)
%PANEL_A  Two-donor schematic. Panel (a).
%
%   Searchers on the ring are split between two donors. Donor 1 absorbs more
%   strongly (larger, darker disk, kbar_1) and carries the cue plume, a tent
%   rising off the ring beside it. Donor 2 absorbs more weakly (kbar_2). The
%   grey arcs from the top of the ring carry the splitting probabilities pi_1
%   and pi_2, the heavier arc toward the donor that wins more searchers.
%
%   The palette and glyphs are those of Fig. 1(a). Geometry is in points on the
%   22 x 19 cm panel canvas, y down, on one axes with equal data aspect that
%   fills the canvas, so assembly at tile = canvas keeps type size.

if nargin < 1 || isempty(outdir), outdir = fullfile('output','figure3'); end
if ~isfolder(outdir), mkdir(outdir); end
P = linger_style();

C.blue  = [0 162 255]/255;                   % searching
C.navy  = [0 77 128]/255;
C.green = [29 177 0]/255;
C.dkgrn = [1 113 0]/255;                     % cue
C.dgrey = [94 94 94]/255;                    % donor
C.mgrey = [146 146 146]/255;
C.ring  = [213 213 213]/255;

cm = 2.54/72;
PAGE = P.FIG/cm;                             % pt
h = P.figure(P);
set(h, 'Color','w', 'InvertHardcopy','off');
ax = axes(h, 'Units','normalized', 'Position',[0 0 1 1], ...
          'XLim',[0 PAGE(1)], 'YLim',[0 PAGE(2)], 'YDir','reverse', ...
          'DataAspectRatio',[1 1 1], 'Visible','off', 'Clipping','off');
hold(ax, 'on');
% white corner marks, so the tight crop of the export keeps the full canvas
for cx = [0 PAGE(1)-0.2], for cy = [0 PAGE(2)-0.2]
    patch(ax, cx + [0 .2 .2 0], cy + [0 0 .2 .2], 'w', 'EdgeColor','none');
end, end

ctr = [PAGE(1)/2 + 22, PAGE(2)/2 + 4];  R = 178;  RLW = 20;
a1 = pi;  a2 = 0;                            % donor angles, y down
d1 = ctr + R*[cos(a1) sin(a1)];
d2 = ctr + R*[cos(a2) sin(a2)];

% cue plume, a tent on the ring beside donor 1
plume(ax, C, ctr, R + RLW/2, a1, 0.62, 92);
ring(ax, C, ctr, R, RLW);

% split arcs from the top of the ring toward each donor
ra = R - 52;  top = -pi/2;
splitArc(ax, C, ctr, ra, top - 0.16, a1 - 2*pi + 0.30, 7.0, 30);
splitArc(ax, C, ctr, ra, top + 0.16, a2 - 0.30, 4.0, 22);
math(ax, ctr + (ra - 60)*[cos(-2.42) sin(-2.42)], '$\pi_1$', P.LSZ, C.dgrey);
math(ax, ctr + (ra - 60)*[cos(-0.72) sin(-0.72)], '$\pi_2$', P.LSZ, C.dgrey);

% searching bees, tangent to the ring in the sense of their motion
SA = 48;
ab = [-pi/2, -2.20, -0.95, 0.55, 1.25, 2.00, 2.62];
sv = [-1 0; 0 1; 1 0.3; 0 -1; 1 0; -1 0; 0 -1];
for k = 1:numel(ab)
    xy = ctr + R*[cos(ab(k)) sin(ab(k))];
    glyph(ax, C, xy, tangent(ctr, xy, sv(k,:)), SA);
end

% donors, drawn last so they sit on the ring
disk(ax, d1, 27, 0.78*C.dgrey, C.mgrey, 3);
disk(ax, d2, 16, C.mgrey, C.ring, 2.5);
math(ax, d1 + [66 52], '$\overline\kappa_1$', P.LSZ, C.dgrey);
math(ax, d2 + [-58 52], '$\overline\kappa_2$', P.LSZ, C.dgrey);
text(ax, d1(1) - 62, d1(2) - 128, 'cue', 'Interpreter','latex', 'FontSize', P.FSZ, ...
     'Color', C.dkgrn, 'HorizontalAlignment','center', 'VerticalAlignment','middle');

set(ax.Children, 'Clipping', 'off');
P.export(P, 'fig2d', outdir, h);
close(h);
out = struct('ctr',ctr,'R',R);
end

% ------------------------------------------------------------------------
function plume(ax, C, ctr, r0, ac, hw, ht)
% Tent of angular half width HW and height HT standing on radius R0 at angle
% AC, shaded dark green at its base under the apex and fading to white at the
% flanks and the tip.
na = 81;  ns = 12;
a = ac + hw*linspace(-1, 1, na);
tent = 1 - abs(a - ac)/hw;
s = linspace(0, 1, ns);
V = zeros(na*ns, 2);  T = zeros(na*ns, 1);
for j = 1:ns
    r = r0 + s(j)*ht*tent;
    idx = (j-1)*na + (1:na);
    V(idx,:) = [ctr(1) + r'.*cos(a'), ctr(2) + r'.*sin(a')];
    T(idx) = 1 - 0.92*tent'.*(1 - 0.75*s(j));
end
F = zeros(2*(na-1)*(ns-1), 3);  m = 0;
for j = 1:ns-1
    for i = 1:na-1
        p = (j-1)*na + i;  q = p + na;
        m = m + 1;  F(m,:) = [p p+1 q+1];
        m = m + 1;  F(m,:) = [p q+1 q];
    end
end
col = (1 - T).*C.dkgrn + T.*[1 1 1];
patch(ax, 'Faces', F, 'Vertices', V, 'FaceVertexCData', col, ...
      'FaceColor','interp', 'EdgeColor','none');
r = r0 + ht*tent;
plot(ax, ctr(1) + r.*cos(a), ctr(2) + r.*sin(a), '-', 'Color', C.dkgrn, 'LineWidth', 1.6);
end

function ring(ax, C, ctr, R, lw)
th = linspace(0, 2*pi, 361);
plot(ax, ctr(1) + R*cos(th), ctr(2) + R*sin(th), '-', 'Color', C.ring, 'LineWidth', lw);
end

function splitArc(ax, C, ctr, r, a0, a1, lw, hl)
% Arc of radius R from angle A0 to A1 with an arrowhead of length HL at A1.
col = C.mgrey;
L = abs(a1 - a0)*r;
ae = a1 - sign(a1 - a0)*hl/r;                % shaft stops at the head base
a = linspace(a0, ae + sign(a1 - a0)*0.3/r, 120);
plot(ax, ctr(1) + r*cos(a), ctr(2) + r*sin(a), '-', 'Color', col, 'LineWidth', lw);
tip  = ctr + r*[cos(a1) sin(a1)];
base = ctr + r*[cos(ae) sin(ae)];
d = (tip - base)/norm(tip - base);  n = [-d(2) d(1)];
hw = 0.5*hl + lw;
patch(ax, [tip(1) base(1)+hw/2*n(1) base(1)-hw/2*n(1)], ...
          [tip(2) base(2)+hw/2*n(2) base(2)-hw/2*n(2)], col, 'EdgeColor','none');
if L <= hl, warning('fig2d:shortArc', 'split arc shorter than its head'); end
end

function disk(ax, xy, rout, fc, ec, lw)
th = linspace(0, 2*pi, 73);  r = rout - lw/2;
patch(ax, xy(1) + r*cos(th), xy(2) + r*sin(th), fc, 'EdgeColor', ec, 'LineWidth', lw);
end

function glyph(ax, C, xy, hd, s)
% Searching bee of body length S, body midpoint at XY, heading along HD.
d = hd/norm(hd);
o = xy + 0.07*s*d;
draw_bee(ax, o(1), o(2), atan2(d(2), d(1)), s, C.blue, 'edge', s/32);
end

function t = tangent(ctr, xy, ref)
% unit tangent to the circle about CTR through XY, on the side of REF
r = xy - ctr;  t = [-r(2) r(1)]/norm(r);
if dot(t, ref) < 0, t = -t; end
end

function t = math(ax, xy, str, fs, col)
t = text(ax, xy(1), xy(2), str, 'Interpreter','latex', 'FontSize', fs, 'Color', col, ...
         'HorizontalAlignment','center', 'VerticalAlignment','middle');
end

% =========================================================================
% Panel (b), fig2e. Averaged splitting probability against reversal rate.
% =========================================================================
function out = panel_b(prm, outdir, bgrid, kbars, nmc, seed)
%PANEL_B  Averaged splitting probability pibar_1 against reversal rate, two
%   donors of unequal absorptivity, no cues. Arc solution and Monte Carlo.
%   Panel (b).
%
%   The arc solution and the Monte Carlo are cached under a stamped key in
%   OUTDIR.

if nargin < 2 || isempty(outdir), outdir = 'output'; end
if nargin < 3 || isempty(bgrid),  bgrid = logspace(-3, 2, 26); end
if nargin < 4 || isempty(kbars),  kbars = [0.06 0.015]; end
if nargin < 5 || isempty(nmc),    nmc = 40000; end
if nargin < 6 || isempty(seed),   seed = 3; end

M = lingermodel();  P = linger_style();
v = prm.v;  eta = tanh(kbars/(2*v));
bal = eta(1)/sum(eta);  dif = 0.5;
bmc = bgrid(round(linspace(1, numel(bgrid), 7)));

R = [];
if isempty(R)
    cfile = fullfile(outdir, 'fig2e_cache.mat');
    % The stamp tag is a fixed part of the cache identity.
    key = M.stamp('fig6b', prm, bgrid, kbars, nmc, seed);
    if isfile(cfile)
        C = load(cfile);
        if isfield(C, 'key') && strcmp(C.key, key), R = C.R; fprintf('fig2e: loaded stamped cache\n'); end
    end
    if isempty(R)
        R = compute(prm, bgrid, kbars, nmc, seed, bmc);
        save(cfile, 'R', 'key');
    end
end

cfv = eta(1)*(pi*R.beta*eta(2) + v)./(2*pi*R.beta*eta(1)*eta(2) + v*eta(1) + v*eta(2));
dev = (R.pi_mc - interp1(R.beta, R.pi1, R.beta_mc))./R.se_mc;
[~, i1] = min(abs(R.beta - 1));
fprintf('fig2e: pibar_1 = %.3f at beta0 = 1 (ballistic %.3f), solver vs closed form %.1e, MC within %.1f SE\n', ...
        R.pi1(i1), bal, max(abs(R.pi1 - cfv)), max(abs(dev)));

% ---- panel --------------------------------------------------------------
h = P.figure(P);  ax = axes('Parent', h);  hold(ax, 'on');
plot(ax, [min(bgrid) max(bgrid)], bal*[1 1], '--', 'Color', P.grey, 'LineWidth', P.ALW);
plot(ax, [min(bgrid) max(bgrid)], dif*[1 1], ':',  'Color', P.grey, 'LineWidth', P.ALW);
plot(ax, R.beta, R.pi1, '-', 'Color', P.greens(4,:), 'LineWidth', P.LW);
errorbar(ax, R.beta_mc, R.pi_mc, 1.96*R.se_mc, 'o', 'Color', P.blues(3,:), ...
         'MarkerFaceColor','w', 'LineWidth', P.ALW, 'MarkerSize', 12, ...
         'CapSize', 10, 'LineStyle','none');
set(ax, 'XScale', 'log');
xlabel(ax, 'reversal rate $\beta$', 'Interpreter','latex', 'FontSize', P.LSZ);
ylabel(ax, 'splitting probability $\bar\pi_1$', 'Interpreter','latex', 'FontSize', P.LSZ);
xlim(ax, [min(bgrid) max(bgrid)]);  ylim(ax, [0.45, bal + 0.05]);
set(ax, 'XTick', 10.^(-3:2));
P.axes(P, ax);
set(ax, 'Position', [0.19 0.17 0.74 0.72]);
P.label(P, ax, 1.4e-3, bal + 0.008, 'ballistic, $\eta_1/(\eta_1+\eta_2)$', ...
        'Color', P.grey, 'VerticalAlignment','bottom');
P.label(P, ax, 1.4e-3, dif + 0.008, 'diffusive, $1/2$', ...
        'Color', P.grey, 'VerticalAlignment','bottom');
P.label(P, ax, 0.9, 0.745, 'arc solution', 'Color', P.greens(4,:), ...
        'HorizontalAlignment','right');
P.label(P, ax, 8, 0.66, 'simulation', 'Color', P.blues(3,:), ...
        'HorizontalAlignment','right');

P.export(P, 'fig2e', outdir, h);
close(h);

out = R;  out.ballistic = bal;  out.diffusive = dif;  out.closedform = cfv;
out.dev = dev;  out.kbar = kbars;  out.eta = eta;
end

% ------------------------------------------------------------------------
function R = compute(prm, bgrid, kbars, nmc, seed, bmc)
A = atwo();  traps = [-pi/2, pi/2];
th = 1 - exp(-kbars/prm.v);
pis = zeros(size(bgrid));
for i = 1:numel(bgrid)
    p = A.split(traps, th, [], prm.sig, bgrid(i), 0);
    pis(i) = p(1);
end
pmc = zeros(size(bmc));  smc = zeros(size(bmc));
for i = 1:numel(bmc)
    who = A.mc(nmc, traps, th, [], prm.sig, bmc(i), 0, prm.sm, prm.v, seed+i);
    p1 = mean(who == 1);
    pmc(i) = p1;  smc(i) = sqrt(p1*(1-p1)/numel(who));
end
R = struct('beta',bgrid,'pi1',pis,'beta_mc',bmc,'pi_mc',pmc,'se_mc',smc);
end

% =========================================================================
% Panel (c), fig2f. One advertised donor.
% =========================================================================
function out = panel_c(prm, outdir, ugrid, umc, nmc, seed)
%PANEL_C  One advertised donor. Arrival rates at the advertised donor and
%   its rival, relative to their no-cue values, against reliance
%   u = gamma/(s_m sigma), with a cue at donor 1 only at zero offset. The
%   conditional first-passage times to either donor sit on the right axis.
%   Panel (c).
%
%   Equal donors at -pi/2 (advertised) and +pi/2 at the representative parameters.
%   Lambda_j = pi_j/Tbar as in atwo. The Monte Carlo rate is the fraction
%   absorbed at j over the mean absorption time, with a delta-method error.

if nargin < 2 || isempty(outdir), outdir = 'output'; end
if nargin < 3 || isempty(ugrid),  ugrid = linspace(0, 1, 41); end
if nargin < 4 || isempty(umc),    umc = [0 0.25 0.5 0.75 1]; end
if nargin < 5 || isempty(nmc),    nmc = 40000; end
if nargin < 6 || isempty(seed),   seed = 61; end

A = atwo();  M = lingermodel();  P = linger_style();
traps = [-pi/2, pi/2];  th = [1 1]*prm.theta;  gs = prm.sm*prm.sig;
gam = min(ugrid*gs, prm.gmax);  gmc = min(umc*gs, prm.gmax);

cfile = fullfile(outdir, 'fig2f_cache.mat');
% The stamp tag is a fixed part of the cache identity.
key = M.stamp('fig6c', prm, ugrid, umc, nmc, seed);
R = [];
if isfile(cfile)
    C = load(cfile);
    if isfield(C, 'key') && strcmp(C.key, key), R = C.R; fprintf('fig2f: loaded stamped cache\n'); end
end
if isempty(R)
    [pi0, T0, Tb0, L0] = A.split(traps, th, [], prm.sig, prm.beta0, 0, prm.sm, prm.v);
    nu = numel(ugrid);
    Pi = zeros(nu, 2);  Tc = zeros(nu, 2);  Tb = zeros(nu, 1);  L = zeros(nu, 2);
    for i = 1:nu
        [Pi(i,:), Tc(i,:), Tb(i), L(i,:)] = A.split(traps, th, traps(1), prm.sig, ...
                                                    prm.beta0, gam(i), prm.sm, prm.v);
    end
    nm = numel(umc);
    Lmc = zeros(nm, 2);  Lse = zeros(nm, 2);  Tmc = zeros(nm, 2);  Tse = zeros(nm, 2);
    for i = 1:nm
        [who, tabs] = A.mc(nmc, traps, th, traps(1), prm.sig, prm.beta0, gmc(i), ...
                           prm.sm, prm.v, seed + i);
        ok = who > 0;  n = nnz(ok);  t = tabs(ok);  w = who(ok);
        b = mean(t);
        for j = 1:2
            I = double(w == j);  a = mean(I);
            cv = cov(I, t);
            vv = (cv(1,1)/b^2 + a^2*cv(2,2)/b^4 - 2*a*cv(1,2)/b^3)/n;
            Lmc(i,j) = a/b;  Lse(i,j) = sqrt(vv);
            Tmc(i,j) = mean(t(w == j));  Tse(i,j) = std(t(w == j))/sqrt(nnz(w == j));
        end
        if n < nmc, warning('fig2f:capped', '%d walkers unabsorbed at u = %g', nmc - n, umc(i)); end
    end
    R = struct('u',ugrid,'pi',Pi,'Tc',Tc,'Tbar',Tb,'Lam',L, ...
               'pi0',pi0,'T0',T0,'Tbar0',Tb0,'Lam0',L0, ...
               'u_mc',umc,'Lam_mc',Lmc,'Lam_se',Lse,'T_mc',Tmc,'T_se',Tse,'n',nmc);
    save(cfile, 'R', 'key');
end

g1 = R.Lam(:,1)/R.Lam0(1);  g2 = R.Lam(:,2)/R.Lam0(2);
m1 = R.Lam_mc(:,1)/R.Lam0(1);  m2 = R.Lam_mc(:,2)/R.Lam0(2);
e1 = 1.96*R.Lam_se(:,1)/R.Lam0(1);  e2 = 1.96*R.Lam_se(:,2)/R.Lam0(2);
dev = [(m1 - interp1(R.u, g1, R.u_mc(:)))./(e1/1.96), ...
       (m2 - interp1(R.u, g2, R.u_mc(:)))./(e2/1.96)];
fprintf('fig2f: no cue, pi = (%.3f, %.3f), T_j = (%.1f, %.1f), Lambda_j = %.5f vs lam0 = %.5f\n', ...
        R.pi0, R.T0, R.Lam0(1), prm.lam0);
fprintf(['fig2f: full reliance, Lambda_1 x %.3f, Lambda_2 x %.3f (%.1f%% loss), total x %.3f\n' ...
         '       pi_2 %.3f -> %.3f, T_2 %.1f -> %.1f, T_1 %.1f -> %.1f, Tbar %.1f -> %.1f\n'], ...
        g1(end), g2(end), 100*(1 - g2(end)), sum(R.Lam(end,:))/sum(R.Lam0), ...
        R.pi0(2), R.pi(end,2), R.T0(2), R.Tc(end,2), R.T0(1), R.Tc(end,1), R.Tbar0, R.Tbar(end));
fprintf('fig2f: Monte Carlo (%d walkers) at full reliance, Lambda ratios %.3f +- %.3f, %.3f +- %.3f; T_2 = %.1f +- %.1f\n', ...
        R.n, m1(end), e1(end), m2(end), e2(end), R.T_mc(end,2), 1.96*R.T_se(end,2));
fprintf('fig2f: Monte Carlo within %.1f SE of the arc solution\n', max(abs(dev(:))));

% ---- panel --------------------------------------------------------------
cL = P.greens(5,:);  cR = P.blues(3,:);
h = P.figure(P);
pos = [0.17 0.17 0.66 0.72];
% Two overlaid axes rather than yyaxis, which the assembly cannot copy.
ax2 = axes('Parent', h, 'Position', pos);  hold(ax2, 'on');
ax  = axes('Parent', h, 'Position', pos);  hold(ax, 'on');
plot(ax, [0 1], [1 1], ':', 'Color', P.grey, 'LineWidth', P.ALW);
plot(ax, R.u, g1, '-', 'Color', cL, 'LineWidth', P.LW);
plot(ax, R.u, g2, '--', 'Color', cL, 'LineWidth', P.LW);
errorbar(ax, R.u_mc, m1, e1, 'o', 'Color', cL, 'MarkerFaceColor','w', ...
         'LineWidth', P.ALW, 'MarkerSize', 12, 'CapSize', 10, 'LineStyle','none');
errorbar(ax, R.u_mc, m2, e2, 'o', 'Color', cL, 'MarkerFaceColor','w', ...
         'LineWidth', P.ALW, 'MarkerSize', 12, 'CapSize', 10, 'LineStyle','none');
xlim(ax, [0 1]);  ylim(ax, [0 4]);
xlabel(ax, 'reliance $\gamma/(s_m\sigma)$', 'Interpreter','latex', 'FontSize', P.LSZ);
ylabel(ax, 'arrival rate $\Lambda_j/\Lambda_j^{0}$', 'Interpreter','latex', 'FontSize', P.LSZ);
P.axes(P, ax);
set(ax, 'Color', 'none', 'YColor', cL, 'Position', pos);
P.label(P, ax, 0.96, 3.62, '$\Lambda_1$, advertised', 'Color', cL, ...
        'HorizontalAlignment','right', 'VerticalAlignment','bottom');
P.label(P, ax, 0.55, 0.62, '$\Lambda_2$, rival', 'Color', cL, ...
        'HorizontalAlignment','center');

plot(ax2, R.u, R.Tc(:,2), '-', 'Color', cR, 'LineWidth', P.LW);
plot(ax2, R.u, R.Tc(:,1), ':', 'Color', cR, 'LineWidth', P.LW);
xlim(ax2, [0 1]);  ylim(ax2, [0 120]);
ylabel(ax2, 'conditional time $T_j$', 'Interpreter','latex', 'FontSize', P.LSZ);
P.axes(P, ax2);
set(ax2, 'YAxisLocation', 'right', 'XColor', 'none', 'YColor', cR, 'Position', pos);
P.label(P, ax2, 0.03, 112, '$T_1\simeq T_2$', 'Color', cR, ...
        'HorizontalAlignment','left', 'VerticalAlignment','middle');

P.export(P, 'fig2f', outdir, h);
close(h);

out = R;  out.gain1 = g1;  out.gain2 = g2;  out.mcdev = dev;
end
