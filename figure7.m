function out = figure7(varargin)
%FIGURE5 Reproduce Figure 7 (search for one donor across parameters).
%
% Usage
%   figure7
%   figure7('fast',true)                  % coarse grids, same parameters
%   figure7('outdir','output/figure7')
%
% Panel map
%   (a) fig5a   attainable gain G(p) over cue offset and donor absorptivity,
%               with the critical offset p*(kbar) and the support bound
%   (b) fig5b   critical offset and gain ceiling against absorptivity at
%               beta0 = 1 and 0.5
%   (c) fig5c   conditional discovery times at two uncued donors of unequal
%               absorptivity, arc solution and Monte Carlo
%
% One gainTable per absorptivity feeds (a) and the beta0 = 1 arm of (b). The
% beta0 = 0.5 arm needs only p* and G(0), so it runs on a three-point offset
% grid, which leaves both unchanged since p* is bisected on [0,1.2] and G(0)
% is the p = 0 node. Both sweeps are cached behind stamped keys.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'src')));

ip = inputParser;
addParameter(ip,'fast',false);
addParameter(ip,'outdir',fullfile('output','figure7'));
addParameter(ip,'figstem',fullfile('output','figs','fig7'));
parse(ip,varargin{:});
fast   = logical(ip.Results.fast);
OUTDIR = ip.Results.outdir;
if ~isfolder(OUTDIR), mkdir(OUTDIR); end
figdir = fileparts(ip.Results.figstem);
if ~isempty(figdir) && ~isfolder(figdir), mkdir(figdir); end
close all

VER = 'fig5-sweep-v1';
A = afpt();  M = lingermodel();  P = linger_style();
prm0 = A.mkprm('kbar',0.0314,'sig',0.8,'beta0',1,'sm',1,'v',1);

% 0.0314 and 0.8 are forced nodes so the quoted values are not interpolated.
if fast
    kbars = unique([logspace(log10(0.025), log10(0.8), 5), 0.0314, 0.8]);
    pgrid = linspace(0, 1.2, 13);
else
    kbars = unique([logspace(log10(0.025), log10(0.8), 11), 0.0314, 0.8]);
    pgrid = linspace(0, 1.2, 25);
end
pgridB = [0 0.6 1.2];
beta0s = [1.0 0.5];

S1 = cachedSweep(A, M, prm0, kbars, beta0s(1), pgrid,  VER, ...
                 fullfile(OUTDIR,'figure7_sweep_b1.mat'));
S2 = cachedSweep(A, M, prm0, kbars, beta0s(2), pgridB, VER, ...
                 fullfile(OUTDIR,'figure7_sweep_b05.mat'));

% Reference values at the representative parameters, as warnings.
iop = find(abs(S1.kbars - prm0.kbar) < 1e-12, 1);
if abs(S1.pstar(iop)/prm0.sig - 0.878) > 5e-3 || abs(S1.G0(iop) - 3.434) > 5e-3
    warning('figure7:reference', 'p*/sigma = %.4f and G(0) = %.4f, references 0.878 and 3.434', ...
            S1.pstar(iop)/prm0.sig, S1.G0(iop));
end

F5a = panel_a(S1, P, prm0, OUTDIR);
F5b = panel_b(S1, S2, P, prm0, OUTDIR);
F5c = panel_c(prm0, P, OUTDIR);

% ---- numbers a caption quotes -------------------------------------------
fprintf('\nFigure 7 caption numbers\n');
for S = [S1 S2]
    for kq = [0.0314 0.8]
        [~, i] = min(abs(S.kbars - kq));
        fprintf('  beta0 = %.2f, kbar = %.4g:  p*/sigma = %.3f,  G(0) = %.3f,  pcap/sigma = %.3f\n', ...
            S.beta0, S.kbars(i), S.pstar(i)/prm0.sig, S.G0(i), S.pcap(i)/prm0.sig);
    end
    fprintf('  beta0 = %.2f over kbar in [%.3g, %.3g]:  p*/sigma %.3f to %.3f,  G(0) %.3f to %.3f\n', ...
        S.beta0, min(S.kbars), max(S.kbars), min(S.pstar)/prm0.sig, ...
        max(S.pstar)/prm0.sig, min(S.G0), max(S.G0));
end
fprintf('  G(p) maximum on the (a) grid %.3f, at p = 0 and kbar = %.4g\n', ...
    max(S1.G(:)), S1.kbars(find(max(S1.G,[],2) == max(S1.G(:)), 1)));

fprintf('  bounded-law transit penalty is at most exp(2 beta0 p/v) = %.2f at p = 0.75\n', ...
    exp(2*prm0.beta0*0.75/prm0.v));

% ---- assembly ------------------------------------------------------------
% One row of three at full textwidth, like Figures 2-3, at the default letter size.
panels = fullfile(OUTDIR, {'fig5a.fig','fig5b.fig','fig5c.fig'});
tiles  = [0 0 22 19; 22 0 22 19; 44 0 22 19];
assemble_figure(ip.Results.figstem, panels, tiles, [66 19]);

out = struct('S1',S1,'S2',S2,'F5a',F5a,'F5b',F5b,'F5c',F5c);
save(fullfile(OUTDIR,'figure7_results.mat'),'-struct','out');
fprintf('Figure 7 of the article (driver figure7) complete. Panel order: fig5a, fig5b, fig5c.\n');
if nargout == 0, clear out; end
end

% ------------------------------------------------------------------------
function S = cachedSweep(A, M, prm0, kbars, beta0, pgrid, ver, cachefile)
key = M.stamp(ver, struct('sig',prm0.sig,'sm',prm0.sm,'v',prm0.v, ...
              'beta0',beta0), kbars, pgrid);
if isfile(cachefile)
    C = load(cachefile);
    if isfield(C,'key') && strcmp(C.key, key)
        S = C.S;
        fprintf('figure7: loaded stamped sweep, beta0 = %.2f\n', beta0);
        return
    end
end
S = sweep(A, M, prm0, kbars, beta0, pgrid);
save(cachefile, 'S', 'key');
end

function S = sweep(A, M, prm0, kbars, beta0, pgrid)
% One gainTable per absorptivity. A point whose arc solve raised a
% singular-matrix warning, or whose p* exceeds the support bound, is dropped
% and reported.
nk = numel(kbars);  np = numel(pgrid);
S = struct('beta0',beta0,'kbars',kbars,'pgrid',pgrid, ...
           'G',nan(nk,np),'Graw',nan(nk,np),'gam',nan(nk,np), ...
           'atcap',false(nk,np),'pstar',nan(1,nk),'pcap',nan(1,nk), ...
           'G0',nan(1,nk),'illcond',false(1,nk),'oversup',false(1,nk));
SING = {'MATLAB:singularMatrix','MATLAB:nearlySingularMatrix', ...
        'MATLAB:illConditionedMatrix'};
ws = warning;
cleanup = onCleanup(@() warning(ws));
for k = 1:numel(SING), warning('off', SING{k}); end
t0 = tic;
for ik = 1:nk
    prm = A.mkprm('kbar',kbars(ik),'sig',prm0.sig,'beta0',beta0, ...
                  'sm',prm0.sm,'v',prm0.v);
    lastwarn('');
    G = M.gainTable(prm, pgrid);
    [~, wid] = lastwarn;
    S.illcond(ik) = any(strcmp(wid, SING));
    S.oversup(ik) = G.pstar > prm0.sm*prm0.sig*(1 + 1e-9);
    if S.illcond(ik) || S.oversup(ik), continue, end
    S.G(ik,:) = G.G;  S.Graw(ik,:) = G.Graw;  S.gam(ik,:) = G.gam;
    S.atcap(ik,:) = G.atcap;
    S.pstar(ik) = G.pstar;  S.pcap(ik) = G.pcap;  S.G0(ik) = G.Graw(1);
end
nbad = nnz(S.illcond | S.oversup);
fprintf('figure7: beta0 = %.2f sweep, %d absorptivities x %d offsets, %.0f s, %d dropped\n', ...
        beta0, nk, np, toc(t0), nbad);
end

% =========================================================================
% Panel (a), fig5a. Attainable gain over offset and absorptivity.
% =========================================================================
function out = panel_a(S, P, prm0, outdir)
%PANEL_A  Attainable gain G(p) over cue offset and donor absorptivity.
%   Panel (a). S is the beta0 = 1 sweep of figure4. The solid curve is the
%   critical offset p*(kbar), bisected in gainTable rather than read off the
%   grid. The dashed line is the support bound p = sigma, and the dotted line
%   and marker are the representative parameters. Dotted contours are levels of G. The
%   dash-dotted curve is pcap, below which the optimal reliance sits on the
%   ceiling gamma = s_m sigma.

ok  = ~(S.illcond | S.oversup);
kb  = S.kbars(ok);
Gm  = S.G(ok,:);
xs  = S.pgrid/prm0.sig;
ps  = S.pstar(ok)/prm0.sig;
kop = prm0.kbar;

% p* is smooth in log kbar, so a pchip through the nodes draws it.
kf  = logspace(log10(min(kb)), log10(max(kb)), 200);
psf = pchip(log(kb), ps, log(kf));
psop = interp1(kb, ps, kop);
pc   = S.pcap(ok)/prm0.sig;
pcf  = pchip(log(kb), pc, log(kf));

h = P.figure(P);  ax = gca;
hs = pcolor(ax, xs, kb, Gm);  set(hs,'EdgeColor','none');
shading(ax,'interp');  colormap(ax, P.plasma(256));
clim(ax, [1, ceil(10*max(Gm(:)))/10]);
set(ax,'YScale','log');
cb = P.cbar(P, 'attainable gain $G(p)$');
hold(ax,'on');

[Cc, hc] = contour(ax, xs, kb, Gm, [1.5 2 2.5 3], 'k:', 'LineWidth', P.ALW);
clabel(Cc, hc, 'FontSize', P.FSZ-6, 'Interpreter','latex', 'LabelSpacing', 1e4);

% Round limits, so the end ticks survive a logspace node that misses 0.025 by rounding.
yl = [0.025 0.8];
plot(ax, [1 1], yl, '--', 'Color','w', 'LineWidth', P.ALW+0.5);
plot(ax, [0 max(xs)], [kop kop], ':', 'Color','w', 'LineWidth', P.ALW+0.5);
plot(ax, pcf, kf, '-.', 'Color','w', 'LineWidth', P.ALW+0.5);
plot(ax, psf, kf, '-', 'Color','w', 'LineWidth', P.LW+2.4);
plot(ax, psf, kf, '-', 'Color',[0.10 0.10 0.10], 'LineWidth', P.LW);
plot(ax, psop, kop, 'o', 'MarkerSize', 14, 'MarkerFaceColor','w', ...
     'MarkerEdgeColor',[0.10 0.10 0.10], 'LineWidth', P.ALW);

P.label(P, ax, 1.03, 0.36, 'support bound', 'Color','w', 'Rotation', 90, ...
        'HorizontalAlignment','center', 'VerticalAlignment','top');
P.label(P, ax, 0.5*(interp1(kf, psf, 0.3) + interp1(kf, pcf, 0.3)), 0.3, ...
        '$p^\ast(\overline\kappa)$', 'Color','w', 'Rotation', 90, ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');
P.label(P, ax, interp1(kf, pcf, 0.07) - 0.02, 0.07, '$\gamma^\ast=s_m\sigma$', 'Color','w', ...
        'Rotation', 90, 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
P.label(P, ax, 0.04, kop*1.06, 'representative', 'Color','w', ...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom');

xlim(ax, [0 max(xs)]);  ylim(ax, yl);
set(ax, 'YTick', [0.025 0.05 0.1 0.2 0.4 0.8], ...
        'YTickLabel', {'0.025','0.05','0.1','0.2','0.4','0.8'}, ...
        'XTick', 0:0.25:1.5);
xlabel(ax, 'cue offset $p/\sigma$', 'Interpreter','latex', 'FontSize', P.LSZ);
ylabel(ax, 'donor absorptivity $\overline\kappa$', 'Interpreter','latex', 'FontSize', P.LSZ);
P.axes(P, ax);
set(ax, 'Position', [0.17 0.16 0.62 0.80]);  cb.Position = [0.81 0.16 0.035 0.80];

P.export(P, 'fig5a', outdir, h);

out = struct('kbars',kb,'poff',xs,'G',Gm,'pstar_norm',ps,'pstar_op',psop,'pcap_norm',pc, ...
             'Gmax',max(Gm(:)));
fprintf('fig5a: G spans %.3f to %.3f; p*/sigma at the representative parameters %.3f\n', ...
        min(Gm(:)), max(Gm(:)), psop);
end

% =========================================================================
% Panel (b), fig5b. Critical offset and gain ceiling.
% =========================================================================
function out = panel_b(S1, S2, P, prm0, outdir)
%PANEL_B  Critical offset p* and gain ceiling G(0) against absorptivity.
%   Panel (b). S1 and S2 are the figure4 sweeps at beta0 = 1 and 0.5.
%   Green on the left axis is p*/sigma, blue on the right axis is G(0). Solid
%   is beta0 = 1, dashed is beta0 = 0.5. Curves are labeled in place. The
%   dotted horizontal line is the support bound and the dotted vertical line
%   the representative parameters, where the beta0 = 1 values are marked.

SS   = [S1 S2];
sty  = {'-','--'};
gcol = [P.greens(6,:); P.greens(4,:)];
bcol = [P.blues(3,:);  P.blues(2,:)];
kop  = prm0.kbar;
kl   = [0.025 0.8];
lgt  = [0.45 0.45 0.45];

% Two overlaid axes rather than yyaxis, since assemble_figure copies axes
% objects and yyaxis axes cannot be copied.
POS = [0.17 0.16 0.64 0.80];
h = P.figure(P);
ax = axes('Parent',h);  hold(ax,'on');
set(ax,'XScale','log');

% ---- left axis, critical offset ------------------------------------------
YL = [0.50 1.02];
plot(ax, kl, [1 1], ':', 'Color', lgt, 'LineWidth', P.ALW);
for ib = 1:2
    S = SS(ib);  ok = ~(S.illcond | S.oversup);
    plot(ax, S.kbars(ok), S.pstar(ok)/prm0.sig, sty{ib}, 'Color', gcol(ib,:), 'LineWidth', P.LW);
end
ps1 = interp1(S1.kbars, S1.pstar/prm0.sig, kop);
plot(ax, kop, ps1, 'o', 'MarkerSize', 12, 'MarkerFaceColor','w', ...
     'MarkerEdgeColor', gcol(1,:), 'LineWidth', P.ALW);
ylim(ax, YL);
set(ax, 'YColor', gcol(1,:), 'YTick', 0.5:0.1:1);
ylabel(ax, 'critical offset $p^\ast/\sigma$', 'Interpreter','latex', 'FontSize', P.LSZ);

P.label(P, ax, 0.037, 0.993, 'support bound $p^\ast=\sigma$', 'Color', lgt, ...
        'VerticalAlignment','top');
P.label(P, ax, 0.2, curveAt(S1, 0.2, prm0) + 0.012, '$p^\ast/\sigma$, $\beta_0=1$', ...
        'Color', gcol(1,:), 'HorizontalAlignment','right', 'VerticalAlignment','bottom');
P.label(P, ax, 0.3, curveAt(S2, 0.3, prm0) - 0.012, '$\beta_0=0.5$', ...
        'Color', gcol(2,:), 'HorizontalAlignment','left', 'VerticalAlignment','top');

% ---- right axis, gain ceiling --------------------------------------------
axL = ax;
xline(axL, kop, ':', 'Color', lgt, 'LineWidth', P.ALW);
xlim(axL, kl);
set(axL, 'XTick', [0.05 0.1 0.2 0.4 0.8], ...
         'XTickLabel', {'0.05','0.1','0.2','0.4','0.8'}, 'XTickLabelRotation', 0);
xlabel(axL, 'donor absorptivity $\overline\kappa$', 'Interpreter','latex', 'FontSize', P.LSZ);
P.axes(P, axL);
set(axL, 'Position', POS);

ax = axes('Parent',h);  hold(ax,'on');
set(ax, 'XScale','log', 'Color','none', 'YAxisLocation','right', ...
        'XColor','none', 'XTick',[]);
for ib = 1:2
    S = SS(ib);  ok = ~(S.illcond | S.oversup);
    plot(ax, S.kbars(ok), S.G0(ok), sty{ib}, 'Color', bcol(ib,:), 'LineWidth', P.LW);
end
G01 = interp1(S1.kbars, S1.G0, kop);
plot(ax, kop, G01, 'o', 'MarkerSize', 12, 'MarkerFaceColor','w', ...
     'MarkerEdgeColor', bcol(1,:), 'LineWidth', P.ALW);
ylim(ax, [1.0 4.8]);
set(ax, 'YColor', bcol(1,:));
ylabel(ax, 'gain ceiling $G(0)$', 'Interpreter','latex', 'FontSize', P.LSZ);

G0at = @(S,k) interp1(log(S.kbars), S.G0, log(k));
P.label(P, ax, 0.25, G0at(S1,0.25) - 0.15, '$G(0)$, $\beta_0=1$', ...
        'Color', bcol(1,:), 'HorizontalAlignment','right', 'VerticalAlignment','top');
P.label(P, ax, 0.2, G0at(S2,0.2) - 0.15, '$\beta_0=0.5$', ...
        'Color', bcol(2,:), 'HorizontalAlignment','left', 'VerticalAlignment','top');

P.label(P, ax, kop*1.05, 1.05, 'representative', 'Color', lgt, ...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom');
xlim(ax, kl);
P.axes(P, ax);
set(ax, 'Position', POS, 'XColor','none');

P.export(P, 'fig5b', outdir, h);

out = struct('kbars',S1.kbars,'pstar_norm',[S1.pstar; S2.pstar]/prm0.sig, ...
             'G0',[S1.G0; S2.G0],'beta0s',[S1.beta0 S2.beta0], ...
             'pstar_op',ps1,'G0_op',G01);
fprintf('fig5b: at kbar = %.4g, p*/sigma = %.3f and G(0) = %.3f (beta0 = 1)\n', ...
        kop, ps1, G01);
end

function y = curveAt(S, k, prm0)
y = interp1(log(S.kbars), S.pstar/prm0.sig, log(k));
end

% =========================================================================
% Panel (c), fig5c. Conditional discovery times at unequal donors.
% =========================================================================
function out = panel_c(prm0, P, outdir, ratios, rmc, nmc, seed)
%PANEL_C  Conditional discovery times at two uncued donors of unequal
%   absorptivity. Panel (c).
%
%   Donors at -pi/2 and +pi/2 with absorptivities kbar*sqrt(r) and kbar/sqrt(r),
%   so their geometric mean stays at the representative value. Curves are the arc
%   solution (atwo) of (T_1 - T_2)/Tbar against r at three reversal rates.
%   Markers are event-driven Monte Carlo at beta0 = 1, bars +-1.96 standard
%   errors. Equal conditional times are what let Eq. (Lambda) split arrivals in
%   proportion to the splitting probabilities. At beta0 = 1 the times agree to
%   a fraction of a percent; they separate only once reversal is frequent.
%
%   The Monte Carlo is cached under outdir behind a stamped key.

if nargin < 4 || isempty(ratios), ratios = logspace(0, log10(16), 25); end
if nargin < 5 || isempty(rmc),    rmc = [2 4 8 16]; end
if nargin < 6 || isempty(nmc),    nmc = 1e5; end
if nargin < 7 || isempty(seed),   seed = 10; end

A2 = atwo();  M = lingermodel();
traps = [-pi/2 pi/2];
betas = [0.1 1 10];
kap   = @(r) prm0.kbar*[sqrt(r) 1/sqrt(r)];
thOf  = @(r) 1 - exp(-kap(r)/prm0.v);

dT = zeros(numel(betas), numel(ratios));
for ib = 1:numel(betas)
    for ir = 1:numel(ratios)
        [~, Tc, Tbar] = A2.split(traps, thOf(ratios(ir)), [], prm0.sig, betas(ib), 0, prm0.sm, prm0.v);
        dT(ib, ir) = 100*(Tc(1) - Tc(2))/Tbar;
    end
end

% ---- Monte Carlo at beta0 = 1 --------------------------------------------
key = M.stamp('fig5c-v1', prm0, rmc, nmc, seed);
cachefile = fullfile(outdir, 'fig5c_mc.mat');
D = [];
if isfile(cachefile)
    C = load(cachefile);
    if isfield(C,'key') && strcmp(C.key, key)
        D = C.D;  fprintf('fig5c: reusing cached Monte Carlo\n');
    end
end
if isempty(D)
    D = struct('d', zeros(size(rmc)), 'se', zeros(size(rmc)), 'pi1', zeros(size(rmc)));
    for i = 1:numel(rmc)
        th = thOf(rmc(i));
        [~, ~, Tbar] = A2.split(traps, th, [], prm0.sig, 1, 0, prm0.sm, prm0.v);
        [who, tabs] = A2.mc(nmc, traps, th, [], prm0.sig, 1, 0, prm0.sm, prm0.v, seed + i);
        a = tabs(who == 1);  b = tabs(who == 2);
        D.d(i)   = 100*(mean(a) - mean(b))/Tbar;
        D.se(i)  = 100*sqrt(var(a)/numel(a) + var(b)/numel(b))/Tbar;
        D.pi1(i) = mean(who == 1);
    end
    save(cachefile, 'D', 'key');
end

% ---- plot ----------------------------------------------------------------
cols = [P.greens(4,:); P.greens(6,:); P.grey];
h = P.figure(P);
ax = axes('Parent', h);  hold(ax, 'on');
yline(ax, 0, ':', 'Color', [0 0 0], 'LineWidth', P.ALW);
for ib = 1:numel(betas)
    plot(ax, ratios, dT(ib,:), '-', 'Color', cols(ib,:), 'LineWidth', P.LW);
end
errorbar(ax, rmc, D.d, 1.96*D.se, 'o', 'Color', P.blues(3,:), 'MarkerFaceColor','w', ...
         'MarkerSize', 12, 'LineWidth', P.ALW, 'CapSize', 10, 'LineStyle','none');
set(ax, 'XScale','log');
xlim(ax, [1 16]);  ylim(ax, [-10 4]);
set(ax, 'XTick', [1 2 4 8 16], 'XTickLabel', {'$1$','$2$','$4$','$8$','$16$'});
xlabel(ax, 'absorptivity ratio $\overline\kappa_1/\overline\kappa_2$', 'Interpreter','latex', 'FontSize', P.LSZ);
ylabel(ax, 'time difference $(T_1-T_2)/\bar T$ (\%)', 'Interpreter','latex', 'FontSize', P.LSZ);
P.label(P, ax, 1.1, 1.6, '$\beta_0=0.1$ and $1$', 'Color', cols(2,:), 'VerticalAlignment','bottom');
P.label(P, ax, 3.0, -6.2, '$\beta_0=10$', 'Color', cols(3,:), 'HorizontalAlignment','left');
P.label(P, ax, 15, 2.6, 'simulation, $\beta_0=1$', 'Color', P.blues(3,:), 'HorizontalAlignment','right');
P.axes(P, ax);
P.export(P, 'fig5c', outdir, h);

out = struct('ratios',ratios,'betas',betas,'dT',dT,'mc',D,'rmc',rmc);
fprintf('fig5c: max |T1-T2|/Tbar over r in [1,16]:  beta0 = 0.1 %.2f%%,  1 %.2f%%,  10 %.2f%%\n', ...
        max(abs(dT),[],2));
fprintf('fig5c: Monte Carlo at beta0 = 1, r = %s:  %s  (+-1.96 SE %s)\n', mat2str(rmc), ...
        mat2str(D.d,3), mat2str(1.96*D.se,3));
end
