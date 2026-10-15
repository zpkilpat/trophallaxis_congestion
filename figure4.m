function out = figure4(varargin)
%FIGURE4 Reproduce the three panels of Figure 4 (recruitment and congestion at
%   one donor) and assemble the lettered figure output/figs/fig4.pdf.
%
% Usage
%   figure4
%   figure4('fast',true)                  % low-cost wiring/smoke test
%   figure4('outdir','output/figure4')
%
% Panel map
%   (a) figL2   stochastic versus substituted-offset recruitment gain
%   (b) figL3   feeding-time landscape versus abandonment and traffic ratio
%   (c) figL4   attainable speedup and lingering-regime boundary
%
% The single-donor gain tables come from linger_setup and are cached with a
% parameter stamp. A normal call uses the paper settings. 'fast' uses coarser
% grids to check that the code runs end to end, and never changes a model
% parameter.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'src')));

ip = inputParser;
addParameter(ip,'fast',false);
addParameter(ip,'outdir',fullfile('output','figure4'));
parse(ip,varargin{:});
fast = logical(ip.Results.fast);
OUTDIR = ip.Results.outdir;
if ~isfolder(OUTDIR), mkdir(OUTDIR); end
close all

S = linger_setup(OUTDIR, fast, 'figure4');
G = S.G;  W = S.W;  psList = S.psList;  ps = S.ps;  cfg = S.cfg;

if fast
    rgrid3 = logspace(log10(0.15),log10(5.0),14);
    rgrid4 = logspace(log10(0.15),log10(5.0),10);
    dgrid3 = logspace(-2.5,1.0,20);
    qgrid4 = linspace(0.10,0.95,8);
    dgrid4 = logspace(-2.5,1.0,14);
else
    rgrid3 = logspace(log10(0.15),log10(5.0),34);
    rgrid4 = logspace(log10(0.15),log10(5.0),61);
    dgrid3 = logspace(-2.5,1.0,44);
    qgrid4 = linspace(0.10,0.95,35);
    dgrid4 = logspace(-2.5,1.0,25);
end

L2 = panel_a(G,psList,OUTDIR);
L3 = panel_b(ps,cfg,OUTDIR,dgrid3,rgrid3,0.5);
L4 = panel_c(G,W,cfg,OUTDIR,rgrid4,qgrid4,cfg.p0,dgrid4);

out = struct('G',G,'W',W,'psList',{psList},'L2',L2,'L3',L3,'L4',L4);
save(fullfile(OUTDIR,'figure4_results.mat'),'-struct','out');

% Lettered figure, one row of three panel canvases at the panel size.
figdir = fullfile('output','figs');
if ~isfolder(figdir), mkdir(figdir); end
T = [22 19];
tiles = [0 0 T; T(1) 0 T; 2*T(1) 0 T];
panels = fullfile(OUTDIR, strcat({'figL2','figL3','figL4'}, '.fig'));
assemble_figure(fullfile(figdir,'fig4'), panels, tiles, [3*T(1) T(2)]);

fprintf('Figure 4 complete. Panel order: figL2, figL3, figL4.\n');
end

% =========================================================================
% Panel (a), figL2. Stochastic versus typical-offset recruitment gain.
% =========================================================================
function out = panel_a(G, psList, outdir)
%PANEL_A  Averaged recruitment gain over scenting population.
%   Psi~(Fbar) from the Poisson-Gaussian average has no threshold and leaves
%   unity with slope g1 - 1 > 0. Substituting the mean offset into G predicts
%   a threshold number of scenters below which recruitment does nothing. The
%   shaded region is what the substitution discards.
%
%   The substitution is pinned at unity for Fbar < (p0/p*)^2, where its
%   argument leaves the cue support. This threshold is ps.Fsub and is printed.
%   It is set by p*, the offset at which the gain reaches unity, not by sigma.
%
%   Line style marks p0 and hue marks the family, keyed by a two-row legend
%   drawn in the axes. The ratio Psi~ / G(p0/sqrt(Fbar)) is computed and
%   printed but not drawn.
%
%   Psi~ uses the reliance chosen before the offset is drawn (ps.psi). The
%   substitution curve uses the ridge table G.G, the reliance optimal at the
%   substituted offset. This is the deterministic shortcut being compared, so
%   keep G.G here rather than ps.gam.
%
%   psList is a cell array of psiTable structs at different p0. Line styles
%   and colors come from P.p0sty and P.p0cols so panel (a), Fig. 5(c), Fig. 5(a) and
%   Fig. 5(b) agree. ps.g1 must exceed 1 for every p0 in the list.
if nargin < 3, outdir = 'output'; end
P = linger_style();
h = P.figure(P);
ax = gca;  hold(ax, 'on')
Fb = logspace(-1, log10(40), 300);
G0 = psList{1}.G0;
np = numel(psList);
ex = nan(np, numel(Fb));  nv = nan(np, numel(Fb));
for i = 1:np
    ps = psList{i};
    ex(i,:) = interp1(ps.F, ps.psi, Fb);
    nv(i,:) = max(interp1(G.p, G.G, min(ps.p0./sqrt(Fb), max(G.p))), 1);
end

yline(ax, 1, '-', 'Color', [0.70 0.72 0.74], 'LineWidth', P.ALW);
for i = 1:np
    plot(ax, Fb, ex(i,:), P.p0sty{i}, 'Color', P.p0cols(i,:), 'LineWidth', P.LW);
end

set(ax,'XScale','log');
xlabel('scenting population $\bar F$','Interpreter','latex','FontSize',P.LSZ);
ylabel('recruitment gain','Interpreter','latex','FontSize',P.LSZ);
xlim([min(Fb) max(Fb)*4.2]);  ylim([1 3]);
P.axes(P);
set(ax, 'Position', [0.17 0.15 0.78 0.73], 'YTick', 1:0.5:3, 'XTick', [0.1 1 10]);

% Curve labels at the right edge, where the three curves are most separated.
% No legend, since line style and hue both encode p0.
for i = 1:np
    P.label(P, ax, Fb(end)*1.12, ex(i,end), sprintf('$p_0=%.2f$', psList{i}.p0), ...
            'FontSize', P.FSZ-2, 'Color', P.p0cols(i,:), ...
            'HorizontalAlignment','left', 'VerticalAlignment','middle');
end

rpk = nan(1,np);  Fpk = nan(1,np);
for i = 1:np
    [rpk(i), j] = max(ex(i,:) ./ nv(i,:));  Fpk(i) = Fb(j);
end

out = struct('g1', cellfun(@(s) s.g1, psList), ...
             'p0', cellfun(@(s) s.p0, psList), ...
             'Fsub', cellfun(@(s) s.Fsub, psList), ...
             'ratiopeak', rpk, 'ratiopeakF', Fpk, ...
             'G0', G0, 'pstar', G.pstar);
fprintf('single-scenter gains g1:');  fprintf(' %.3f', out.g1);  fprintf('\n');
fprintf('substitution pinned at unity below Fbar =');
fprintf(' %.2f', out.Fsub);  fprintf('  (p* = %.4f)\n', G.pstar);
fprintf('peak ratio Psi~/G(p0/sqrt(Fbar)):');
for i = 1:np, fprintf('  %.2f at Fbar = %.2f', rpk(i), Fpk(i)); end
fprintf('\n');
if any(out.g1 <= 1)
    warning('figL2:g1NotAboveOne', ...
        ['g1 <= 1 at p0 = %s; the no-threshold argument fails there, since ' ...
         'it rests on Psi~''(0) = g1 - 1 > 0'], ...
        mat2str(out.p0(out.g1 <= 1), 3));
end
if abs(out.g1(1) - 1.644) < 5e-3
    warning('figL2:staleReliance', ...
        ['g1 = %.3f at the reference p0, which is the pre-correction value. ' ...
         'psiTable is averaging along the ridge again, i.e. crediting a ' ...
         'searcher with the realized peak.  Expect 1.552.  Check the cache ' ...
         'stamp reads v2-reliance-outside-average.'], out.g1(1));
end
% Inset, a scenting cluster at the donor. The solid tent is one realization
% of the aggregate cue and the dashed outlines are others, so the peak offset
% is drawn as a random quantity.
inset_cluster(h, P);
P.export(P, 'figL2', outdir, h);
end

% ------------------------------------------------------------------------
function inset_cluster(h, P)
%INSET_CLUSTER  Scenting cluster sketch, in the empty region under the p0 labels.
K = inset_kit();
ax = K.axes(h, P, [4.2 5.8 5.6 4.0]);
y0 = 1.0;  xd = 3.0;
for xg = [1.95 3.4], K.ghost(ax, xg, 1.4, y0, 2.3); end
K.tent(ax, 2.55, 1.4, y0, 2.3, 0.30, 0.93);
K.track(ax, 0.1, 5.5, y0);
K.donor(ax, xd, y0, 0.25);
B = [1.8 y0+0.3 -0.25; 2.3 y0-0.32 0.35; 3.75 y0+0.32 pi+0.3; 4.2 y0-0.28 pi-0.25];
for k = 1:size(B,1), K.bee(ax, B(k,1), B(k,2), B(k,3), 0.8, 'scent'); end
end

% =========================================================================
% Panel (b), figL3. Feeding-time landscape over abandonment and traffic.
% =========================================================================
function out = panel_b(ps, cfg, outdir, dgrid, rgrid, q)
%PANEL_B  Feeding time relative to no lingering, over abandonment rate and
%   traffic ratio. Diverging map centered on unity. Blue is faster than never
%   lingering, red is slower. The solid curve is delta*(rho). The dashed
%   contour is the break-even locus, above which lingering never pays.
%
%   delta*(rho) is a trough, the minimum of the ratio at each rho. It is
%   refined off the grid with refineDelta so the curve does not step.
%
%   rgrid is shared with panel (c) so both panels cover the same decades of rho.
%
%   rho_detach here is where any benefit disappears. Panel (c)'s rho*(q) is the
%   one percent level. They are close at q = 0.5 but are different quantities,
%   and out.detach_threshold records which is which.
%
%   Detachment is extrapolated from the paying side. The benefit 1 - Rmin goes
%   continuously to zero, so the last paying grid point understates the
%   boundary by up to one spacing. The first non-paying point has dstar = NaN
%   and cannot serve as a node. Extrapolate from the last two resolved rho and
%   clamp to the bracket. Quote the bracket with the number.
%
%   The delta -> infinity limit depends on rho. At delta = inf no bee scents
%   and every service event is productive. A(0) = lam0*N = rho*mu, so for
%   rho <= 1 the queue stays empty and Tinf = L_q/lam0 exactly. Above rho = 1
%   the queue binds and Tinf approaches the service bound q*N/mu.
%
%   Tinf must equal log(1/(1-q))/lam0 for rho <= 1 and must not fall below
%   max(log(1/(1-q))/lam0, q*N/mu) at any rho. Test Tinf itself, not the last
%   column of dgrid, where some lingering remains and Tq sits below Tinf.

if nargin < 3, outdir = 'output'; end
if nargin < 4, dgrid = logspace(-2.5, 1.0, 44); end
if nargin < 5, rgrid = logspace(log10(0.15), log10(5.0), 34); end
if nargin < 6, q = 0.5; end

M = lingermodel();  P = linger_style();
R = zeros(numel(rgrid), numel(dgrid));
Tinf = zeros(size(rgrid));  dstar = nan(size(rgrid));  Rmin = nan(size(rgrid));
for i = 1:numel(rgrid)
    c = cfg;  c.lam0 = rgrid(i)*cfg.mu/cfg.N;
    f = @(d) M.Tq(d, ps, c, q);
    Tinf(i) = f(1e4);
    for j = 1:numel(dgrid)
        R(i,j) = f(dgrid(j))/Tinf(i);
    end
    [dopt, Topt] = M.refineDelta(f, dgrid, R(i,:)*Tinf(i));
    if Topt < Tinf(i)*(1 - 1e-6)
        dstar(i) = dopt;  Rmin(i) = Topt/Tinf(i);
    end
end

rdetach = detachRho(rgrid, dstar, Rmin);

% The drawn delta*(rho) is continued past the last grid row on a finer rho
% bracket and a wider delta range, then clipped at the right edge of the map,
% so it leaves the panel where the optimum passes delta = max(dgrid) instead
% of stopping at the grid spacing. Printed quantities use the grid values.
[dcurve, rcurve] = continueRidge(M, ps, cfg, q, dgrid, rgrid, dstar);

h = P.figure(P);
ax = gca;
% pcolor with interpolated shading, since imagesc cells are visible at print
% size. Contours are drawn on the swept grid.
hs = pcolor(ax, log10(dgrid), log10(rgrid), R);
set(hs, 'EdgeColor','none');  shading(ax, 'interp');
colormap(ax, P.divmap(256));
m = max(abs(R(:) - 1));  clim(ax, [1-m, 1+m]);
cb = P.cbar(P, 'feeding time ratio $T_q(\delta)/T_q(\infty)$');
cb.Label.FontSize = P.LSZ-4;
hold(ax, 'on')
plot(ax, log10(dcurve), log10(rcurve), '-', 'Color','w','LineWidth', P.LW+1.5);
plot(ax, log10(dcurve), log10(rcurve), '-', 'Color',P.greens(6,:),'LineWidth', P.LW);
contour(ax, log10(dgrid), log10(rgrid), R, [1 1], 'k--','LineWidth',P.ALW);

xlabel('abandonment rate magnitude $\log_{10}\delta$','Interpreter','latex','FontSize',P.LSZ);
ylabel('traffic magnitude $\log_{10}\rho$','Interpreter','latex','FontSize',P.LSZ);
P.axes(P);
set(ax, 'YTick', -0.8:0.2:0.6, 'Position', [0.17 0.15 0.62 0.73]);  cb.Position = [0.81 0.15 0.035 0.73];
P.label(P, ax, 0.92, 0.17, 'break even', 'FontSize', P.FSZ+4, ...
        'HorizontalAlignment','right', 'VerticalAlignment','bottom');
P.label(P, ax, -0.30, -0.21, 'optimal $\delta^\ast$', 'FontSize', P.FSZ+4, ...
        'Color', P.greens(6,:), 'HorizontalAlignment','left');

Lq    = log(1/(1-q));
Tserv = q*cfg.N/cfg.mu;

out = struct('dgrid',dgrid,'rgrid',rgrid,'R',R,'dstar',dstar,'Rmin',Rmin, ...
             'Tinf',Tinf,'q',q,'rho_detach',rdetach,'detach_threshold',0, ...
             'Tserv',Tserv,'best_speedup',100*(1-min(Rmin)));
k = find(~isnan(dstar), 1, 'last');
out.detach_bracket = [rgrid(k), rgrid(min(k+1,numel(rgrid)))];
fprintf(['best speedup %.1f%%, detaches above rho = %.2f ' ...
         '(any benefit, bracketed by [%.2f, %.2f]; figL4 uses 1%%)\n'], ...
        out.best_speedup, out.rho_detach, out.detach_bracket(1), ...
        out.detach_bracket(2));

% ---- the delta -> infinity limit, tested where each statement holds ------
lam   = rgrid*cfg.mu/cfg.N;
Tsrch = Lq./lam;                        % search limited, exact for rho <= 1
sub   = rgrid <= 1;
if any(sub), chk = max(abs(Tinf(sub).*lam(sub)/Lq - 1)); else, chk = 0; end
if chk > 1e-3
    warning('figL3:TinfNotExact', ...
        ['Tq(inf) departs from log(1/(1-q))/lam0 by %.1e at rho <= 1, where ' ...
         'the queue never binds and the limit is exact'], chk);
end

% Both are lower bounds at every rho. Queueing only slows a run, and no bee
% is fed faster than mu. A violation is an error in the reduced system.
gap = Tinf - max(Tsrch, Tserv);
if min(gap) < -1e-3*Tserv
    [~, kb] = min(gap);
    warning('figL3:TinfBelowBound', ...
        ['Tq(inf) falls %.2f below max(L_q/lam0, q*N/mu) at rho = %.2f; the ' ...
         'reduced system is feeding faster than either limit allows'], ...
        -gap(kb), rgrid(kb));
end
out.Tinf_dev_sub = chk;
out.Tinf_end     = Tinf(end);
fprintf(['Tq(inf): exact to %.1e for rho <= 1;  at rho = %.2f it is %.1f ' ...
         'against the service bound %.1f\n'], chk, rgrid(end), Tinf(end), Tserv);

P.export(P, 'figL3', outdir, h);
end

% ------------------------------------------------------------------------
function [dc, rc] = continueRidge(M, ps, cfg, q, dgrid, rgrid, dstar)
%CONTINUERIDGE  delta*(rho) on the grid rows, extended between the last row
%   with an optimum and the next one, and clipped at delta = max(dgrid).
k  = find(~isnan(dstar), 1, 'last');
dc = dstar(1:k);  rc = rgrid(1:k);
if isempty(k) || k == numel(rgrid), return; end
rf   = logspace(log10(rgrid(k)), log10(rgrid(k+1)), 14);  rf = rf(2:end-1);
dext = logspace(log10(dgrid(1)), 3, 60);
for i = 1:numel(rf)
    c = cfg;  c.lam0 = rf(i)*cfg.mu/cfg.N;
    f = @(d) M.Tq(d, ps, c, q);
    Ti = f(1e4);
    [dopt, Topt] = M.refineDelta(f, dext, arrayfun(f, dext));
    if Topt >= Ti*(1 - 1e-6), break, end
    dc(end+1) = dopt;  rc(end+1) = rf(i); %#ok<AGROW>
end
j = find(dc > dgrid(end), 1);
if ~isempty(j) && j > 1
    u = interp1(log10(dc([j-1 j])), log10(rc([j-1 j])), log10(dgrid(end)));
    dc = [dc(1:j-1) dgrid(end)];  rc = [rc(1:j-1) 10^u];
end
end

function r = detachRho(rgrid, dstar, Rmin)
%DETACHRHO  Where the attainable benefit reaches zero, extrapolated from the
%   last two resolved rho and clamped to the grid bracket. Returns NaN if the
%   boundary is not on the grid.
idx = find(~isnan(dstar));
if numel(idx) < 2 || idx(end) == numel(rgrid), r = NaN; return; end
k = idx(end);  km = idx(end-1);
ben = 1 - Rmin([km k]);                 % benefit, positive and decreasing
u   = log10(rgrid([km k]));
if ben(1) > ben(2) && ben(2) > 0
    r = 10^interp1(ben, u, 0, 'linear', 'extrap');
else
    r = rgrid(k);                       % not resolvable, fall back to the grid
end
r = min(max(r, rgrid(k)), rgrid(k+1));  % the grid cannot support more than this
end

% =========================================================================
% Panel (c), figL4. Attainable speedup and the lingering boundary.
% =========================================================================
function out = panel_c(G, W, cfg, outdir, rgrid, qgrid, p0, dgrid)
%PANEL_C  Phase diagram of the attainable speedup over traffic ratio and
%   fed fraction. The solid contour is the one percent level rho*(q), the
%   practical boundary of the lingering regime. Feeding a larger share of the
%   colony keeps lingering worthwhile at higher load, because the cue comes
%   from the fed population and its value accrues over the run.
%
%   Takes W because psiTable selects one reliance before the offset is drawn
%   and needs the gain surface as well as the ridge table.
%
%   Cue quality is not an axis. The existence criterion
%   (g1-1)*q*N*(mu-A_q) > A_q has slack of about two orders of magnitude at
%   these parameters. Fig. 5(c) tests this by sweeping q at three values of p_0.
%   The speedup at p_0 = 1.52, q = 0.1 is about 1.1 percent, close to the one
%   percent contour. If it falls below one percent, rho*(q) depends on p_0.
%
%   The speedup at each cell uses the minimum of Tq over delta, so a coarse grid
%   biases it low and moves rho*(q) left. The grid is scanned and then refined
%   with refineDelta.
%
%   P.CLIM_SPEEDUP is a floor for the color limit. The upper limit is the data
%   maximum rounded up to a multiple of five.
%
%   The sweep is cached under outdir and reused when the stamp matches. The
%   stamp carries a version tag so a convention change without a parameter
%   change forces a rebuild.

if nargin < 4 || isempty(outdir), outdir = 'output'; end
if nargin < 5 || isempty(rgrid), rgrid = logspace(log10(0.15), log10(5.0), 61); end
if nargin < 6 || isempty(qgrid), qgrid = linspace(0.10, 0.95, 35); end
if nargin < 7 || isempty(p0),    p0 = cfg.p0; end
if nargin < 8 || isempty(dgrid), dgrid = logspace(-2.5, 1.0, 25); end
VER = 'v2-reliance-outside-average';

M = lingermodel();  P = linger_style();
ps  = M.psiTable(G, p0, W);
key = M.stamp(VER, cfg, G.prm, rgrid, qgrid, dgrid, p0);
cachefile = fullfile(outdir, 'figL4_cache.mat');

Sp = [];
if isfile(cachefile)
    C = load(cachefile);
    if isfield(C,'key') && strcmp(C.key, key)
        Sp = C.Sp;  Ds = C.Ds;
        fprintf('figL4: reusing cached sweep\n');
    else
        fprintf('figL4: cache stamp mismatch, recomputing\n');
    end
end
if isempty(Sp)
    Sp = zeros(numel(qgrid), numel(rgrid));
    Ds = nan(size(Sp));
    for b = 1:numel(rgrid)
        c = cfg;  c.lam0 = rgrid(b)*cfg.mu/cfg.N;
        for a = 1:numel(qgrid)
            [dopt, Topt, Tinf] = M.optDelta(ps, c, qgrid(a), dgrid);
            Sp(a,b) = 100*max(0, 1 - Topt/Tinf);
            if Sp(a,b) > 1, Ds(a,b) = dopt; end
        end
        fprintf('  rho = %.2f done\n', rgrid(b));
    end
    save(cachefile, 'Sp', 'Ds', 'key');
end

h = P.figure(P);
ax = gca;
% pcolor with interpolated shading, since imagesc cells are visible at print
% size. Contours are computed on the swept grid.
hs = pcolor(ax, log10(rgrid), qgrid, Sp);
set(hs, 'EdgeColor', 'none');  shading(ax, 'interp');
colormap(ax, flipud(P.plasma(256)));
cl = [0, max(P.CLIM_SPEEDUP(2), 5*ceil(max(Sp(:))/5))];
clim(ax, cl);
cb = P.cbar(P, 'speedup (\%)');
cb.Label.FontSize = P.LSZ-4;
hold(ax, 'on')
[C1, h1] = contour(ax, log10(rgrid), qgrid, Sp, [1 1],  'k-','LineWidth',P.LW-1);
[Cc, hc] = contour(ax, log10(rgrid), qgrid, Sp, [10 25 40], 'k:','LineWidth',P.ALW);
% Large LabelSpacing gives one label per contour level.
clabel(Cc, hc, 'FontSize', P.FSZ-6, 'Interpreter','latex', 'LabelSpacing', 1e4);
clabel(C1, h1, 'FontSize', P.FSZ-6, 'Interpreter','latex', 'LabelSpacing', 1e4);

% Region labels. Light text on the dark upper left, dark text on the yellow right.
xl = xlim(ax);  yl = ylim(ax);
text(ax, xl(1)+0.04*diff(xl), yl(1)+0.93*diff(yl), {'lingering','speeds feeding'}, ...
     'Interpreter','latex', 'FontSize', P.FSZ, 'Color', [1 1 1], ...
     'HorizontalAlignment','left', 'VerticalAlignment','top');
text(ax, xl(1)+0.90*diff(xl), yl(1)+0.12*diff(yl), {'no net','benefit'}, ...
     'Interpreter','latex', 'FontSize', P.FSZ, 'Color', [0.1 0.1 0.1], ...
     'HorizontalAlignment','right', 'VerticalAlignment','bottom');
xlabel('traffic magnitude $\log_{10}\rho$','Interpreter','latex','FontSize',P.LSZ);
ylabel('fed fraction $q$','Interpreter','latex','FontSize',P.LSZ);
P.axes(P);
set(ax, 'Position', [0.17 0.15 0.62 0.73]);  cb.Position = [0.81 0.15 0.035 0.73];

rstar = nan(size(qgrid));
for a = 1:numel(qgrid)
    k = find(Sp(a,:) > 1, 1, 'last');
    if isempty(k)
        rstar(a) = -Inf;                       % never pays anywhere on the grid
    elseif k == numel(rgrid)
        rstar(a) = Inf;                        % boundary lies past the grid top
    else                                       % interpolate the crossing in
        u = log10(rgrid(k:k+1));               % log rho rather than snapping
        rstar(a) = 10^interp1(Sp(a,k:k+1), u, 1);
    end
end
out = struct('rgrid',rgrid,'qgrid',qgrid,'speedup',Sp,'dstar',Ds, ...
             'rstar',rstar,'p0',p0,'g1',ps.g1,'clim',cl);
if any(~isfinite(rstar))
    fprintf(['WARNING: rho*(q) not bracketed at %d of %d values of q; ' ...
             'widen rgrid\n'], sum(~isfinite(rstar)), numel(qgrid));
end
fprintf('rho*(q) runs %.2f at q=%.2f to %.2f at q=%.2f;  max speedup %.1f%% (clim %.0f)\n', ...
        rstar(1), qgrid(1), rstar(end), qgrid(end), max(Sp(:)), cl(2));
P.export(P, 'figL4', outdir, h);
end
