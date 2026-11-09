function out = figure6(varargin)
%FIGURE6 Reproduce the three panels of Figure 6 (allocation sets the optimal
%   abandonment rate) and assemble the lettered figure output/figs/fig6.pdf.
%
% Usage
%   figure6
%   figure6('fast',true)                  % low-cost wiring/smoke test
%   figure6('outdir','output/figure6')
%
% Panel map
%   (a) figT4   gain elasticity and fed populations at the optimum
%   (b) figT5   reduction/asymptotics versus individual-based simulation
%   (c) figT7   stochastic donor asymmetry at identical donors
%
% The exact two-donor arrival table and the reduced sweep come from
% twodonor_setup, the same table as Fig. 5, with the single-donor tables at
% the same parameters. The large-delta expansion is twoasym and the full
% individual-based model is abm2. A normal call uses the production grids and
% replicate counts. 'fast' uses coarse tables and small ensembles as an
% end-to-end smoke test.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'src')));

ip = inputParser;
addParameter(ip,'fast',false);
addParameter(ip,'outdir',fullfile('output','figure6'));
parse(ip,varargin{:});
fast = logical(ip.Results.fast);
OUTDIR = ip.Results.outdir;
if ~isfolder(OUTDIR), mkdir(OUTDIR); end
close all

S = twodonor_setup(OUTDIR, fast, 'figure6', true);
G2 = S.G2;  res = S.res;  cfg0 = S.cfg0;  G = S.G;  W = S.W;  ps = S.ps;
p0 = S.prm.p0;

if fast
    T4 = panel_a(res,G,W,p0,OUTDIR,12);
    T5 = panel_b(G2,cfg0,ps.g1,OUTDIR,4,[0.05 0.2 1 5],12,0.5,logspace(-2,1,20));
    T7 = panel_c(G2,cfg0,OUTDIR,0.2,16,0.9,1.0,0,10);
else
    T4 = panel_a(res,G,W,p0,OUTDIR);
    T5 = panel_b(G2,cfg0,ps.g1,OUTDIR);
    T7 = panel_c(G2,cfg0,OUTDIR,[],[],0.9);
end

out = struct('G2',G2,'res',res,'G',G,'W',W,'T4',T4,'T5',T5,'T7',T7);
save(fullfile(OUTDIR,'figure6_results.mat'),'-struct','out');

% Lettered figure, one row of three panel canvases at the panel size.
figdir = fullfile('output','figs');
if ~isfolder(figdir), mkdir(figdir); end
T = [22 19];
tiles = [0 0 T; T(1) 0 T; 2*T(1) 0 T];
panels = fullfile(OUTDIR, strcat({'figT4','figT5','figT7'}, '.fig'));
assemble_figure(fullfile(figdir,'fig6'), panels, tiles, [3*T(1) T(2)]);

fprintf('Figure 6 complete. Panel order: figT4, figT5, figT7.\n');
end

% =========================================================================
% Panel (a), figT4. Gain elasticity at the optimum.
% =========================================================================
function out = panel_a(res, G, W, p0, outdir, Fmax)
%PANEL_A  Gain elasticity e(F) = F Gtilde'(F)/Gtilde(F) with the two fed
%   populations overlaid. The contrast condition e(F_1) = e(F_2) places them on
%   opposite flanks of the elasticity peak, and the pairs separate about the
%   peak as the asymmetry grows.
%
%   Takes W because psiTable selects one reliance before the offset is drawn
%   and needs the gain surface as well as the ridge table.
%
%   The elasticity is the single-donor one. The product form reproduces
%   Lambda_1/Lambda_2 to within six percent though it misstates the total by
%   twenty. The contrast condition depends on the ratio alone, so this
%   substitution holds for this panel but not for the stability threshold.
%
%   The markers are the fed populations at the feeding threshold, the
%   representative parameters where decompose stops.
%
%   Circles are the slower donor and squares the faster one. Both are labeled
%   at the widest pair. If the pairs move, check the labels clear the curve
%   and the colorbar.
%
%   At r = 1 the two markers must coincide.

if nargin < 5 || isempty(outdir), outdir = 'output'; end
if nargin < 6 || isempty(Fmax),   Fmax   = 12; end
M = lingermodel();  P = linger_style();
ps = M.psiTable(G, p0, W);

F  = linspace(0.05, Fmax, 400);
el = F .* M.dlogPsi(ps, F);
[~, kp] = max(el);

h = P.figure(P, true);
ax = gca;  hold(ax, 'on');
plot(ax, F, el, '-', 'Color', [0.35 0.35 0.35], 'LineWidth', P.LW);
xline(F(kp), ':', 'Color', [0.35 0.35 0.35], 'LineWidth', P.ALW);

over = res.F1 > Fmax | res.F2 > Fmax;
if any(over)
    fprintf(['figT4: %d of %d pairs fall outside Fmax = %g and are clipped ' ...
             'onto the axis edge; raise Fmax to show them\n'], ...
            sum(over), numel(over), Fmax);
end
ok = ~over;
for i = find(ok)
    e1 = interp1(F, el, min(res.F1(i), Fmax));
    e2 = interp1(F, el, max(res.F2(i), F(1)));
    plot(ax, [res.F2(i) res.F1(i)], [e2 e1], '-', 'Color', [0.75 0.75 0.75], ...
         'LineWidth', P.ALW-0.5);
end
scatter(ax, res.F2(ok), interp1(F, el, max(res.F2(ok), F(1))), 140, ...
        log10(res.r(ok)), 'filled', 'MarkerEdgeColor','k', 'LineWidth', 1.2);
scatter(ax, res.F1(ok), interp1(F, el, min(res.F1(ok), Fmax)), 140, ...
        log10(res.r(ok)), 'filled', 's', 'MarkerEdgeColor','k', 'LineWidth', 1.2);
colormap(ax, P.plasma(256));
cb = P.cbar(P, 'asymmetry $\log_{10} r$');
cb.Label.FontSize = P.FSZ;

xlabel('scenting population $\bar F$','Interpreter','latex','FontSize',P.LSZ);
ylabel('gain elasticity $e(\bar F)$','Interpreter','latex','FontSize',P.LSZ);
xlim([0 Fmax]);  ylim([0 0.05*ceil(1.12*max(el)/0.05)]);
set(ax, 'YTick', 0:0.05:1);
cb.Ticks = 0:0.2:log10(max(res.r));
P.axes(P);
set(ax, 'PositionConstraint','innerposition', 'Position', [0.16 0.17 0.60 0.73]);
set(cb, 'Position', [0.79 0.17 0.03 0.73]);

% ---- inline annotation, no legend box -----------------------------------
% The peak label sits right of the dotted rule at the top edge. The marker
% labels go on the widest pair, the last unclipped r, faster donor above and
% right of its square, slower donor below and right of its circle.
yl = ylim(ax);
text(ax, F(kp), yl(2), {' elasticity',' peak'}, 'Interpreter','latex', ...
     'FontSize', P.FSZ, 'Color','k', ...
     'HorizontalAlignment','left', 'VerticalAlignment','middle');
iw = find(ok, 1, 'last');
if ~isempty(iw)
    text(ax, res.F2(iw) + 0.05*Fmax, interp1(F, el, max(res.F2(iw), F(1))), ...
         {'slower','donor'}, 'Interpreter','latex', 'FontSize', P.FSZ, ...
         'HorizontalAlignment','left', 'VerticalAlignment','top');
    text(ax, res.F1(iw) - 0.02*Fmax, interp1(F, el, min(res.F1(iw), Fmax)) ...
         + 0.06*diff(yl), {'faster','donor'}, 'Interpreter','latex', ...
         'FontSize', P.FSZ, 'HorizontalAlignment','left', ...
         'VerticalAlignment','bottom');
end

d1 = abs(res.F1(1) - res.F2(1));
if d1 > 1e-3
    warning('figT4:asymmetricAtUnity', ...
        'the r = 1 pair differs by %.3g; identical donors must give F1 = F2', d1);
end

out = struct('F',F,'elasticity',el,'Fpeak',F(kp),'emax',max(el), ...
             'F1',res.F1,'F2',res.F2,'r',res.r,'clipped',over);
fprintf('figT4: elasticity peaks at F = %.2f with e = %.3f\n', F(kp), max(el));
fprintf('       at r = 1 the pair is (%.2f, %.2f); at r = %g it is (%.2f, %.2f)\n', ...
        res.F2(1), res.F1(1), res.r(end), res.F2(end), res.F1(end));
P.export(P, 'figT4', outdir, h);
end

% =========================================================================
% Panel (b), figT5. Reduction and asymptotics against simulation.
% =========================================================================
function out = panel_b(G2, cfg0, g1, outdir, r, dcheck, nrep, q, dgrid)
%PANEL_B  Simple asymptotics against the full two-donor simulation.
%
%   Three objects on one pair of axes, as in Fig. 4(a). The reduction
%   is the two-donor ODE with the tabulated arrival rates. The expansion is
%   the closed form of twoasym, a straight line in 1/delta anchored on the
%   exact no-lingering limit, with no fitted constant. The markers are the
%   full individual-based model, in which the plume peak is the mean of the
%   scenting marks present.
%
%   The expansion derives from the reduction, so their agreement at large
%   delta tests only the algebra. The ABM tests two predictions. First, the
%   expansion is asymptotically exact, so the ABM joins it as delta grows.
%   Second, by app:gap the two-plume spectral gap 0.274 is a factor 1.64
%   below the single-donor 0.449, so the reduction over-credits at an
%   abandonment rate about 1.6 times smaller than at one donor.
%
%   The realized fraction (Tinf - ABM)/(Tinf - reduction) is printed at every
%   delta. It is near 1 inside the window and falls as the band is entered.
%
%   The expansion is drawn only where it is positive. Below delta of order
%   slope/Tinf it goes negative and is masked. It is not masked for leaving
%   the frame, since its approach to the reduction happens below the axis
%   floor. Both thresholds are printed.
%
%   The ABM must not fall below the reduction anywhere, and at the largest
%   delta all three must agree with Tinf, since no cue forms there. Both are
%   checked below.

if nargin < 4 || isempty(outdir), outdir = 'output'; end
if nargin < 5 || isempty(r),      r = 4; end   % inside twoasym's validity
if nargin < 6 || isempty(dcheck), dcheck = [0.02 0.05 0.1 0.2 0.5 1.0 2.0 5.0]; end
if nargin < 7 || isempty(nrep),   nrep = 200; end
if nargin < 8 || isempty(q),      q = 0.5; end
if nargin < 9 || isempty(dgrid),  dgrid = logspace(-2, 1.0, 40); end
VER = 'v1-two-donor-asymptotics';

M2 = twodonor();  AS = twoasym();  P = linger_style();
cfg = cfg0;  cfg.mu2 = cfg0.mutot/(1+r);  cfg.mu1 = r*cfg.mu2;
[GAP1, GAP2] = spectralgap(G2.prm); % app:gap, 0.449 one plume and 0.274 two

% The expansion needs A_q < mu_j at both donors (see twoasym). The slow donor
% sets the limit, so it is checked here.
rmax = cfg0.mutot/(G2.lam0*(1-q)*cfg0.N) - 1;
if r > rmax
    warning('figT5:outsideExpansion', ...
        ['r = %g is above %.2f, where the slow donor is at or beyond its own ' ...
         'capacity at the end of the run and the large-delta expansion is ' ...
         'undefined.  The reduction and the ABM are still valid; only the ' ...
         'dashed expansion line is not.  Rerun at r <= %.1f.'], r, rmax, floor(rmax));
end

[Tinf, Aq] = AS.Tinf(G2, cfg, q);
s = AS.slope(G2, cfg, q, g1);

% ---- reduction over a fine grid, expansion on the same grid -------------
Tred = arrayfun(@(d) M2.Tq(d, G2, cfg, 'full', q), dgrid);
Texp = AS.Tlarge(dgrid, G2, cfg, q, g1, Tinf);

% ---- the full simulation, cached ---------------------------------------
key = sprintf('%s|r=%g|q=%g|n=%d|d=%s|%s', VER, r, q, nrep, ...
              mat2str(dcheck,4), G2.key);
cachefile = fullfile(outdir, 'figT5_abm.mat');
mu_ = [];
if isfile(cachefile)
    C = load(cachefile);
    if isfield(C,'key') && strcmp(C.key, key)
        mu_ = C.mu_;  se = C.se;
        fprintf('figT5: reusing %d cached realizations per delta\n', nrep);
    else
        fprintf('figT5: ABM cache key mismatch, rerunning\n');
    end
end
if isempty(mu_)
    mu_ = zeros(size(dcheck));  se = zeros(size(dcheck));
    for j = 1:numel(dcheck)
        v = zeros(1,nrep);
        for k = 1:nrep
            v(k) = abm2(dcheck(j), G2, cfg, q, k);
        end
        nbad = sum(isnan(v));
        if nbad > 0
            fprintf('  %d of %d runs hit the time cap at delta = %.3f\n', ...
                    nbad, nrep, dcheck(j));
        end
        v = v(~isnan(v));
        mu_(j) = mean(v);  se(j) = std(v)/sqrt(numel(v));
        fprintf('  delta = %.3f done, T = %.1f +- %.1f\n', ...
                dcheck(j), mu_(j), se(j));
    end
    save(cachefile, 'mu_', 'se', 'key');
end

% ---- immediate departure, simulated -------------------------------------
% With delta = Inf no bee scents, so this run measures the unaided two-donor
% search of the individual-based model. The offset of the large-delta
% simulations from Tinf should match its offset; cached separately so the
% sweep cache above is untouched.
keyInf = sprintf('%s|inf|r=%g|q=%g|n=%d|%s', VER, r, q, nrep, G2.key);
cacheInf = fullfile(outdir, 'figT5_abm_inf.mat');
muInf = [];
if isfile(cacheInf)
    C = load(cacheInf);
    if isfield(C,'keyInf') && strcmp(C.keyInf, keyInf)
        muInf = C.muInf;  seInf = C.seInf;
    end
end
if isempty(muInf)
    v = arrayfun(@(k) abm2(Inf, G2, cfg, q, k), 1:nrep);
    v = v(~isnan(v));
    muInf = mean(v);  seInf = std(v)/sqrt(numel(v));
    save(cacheInf, 'muInf', 'seInf', 'keyInf');
end
fprintf('figT5: immediate departure, ABM %.1f +- %.1f against Tinf %.1f (%+.1f%%)\n', ...
        muInf, seInf, Tinf, 100*(muInf/Tinf - 1));

% ---- diagnostics --------------------------------------------------------
Tr = interp1(dgrid, Tred, dcheck);
Te = AS.Tlarge(dcheck, G2, cfg, q, g1, Tinf);
real_ = (Tinf - mu_)./(Tinf - Tr);
fprintf('\n  delta   ABM        reduction  expansion   realized   dev(exp)\n');
for j = 1:numel(dcheck)
    fprintf('%7.3f %8.1f+-%4.1f %9.1f %10.1f %9.2f %8.1f%%\n', ...
        dcheck(j), mu_(j), se(j), Tr(j), Te(j), real_(j), ...
        100*(mu_(j)/Te(j) - 1));
end
fprintf(['       the dev(exp) column is only meaningful where the expansion ' ...
         'is; see the validity threshold below\n']);

below = (mu_ + 1.96*se) < Tr;
if any(below)
    warning('figT5:abmBelowReduction', ...
        ['the ABM is significantly below the reduction at delta = %s; a bee ' ...
         'is being returned to the searching pool somewhere'], ...
        mat2str(dcheck(below), 3));
end
k = numel(dcheck);
if abs(mu_(k) - muInf) > 1.96*hypot(se(k), seInf) && dcheck(k) > 5*GAP2
    warning('figT5:tailNotFlat', ...
        ['at the largest delta the ABM differs from its own immediate-departure ' ...
         'run by more than 1.96 SE (%.1f against %.1f); a cue may still be ' ...
         'forming at delta = %.3g'], mu_(k), muInf, dcheck(k));
end

% ---- the panel ----------------------------------------------------------
h = P.figure(P);
ax = gca;  hold(ax, 'on');
yl = [min([mu_ - 2*se, Tred]) - 4, max([Tinf, mu_ + 2*se, muInf + 2*seInf]) + 6];

% Mask the expansion to where it is positive, not to the axis floor. It runs
% off the bottom edge where the asymptotics stop working. dvalid is the
% smallest delta on the grid past which the expansion is positive.
inframe = Texp > 0;
kv = find(~inframe, 1, 'last');
if isempty(kv), kv = 0; end
mask = false(size(dgrid));  mask(kv+1:end) = true;
if kv > 0
    dvalid = dgrid(kv+1);
    fprintf(['       expansion is negative below delta = %.3g (it is %.0f at ' ...
             'the left edge of the grid) and is masked there; it enters the ' ...
             'frame at about delta = %.3g\n'], ...
            dvalid, Texp(1), local_enters(dgrid, Texp, yl(1)));
else
    dvalid = dgrid(1);
end

plot(ax, dgrid, Tred, '-',  'Color', P.greens(4,:), 'LineWidth', P.LW);
plot(ax, dgrid(mask), Texp(mask), '--', 'Color', [0.25 0.25 0.25], ...
     'LineWidth', P.ALW);
yline(Tinf, ':', 'Color', [0 0 0], 'LineWidth', P.ALW);
errorbar(ax, dcheck, mu_, 1.96*se, 'o', 'Color', P.blues(3,:), ...
         'MarkerFaceColor','w', 'LineWidth', P.ALW, 'MarkerSize', 12, ...
         'CapSize', 10, 'LineStyle','none');
errorbar(ax, max(dgrid), muInf, 1.96*seInf, 's', 'Color', P.blues(3,:), ...
         'MarkerFaceColor', P.blues(3,:), 'LineWidth', P.ALW, 'MarkerSize', 12, ...
         'CapSize', 10, 'LineStyle','none', 'Clipping','off');

set(ax, 'XScale', 'log');
xlabel('abandonment rate $\delta$','Interpreter','latex','FontSize',P.LSZ);
ylabel('half-colony feeding time $T_{1/2}$','Interpreter','latex','FontSize',P.LSZ);
xlim([min(dgrid) max(dgrid)]);  ylim(yl);
set(ax, 'XTick', [1e-2 1e-1 1e0 1e1]);
P.axes(P);
set(ax, 'PositionConstraint','innerposition', 'Position', [0.18 0.17 0.78 0.73]);

% ---- inline labels, no legend box ---------------------------------------
dk = 0.4;                                   % gray for the reference lines
text(ax, 1.25*dcheck(1), mu_(1), 'simulation', 'Interpreter','latex', ...
     'FontSize', P.FSZ, 'Color', P.blues(3,:), ...
     'HorizontalAlignment','left', 'VerticalAlignment','middle');
text(ax, 1.7e-2, min(Tred) + 0.04*diff(yl), 'reduction', ...
     'Interpreter','latex', 'FontSize', P.FSZ, 'Color', P.greens(4,:), ...
     'HorizontalAlignment','left', 'VerticalAlignment','middle');
text(ax, 1.4*max(dgrid), Tinf - 0.035*diff(yl), 'no lingering', ...
     'Interpreter','latex', 'FontSize', P.FSZ, 'Color', dk*[1 1 1], ...
     'HorizontalAlignment','right', 'VerticalAlignment','top');
text(ax, 1.05*max(dgrid), muInf + 1.96*seInf + 0.03*diff(yl), 'immediate departure', ...
     'Interpreter','latex', 'FontSize', P.FSZ - 4, 'Color', P.blues(3,:), ...
     'HorizontalAlignment','right', 'VerticalAlignment','bottom');
de = local_enters(dgrid, Texp, yl(1));
if isnan(de), de = max(dgrid); end
text(ax, 0.36*de, yl(1) + 0.01*diff(yl), {'large $\delta$','expansion'}, ...
     'Interpreter','latex', 'FontSize', P.FSZ, 'Color', dk*[1 1 1], ...
     'HorizontalAlignment','center', 'VerticalAlignment','bottom');

out = struct('delta',dcheck,'abm',mu_,'se',se,'red',Tr,'exp',Te, ...
             'realized',real_,'dgrid',dgrid,'Tred',Tred,'Texp',Texp, ...
             'expmask',mask,'dvalid',dvalid, ...
             'Tinf',Tinf,'slope',s,'Aend',Aq,'r',r,'q',q,'nrep',nrep, ...
             'gap2',GAP2,'gap_ratio',GAP1/GAP2,'below',below, ...
             'tail_dev',mu_(k)/Tinf - 1,'abmInf',muInf,'seInf',seInf);
fprintf(['\nexpansion slope %.1f (cost %.1f, benefit %.1f), Tinf %.1f\n' ...
         'two-plume gap %.3f against single-plume %.3f, a factor %.2f\n'], ...
        s, s + (g1-1)*q*cfg.N, (g1-1)*q*cfg.N, Tinf, GAP2, GAP1, GAP1/GAP2);
P.export(P, 'figT5', outdir, h);
end

% ------------------------------------------------------------------------
function d = local_enters(dgrid, Texp, floorval)
%LOCAL_ENTERS  Smallest delta past which the expansion stays inside the frame.
k = find(Texp <= floorval, 1, 'last');
if isempty(k) || k >= numel(dgrid)
    d = NaN;
else
    d = dgrid(k+1);
end
end

% =========================================================================
% Panel (c), figT7. Stochastic asymmetry at identical donors.
% =========================================================================
function out = panel_c(G2, cfg0, outdir, delta, nrep, q, thresh, nshow, nbin)
%PANEL_C  Does the symmetric two-donor state actually hold up.
%
%   The symmetry appendix argues that recruitment cannot break the
%   symmetry between identical donors, because the antisymmetric mode grows
%   only if the elasticity difference exceeds one and the self elasticity is
%   bounded by max_F e(F), about a quarter. That argument linearizes a mean
%   field system. This panel tests whether a finite stochastic colony, in
%   which one donor can get several bees ahead by chance, recovers.
%
%   Under the theory |F_1 - F_2| is a fluctuation that scales as the square
%   root of the fed population and does not grow. Under a runaway it grows
%   and the donors split into a served one and a starved one. The panel shows
%   the ensemble, since single runs cannot tell these apart.
%
%   The comparison line is the no-feedback fluctuation scale. For independent
%   Poisson occupancies at mean F, E|F_1 - F_2| = sqrt(4F/pi) at large F.
%   Lying on the line means feedback does not affect the symmetry. Lying
%   above it means feedback amplifies the asymmetry.
%
%   The ensemble is shown as one band, mean plus or minus one standard
%   deviation over all nrep runs. No traces are drawn. Set nshow above zero
%   to draw traces for checking.
%
%   Statistics are binned, not interpolated. |F_1 - F_2| is a jump process,
%   so each run is first averaged within bins of fed fraction, and the mean
%   and standard deviation are taken across runs. nbin sets the bin count.
%   Fewer, wider bins are smoother.
%
%   The band is clamped at zero, since |F_1 - F_2| is nonnegative and right
%   skewed. The clipping is reported.
%
%   The band is the spread over runs, about ten times the standard error of
%   the mean. Wiggle in the mean line is sampling noise. Quote the precision
%   of the mean in the text and do not add error bars to the band.
%
%   mu1 must equal mu2, or the panel measures an imposed asymmetry rather
%   than a spontaneous one.

if nargin < 3 || isempty(outdir), outdir = 'output'; end
if nargin < 4 || isempty(delta),  delta = 0.2; end
if nargin < 5 || isempty(nrep),   nrep = 120; end
if nargin < 6 || isempty(q),      q = 0.9; end
if nargin < 7 || isempty(thresh), thresh = 1.0; end
if nargin < 8 || isempty(nshow),  nshow = 0; end
if nargin < 9 || isempty(nbin),   nbin = 16; end
VER = 'v1-symmetry-stochastic';

P = linger_style();
cfg = cfg0;  cfg.mu1 = cfg0.mutot/2;  cfg.mu2 = cfg.mu1;
if abs(cfg.mu1 - cfg.mu2) > 0
    error('figT7:notSymmetric', 'the donors must be identical for this test');
end

key = sprintf('%s|d=%g|q=%g|n=%d|%s', VER, delta, q, nrep, G2.key);
cachefile = fullfile(outdir, 'figT7_abm.mat');
Tr = [];
if isfile(cachefile)
    C = load(cachefile);
    if isfield(C,'key') && strcmp(C.key, key)
        Tr = C.Tr;  fprintf('figT7: reusing %d cached realizations\n', nrep);
    else
        fprintf('figT7: cache key mismatch, rerunning\n');
    end
end
if isempty(Tr)
    Tr = cell(1, nrep);
    for k = 1:nrep
        [~, tr] = abm2(delta, G2, cfg, q, 1000+k, struct('trace',true));
        Tr{k} = tr;
        if mod(k, 20) == 0, fprintf('  %d of %d runs done\n', k, nrep); end
    end
    save(cachefile, 'Tr', 'key');
end

% ---- resample every run onto a common fed-count axis --------------------
% The abscissa is progress through the run, not time, so runs of different
% speed are averaged at the same stage.
edges = linspace(0.05, q, nbin+1);
fg    = 0.5*(edges(1:end-1) + edges(2:end));   % bin centres
D  = nan(nrep, nbin);            % per-run mean |F1 - F2| in each bin
Fb = nan(nrep, nbin);            % per-run mean fed population per donor
nsamp = zeros(1, nbin);          % total samples contributing, for reporting
for k = 1:nrep
    tr = Tr{k};
    if isempty(tr), continue; end
    fed = (tr(:,4) + tr(:,6) + (cfg.N - tr(:,2) - tr(:,3) - tr(:,4) ...
           - tr(:,5) - tr(:,6)))/cfg.N;
    d = abs(tr(:,4) - tr(:,6));
    f = 0.5*(tr(:,4) + tr(:,6));
    ok = fed > 0 & fed >= edges(1) & fed <= edges(end);
    if sum(ok) < 3, continue; end
    b = discretize(fed(ok), edges);
    dk = d(ok);  fk = f(ok);
    for j = 1:nbin
        m = b == j;
        if any(m)
            D(k,j)  = mean(dk(m));
            Fb(k,j) = mean(fk(m));
            nsamp(j) = nsamp(j) + sum(m);
        end
    end
end

nk = sum(~isnan(D), 1);
fprintf(['figT7: %d bins over fed fraction, %d to %d runs per bin, %.0f ' ...
         'samples per bin on average\n'], nbin, min(nk), max(nk), mean(nsamp));
if min(nk) < 0.8*nrep
    warning('figT7:thinGrid', ...
        ['only %d of %d runs reach some points of the fed grid, so the mean ' ...
         'there is over a biased subset of the faster runs'], min(nk), nrep);
end

mD  = mean(D, 1, 'omitnan');
sdD = std(D, 0, 1, 'omitnan');            % spread over runs
seD = sdD./sqrt(max(nk,1));               % precision of the mean
lo  = max(mD - sdD, 0);                   % clamped: the quantity is >= 0
hi  = mD + sdD;
nclip = sum(mD - sdD < 0);
mF  = mean(Fb, 1, 'omitnan');
pred = sqrt(4*mF/pi);                     % independent-Poisson fluctuation

h = P.figure(P);
ax = gca;  hold(ax, 'on');

% ensemble spread, behind everything
fill(ax, [fg fliplr(fg)], [lo fliplr(hi)], P.blues(1,:), ...
     'FaceAlpha', 0.30, 'EdgeColor', 'none');

% single runs are off by default; nshow > 0 is for checking, not for the paper
if nshow > 0
    sel = unique(round(linspace(1, nrep, min(nshow, nrep))));
    for k = sel
        plot(ax, fg, D(k,:), '-', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.9);
    end
else
    sel = [];
end

plot(ax, fg, pred, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', P.ALW);
plot(ax, fg, mD, 'o-', 'Color', P.blues(3,:), ...
     'MarkerFaceColor', 'w', 'LineWidth', P.ALW, 'MarkerSize', 10);

xlabel('fed fraction $q$','Interpreter','latex','FontSize',P.LSZ);
ylabel('donor asymmetry $|F_1-F_2|$','Interpreter','latex','FontSize',P.LSZ);
xlim([min(fg) max(fg)]);
ytop = 1.12*max([hi(:); pred(:)], [], 'omitnan');
ylim(ax, [0 ytop]);
set(ax, 'YTick', 0:0.5:ytop);
P.axes(P);
set(ax, 'PositionConstraint','innerposition', 'Position', [0.18 0.17 0.78 0.73]);

% ---- inline labels, no legend box ---------------------------------------
% The fluctuation label clears both curves from above over the span it
% covers, the ensemble label from below. The band label goes in the upper
% right of the band.
yl = ylim(ax);
xs = [0.18 0.63];  js = fg >= xs(1) & fg <= xs(2);
text(ax, xs(1), max([pred(js) mD(js)]) + 0.02*diff(yl), ...
     'no feedback fluctuations', 'Interpreter','latex', ...
     'FontSize', P.FSZ, 'Color', [0.25 0.25 0.25], ...
     'HorizontalAlignment','left', 'VerticalAlignment','bottom');
xs = [0.13 0.37];  js = fg >= xs(1) & fg <= xs(2);
text(ax, xs(1), min([pred(js) mD(js)]) - 0.02*diff(yl), 'ensemble mean', ...
     'Interpreter','latex', 'FontSize', P.FSZ, 'Color', P.blues(3,:), ...
     'HorizontalAlignment','left', 'VerticalAlignment','top');
text(ax, 0.6, 0.81*yl(2), 'mean $\pm$ std', 'Interpreter','latex', ...
     'FontSize', P.FSZ, 'Color', [0.55 0.72 0.90], ...
     'HorizontalAlignment','center', 'VerticalAlignment','middle');

ratio = mD./max(pred, eps);
fprintf(['figT7: asymmetry over the no-feedback scale runs %.2f at the start ' ...
         'to %.2f at the end (max %.2f)\n'], ratio(1), ratio(end), max(ratio));
fprintf(['       %d realizations, %d drawn; spread (s.d.) %.2f at its widest ' ...
         'against a mean se of %.3f, a factor of %.0f\n'], ...
        nrep, numel(sel), max(sdD), max(seD), max(sdD)/max(seD));
if nclip > 0
    fprintf(['       the lower band edge is clamped at zero at %d of %d grid ' ...
             'points, where the mean is within one s.d. of the floor\n'], ...
            nclip, numel(fg));
end
if max(ratio) > thresh + 0.5
    warning('figT7:amplified', ...
        ['the asymmetry exceeds the no-feedback fluctuation scale by %.2f, ' ...
         'so the feedback is amplifying it; the linearization in app:symmetry ' ...
         'says it should not'], max(ratio));
end
grosslev = 0.5*cfg0.N*0.1;                % half the colony split ten to one
frac = mean(D(:,end) > grosslev, 'omitnan');
fprintf('       %.1f%% of runs end above |F1-F2| = %g\n', 100*frac, grosslev);

out = struct('fedgrid',fg,'edges',edges,'nbin',nbin,'nsamp',nsamp, ...
             'asym',D,'meanasym',mD,'sd',sdD,'se',seD, ...
             'lo',lo,'hi',hi,'nclipped',nclip,'meanF',mF,'nvalid',nk, ...
             'shown',sel,'pred',pred,'ratio',ratio,'delta',delta,'q',q, ...
             'nrep',nrep,'gross_level',grosslev,'gross_fraction',frac);
% Inset, two identical donors, each with a scenting pair, in the clear strip
% under the band.
inset_pair(h, P);
P.export(P, 'figT7', outdir, h);
end

% ------------------------------------------------------------------------
function inset_pair(h, P)
%INSET_PAIR  Identical-donor sketch, between the band and the x axis.
K = inset_kit();
ax = K.axes(h, P, [9.6 3.4 6.8 1.65]);
y0 = 0.55;  xd = [1.7 5.1];
for k = 1:2, K.tent(ax, xd(k) - 0.1, 0.85, y0, 0.95, 0.45, 0.94); end
K.track(ax, 0.1, 6.7, y0);
for k = 1:2
    K.donor(ax, xd(k), y0, 0.2);
    K.bee(ax, xd(k) - 0.62, y0 + 0.22, -0.3, 0.68, 'scent');
    K.bee(ax, xd(k) + 0.62, y0 - 0.16, pi - 0.3, 0.68, 'scent');
end
end
