function out = figure4(varargin)
%FIGURE5 Reproduce the three panels of Figure 4 (feeding is fastest at an
%   intermediate abandonment rate) and assemble the lettered figure
%   output/figs/fig4.pdf.
%
% Usage
%   figure4
%   figure4('fast',true)                  % low-cost wiring/smoke test
%   figure4('outdir','output/figure4')
%
% Panel map
%   (a) figL6   reduction versus individual-based model
%   (b) figL7   fed population and total donor occupancy over time
%   (c) figL5   optimal abandonment versus fed fraction
%
% The single-donor gain tables come from linger_setup, the same tables as
% Fig. 3, and are cached with a parameter stamp. The stochastic panel data are
% cached in outdir by each panel routine. A normal call uses the paper
% settings. 'fast' uses coarser grids and fewer realizations to check that the
% code runs end to end, and never changes a model parameter.

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
G = S.G;  psList = S.psList;  cfg = S.cfg;

if fast
    qgrid5 = linspace(0.10,0.95,10);
    dgrid5 = logspace(-2.5,1.0,20);
    dcheck = [0.03 0.10 0.20 0.50 1.0];
    nrep = 12;
else
    qgrid5 = linspace(0.10,0.95,22);
    dgrid5 = logspace(-2.5,1.0,44);
    dcheck = [0.01 0.03 0.05 0.1 0.15 0.2 0.3 0.5 1.0 3.0];
    nrep = 200;
end

L5 = panel_c(psList,cfg,OUTDIR,qgrid5,dgrid5);
L6 = panel_a(psList,G,cfg,OUTDIR,dcheck,nrep,0.5,1);
L7 = panel_b(psList,G,cfg,OUTDIR,[0.20 0.248 1.00],0.5,1,false,0);

out = struct('G',G,'psList',{psList},'L5',L5,'L6',L6,'L7',L7);
save(fullfile(OUTDIR,'figure4_results.mat'),'-struct','out');

% Lettered figure, one row of three panel canvases at the panel size.
figdir = fullfile('output','figs');
if ~isfolder(figdir), mkdir(figdir); end
T = [22 19];
tiles = [0 0 T; T(1) 0 T; 2*T(1) 0 T];
panels = fullfile(OUTDIR, strcat({'figL6','figL7','figL5'}, '.fig'));
assemble_figure(fullfile(figdir,'fig4'), panels, tiles, [3*T(1) T(2)]);

fprintf('Figure 4 complete. Panel order: figL6, figL7, figL5.\n');
end

% =========================================================================
% Panel (a), figL6. Reduction against the individual-based model.
% =========================================================================
function out = panel_a(psList, G, cfg, outdir, dcheck, nrep, q, abmIdx)
%PANEL_A  Reduced system against the full individual-based model.
%   The reduction credits the recruitment gain Psi~(Fbar) in full, but a
%   searcher can exploit a plume only if it persists long enough for its
%   phase profile on the ring to build. That profile relaxes on the lap time
%   L = 2*pi/v, so the closure holds in the window
%
%       2*pi/v  <<  1/delta  <<  1/lam0,
%
%   and over-credits when the scenting population turns over faster than a
%   circuit. The printed realized fraction (Tinf - ABM)/(Tinf - reduction) is
%   near 1 inside the window and falls to 0 as delta*L passes 1, where the
%   ABM reaches the no-lingering line.
%
%   Only the reference p_0 is drawn. The closure turns on delta*L, not p_0,
%   and Fig. 3(a) and panel (c) show the p_0 family. The gap between curve and markers
%   is shaded, since it is the realized fraction quoted in the text. Set
%   SHOWALL true to draw all three reduction curves.
%
%   T_{1/2} is the time to feed half the colony, not a median. Markers are
%   the mean over nrep realizations and bars are 1.96 standard errors.
%
%   Crowd occupancy is not the relevant scale. The crowd is largest where the
%   agreement is best. Panel (b) shows the occupancy.
%
%   The ABM must not fall below the reduction anywhere. If it does, a bee is
%   being returned to the searching pool. This is checked below.
%
%   The ABM runs are cached under outdir and reused when the stamp matches.
%   The stamp carries a version tag so a convention change without a
%   parameter change forces a rerun.

if nargin < 4, outdir = 'output'; end
if nargin < 5, dcheck = [0.01 0.03 0.05 0.1 0.15 0.2 0.3 0.5 1.0 3.0]; end
if nargin < 6, nrep = 200; end
if nargin < 7, q = 0.5; end
if nargin < 8, abmIdx = 1; end          % which p_0 gets simulated markers
SHOWALL = false;
VER = 'v2-reliance-outside-average';

M = lingermodel();  P = linger_style();
Lap = 2*pi/G.prm.v;
dg  = logspace(-2.5, 1.0, 44);
np  = numel(psList);

Tred = zeros(np, numel(dg));
for i = 1:np
    Tred(i,:) = arrayfun(@(d) M.Tq(d, psList{i}, cfg, q), dg);
end
ps   = psList{abmIdx};                   % the curve the markers belong to
Tinf = M.Tq(1e4, ps, cfg, q);            % independent of p_0, so any will do

key = M.stamp(VER, cfg, G.prm, dcheck, nrep, q, ps.p0);
cachefile = fullfile(outdir, 'figL6_abm.mat');
mu_ = [];
if isfile(cachefile)
    C = load(cachefile);
    if isfield(C,'key') && strcmp(C.key, key)
        mu_ = C.mu_;  se = C.se;
        fprintf('figL6: reusing %d cached ABM realizations per delta\n', nrep);
    else
        fprintf('figL6: ABM cache stamp mismatch, rerunning\n');
    end
end
if isempty(mu_)
    mu_ = zeros(size(dcheck));  se = zeros(size(dcheck));
    for j = 1:numel(dcheck)
        r = arrayfun(@(k) M.abm(dcheck(j), G, ps, cfg, q, k), 1:nrep);
        if any(isnan(r))
            fprintf('  %d of %d runs hit the time cap at delta = %.3f\n', ...
                    sum(isnan(r)), nrep, dcheck(j));
        end
        r = r(~isnan(r));
        mu_(j) = mean(r);  se(j) = std(r)/sqrt(numel(r));
    end
    save(cachefile, 'mu_', 'se', 'key');
end

Tr = zeros(size(dcheck));  real_ = zeros(size(dcheck));
for j = 1:numel(dcheck)
    Tr(j)    = M.Tq(dcheck(j), ps, cfg, q);
    real_(j) = (Tinf - mu_(j))/(Tinf - Tr(j));
    fprintf(['delta=%6.3f  ABM %7.1f +- %4.1f  reduction %7.1f  ' ...
             'gap %+5.1f%%  delta*L = %5.2f  realized %.2f\n'], ...
        dcheck(j), mu_(j), se(j), Tr(j), 100*(mu_(j)/Tr(j)-1), ...
        dcheck(j)*Lap, real_(j));
end

below = (mu_ + 1.96*se) < Tr;
if any(below)
    warning('figL6:abmBelowReduction', ...
        ['ABM is significantly BELOW the reduction at delta = %s; a bee is ' ...
         'being returned to the searching pool somewhere'], ...
        mat2str(dcheck(below), 3));
end

% optimum located from the simulation rather than from the reduction
[~, kA] = min(mu_);
near = abs(mu_ - mu_(kA)) < 1.96*sqrt(se.^2 + se(kA)^2);
[dR, TR] = M.refineDelta(@(d) M.Tq(d, ps, cfg, q), dg, Tred(abmIdx,:));

h = P.figure(P);
ax = gca;  hold(ax, 'on');
yl = [min([mu_ - 2*se, Tred(abmIdx,:)]) - 4, Tinf + 6];
hb = patch([1/Lap max(dg) max(dg) 1/Lap], [yl(1) yl(1) yl(2) yl(2)], ...
           [0.90 0.90 0.90], 'EdgeColor','none');       % delta*L > 1

% Shaded gap between reduction and simulation, the credited benefit that is
% not realized. Drawn on the marker abscissae, where both are known.
Tri = interp1(dg, Tred(abmIdx,:), dcheck);
patch(ax, [dcheck fliplr(dcheck)], [Tri fliplr(mu_)], P.p0cols(1,:), ...
      'FaceAlpha', 0.16, 'EdgeColor','none');

if SHOWALL
    ord = [setdiff(1:np, abmIdx), abmIdx];   % reference drawn last, on top
    for j = ord
        plot(ax, dg, Tred(j,:), P.p0sty{min(j,3)}, 'Color', P.p0cols(min(j,3),:), ...
             'LineWidth', P.LW - 1.0*(j ~= abmIdx));
    end
else
    plot(ax, dg, Tred(abmIdx,:), '-', 'Color', P.p0cols(1,:), 'LineWidth', P.LW);
end

yline(Tinf, ':', 'Color', [0 0 0], 'LineWidth', P.ALW);
errorbar(ax, dcheck, mu_, 1.96*se, 'o', 'Color', P.blues(3,:), ...
         'MarkerFaceColor','w', 'LineWidth', P.ALW, 'MarkerSize', 12, ...
         'CapSize', 10, 'LineStyle','none');
uistack(hb, 'bottom');
set(ax,'XScale','log');
xlabel('abandonment rate $\delta$','Interpreter','latex','FontSize',P.LSZ);
ylabel('half-colony feeding time $T_{1/2}$','Interpreter','latex','FontSize',P.LSZ);
xlim([min(dg) max(dg)]);  ylim(yl);
% Label every decade so the optimum is easy to place.
set(ax, 'XTick', [1e-2 1e-1 1e0 1e1]);
P.axes(P);
set(ax, 'Position', [0.17 0.15 0.78 0.73]);

x0 = 1.2*min(dg);
P.label(P, ax, x0, Tinf + 0.8, 'no lingering', 'FontSize', P.FSZ+4, ...
        'VerticalAlignment','bottom');
P.label(P, ax, x0, Tinf - 0.8, sprintf('$T_{1/2}(\\infty)=%.0f$', Tinf), ...
        'FontSize', P.FSZ+4, 'VerticalAlignment','top');
P.label(P, ax, 1.15/Lap, Tinf - 9, '$\delta L>1$', 'FontSize', P.FSZ+4, ...
        'VerticalAlignment','top');
% Curve labels along the descending branch, rotated to its screen slope.
xs = [0.06 0.017];  ys = interp1(dg, Tred(abmIdx,:), xs);
P.label(P, ax, xs(1), ys(1) + 12, 'simulation', 'FontSize', P.FSZ+4, ...
        'Color', P.blues(3,:), 'Rotation', -52, 'HorizontalAlignment','center');
P.label(P, ax, xs(2), ys(2) - 7, 'theory', 'FontSize', P.FSZ+4, ...
        'Color', P.p0cols(1,:), 'Rotation', -52, 'HorizontalAlignment','center');

fprintf(['\nreduction optimum  delta = %.3f, T = %.1f, speedup %.1f%%\n' ...
         'simulated optimum  delta = %.3f, T = %.1f +- %.1f, speedup %.1f%%\n' ...
         'indistinguishable from it: delta in [%.3f, %.3f]\n'], ...
    dR, TR, 100*(1-TR/Tinf), ...
    dcheck(kA), mu_(kA), se(kA), 100*(1-mu_(kA)/Tinf), ...
    min(dcheck(near)), max(dcheck(near)));

inwin = dcheck*Lap < 1;
fprintf('realized fraction inside the window: %.2f to %.2f\n', ...
        min(real_(inwin)), max(real_(inwin)));

out = struct('delta',dcheck,'abm',mu_,'se',se,'red',Tr,'realized',real_, ...
             'lapd',dcheck*Lap,'dg',dg,'Tred',Tred,'Tinf',Tinf, ...
             'dstar_abm',dcheck(kA),'dstar_red',dR,'near',near, ...
             'nrep',nrep,'q',q,'below',below, ...
             'p0',cellfun(@(z) z.p0, psList),'abmIdx',abmIdx);
P.export(P, 'figL6', outdir, h);
end

% =========================================================================
% Panel (b), figL7. Crowd composition at the donor over time.
% =========================================================================
function out = panel_b(psList, G, cfg, outdir, deltas, q, idx, showload, nover)
%PANEL_B  Crowd composition at the donor.
%
%   The no-imposed-capacity claim is scoped to the scenting population F, not
%   the crowd U+F. F stays near three bees while U+F reaches about fifteen at
%   the optimum, so a capacity on scenters would be inactive and a capacity
%   on the crowd would bind.
%
%   The focal abandonment rate deltas(1) is drawn as a stacked composition, U
%   filled on top of F. The other two rates are thinner curves. Hue marks the
%   family and line style marks delta. Peaks of the focal pair are marked,
%   since the text quotes maxima.
%
%   G is the second argument, matching panel (a), because the overlay calls the
%   ABM and the ABM needs G.prm.
%
%   Optional overlay of nover realizations of F(t) at the focal rate, drawn as
%   thin stairs over the fill, off by default. The realizations are small
%   integers with large spread, the smallness the substitution argument of
%   panel (a) turns on. Seeds are 1:nover, the first members of the ensemble
%   panel (a) averages at this delta.
%
%   Stairs are used because the count is piecewise constant between events.
%   Each realization stops at its own T, since the ABM halts at nFed >= q*N.
%   Requires the two-output abm in lingermodel.m.
%
%   Measured values
%       delta      max(U+F)    max(F)   max(U)   max(load)
%       0.200        15.0       3.75     11.3      1.01     ABM optimum
%       0.248        12.4       3.04      9.4      0.97     reduced optimum
%       1.000         small     small    small     small
%
%   The load is near unity at the optimum, so the quasi-stationary relations
%   sit at the edge of their range at delta*. The load is returned but not
%   drawn by default. Pass showload = true to draw it on a second y axis.
%
%   q is the stopping rule, so the panel ends when a fraction q is fed. The
%   horizon scales as log(1/(1-q))/lam0, 86 at q = 0.5. The rise and peak
%   occupy the first fifty time units. q = 0.5 matches the other panels.

if nargin < 4 || isempty(outdir), outdir = 'output'; end
if nargin < 5 || isempty(deltas), deltas = [0.20 0.248 1.00]; end  % ABM opt, reduced, fast
if nargin < 6 || isempty(q),      q = 0.5; end
if nargin < 7 || isempty(idx),    idx = 1; end                     % reference p_0
if nargin < 8 || isempty(showload), showload = false; end
if nargin < 9 || isempty(nover),  nover = 0; end

M = lingermodel();  P = linger_style();
ps = psList{idx};
sty = {'-','--',':'};
kf = 1;                                          % focal delta, the ABM optimum

nd = numel(deltas);
out = struct('delta',deltas,'q',q,'p0',ps.p0,'focal',deltas(kf),'kf',kf, ...
             'maxCrowd',nan(1,nd),'maxF',nan(1,nd),'maxU',nan(1,nd), ...
             'maxLoad',nan(1,nd),'Tq',nan(1,nd), ...
             'nover',nover,'overMaxF',[],'overT',[]);

h = P.figure(P);
ax = gca;  hold(ax,'on');
tr = cell(1,nd);
for k = 1:nd
    tr{k} = M.traj(deltas(k), ps, cfg, q);
    out.maxCrowd(k) = max(tr{k}.crowd);
    out.maxF(k)     = max(tr{k}.F);
    out.maxU(k)     = max(tr{k}.U);
    out.maxLoad(k)  = max(tr{k}.load);
    out.Tq(k)       = tr{k}.T;
end

if showload, yyaxis(ax,'left'); set(ax,'YColor',[0.15 0.15 0.15]); end

% Focal rate as a composition. Scenters at the bottom in the F green, unfed
% stacked above in the crowd blue, matching the line colors.
tf = tr{kf}.t(:)';  Ff = tr{kf}.F(:)';  Cf = tr{kf}.crowd(:)';
fill(ax, [tf fliplr(tf)], [zeros(size(tf)) fliplr(Ff)], P.p0cols(1,:), ...
     'FaceAlpha', 0.85, 'EdgeColor','none');
fill(ax, [tf fliplr(tf)], [Ff fliplr(Cf)], P.blues(3,:), ...
     'FaceAlpha', 0.30, 'EdgeColor','none');

% Realizations of F at the focal rate, over the fill and under the curves.
if nover > 0
    oM = nan(1,nover);  oT = nan(1,nover);
    for k = 1:nover
        [~, o] = M.abm(deltas(kf), G, ps, cfg, q, k);
        stairs(ax, o.t, o.F, '-', 'Color', [P.greens(2,:) 0.55], ...
               'LineWidth', max(0.8, P.LW-2.3));
        oM(k) = max(o.F);  oT(k) = o.T;
    end
    out.overMaxF = oM;  out.overT = oT;
end

% Crowd above, scenters below, one hue each, line style keyed to delta. The
% focal rate keeps full weight and the others are thinned.
for k = 1:nd
    lw = P.LW;  if k ~= kf, lw = P.LW - 0.9; end
    plot(ax, tr{k}.t, tr{k}.crowd, sty{k}, 'Color', P.blues(3,:), 'LineWidth', lw);
end
for k = 1:nd
    lw = P.LW;  if k ~= kf, lw = P.LW - 0.9; end
    plot(ax, tr{k}.t, tr{k}.F, sty{k}, 'Color', P.p0cols(1,:), 'LineWidth', lw);
end

% The two maxima quoted in the text are located and printed, not drawn.
[~, iC] = max(Cf);  [~, iF] = max(Ff);

xlabel(ax,'time $t$','Interpreter','latex','FontSize',P.LSZ);
ylabel(ax,'bees at the donor','Interpreter','latex','FontSize',P.LSZ);
% The x axis ends at the focal T_q, not the largest, to avoid a long flat tail.
xlim(ax,[0 out.Tq(kf)]);
ytop = 1.12*max(out.maxCrowd);
if nover > 0, ytop = max(ytop, 1.08*max(out.overMaxF)); end
ylim(ax,[0 ytop]);

if showload
    yyaxis(ax,'right');
    set(ax,'YColor',P.blues(1,:));
    yline(ax, 1, ':', 'Color', [0.35 0.35 0.35], 'LineWidth', P.ALW);
    for k = 1:nd
        lw = P.LW - 0.9;  if k == kf, lw = P.LW - 0.4; end
        plot(ax, tr{k}.t, tr{k}.load, sty{k}, 'Color', P.blues(1,:), 'LineWidth', lw);
    end
    ylim(ax,[0 max(1.45, 1.2*max(out.maxLoad))]);
    ylabel(ax,'load $A/\mu$','Interpreter','latex','FontSize',P.LSZ);
    yyaxis(ax,'left');
end

P.axes(P);
set(ax, 'Position', [0.17 0.15 0.78 0.73]);

% Line style key for delta, drawn in the axes, and family labels at the peaks.
yk = ytop*[0.95 0.885 0.82];  xk = out.Tq(kf)*[0.02 0.11];
for k = 1:nd
    plot(ax, xk, yk(k)*[1 1], sty{k}, 'Color','k', 'LineWidth', P.LW-0.5);
    dl = sprintf('%.3f', deltas(k));
    if deltas(k) == round(deltas(k)), dl = sprintf('%d', deltas(k)); end
    P.label(P, ax, xk(2) + 0.02*out.Tq(kf), yk(k), sprintf('$\\delta=%s$', dl), ...
            'FontSize', P.FSZ+4, 'VerticalAlignment','middle');
end
P.label(P, ax, tf(iC), Cf(iC) + 0.04*ytop, '$U+F$', 'FontSize', P.FSZ+4, ...
        'Color', P.blues(3,:), 'HorizontalAlignment','center', 'VerticalAlignment','bottom');

% Region labels inside the focal composition, placed between the focal and
% dashed curves so no line crosses them.
xr = 0.60*out.Tq(kf);  Fr = interp1(tf, Ff, xr);  Cr = interp1(tf, Cf, xr);
P.label(P, ax, xr, Fr + 0.35*(Cr - Fr), 'waiting, $U$', 'FontSize', P.FSZ+2, ...
        'Color', P.blues(3,:), 'HorizontalAlignment','center', 'VerticalAlignment','middle');
P.label(P, ax, xr, 0.45*Fr, 'lingering, $F$', 'FontSize', P.FSZ+2, ...
        'Color', 'w', 'HorizontalAlignment','center', 'VerticalAlignment','middle');

fprintf('\nfigL7: crowd composition at p0 = %.2f, q = %.2f\n', ps.p0, q);
for k = 1:nd
    mark = '  ';  if k == kf, mark = ' *'; end
    fprintf(['%s delta=%6.3f  T_q=%7.2f  max(U+F)=%6.1f  max(F)=%5.2f  ' ...
             'max(U)=%6.1f  max(load)=%5.2f\n'], mark, ...
        deltas(k), out.Tq(k), out.maxCrowd(k), out.maxF(k), ...
        out.maxU(k), out.maxLoad(k));
end
fprintf('  (* focal, drawn as the stacked composition)\n');
if out.maxLoad(kf) > 1
    fprintf(['  note: at the focal delta the arrival rate reaches %.2f mu, ' ...
             'so eq:qss_crowd sits at the edge of its range; App. J should ' ...
             'say grazes, not exceeds\n'], out.maxLoad(kf));
end
if any(out.maxLoad > 1.10)
    fprintf(['  note: load exceeds 1.10 at %d of %d values of delta, ' ...
             'so eq:qss_crowd is outside its own range there\n'], ...
            sum(out.maxLoad > 1.10), nd);
end
if nover > 0
    fprintf(['  overlay: %d realizations, peak F per run %d to %d against a ' ...
             'reduction peak of %.2f;  T_q %.1f to %.1f\n'], ...
        nover, min(out.overMaxF), max(out.overMaxF), out.maxF(kf), ...
        min(out.overT), max(out.overT));
    if max(out.overMaxF) < 1.5*out.maxF(kf)
        fprintf(['  the realizations track the mean closely, so the overlay ' ...
                 'is not earning its ink; drop it rather than keeping it ' ...
                 'for decoration\n']);
    end
end
if out.maxF(kf) > 5
    warning('figL7:scentersLarge', ...
        ['max(F) = %.1f at the focal delta; the no-capacity claim is ' ...
         'scoped to F, so check Sec. 4.2 still reads correctly'], out.maxF(kf));
end
% Inset, the crowd at the donor, with scenting and waiting bees mixed.
inset_crowd(h, P);
P.export(P, 'figL7', outdir, h);
end

% ------------------------------------------------------------------------
function inset_crowd(h, P)
%INSET_CROWD  Donor crowd icon, in the upper right beside the U+F label.
K = inset_kit();
ax = K.axes(h, P, [17.0 15.0 4.2 2.4]);
xd = 2.1;  yd = 1.2;
K.donor(ax, xd, yd, 0.24);
% Scenting and waiting bees are mixed around the donor, since random-order
% service does not distinguish them by position.
a = [30 90 150 210 270 330]*pi/180;
kind = {'wait','scent','wait','wait','scent','wait'};
for k = 1:numel(a)
    K.bee(ax, xd + 0.72*cos(a(k)), yd + 0.62*sin(a(k)), a(k) + pi, 0.7, kind{k});
end
end

% =========================================================================
% Panel (c), figL5. Optimal abandonment against fed fraction.
% =========================================================================
function out = panel_c(psList, cfg, outdir, qgrid, dgrid)
%PANEL_C  How the optimum moves with cue quality.
%   delta*(q) at three baseline offsets p_0, with the same values, line styles
%   and colors (P.p0sty, P.p0cols) as Fig. 3(a), so no legend is needed.
%
%   The panel tests whether cue quality sets how much lingering is worth
%   rather than whether it is worth anything. If so, no curve terminates
%   partway across q. A terminating curve means Fig. 3(c) needs p_0 as an axis.
%   The speedup at p_0 = 1.52, q = 0.1 is about one percent, the level Fig. 3(c)
%   contours. Report the margin along with the verdict.
%
%   The ordering of delta* with p_0 reverses in q. A noisier cue is abandoned
%   sooner when a tenth of the colony must be fed and held longer when nearly
%   all must be. The large-delta expansion predicts this, since the benefit
%   coefficient carries a factor q*N and the congestion cost does not. The
%   payoff falls with p_0 at every q. State the reversal in the caption.
%
%   Speedup is not drawn here. Fig. 3(c) contours it, and out.speedup holds it.
%   P.CLIM_SPEEDUP applies to Fig. 3(c) only.
%
%   Every curve is at the run's traffic ratio rho = lam0*N/mu, printed and
%   returned in out.rho. State it in the caption.

if nargin < 3, outdir = 'output'; end
if nargin < 4, qgrid = linspace(0.1, 0.95, 22); end
if nargin < 5, dgrid = logspace(-2.5, 1.0, 44); end

M = lingermodel();  P = linger_style();
rho = cfg.lam0*cfg.N/cfg.mu;
np  = numel(psList);

dstar = nan(np, numel(qgrid));
spd   = zeros(np, numel(qgrid));
chk   = zeros(np, numel(qgrid));
for i = 1:np
    for k = 1:numel(qgrid)
        [dopt, Topt, Tinf] = M.optDelta(psList{i}, cfg, qgrid(k), dgrid);
        dstar(i,k) = dopt;
        spd(i,k)   = 100*max(0, 1 - Topt/Tinf);
        chk(i,k)   = Tinf/(log(1/(1-qgrid(k)))/cfg.lam0);
    end
    fprintf('  p0 = %.2f done\n', psList{i}.p0);
end

h = P.figure(P);
ax = gca;  hold(ax, 'on')
for i = 1:np
    ci = P.p0cols(min(i,3),:);
    plot(ax, qgrid, dstar(i,:), P.p0sty{min(i,3)}, 'Color', ci, 'LineWidth', P.LW);
    scatter(ax, qgrid(2:4:end), dstar(i,2:4:end), 60, ci, 'filled', ...
            'MarkerEdgeColor','k', 'LineWidth', 0.9);
end
set(ax, 'YScale', 'log');
% Curve labels at the right edge, where the fan is widest. No legend, since
% the labels in Fig. 3(a) serve as the key.
ylab = dstar(:,end);   % spread the end labels so they do not touch
if np == 3
    ylab(1) = ylab(1)*1.4;
    ylab(3) = ylab(3)/1.3;
end
for i = 1:np
    text(ax, max(qgrid)+0.015, ylab(i), sprintf('$p_0=%.2f$', psList{i}.p0), ...
         'Interpreter','latex', 'FontSize', P.FSZ+2, 'Color', P.p0cols(min(i,3),:), ...
         'HorizontalAlignment','left', 'VerticalAlignment','middle');
end
xlabel('fed fraction $q$','Interpreter','latex','FontSize',P.LSZ);
ylabel('optimal abandonment rate $\delta^\ast$','Interpreter','latex','FontSize',P.LSZ);
xlim([min(qgrid)-0.02, max(qgrid)+0.34]);   % right margin holds the p_0 labels
ylim([0.45*min(dstar(:)), 1.6*max(dstar(:))]);
P.axes(P);
set(ax, 'Position', [0.17 0.15 0.78 0.73], 'XTick', 0.2:0.2:1);

out = struct('q',qgrid,'dstar',dstar,'speedup',spd,'Tinf_check',chk, ...
             'rho',rho,'p0',cellfun(@(s) s.p0, psList), ...
             'g1',cellfun(@(s) s.g1, psList));
fprintf('rho for this panel: %.3f\n', rho);
fprintf('Tq(inf) check, max deviation: %.2e\n', max(abs(chk(:)-1)));
for i = 1:np
    fprintf('p0 = %.2f (g1 = %.3f): delta* %.3f -> %.3f, speedup %.1f%% -> %.1f%%\n', ...
        psList{i}.p0, psList{i}.g1, dstar(i,1), dstar(i,end), spd(i,1), spd(i,end));
end
if any(~isfinite(dstar(:)))
    fprintf(['WARNING: delta* undefined at %d of %d points; a curve dies ' ...
             'partway across q, so cue quality CAN extinguish the optimum ' ...
             'and figL4 needs p_0 as an axis\n'], ...
            sum(~isfinite(dstar(:))), numel(dstar));
end
out.margin = min(spd(:));
fprintf('smallest speedup anywhere on the sweep: %.1f%% (figL4 contours 1%%)\n', out.margin);
if out.margin < 2
    fprintf(['  the p_0-independence of rho*(q) holds to within about ' ...
             'a percent\n']);
end
ordD = all(diff(dstar,1,1) > 0, 'all');      % worse cue -> abandon sooner
ordS = all(diff(spd,1,1)  < 0, 'all');       % worse cue -> smaller payoff
fprintf('ordering with p_0: delta* monotone %d, speedup decreasing %d\n', ordD, ordS);
if ~ordD && ordS
    kx = find(any(diff(sign(diff(dstar,1,1)),1,2),1), 1);
    fprintf(['the delta* ordering reverses in q: a noisier cue is ' ...
             'abandoned sooner at small q and held longer at large q, ' ...
             'because the benefit coefficient in the large-delta expansion ' ...
             'carries a factor q*N and the congestion cost does not. Only ' ...
             'the payoff ordering is monotone in p_0\n']);
    if ~isempty(kx)
        fprintf('  the curves cross near q = %.2f\n', qgrid(kx));
    end
elseif ~ordS
    fprintf(['the speedup ordering is not monotone in p_0, which the ' ...
             'expansion does not predict\n']);
end
P.export(P, 'figL5', outdir, h);
end
