function out = figure5(varargin)
%FIGURE6 Reproduce the three panels of Figure 5 (competing donors) and assemble
%   the lettered figure output/figs/fig5.pdf.
%
% Usage
%   figure5
%   figure5('fast',true)                  % low-cost wiring/smoke test
%   figure5('outdir','output/figure5')
%
% Panel map
%   (a) figT3   total recruitment benefit versus service asymmetry
%   (b) figT1   feeding time and cue contrast versus abandonment rate
%   (c) figT2   optimal abandonment versus service asymmetry
%
% The exact two-donor arrival table is built by twogain (on atwo) in
% twodonor_setup and the reduced system by twodonor. A normal call uses the
% production grids. 'fast' uses coarse tables as an end-to-end smoke test.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'src')));

ip = inputParser;
addParameter(ip,'fast',false);
addParameter(ip,'outdir',fullfile('output','figure5'));
parse(ip,varargin{:});
fast = logical(ip.Results.fast);
OUTDIR = ip.Results.outdir;
if ~isfolder(OUTDIR), mkdir(OUTDIR); end
close all

S = twodonor_setup(OUTDIR, fast, 'figure5', false);
G2 = S.G2;  res = S.res;  cfg0 = S.cfg0;

if fast
    T1 = panel_b(G2,cfg0,OUTDIR,20,logspace(-2.5,1.0,24),0.5,0.35);
    T2 = panel_c(res,OUTDIR);
    T3 = panel_a(res,OUTDIR,true);
else
    T1 = panel_b(G2,cfg0,OUTDIR);
    T2 = panel_c(res,OUTDIR);
    T3 = panel_a(res,OUTDIR);
end

out = struct('G2',G2,'res',res,'T1',T1,'T2',T2,'T3',T3);
save(fullfile(OUTDIR,'figure5_results.mat'),'-struct','out');

% Lettered figure, one row of three panel canvases at the panel size.
figdir = fullfile('output','figs');
if ~isfolder(figdir), mkdir(figdir); end
T = [22 19];
tiles = [0 0 T; T(1) 0 T; 2*T(1) 0 T];
panels = fullfile(OUTDIR, strcat({'figT3','figT1','figT2'}, '.fig'));
assemble_figure(fullfile(figdir,'fig5'), panels, tiles, [3*T(1) T(2)]);

fprintf('Figure 5 complete. Panel order: figT3, figT1, figT2.\n');
end

% =========================================================================
% Panel (a), figT3. Total recruitment benefit.
% =========================================================================
function out = panel_a(res, outdir, showbase)
%PANEL_A  Total recruitment benefit against service asymmetry.
%   Benefit is 100*(Tnone-Tfull)/Tnone, with each donor-specific cue run at
%   its optimal abandonment rate. The unaided feeding time is shown on the
%   right axis because it also changes with asymmetry.
%   The minimum is interpolated within its neighboring log-r grid points.

if nargin < 2 || isempty(outdir),   outdir = 'output'; end
if nargin < 3 || isempty(showbase), showbase = true;   end
P = linger_style();

h = P.figure(P);
% The baseline goes on a second, transparent axes on top rather than yyaxis,
% which copyobj cannot copy at assembly.
POS = [0.17 0.17 0.66 0.73];
ax  = axes(h, 'Units','normalized', 'Position', POS, ...
           'PositionConstraint','innerposition');
plot(ax, res.r, res.total, '-', 'Color', P.greens(4,:), 'LineWidth', P.LW);
hold(ax, 'on');
set(ax, 'XScale', 'log');

% ---- minimum of the total ---------------------------------------------
% The sweep samples r finely enough near the dip that the marker sits on the
% sampled minimum, so it lies on the plotted curve.
[tg, kg] = min(res.total);
rmin = res.r(kg);  tmin = tg;  bracket = res.r(max(kg-1,1):min(kg+1,numel(res.r)));
bracket = [bracket(1) bracket(end)];
if kg == 1 || kg == numel(res.r)
    warning('figT3:minAtEdge', ...
        ['the total is smallest at the edge of rlist (r = %g), so the dip ' ...
         'is not bracketed; widen the sweep'], res.r(kg));
end
plot(ax, [rmin rmin], [0 tmin], ':', 'Color', [0.45 0.45 0.45], 'LineWidth', P.ALW-0.5);
plot(ax, rmin, tmin, 'o', 'Color','k','MarkerFaceColor','w', ...
     'MarkerSize', 14, 'LineWidth', P.ALW);

xlabel(ax, 'service asymmetry $r=\mu_1/\mu_2$','Interpreter','latex','FontSize',P.LSZ);
ylabel(ax, 'benefit (\% of $T_{1/2}$ unaided)','Interpreter','latex','FontSize',P.LSZ);
xlim(ax, [min(res.r) max(res.r)]);  ylim(ax, [0 1.1*max(res.total)]);
P.axes(P, ax);

if showbase
    axr = axes(h, 'Units','normalized', 'Position', POS, 'Color','none', ...
               'YAxisLocation','right', 'XScale','log');
    hold(axr, 'on');
    plot(axr, res.r, res.Tnone, '-', 'Color', [0.72 0.72 0.72], ...
         'LineWidth', P.ALW);
    ylabel(axr, '$T_{1/2}$ unaided','Interpreter','latex','FontSize',P.LSZ);
    xlim(axr, xlim(ax));
    ylim(axr, 50*[floor(0.8*min(res.Tnone)/50), ceil(1.2*max(res.Tnone)/50)]);
    P.axes(P, axr);
    set(axr, 'YColor', [0.55 0.55 0.55], 'XColor','none', 'Color','none');
    set([ax axr], 'PositionConstraint','innerposition', 'Position', POS);
    text(axr, 17, 103, {'unaided','feeding time'}, 'Interpreter','latex', ...
        'FontSize',P.FSZ,'Color',[0.55 0.55 0.55], ...
        'HorizontalAlignment','center','VerticalAlignment','middle');
end

text(ax, 12, 57, {'recruitment','benefit'}, 'Interpreter','latex', ...
    'FontSize',P.FSZ,'Color',P.greens(4,:), ...
    'HorizontalAlignment','center','VerticalAlignment','middle');

out = struct('r',res.r,'total',res.total,'Tnone',res.Tnone, ...
             'rmin',rmin,'tmin',tmin,'rmin_bracket',bracket, ...
             'base_growth',res.Tnone(end)/res.Tnone(1));
fprintf(['figT3: benefit %.1f%% at r=1, minimum %.1f%% at r=%.2f ' ...
         '(bracketed by %g and %g), %.1f%% at r=%g\n'], ...
        res.total(1),tmin,rmin,bracket(1),bracket(2),res.total(end),res.r(end));
fprintf('       unaided feeding time grows by a factor %.2f\n',out.base_growth);
% Inset, a fast and a slow donor. The fast donor carries the scenters and the
% strong cue, and searchers stream toward it.
inset_twodonor(h, P);
P.export(P, 'figT3', outdir, h);
end

% ------------------------------------------------------------------------
function inset_twodonor(h, P)
%INSET_TWODONOR  Two-donor sketch, in the empty upper left above the total.
K = inset_kit();
ax = K.axes(h, P, [4.3 9.6 7.0 6.0]);
y0 = 1.4;  x1 = 4.9;  x2 = 1.2;
K.tent(ax, x1 - 0.15, 1.5, y0, 3.6, 0.22, 0.90);
K.tent(ax, x2 + 0.1, 0.75, y0, 1.1, 0.68, 0.96);
K.track(ax, 0.1, 7.0, y0);
K.donor(ax, x1, y0, 0.4);
K.donor(ax, x2, y0, 0.2);
B = [x1-1.1 y0+0.34 -0.25; x1+1.0 y0+0.36 pi+0.25; x1-0.95 y0-0.38 0.3; x1+0.95 y0-0.38 pi-0.3];
for k = 1:size(B,1), K.bee(ax, B(k,1), B(k,2), B(k,3), 0.75, 'scent'); end
K.bee(ax, x2 + 0.55, y0 + 0.3, pi + 0.3, 0.6, 'scent');
S = [2.45 y0 0; 3.2 y0 0; 6.65 y0 pi];
for k = 1:size(S,1), K.bee(ax, S(k,1), S(k,2), S(k,3), 0.7, 'search'); end
K.label(ax, x1, y0 - 0.5, '$\mu_1$', P, 'top');
K.label(ax, x2, y0 - 0.3, '$\mu_2$', P, 'top');
end

% =========================================================================
% Panel (b), figT1. Feeding time and cue contrast against abandonment.
% =========================================================================
function out = panel_b(G2, cfg0, outdir, r, dgrid, q, tol)
%PANEL_B  Cue contrast and feeding time against abandonment rate, at fixed
%   service asymmetry. The contrast Lambda_1/Lambda_2 peaks close to the delta
%   that minimizes T_q. delta^ast is set by where the cue discriminates, not
%   by where it is sharpest.
%
%   Both optima are refined off the grid. One step of dgrid is a factor of
%   about 1.2, so grid argmins in adjacent cells differ by about 21 percent
%   from quantization alone.
%
%   eq:contrast_condition locates where the contrast is stationary at a given
%   representative parameters, not where T_q is least. The panel shows that the two fall
%   together to within the resolution of the sweep. out.sep_steps gives the
%   separation in grid spacings and out.sep_rel as a fraction. Quote the
%   tolerance with either.
%
%   sep_rel must stay under tol with both optima refined. If it does not,
%   refine dgrid near the optimum before treating the separation as real,
%   since a coarse grid biases each refinement toward its bracket.

if nargin < 3 || isempty(outdir), outdir = 'output'; end
if nargin < 4 || isempty(r),      r = 20; end
if nargin < 5 || isempty(dgrid),  dgrid = logspace(-2.5, 1.0, 44); end
if nargin < 6 || isempty(q),      q = 0.5; end
if nargin < 7 || isempty(tol),    tol = 0.25; end   % one grid step is 0.21

M2 = twodonor();  M = lingermodel();  P = linger_style();
cfg = cfg0;  cfg.mu2 = cfg0.mutot/(1+r);  cfg.mu1 = r*cfg.mu2;

T = zeros(size(dgrid));  C = zeros(size(dgrid));
F1 = zeros(size(dgrid)); F2 = zeros(size(dgrid));
for j = 1:numel(dgrid)
    [T(j), C(j), F1(j), F2(j)] = sweepPoint(dgrid(j), G2, cfg, q, M2);
end

% ---- refine both optima off the grid ------------------------------------
% refineDelta minimizes, so the contrast is passed negated. Each call costs a
% few ODE solves.
fT = @(d) M2.Tq(d, G2, cfg, 'full', q);
fC = @(d) -contrastAt(d, G2, cfg, q, M2);
[dT, Topt] = M.refineDelta(fT, dgrid, T);
[dC, mCopt] = M.refineDelta(fC, dgrid, -C);
Copt = -mCopt;

step     = dgrid(2)/dgrid(1);
sep_rel  = abs(dT/dC - 1);
sep_step = abs(log(dT/dC))/log(step);

[~, kT] = min(T);  [~, kC] = max(C);        % grid argmins, for the markers
[~, jT] = min(abs(dgrid - dT));
F1o = interp1(dgrid, F1, dT);  F2o = interp1(dgrid, F2, dT);

% ---- panel ---------------------------------------------------------------
% Only the contrast peak is marked. The dotted rule sits at the refined
% feeding-time optimum, so the gap between marker and rule is sep_rel.
h = P.figure(P);                            % handle captured, not left to gcf
% Two overlaid axes rather than yyaxis, which copyobj cannot copy at assembly.
POS = [0.17 0.17 0.66 0.73];
ax  = axes(h, 'Units','normalized', 'Position', POS, ...
           'PositionConstraint','innerposition');
hold(ax, 'on');
semilogx(ax, dgrid, T, '-', 'Color', P.greens(4,:), 'LineWidth', P.LW);
xline(ax, dT, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', P.ALW);
set(ax, 'XScale','log');
xlabel(ax, 'abandonment rate $\delta$','Interpreter','latex','FontSize',P.LSZ);
ylabel(ax, 'half-colony feeding time $T_{1/2}$','Interpreter','latex', ...
       'FontSize',P.LSZ);
xlim(ax, [min(dgrid) max(dgrid)]);
ylim(ax, 10*[floor(0.98*min(T)/10), ceil(1.02*max(T)/10)]);   % room for the line width at the minimum
P.axes(P, ax);
set(ax, 'YColor', P.greens(4,:));

axr = axes(h, 'Units','normalized', 'Position', POS, 'Color','none', ...
           'YAxisLocation','right', 'XScale','log');
hold(axr, 'on');
plot(axr, dgrid, C, '-', 'Color', P.blues(3,:), 'LineWidth', P.LW);
plot(axr, dC, Copt, 'o', 'Color', P.blues(3,:), ...
     'MarkerFaceColor','w', 'MarkerSize', 14, 'LineWidth', P.ALW);
ylabel(axr, 'cue contrast $\widetilde\Lambda_1/\widetilde\Lambda_2$', ...
       'Interpreter','latex','FontSize',P.LSZ);
xlim(axr, xlim(ax));
P.axes(P, axr);
set(axr, 'YColor', P.blues(3,:), 'XColor','none', 'Color','none');
set([ax axr], 'PositionConstraint','innerposition', 'Position', POS);

yl = ylim(axr);
text(axr, dT, yl(2), '$\delta^\ast$', 'Interpreter','latex', ...
     'FontSize', P.LSZ, 'Color', [0.6 0.6 0.6], ...
     'HorizontalAlignment','left', 'VerticalAlignment','bottom');

out = struct('r',r,'delta',dgrid,'T',T,'contrast',C,'F1',F1,'F2',F2, ...
             'dT',dT,'dC',dC,'Topt',Topt,'Copt',Copt, ...
             'dT_grid',dgrid(kT),'dC_grid',dgrid(kC), ...
             'sep_rel',sep_rel,'sep_steps',sep_step,'tol',tol, ...
             'step',step,'F1opt',F1o,'F2opt',F2o,'q',q);

fprintf(['figT1: r = %g, T minimized at delta = %.4f, contrast maximized at ' ...
         '%.4f\n'], r, dT, dC);
fprintf(['       separation %.1f%% of the optimum, %.2f grid steps ' ...
         '(one step = %.1f%%)  ->  %s\n'], 100*sep_rel, sep_step, ...
        100*(step-1), string(sep_rel <= tol));
fprintf('       on the grid alone they sat at %.4f and %.4f, %d cells apart\n', ...
        dgrid(kT), dgrid(kC), abs(kT-kC));
fprintf('       at the optimum F1 = %.2f, F2 = %.2f\n', F1o, F2o);

if sep_rel > tol
    warning('figT1:optimaSeparate', ...
        ['the contrast maximum and the feeding-time minimum differ by %.1f%%, ' ...
         'above the %.0f%% tolerance, with both refined off the grid.  The ' ...
         'caption claims they fall together; refine dgrid near the optimum ' ...
         'before treating the separation as real'], 100*sep_rel, 100*tol);
end
if sep_step < 0.5 && abs(kT-kC) >= 1
    fprintf(['       note: the refinement moved them into the same cell, so ' ...
             'the grid argmins were quantization and not a result\n']);
end

P.export(P, 'figT1', outdir, h);
end

% ------------------------------------------------------------------------
function [T, C, F1, F2] = sweepPoint(d, G2, cfg, q, M2)
[T, sol] = M2.Tq(d, G2, cfg, 'full', q);
y = sol.y(:,end);  F1 = y(3);  F2 = y(5);
[L1, L2] = M2.rates(G2, F1, F2, 'full');
C = L1/L2;
end

function C = contrastAt(d, G2, cfg, q, M2)
%CONTRASTAT  Cue contrast at the feeding threshold, as a function of delta.
%   The fed populations are read at the stopping event, matching
%   twodonor.decompose, so the contrast is a representative-parameter quantity, not a
%   run average.
[~, C] = sweepPoint(d, G2, cfg, q, M2);
end

% =========================================================================
% Panel (c), figT2. Optimal abandonment versus service asymmetry.
% =========================================================================
function out = panel_c(res, outdir)
%PANEL_C  Optimal abandonment rate against service asymmetry.
%   Uses the donor-specific optimum saved by M2.sweep; no comparison cue.
if nargin < 2 || isempty(outdir), outdir = 'output'; end
P = linger_style();
h = P.figure(P);
ax = gca;
loglog(ax, res.r, res.dfull, '-', 'Color', P.greens(4,:), 'LineWidth', P.LW);
hold(ax, 'on');
plot(ax, res.r, res.dfull, 'o', 'Color', P.greens(4,:), 'MarkerFaceColor','w', ...
     'MarkerSize', 10, 'LineWidth', P.ALW-0.5);
xlabel('service asymmetry $r=\mu_1/\mu_2$','Interpreter','latex','FontSize',P.LSZ);
ylabel('optimal abandonment rate $\delta^\ast$','Interpreter','latex','FontSize',P.LSZ);
xlim([min(res.r) max(res.r)]);
ylim([0.8*min(res.dfull) 1.25*max(res.dfull)]);
P.axes(P);
yt = [0.01 0.02 0.05 0.1 0.2 0.3 0.5 1]; yl = ylim(ax);
yt = yt(yt >= yl(1) & yt <= yl(2));
set(ax, 'YTick', yt, 'YTickLabel', arrayfun(@(y) sprintf('$%g$', y), yt, 'uni', 0));
set(ax, 'PositionConstraint','innerposition', 'Position', [0.18 0.17 0.78 0.73]);
% Dotted line at the r of minimum benefit in panel (a), where the slower donor
% first receives more searchers than it can serve.
[~, kmin] = min(res.total);  rmin = res.r(kmin);
plot(ax, [rmin rmin], ylim(ax), ':', 'Color', [0.45 0.45 0.45], 'LineWidth', P.ALW-0.5);
uistack(findobj(ax,'Type','line','Color',P.greens(4,:)),'top');
d = res.dfull(:).';
q.nrep = sum(abs(diff(d))./max(abs(d(1:end-1)),eps) < 1e-9);
q.nlev = numel(uniquetol(d,1e-9,'DataScale',max(abs(d))));
if q.nrep > 0
    warning('figT2:quantized','delta* repeats at %d consecutive r; check optimum refinement.',q.nrep);
end
out = struct('r',res.r,'dfull',res.dfull,'quant_full',q, ...
             'range_factor',max(d)/min(d));
fprintf('figT2: optimal abandonment varies by a factor %.2f over r=%g to %g\n', ...
        out.range_factor,res.r(1),res.r(end));
P.export(P, 'figT2', outdir, h);
end
