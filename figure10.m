function out = figure10(varargin)
%FIGURE7 Reproduce the three panels of Figure 10 (limits and robustness)
%   and assemble the lettered figure.
%
% Usage
%   figure10
%   figure10('fast',true)                  % coarse grids, few realizations
%   figure10('outdir','output/figure10')
%
% Panel map
%   (a) fig7a  T_q(delta)/T_q(inf) with the large-delta expansion and the
%              delta = 0 values, q = 0.25, 0.5, 0.9
%   (b) fig7b  reduction error against delta scaled by the spectral gap
%   (c) fig7c  half-colony feeding time under four service laws
%
% Caches. The single-donor gain tables and the service-law and reduction-error
% ensembles are cached in outdir under stamped keys. The
% two-donor table is read from output/figure4/figure4_twogain_cache.mat and
% the two-donor ensemble from output/figure4/figT5_abm.mat. Neither is
% written. A fast run writes its own caches with a _fast suffix.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'src')));

ip = inputParser;
addParameter(ip,'fast',false);
addParameter(ip,'outdir',fullfile('output','figure10'));
addParameter(ip,'figdir',fullfile('output','figs'));
parse(ip,varargin{:});
fast   = logical(ip.Results.fast);
OUTDIR = ip.Results.outdir;
FIGDIR = ip.Results.figdir;
if ~isfolder(OUTDIR), mkdir(OUTDIR); end
if ~isfolder(FIGDIR), mkdir(FIGDIR); end
close all

VER  = 'vS1-si-panels';            % cache version tag
F6   = fullfile('output','figure6');
F7   = fullfile('output','figure7');
if fast, sfx = '_fast'; else, sfx = ''; end

A = afpt();
M = lingermodel();
P = linger_style();
prm = A.mkprm('kbar',0.0314,'sig',0.8,'beta0',1,'sm',1,'v',1);
cfg = M.mkcfg('N',120,'mu',1,'lam0',prm.lam0,'p0',0.75);
fprintf('figure10: lambda0 = %.5f, rho = %.3f, q0 = %.5f\n', ...
        prm.lam0, prm.lam0*cfg.N/cfg.mu, prm.theta);

% ---- single-donor tables ------------------------------------------------
if fast
    pgrid = linspace(0,1.2,31);   ggrid = linspace(0,prm.gmax,25);
else
    pgrid = linspace(0,1.2,61);   ggrid = linspace(0,prm.gmax,51);
end
cachefile = fullfile(OUTDIR, ['opgain_cache' sfx '.mat']);
key = M.stamp(VER, prm, cfg, pgrid, ggrid);
G = [];
if isfile(cachefile)
    S = load(cachefile);
    if isfield(S,'key') && strcmp(S.key,key)
        G = S.G;  W = S.W;  ps = S.ps;
        fprintf('figure10: loaded stamped gain cache\n');
    end
end
if isempty(G)
    G  = M.gainTable(prm,pgrid);
    W  = M.gainSurface(prm,pgrid,ggrid);
    ps = M.psiTable(G,cfg.p0,W);
    save(cachefile,'G','W','ps','key');
end
fprintf('figure10 checks: p*/sigma = %.3f, G(0) = %.3f, g1 = %.3f\n', ...
        G.pstar/prm.sig, G.Graw(1), ps.g1);

% ---- two-donor table, production grids, read only ----------------------
prm2 = struct('kbar',0.0314,'sig',0.8,'beta0',1,'sm',1,'v',1,'p0',0.75);
prm2.gmax = prm2.sm*prm2.sig - 1e-6;
PGRID = linspace(-prm2.sig,prm2.sig,25);  KMAX = 25;  GGRID = linspace(0,prm2.gmax,13);
want = twogain_key(prm2,PGRID,KMAX,GGRID);
G2 = [];
for f = {fullfile(F6,'figure6_twogain_cache.mat'), fullfile(OUTDIR,'figure10_twogain_cache.mat')}
    if isfile(f{1})
        C = load(f{1},'G2');
        if isfield(C,'G2') && isfield(C.G2,'key') && strcmp(C.G2.key,want)
            G2 = C.G2;
            fprintf('figure10: two-donor table read from %s\n', f{1});
            break
        end
    end
end
if isempty(G2)
    G2 = twogain(prm2,PGRID,KMAX,GGRID);
    save(fullfile(OUTDIR,'figure10_twogain_cache.mat'),'G2');
end

% ---- panels -------------------------------------------------------------
two = struct('r',4,'nrep',200,'deltas',[0.02 0.05 0.1 0.2 0.5 1 2 5], ...
             'cachefile',fullfile(F7,'figT5_abm.mat'), ...
             'fallback',fullfile(OUTDIR,'fig7b_abm2.mat'));
if fast
    Ra = panel_a(M,P,ps,cfg,OUTDIR,logspace(-3,2,20));
    Rb = panel_b(M,P,G,ps,cfg,G2,OUTDIR,[0.05 0.2 1 5],12, ...
                    fullfile(OUTDIR,'fig7b_abm_fast.mat'),two);
    Rc = panel_c(M,P,G,ps,cfg,OUTDIR,logspace(log10(0.05),log10(3),5),12, ...
                    fullfile(OUTDIR,'fig7c_abm_fast.mat'));
else
    Ra = panel_a(M,P,ps,cfg,OUTDIR);
    Rb = panel_b(M,P,G,ps,cfg,G2,OUTDIR,[],[],fullfile(OUTDIR,'fig7b_abm.mat'),two);
    Rc = panel_c(M,P,G,ps,cfg,OUTDIR,[],[],fullfile(OUTDIR,'fig7c_abm.mat'));
end

% ---- checks quoted in the text, printout only ----------------------------
Klim = check_linger_limit(M, cfg);
Ksym = check_symmetry(M, ps, G2);

out = struct('prm',prm,'cfg',cfg,'g1',ps.g1,'A',Ra,'B',Rb,'C',Rc, ...
             'limit',Klim,'symmetry',Ksym);
save(fullfile(OUTDIR,['figure10_results' sfx '.mat']),'-struct','out');

% ---- assembly -----------------------------------------------------------
names  = {'fig7a','fig7b','fig7c'};
panels = cellfun(@(n) fullfile(OUTDIR,[n '.fig']), names, 'uni', 0);
tiles  = [0 0 22 19; 22 0 22 19; 44 0 22 19];
assemble_figure(fullfile(FIGDIR,'fig10'), panels, tiles, [66 19], ...
                'letterxy', [0.2 -0.2]);
fprintf('figure10 complete. Panels in %s, assembled figure %s\n', ...
        OUTDIR, fullfile(FIGDIR,'fig10.pdf'));
end

% =========================================================================
% Panel (a), fig7a. The two abandonment limits.
% =========================================================================
function out = panel_a(M, P, ps, cfg, outdir, dgrid, qlist)
%PANEL_A  The two abandonment limits and the interior optimum. Panel (a).
%
%   T_q(delta)/T_q(inf) from the reduced dynamics at three fed fractions, with
%   the large-delta expansion (main Eq. large_delta)
%
%       T_q = L_q/lam0 + (1/delta)*[A_q/(mu - A_q) - (g1-1)*q*N],
%
%   dashed at large delta, and two delta = 0 values at the left edge. The
%   filled marker integrates the reduced dynamics at delta = 0. The open
%   marker is the quasi-stationary integral of appendix Eq. T0,
%
%       T_q(0) = int_0^{qN} dF / A(F),
%       A = [(mu + Gam N) - sqrt((mu + Gam N)^2 - 4 Gam mu (N - F))]/2,
%
%   with Gam = lam0*Psi~(F). The two differ because the aggregation does not
%   relax quickly at delta = 0, which the printout quantifies.

if nargin < 6 || isempty(dgrid), dgrid = logspace(-3, 2, 60); end
if nargin < 7 || isempty(qlist), qlist = [0.25 0.5 0.9]; end

N = cfg.N;  mu = cfg.mu;  lam0 = cfg.lam0;  g1 = ps.g1;
nq = numel(qlist);  nd = numel(dgrid);
R = zeros(nq, nd);  Rexp = zeros(nq, nd);
T0ode = zeros(1,nq);  T0int = zeros(1,nq);  Tinf = zeros(1,nq);
slope = zeros(1,nq);  dopt = zeros(1,nq);  Ropt = zeros(1,nq);

fprintf('fig7a: g1 = %.4f, lam0 = %.5f, rho = %.3f\n', g1, lam0, lam0*N/mu);
fprintf(['fig7a: first chord of the Psi table, (Psi(F_2)-1)/F_2 = %.4f at F_2 = %.3g, ' ...
         'against g1 - 1 = %.4f\n'], (ps.psi(2) - 1)/ps.F(2), ps.F(2), g1 - 1);
for i = 1:nq
    q = qlist(i);
    f = @(d) M.Tq(d, ps, cfg, q);
    Tinf(i) = f(1e4);
    Lq = log(1/(1-q))/lam0;
    T  = arrayfun(f, dgrid);
    R(i,:) = T/Tinf(i);
    [dopt(i), Topt] = M.refineDelta(f, dgrid, T);
    Ropt(i) = Topt/Tinf(i);

    Aq = lam0*(1-q)*N;
    slope(i) = Aq/(mu - Aq) - (g1 - 1)*q*N;
    Rexp(i,:) = (Lq + slope(i)./dgrid)/Tinf(i);

    T0ode(i) = f(0);
    Gam = @(F) lam0*M.psiAt(ps, F);
    Af  = @(F) 0.5*((mu + Gam(F)*N) - sqrt((mu + Gam(F)*N).^2 - 4*Gam(F)*mu.*(N - F)));
    T0int(i) = integral(@(F) 1./arrayfun(Af, F), 0, q*N, 'RelTol', 1e-8);

    % bracket of the expansion, extrapolated from the reduced dynamics in 1/delta
    dd = [300 1000];
    b  = arrayfun(@(d) d*(f(d) - Lq), dd);
    bx = b(2) + (b(2) - b(1))*(1/dd(2))/(1/dd(1) - 1/dd(2));
    fprintf(['fig7a: q = %.2f  T(inf) = %.2f (L_q/lam0 = %.2f)  optimum delta* = %.3f, ' ...
             'T*/T(inf) = %.3f\n'], q, Tinf(i), Lq, dopt(i), Ropt(i));
    fprintf(['        bracket: expansion %.3f, extrapolated %.3f   ' ...
             'T(0): reduced dynamics %.2f (ratio %.3f), Eq. T0 %.2f (ratio %.3f), ' ...
             'Eq. T0 off by %+.1f%%\n'], slope(i), bx, T0ode(i), T0ode(i)/Tinf(i), ...
            T0int(i), T0int(i)/Tinf(i), 100*(T0int(i)/T0ode(i) - 1));
    fprintf('        reduced ratio at delta = %.0e is %.4f\n', dgrid(1), R(i,1));
end
if any(T0ode < Tinf)
    fprintf(['fig7a: with recruitment T_q(0) < T_q(inf) at q = %s, so delta = 0 is ' ...
             'faster than immediate departure there\n'], mat2str(qlist(T0ode < Tinf)));
end

% ---- panel --------------------------------------------------------------
cols = [P.greens(3,:); P.greens(5,:); P.blues(3,:)];
h  = P.figure(P);
ax = axes('Parent', h);  hold(ax, 'on');
yl = [0.35 1.1];
yline(ax, 1, ':', 'Color', [0 0 0], 'LineWidth', P.ALW);
for i = 1:nq
    ok = dgrid >= 3*abs(slope(i))/Tinf(i);        % correction at most a third
    plot(ax, dgrid(ok), Rexp(i,ok), '--', 'Color', cols(i,:), 'LineWidth', P.ALW);
end
for i = 1:nq
    plot(ax, dgrid, R(i,:), '-', 'Color', cols(i,:), 'LineWidth', P.LW);
    plot(ax, dopt(i), Ropt(i), 'o', 'Color', cols(i,:), 'MarkerFaceColor', 'w', ...
         'MarkerSize', 11, 'LineWidth', P.ALW);
end
% delta = 0 values on the left spine
for i = 1:nq
    plot(ax, dgrid(1), T0ode(i)/Tinf(i), 'o', 'Color', cols(i,:), ...
         'MarkerFaceColor', cols(i,:), 'MarkerSize', 13, 'Clipping', 'off');
    plot(ax, dgrid(1), T0int(i)/Tinf(i), 'd', 'Color', cols(i,:), ...
         'MarkerFaceColor', 'w', 'MarkerSize', 13, 'LineWidth', P.ALW, 'Clipping', 'off');
end
set(ax, 'XScale', 'log');
xlim(ax, [dgrid(1) dgrid(end)]);  ylim(ax, yl);
set(ax, 'XTick', 10.^(-3:2));
xlabel(ax, 'abandonment rate $\delta$', 'Interpreter','latex', 'FontSize', P.LSZ);
ylabel(ax, 'feeding time ratio $T_q(\delta)/T_q(\infty)$', 'Interpreter','latex', ...
       'FontSize', P.LSZ);
P.axes(P, ax);
set(ax, 'Position', [0.17 0.16 0.78 0.74]);

% in-place labels, just right of the left spine where the curves separate
for i = 1:nq
    P.label(P, ax, 2.2e-3, R(i, find(dgrid >= 2.2e-3, 1)) + 0.012, ...
            sprintf('$q=%.2g$', qlist(i)), 'Color', cols(i,:), ...
            'HorizontalAlignment','left', 'VerticalAlignment','bottom');
end
P.label(P, ax, 90, 1.012, 'no lingering', 'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom', 'FontSize', P.FSZ-4);
% key to the delta = 0 markers, right of the expansions
P.label(P, ax, 2.0, 0.57, 'at $\delta=0$', 'HorizontalAlignment','left', ...
        'VerticalAlignment','middle', 'FontSize', P.FSZ-4);
plot(ax, 2.6, 0.50, 'o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'MarkerSize', 11);
plot(ax, 2.6, 0.43, 'd', 'Color', [0 0 0], 'MarkerFaceColor', 'w', 'MarkerSize', 11, ...
     'LineWidth', P.ALW);
P.label(P, ax, 4.0, 0.50, 'reduced', 'HorizontalAlignment','left', ...
        'VerticalAlignment','middle', 'FontSize', P.FSZ-4);
P.label(P, ax, 4.0, 0.43, 'quasi-stationary', 'HorizontalAlignment','left', ...
        'VerticalAlignment','middle', 'FontSize', P.FSZ-4);
P.label(P, ax, 1.6, 0.72, 'large $\delta$ expansion', 'Color', [0.25 0.25 0.25], ...
        'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize', P.FSZ-4);

P.export(P, 'fig7a', outdir, h);
out = struct('q',qlist,'dgrid',dgrid,'R',R,'Rexp',Rexp,'Tinf',Tinf,'slope',slope, ...
             'T0ode',T0ode,'T0int',T0int,'dopt',dopt,'Ropt',Ropt,'g1',g1);
end

% =========================================================================
% Panel (b), fig7b. Reduction error against the spectral gap.
% =========================================================================
function out = panel_b(M, P, G, ps, cfg, G2, outdir, deltas, nrep, cachefile, two)
%PANEL_B  Reduction error against abandonment rate scaled by the spectral
%   gap, one donor and two. Panel (b).
%
%       E(delta) = [T_q^sim - T_q^red] / T_q^red
%
%   One donor. abm ensemble cached under cachefile with a stamped key.
%   Two donors at r = 4. abm2 ensemble read from two.cachefile, Fig. 7(b)'s cache,
%   when its key matches the table G2, otherwise run into two.fallback. The
%   reduction is recomputed here from twodonor on G2.
%
%   Gaps from spectralgap at full reliance and p = 0.3, 0.449 for one plume
%   and 0.274 for two.

if nargin < 8 || isempty(deltas), deltas = logspace(log10(0.05),log10(6),11); end
if nargin < 9 || isempty(nrep),   nrep   = 150; end

q    = 0.5;
[GAP1, GAP2] = spectralgap(G.prm);
L    = 2*pi/G.prm.v;
Tinf = log(1/(1-q))/cfg.lam0;

% ---- one donor ----------------------------------------------------------
key = sprintf('vS4b|n=%d|d=%s|p0=%g', nrep, mat2str(deltas,4), ps.p0);
Tsim = [];
if isfile(cachefile)
    C = load(cachefile);
    if isfield(C,'key') && strcmp(C.key,key)
        Tsim = C.Tsim;  SEs = C.SEs;
        fprintf('fig7b: reusing %d cached one-donor realizations per delta\n', nrep);
    end
end
if isempty(Tsim)
    nd = numel(deltas);  Tsim = nan(1,nd);  SEs = nan(1,nd);
    for id = 1:nd
        v = arrayfun(@(k) M.abm(deltas(id), G, ps, cfg, q, k), 1:nrep);
        v = v(isfinite(v));
        Tsim(id) = mean(v);  SEs(id) = std(v)/sqrt(numel(v));
        fprintf('  delta = %.3f done, T = %.1f +- %.1f\n', deltas(id), Tsim(id), SEs(id));
    end
    save(cachefile, 'Tsim', 'SEs', 'key');
end
Tred = arrayfun(@(d) M.Tq(d, ps, cfg, q), deltas);
E1   = (Tsim - Tred)./Tred;
SE1  = SEs./Tred;

% ---- two donors ---------------------------------------------------------
M2  = twodonor();
c2  = struct('N', cfg.N, 'mutot', 2.0);
c2.mu2 = c2.mutot/(1 + two.r);  c2.mu1 = two.r*c2.mu2;
key2 = sprintf('%s|r=%g|q=%g|n=%d|d=%s|%s', 'v1-two-donor-asymptotics', two.r, q, ...
               two.nrep, mat2str(two.deltas,4), G2.key);
mu2 = [];
for f = {two.cachefile, two.fallback}
    if isfile(f{1})
        C = load(f{1});
        if isfield(C,'key') && strcmp(C.key, key2)
            mu2 = C.mu_;  se2 = C.se;
            fprintf('fig7b: two-donor ensemble read from %s (%d realizations)\n', f{1}, two.nrep);
            break
        end
    end
end
if isempty(mu2)
    fprintf('fig7b: no matching two-donor ensemble, running abm2 into %s\n', two.fallback);
    nd = numel(two.deltas);  mu2 = zeros(1,nd);  se2 = zeros(1,nd);
    for j = 1:nd
        v = arrayfun(@(k) abm2(two.deltas(j), G2, c2, q, k), 1:two.nrep);
        v = v(isfinite(v));
        mu2(j) = mean(v);  se2(j) = std(v)/sqrt(numel(v));
    end
    mu_ = mu2;  se = se2;  key = key2;                               %#ok<NASGU>
    save(two.fallback, 'mu_', 'se', 'key');
end
d2   = two.deltas;
Tr2  = arrayfun(@(d) M2.Tq(d, G2, c2, 'full', q), d2);
E2   = (mu2 - Tr2)./Tr2;
SE2  = se2./Tr2;

% ---- printout -----------------------------------------------------------
fprintf('\n  one donor                                  two donors, r = %g\n', two.r);
fprintf('  delta  delta/gap   error              delta  delta/gap   error\n');
for k = 1:max(numel(deltas), numel(d2))
    s1 = repmat(' ', 1, 43);  s2 = '';
    if k <= numel(deltas)
        s1 = sprintf('  %5.3f  %7.3f   %+5.1f +- %3.1f%%     ', deltas(k), deltas(k)/GAP1, ...
                     100*E1(k), 100*1.96*SE1(k));
    end
    if k <= numel(d2)
        s2 = sprintf('%5.3f  %7.3f   %+5.1f +- %3.1f%%', d2(k), d2(k)/GAP2, 100*E2(k), 100*1.96*SE2(k));
    end
    fprintf('%s%s\n', s1, s2);
end
lo1 = deltas/GAP1 < 1;  lo2 = d2/GAP2 < 1;
[pk1, i1] = max(E1);
fprintf(['fig7b: one donor, mean error %.1f%% below the gap, peak %.1f%% at delta/gap = %.2f ' ...
         '(delta L = %.2f)\n'], 100*mean(E1(lo1)), 100*pk1, deltas(i1)/GAP1, deltas(i1)*L);
fprintf('fig7b: two donors, mean error %.1f%% below the gap, %.1f%% above, range [%.1f, %.1f]%%\n', ...
        100*mean(E2(lo2)), 100*mean(E2(~lo2)), 100*min(E2), 100*max(E2));
fprintf('fig7b: gaps %.3f and %.3f, ratio %.2f;  T_1/2(inf) one donor %.1f\n', ...
        GAP1, GAP2, GAP1/GAP2, Tinf);

% ---- panel --------------------------------------------------------------
c1 = P.greens(5,:);  cb = P.blues(3,:);
h  = P.figure(P);
ax = axes('Parent', h);  hold(ax, 'on');
x1 = deltas/GAP1;  x2 = d2/GAP2;
xl = [0.8*min([x1 x2]) 1.25*max([x1 x2])];
yl = [-5, 100*max([E1 + 1.96*SE1, E2 + 1.96*SE2]) + 3];
hb = patch(ax, [1 xl(2) xl(2) 1], [yl(1) yl(1) yl(2) yl(2)], [0.92 0.92 0.92], ...
           'EdgeColor', 'none');
yline(ax, 0, ':', 'Color', [0 0 0], 'LineWidth', P.ALW);
xline(ax, 1, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', P.ALW);
plot(ax, x1, 100*E1, '-', 'Color', P.greens(2,:), 'LineWidth', P.ALW);
plot(ax, x2, 100*E2, '-', 'Color', P.blues(1,:), 'LineWidth', P.ALW);
errorbar(ax, x1, 100*E1, 100*1.96*SE1, 'o', 'Color', c1, 'MarkerFaceColor', 'w', ...
         'LineWidth', P.ALW, 'MarkerSize', 12, 'CapSize', 10, 'LineStyle', 'none');
errorbar(ax, x2, 100*E2, 100*1.96*SE2, 's', 'Color', cb, 'MarkerFaceColor', 'w', ...
         'LineWidth', P.ALW, 'MarkerSize', 12, 'CapSize', 10, 'LineStyle', 'none');
uistack(hb, 'bottom');
set(ax, 'XScale', 'log', 'Layer', 'top');
xlim(ax, xl);  ylim(ax, yl);
set(ax, 'XTick', [0.1 1 10], 'XTickLabel', {'$0.1$','$1$','$10$'});
xlabel(ax, 'scaled abandonment rate $\delta/(\vartheta_2-\vartheta_1)$', ...
       'Interpreter','latex', 'FontSize', P.LSZ);
ylabel(ax, 'reduction error (\%)', 'Interpreter','latex', 'FontSize', P.LSZ);
P.axes(P, ax);
set(ax, 'Position', [0.17 0.16 0.78 0.74]);

P.label(P, ax, 1.55, 100*pk1, 'one donor', 'Color', c1, ...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom');
P.label(P, ax, x2(1), 100*max(E2(1:3) + 1.96*SE2(1:3)) + 0.8, ...
        {'two donors,', sprintf('$r=%g$', two.r)}, 'Color', cb, ...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom');
P.label(P, ax, 1.1, yl(1) + 0.8, '$\delta=\vartheta_2-\vartheta_1$', ...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom');

P.export(P, 'fig7b', outdir, h);
out = struct('deltas',deltas,'E1',E1,'SE1',SE1,'Tsim',Tsim,'Tred',Tred,'d2',d2,'E2',E2, ...
             'SE2',SE2,'abm2',mu2,'red2',Tr2,'gap1',GAP1,'gap2',GAP2,'nrep',nrep, ...
             'nrep2',two.nrep,'r',two.r);
end

% =========================================================================
% Panel (c), fig7c. Non-exponential service.
% =========================================================================
function out = panel_c(M, P, G, ps, cfg, outdir, deltas, nrep, cachefile)
%PANEL_C  Interior optimum under non-exponential service. Panel (c),
%   with the four service-time survival functions inset.
%
%   abm_svc runs with the service interval drawn from laws of equal mean 1/mu.
%       exponential      CV = 1      the control
%       deterministic    CV = 0
%       lognormal        CV = 1.99
%       Pareto, censored alpha = 1.21, censored at 179/7.9 mean durations
%   The ensemble is cached under cachefile with a stamped key.

if nargin < 7 || isempty(deltas), deltas = logspace(log10(0.05),log10(3),9); end
if nargin < 8 || isempty(nrep),   nrep   = 80; end

q    = 0.5;
laws = local_laws(cfg);
nl   = numel(laws);  nd = numel(deltas);

key = sprintf('vS3|n=%d|d=%s|laws=%s|p0=%g', nrep, mat2str(deltas,4), ...
              strjoin({laws.name},','), ps.p0);
T = [];
if isfile(cachefile)
    C = load(cachefile);
    if isfield(C,'key') && strcmp(C.key,key)
        T = C.T;  SE = C.SE;
        fprintf('fig7c: reusing %d cached realizations per (law, delta)\n', nrep);
    end
end
if isempty(T)
    T = nan(nl,nd);  SE = nan(nl,nd);
    for il = 1:nl
        for id = 1:nd
            r = arrayfun(@(k) abm_svc(M, deltas(id), G, ps, cfg, q, k, laws(il).draw), 1:nrep);
            r = r(isfinite(r));
            T(il,id) = mean(r);  SE(il,id) = std(r)/sqrt(numel(r));
        end
        fprintf('fig7c: %-18s done\n', laws(il).name);
    end
    save(cachefile, 'T', 'SE', 'key');
end

dbest = nan(1,nl);  Tbest = nan(1,nl);  SEbest = nan(1,nl);
for il = 1:nl
    [~, ib] = min(T(il,:));
    dbest(il) = deltas(ib);  Tbest(il) = T(il,ib);  SEbest(il) = SE(il,ib);
    fprintf('fig7c: %-18s CV = %.2f  optimum delta = %.3f, T_1/2 = %.1f +- %.1f (1.96 SE %.1f)\n', ...
            laws(il).name, laws(il).cv, deltas(ib), T(il,ib), SE(il,ib), 1.96*SE(il,ib));
end
Tinf = log(1/(1-q))/cfg.lam0;
z = abs(Tbest(1) - 90.2)/sqrt(SEbest(1)^2 + 0.7^2);
fprintf('fig7c: exponential control %.1f +- %.1f against abm 90.2 +- 0.7 (%.1f sigma)\n', ...
        Tbest(1), SEbest(1), z);
step  = exp(mean(diff(log(deltas))));
cells = log(max(dbest)/min(dbest))/log(step);
fprintf('fig7c: optima span a factor %.2f in delta, %.2f grid cells (one cell = %.0f%%)\n', ...
        max(dbest)/min(dbest), cells, 100*(step - 1));
fprintf('fig7c: speedup at each optimum over T_1/2(inf) = %.1f:%s\n', Tinf, ...
        sprintf(' %.0f%%', 100*(1 - Tbest/Tinf)));

% ---- panel --------------------------------------------------------------
cols = [P.blues(3,:); P.greens(6,:); P.greens(4,:); [0.55 0.55 0.55]];
mk   = {'o','s','^','d'};
h  = P.figure(P);
ax = axes('Parent', h);  hold(ax, 'on');
yline(ax, Tinf, ':', 'Color', [0 0 0], 'LineWidth', P.ALW);
for il = 1:nl
    plot(ax, deltas, T(il,:), '-', 'Color', cols(il,:), 'LineWidth', P.ALW);
end
for il = 1:nl
    errorbar(ax, deltas, T(il,:), 1.96*SE(il,:), mk{il}, 'Color', cols(il,:), ...
             'MarkerFaceColor', 'w', 'LineWidth', P.ALW, 'MarkerSize', 11, ...
             'CapSize', 8, 'LineStyle', 'none');
    plot(ax, dbest(il), Tbest(il), mk{il}, 'Color', cols(il,:), ...
         'MarkerFaceColor', cols(il,:), 'MarkerSize', 11, 'LineWidth', P.ALW);
end
set(ax, 'XScale', 'log');
xl = [0.04 4];  yl = [82 210];
xlim(ax, xl);  ylim(ax, yl);
set(ax, 'XTick', [0.1 1], 'YTick', 100:20:140);
xlabel(ax, 'abandonment rate $\delta$', 'Interpreter','latex', 'FontSize', P.LSZ);
ylabel(ax, 'half-colony feeding time $T_{1/2}$', 'Interpreter','latex', 'FontSize', P.LSZ);
P.axes(P, ax);
set(ax, 'Position', [0.17 0.16 0.78 0.74]);
P.label(P, ax, 0.045, Tinf + 1, 'no lingering', 'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', 'FontSize', P.FSZ-4);

% key, marker and name, between the curves and the no-lingering line
names = {'exponential', 'deterministic', 'lognormal', 'censored Pareto'};
for il = 1:nl
    yk = 111 + 7.5*(il - 1);
    plot(ax, 0.052, yk, mk{il}, 'Color', cols(il,:), 'MarkerFaceColor', 'w', ...
         'MarkerSize', 11, 'LineWidth', P.ALW);
    P.label(P, ax, 0.06, yk, names{il}, 'Color', cols(il,:), 'HorizontalAlignment','left', ...
            'VerticalAlignment','middle', 'FontSize', P.FSZ-4);
end

% inset, survival functions of the four laws on log axes, above the line
p  = get(ax, 'Position');
ai = axes('Parent', h, 'Position', [p(1)+0.20*p(3), p(2)+0.68*p(4), 0.50*p(3), 0.31*p(4)]);
hold(ai, 'on');
tt = logspace(-1, 1.5, 600);
for il = 1:nl
    [tx, Sx] = laws(il).surv(tt);
    plot(ai, tx, Sx, '-', 'Color', cols(il,:), 'LineWidth', P.LW-1);
end
set(ai, 'XScale','log', 'YScale','log', 'FontSize', P.FSZ-4, 'LineWidth', P.ALW-0.5, ...
    'TickDir','out', 'Box','off', 'TickLabelInterpreter','latex', 'Color','none', ...
    'XTick', [1e-1 1 1e1], 'YTick', [1e-3 1e-2 1e-1 1], 'XMinorTick','off', 'YMinorTick','off', ...
    'Layer', 'top');
xlim(ai, [1e-1 10^1.5]);  ylim(ai, [1e-3 1.5]);
xlabel(ai, 'service time $\times\,\mu$', 'Interpreter','latex', 'FontSize', P.FSZ-4);
ylabel(ai, 'survival', 'Interpreter','latex', 'FontSize', P.FSZ-4);

P.export(P, 'fig7c', outdir, h);
out = struct('deltas',deltas,'T',T,'SE',SE,'laws',{{laws.name}},'cv',[laws.cv], ...
             'dbest',dbest,'Tbest',Tbest,'SEbest',SEbest,'nrep',nrep,'q',q,'Tinf',Tinf);
end

% ------------------------------------------------------------------------
function laws = local_laws(cfg)
%LOCAL_LAWS  Service-interval samplers, each with its survival
%   function, every one rescaled to mean 1/mu.
m = 1/cfg.mu;

laws(1).name = 'exponential';
laws(1).cv   = 1.00;
laws(1).draw = @() -m*log(rand);
laws(1).surv = @(t) deal(t, exp(-t/m));

laws(2).name = 'deterministic';
laws(2).cv   = 0.00;
laws(2).draw = @() m;
laws(2).surv = @(t) deal([t(t < m), m, m], [ones(1, nnz(t < m)), 1, 1e-6]);

cv    = 1.99;
s2    = log(1 + cv^2);
mu_ln = log(m) - s2/2;
laws(3).name = 'lognormal';
laws(3).cv   = cv;
laws(3).draw = @() exp(mu_ln + sqrt(s2)*randn);
laws(3).surv = @(t) deal(t, 0.5*erfc((log(t) - mu_ln)/sqrt(2*s2)));

al = 1.21;
b  = 179/7.9;
mp = integral(@(x) x.*(al*x.^(-al-1))./(1 - b^(-al)), 1, b);
laws(4).name = 'Pareto (censored)';
laws(4).cv   = local_tp_cv(al, b);
laws(4).draw = @() (m/mp)*local_tp(al, b);
laws(4).surv = @(t) local_tp_surv(t, al, b, m/mp);
end

function x = local_tp(al, b)
u = rand;
x = (1 - u*(1 - b^(-al)))^(-1/al);
end

function c = local_tp_cv(al, b)
Z  = 1 - b^(-al);
m1 = integral(@(x) x   .*(al*x.^(-al-1))/Z, 1, b);
m2 = integral(@(x) x.^2.*(al*x.^(-al-1))/Z, 1, b);
c  = sqrt(max(m2 - m1^2,0))/m1;
end

function [t, S] = local_tp_surv(t, al, b, sc)
%LOCAL_TP_SURV  Survival of sc*Y with Y Pareto(alpha) on [1,b].
y = t/sc;
S = ones(size(y));
in = y >= 1 & y < b;
S(in) = (y(in).^(-al) - b^(-al))/(1 - b^(-al));
S(y >= b) = 1e-6;
k = find(y < b, 1, 'last');
t = [t(1:k), sc*b, sc*b];  S = [S(1:k), (min(y(k),b)^(-al) - b^(-al))/(1 - b^(-al)), 1e-6];
end
