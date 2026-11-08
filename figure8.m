function out = figure8(varargin)
%FIGURE6 Reproduce the three panels of Figure 8 (the averaged recruitment gain)
%   and assemble the lettered figure.
%
% Usage
%   figure8
%   figure8('fast',true)                  % coarse grids, wiring check
%   figure8('outdir','output/figure8')
%
% Panel map
%   (a) fig6b  first-scenter value g_1 against p_0 with its two-sided bound
%   (b) fig6c  gain elasticity with its small- and large-Fbar asymptotes
%   (c) fig6d  variance ratio of the aggregate plume peak against scenter count
%   Panel files keep their stems. The averaged-gain surface that was a fourth
%   panel duplicates the slices drawn in Figure 3(a) of the article and is retired.
%
% The gain surface is extended past the cue support to p = 3, so the Gaussian
% average at large p_0 samples the flat far field before offsets are clamped
% to max(W.p). The representative-parameter tables on [0, 1.2] come
% from the representative-parameter cache in outdir. The extension and the
% kbar = 0.8 surface are cached here behind a stamp.
%
% Every a_k is the Gaussian smoothing of W at spread p_0/sqrt(k), on the
% quadrature grid psiTable uses. The averaged gain and its derivative come from
% mixGain, which takes the derivative by the envelope theorem.
%
% Panel (d) draws on its own peak ensemble, cached in outdir under a stamped
% key. A fast run writes its panels and caches to outdir/fast so a production
% cache is never replaced. The lettered figure is written to output/figs/fig8 (Fig. 8 of the article).

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'src')));

ip = inputParser;
addParameter(ip,'fast',false);
addParameter(ip,'outdir',fullfile('output','figure8'));
addParameter(ip,'figdir',fullfile('output','figs'));
parse(ip,varargin{:});
fast   = logical(ip.Results.fast);
OUTDIR = ip.Results.outdir;
FIGDIR = ip.Results.figdir;
if fast, OUTDIR = fullfile(OUTDIR, 'fast');  FIGDIR = OUTDIR; end
if ~isfolder(OUTDIR), mkdir(OUTDIR); end
if ~isfolder(FIGDIR), mkdir(FIGDIR); end
close all

VER    = 'v7-averaged-gain';
P0LIST = [0.75 1.13 1.52];
KBARS  = [0.0314 0.8];

A = afpt();
M = lingermodel();
P = linger_style();
prm = A.mkprm('kbar',0.0314,'sig',0.8,'beta0',1,'sm',1,'v',1);
cfg = M.mkcfg('N',120,'mu',1,'lam0',prm.lam0,'p0',0.75);
fprintf('figure8: lambda0 = %.5f, rho = %.3f, gmax = %.4f\n', ...
        prm.lam0, prm.lam0*cfg.N/cfg.mu, prm.gmax);

% Grids are part of both cache keys, so fix them before stamping.
if fast
    pbase = linspace(0,1.2,31);  pext = 1.4:0.2:3.0;  ggrid = linspace(0,prm.gmax,25);
    Fb  = logspace(-1, log10(50), 60);
    p0B = logspace(log10(0.3), log10(5), 16);
else
    pbase = linspace(0,1.2,61);  pext = 1.25:0.05:3.0; ggrid = linspace(0,prm.gmax,51);
    Fb  = logspace(-1, log10(50), 160);
    p0B = logspace(log10(0.3), log10(5), 48);
end
pgrid = [pbase pext];

% ---- representative-parameter tables on [0, 1.2] ---------------------------------
opfile = fullfile(OUTDIR,'opgain_cache.mat');
keyS = M.stamp('vS1-si-panels', prm, cfg, pbase, ggrid);
G = [];  Wb = [];
if isfile(opfile)
    S = load(opfile);
    if isfield(S,'key') && strcmp(S.key,keyS)
        G = S.G;  Wb = S.W;
        fprintf('figure8: loaded representative-parameter tables from cache\n');
    end
end
if isempty(G)
    fprintf('figure8: building representative-parameter tables on %d x %d\n', numel(pbase), numel(ggrid));
    G  = M.gainTable(prm, pbase);
    Wb = M.gainSurface(prm, pbase, ggrid);
end

% ---- extended surfaces, both absorptivities -----------------------------
cachefile = fullfile(OUTDIR,'figure8_gain_cache.mat');
key = M.stamp(VER, prm, cfg, pbase, pext, ggrid, KBARS);
T = [];
if isfile(cachefile)
    S = load(cachefile);
    if isfield(S,'key') && strcmp(S.key,key)
        T = S.T;
        fprintf('figure8: loaded stamped surface cache (%d x %d)\n', numel(pgrid), numel(ggrid));
    end
end
if isempty(T)
    T = struct('kbar',KBARS,'W',{cell(1,2)},'G',{cell(1,2)},'D0',{cell(1,2)}, ...
               'pgrid',pgrid,'ggrid',ggrid);
    dg = prm.gmax/500;
    for ik = 1:numel(KBARS)
        pk = A.mkprm('kbar',KBARS(ik),'sig',prm.sig,'beta0',prm.beta0,'sm',prm.sm,'v',prm.v);
        fprintf('figure8: kbar = %.4g, surface on %d x %d\n', KBARS(ik), numel(pgrid), numel(ggrid));
        if ik == 1
            Gk = G;  Wk = Wb;
        else
            Gk = M.gainTable(pk, pbase);
            Wk = M.gainSurface(pk, pbase, ggrid);
        end
        We = M.gainSurface(pk, pext, ggrid);
        if abs(We.T0 - Wk.T0) > 1e-8*Wk.T0
            warning('figure8:T0mismatch', 'unaided MFPT differs between surface pieces');
        end
        Wk.p = pgrid;  Wk.G = [Wk.G We.G];
        % dW/dgamma at gamma = 0, one-sided second order.
        Wd = M.gainSurface(pk, pgrid, [0 dg 2*dg]);
        T.D0{ik} = (-3*Wd.G(1,:) + 4*Wd.G(2,:) - Wd.G(3,:))/(2*dg);
        T.W{ik} = Wk;  T.G{ik} = Gk;
    end
    % One-sided slope of W at the ceiling at the peak. W has a cusp at p = 0.
    h1 = 1e-4;
    T.c1 = (T.W{1}.T0/A.moments(0, prm.sig, prm.beta0, prm.gmax, prm.sm, prm.theta, prm.v) - ...
            T.W{1}.T0/A.moments(h1, prm.sig, prm.beta0, prm.gmax, prm.sm, prm.theta, prm.v))/h1;
    save(cachefile,'T','key');
end

W = T.W{1};
psList = arrayfun(@(p0) M.psiTable(G, p0, W), P0LIST, 'uni', 0);
G0 = max(W.G(:,1));
fprintf('figure8 checks: p*/sigma = %.3f, G(0) = %.3f, g1 =', G.pstar/prm.sig, G0);
fprintf(' %.3f', cellfun(@(s) s.g1, psList));
fprintf(', F_sub ='); fprintf(' %.2f', cellfun(@(s) s.Fsub, psList)); fprintf('\n');
if abs(G0 - 3.434) > 5e-3
    warning('figure8:G0Moved', 'G(0) = %.4f differs from the reference 3.434', G0);
end
if any(abs(cellfun(@(s) s.g1, psList) - [1.552 1.298 1.159]) > 8e-3)
    warning('figure8:g1Moved', 'g_1 differs from the reference 1.552 / 1.298 / 1.159');
end
if any(abs(cellfun(@(s) s.Fsub, psList) - [1.14 2.59 4.68]) > 1e-2)
    warning('figure8:FsubMoved', 'F_sub differs from the reference 1.14 / 2.59 / 4.68');
end

Q = struct('ak', @akTable, 'akv', @akTableV, 'mix', @mixGain, 'G0', G0, 'p0list', P0LIST, ...
           'prm', prm, 'psList', {psList});

out = struct('T',T,'G',G,'psList',{psList});
out.b = panel_a(Q, T, P, OUTDIR, p0B);
if fast
    out.d = panel_c(prm, cfg.p0, OUTDIR, [1 2 3 5 10 20 30], 4000, 5e-3, 40000);
else
    out.d = panel_c(prm, cfg.p0, OUTDIR);
end
Q.vk = struct('k', out.d.k, 'ratio', out.d.ratio, 'vinf', out.d.vinf);   % measured k Var(peak)/p0^2 of the summed tents
out.c = panel_b(Q, T, P, OUTDIR, Fb);
save(fullfile(OUTDIR,'figure8_results.mat'),'-struct','out');

% One row of three at the panel canvas size, 66 cm, included at full textwidth
% like the other wide figures.
names  = {'fig6b','fig6c','fig6d'};
panels = cellfun(@(n) fullfile(OUTDIR,[n '.fig']), names, 'uni', 0);
tiles  = [0 0 22 19; 22 0 22 19; 44 0 22 19];
assemble_figure(fullfile(FIGDIR,'fig8'), panels, tiles, [66 19]);
fprintf('figure8 complete. Assembled figure written to %s\n', fullfile(FIGDIR,'fig8'));
end

% ------------------------------------------------------------------------
function rk = akTable(p, Gm, p0, K)
%AKTABLE  a_k for k = 1..K, the smoothing of each row of Gm (rows gamma,
%   columns p) by N(0, p0^2/k), offsets clamped to max(p) as in psiTable.
z  = linspace(-6, 6, 1201).';
w  = exp(-z.^2/2);  w = w/trapz(z, w);
Gt = Gm.';
pmax = max(p);
rk = zeros(size(Gm,1), K);
for k = 1:K
    pp = min(abs(z*p0/sqrt(k)), pmax);
    rk(:,k) = trapz(z, w.*interp1(p(:), Gt, pp), 1).';
end
end

% ------------------------------------------------------------------------
function rk = akTableV(p, Gm, p0, K, kk, rr, vinf)
%AKTABLEV  As akTable, but with the peak variance of a literal sum of tents.
%   The scaled variance rr at the simulated counts kk is interpolated in k and
%   held at its last value beyond them, so only F with Poisson mass
%   mostly below k = kk(end) should be read, so the peak spread is p0 sqrt(r(k)/k).
z  = linspace(-6, 6, 1201).';
w  = exp(-z.^2/2);  w = w/trapz(z, w);
Gt = Gm.';
pmax = max(p);
r = interp1(kk, rr, 1:K, 'pchip', rr(end));   % measured counts only; held beyond k = kk(end)
rk = zeros(size(Gm,1), K);
for k = 1:K
    pp = min(abs(z*p0*sqrt(r(k)/k)), pmax);
    rk(:,k) = trapz(z, w.*interp1(p(:), Gt, pp), 1).';
end
end

% ------------------------------------------------------------------------
function [vals, jsel, dG] = mixGain(rk, F)
%MIXGAIN  Gtilde(F) = max_gamma [e^-F + sum_k pi_k a_k], the selected row, and
%   Gtilde'(F) = sum_{k>=0} pi_k (a_{k+1} - a_k) at that row, a_0 = 1. At F = 0
%   the maximizer is the whole interval and the right derivative is g_1 - 1.
K  = size(rk,2);  ks = 1:K;
inc = [rk(:,1) - 1, diff(rk,1,2)];        % a_{k+1} - a_k for k = 0..K-1
vals = ones(size(F));  jsel = ones(size(F));  dG = zeros(size(F));
for i = 1:numel(F)
    f = F(i);
    if f <= 0
        [g1, j] = max(rk(:,1));  jsel(i) = j;  dG(i) = g1 - 1;
        continue
    end
    pk = exp(-f + ks*log(f) - gammaln(ks+1));
    [vals(i), j] = max(exp(-f) + rk*pk(:));
    jsel(i) = j;
    dG(i) = [exp(-f), pk(1:K-1)]*inc(j,:).';
end
end

% =========================================================================
% Panel (a), fig6b. Value of the first scenter.
% =========================================================================
function out = panel_a(Q, T, P, outdir, p0s)
%PANEL_A  Value of the first scenter, g_1, against the single-scenter
%   offset p_0, with the two-sided bound. Panel (a),
%   with in-panel labels.
%
%       g_1(p_0) = max_gamma  Int W(gamma,|p|) N(p; 0, p_0^2) dp
%
%       1 - (1 - Wmin)(1 - Phi*)  <=  a_1(gamma)  <=  1 + (G(0) - 1) Phi*,
%       Phi* = 2 Phi(p*/p_0) - 1
%
%   Solid is kbar = 0.0314, dashed kbar = 0.8. The crossing where g_1 returns
%   to one is bisected on the reliance derivative of a_1 at gamma = 0,
%   s(p_0) = Int dW/dgamma|_0 N(p; 0, p_0^2) dp.

nk = numel(T.kbar);  npo = numel(p0s);
g1 = nan(nk,npo);  gstar = nan(nk,npo);  ub = g1;  lb = g1;
p0crit = nan(nk,1);  p0need = nan(nk,1);  pstars = nan(nk,1);
G0s = nan(nk,1);  Wmins = nan(nk,1);  g1ref = nan(nk,3);

z = linspace(-6,6,1201);
w = exp(-z.^2/2);  w = w/trapz(z,w);
for ik = 1:nk
    Wk = T.W{ik};  pmax = max(Wk.p);
    Wmins(ik)  = min(Wk.G(:));
    G0s(ik)    = max(Wk.G(:,1));
    pstars(ik) = T.G{ik}.pstar;
    for j = 1:npo
        a1 = Q.ak(Wk.p, Wk.G, p0s(j), 1);
        [g1(ik,j), imax] = max(a1);
        gstar(ik,j) = Wk.gam(imax);
        Phi = erf(pstars(ik)/p0s(j)/sqrt(2));
        ub(ik,j) = 1 + (G0s(ik) - 1)*Phi;
        lb(ik,j) = 1 - (1 - Wmins(ik))*(1 - Phi);
    end
    g1ref(ik,:) = arrayfun(@(q) max(Q.ak(Wk.p, Wk.G, q, 1)), Q.p0list);

    sfun = @(qq) trapz(z, w.*interp1(Wk.p, T.D0{ik}, min(abs(z*qq), pmax), 'linear'));
    plo = min(p0s);  phi = max(p0s);
    if sfun(plo) > 0 && sfun(phi) <= 0
        for it = 1:60
            pm = 0.5*(plo + phi);
            if sfun(pm) > 0, plo = pm; else, phi = pm; end
            if phi - plo < 1e-8*phi, break, end
        end
        p0crit(ik) = 0.5*(plo + phi);
    else
        warning('fig6b:noCrossing', 'at kbar = %.4g the crossing is not bracketed by the sweep', T.kbar(ik));
    end
    Phineed = (1 - Wmins(ik))/(G0s(ik) - Wmins(ik));
    p0need(ik) = pstars(ik)/max(sqrt(2)*erfinv(Phineed), eps);
    fprintf(['fig6b: kbar = %.4g, G(0) = %.3f, Wmin = %.3f, p* = %.4f, ' ...
             'g_1 returns to 1 at p_0 = %.3f (%.2f sigma), necessary bound %.2f\n'], ...
            T.kbar(ik), G0s(ik), Wmins(ik), pstars(ik), p0crit(ik), ...
            p0crit(ik)/Q.prm.sig, p0need(ik));
    fprintf('fig6b: kbar = %.4g, g_1 at p0 = 0.75, 1.13, 1.52: %.4f %.4f %.4f\n', ...
            T.kbar(ik), g1ref(ik,:));
end
if any(g1(:) > ub(:) + 1e-6) || any(g1(:) < lb(:) - 1e-6)
    warning('fig6b:outsideBracket', 'g_1 fell outside the two-sided bound somewhere');
end

% ---- panel --------------------------------------------------------------
h  = P.figure(P);
ax = axes('Parent',h);  hold(ax,'on');
cols = [P.greens(6,:); P.blues(3,:)];
for ik = 1:nk
    fill(ax, [p0s fliplr(p0s)], [lb(ik,:) fliplr(ub(ik,:))], cols(ik,:), ...
         'FaceAlpha', 0.09, 'EdgeColor','none');
end
yline(ax, 1, '-', 'Color', [0.82 0.84 0.86], 'LineWidth', P.ALW);
for ik = 1:nk
    plot(ax, p0s, g1(ik,:), P.p0sty{ik}, 'Color', cols(ik,:), 'LineWidth', P.LW);
    if ~isnan(p0crit(ik))
        plot(ax, p0crit(ik), 1, 'o', 'Color', cols(ik,:), 'MarkerFaceColor','w', ...
             'MarkerSize', 13, 'LineWidth', P.ALW);
    end
end
g1op = g1ref(1,1);
plot(ax, 0.75, g1op, 'v', 'Color', cols(1,:), 'MarkerFaceColor', cols(1,:), 'MarkerSize', 15);
for i = 2:3
    plot(ax, Q.p0list(i), g1ref(1,i), 'o', 'MarkerSize', 7, ...
         'MarkerFaceColor', P.p0cols(i,:), 'MarkerEdgeColor','none');
end
set(ax,'XScale','log');
xlim(ax, [min(p0s) max(p0s)]);
ylim(ax, [0.5 max(ub(:))*1.02]);
set(ax, 'XTick', [0.3 0.5 1 2 5]);

fs = P.FSZ - 2;
text(ax, 0.8, 1.78, 'representative', 'Interpreter','latex', 'FontSize', fs, ...
     'Color', cols(1,:), 'HorizontalAlignment','left');
text(ax, 0.32, 3.12, 'two-sided bound', 'Interpreter','latex', 'FontSize', fs, ...
     'Color', P.grey, 'HorizontalAlignment','left');
text(ax, 0.32, 2.52, '$\overline\kappa=0.0314$', 'Interpreter','latex', ...
     'FontSize', fs, 'Color', cols(1,:));
text(ax, 0.32, 1.45, '$\overline\kappa=0.8$', 'Interpreter','latex', ...
     'FontSize', fs, 'Color', cols(2,:));
text(ax, p0crit(1), 0.9, {'$g_1$ returns', 'to one'}, 'Interpreter','latex', ...
     'FontSize', fs, 'Color', 'k', 'HorizontalAlignment','center', 'VerticalAlignment','top');

xlabel(ax, 'single-scenter offset $p_0$', 'Interpreter','latex', 'FontSize', P.LSZ);
ylabel(ax, 'first-scenter value $g_1$', 'Interpreter','latex', 'FontSize', P.LSZ);
P.axes(P, ax);
set(ax, 'Position', [0.17 0.16 0.77 0.72]);
P.export(P, 'fig6b', outdir, h);

out = struct('p0s',p0s,'g1',g1,'gstar',gstar,'ub',ub,'lb',lb,'p0crit',p0crit, ...
             'p0need',p0need,'pstar',pstars,'G0',G0s,'Wmin',Wmins,'g1ref',g1ref);
end

% =========================================================================
% Panel (b), fig6c. Gain elasticity and its limits.
% =========================================================================
function out = panel_b(Q, T, P, outdir, Fb)
%PANEL_B  Gain elasticity e(Fbar) = Fbar Gtilde'/Gtilde for the three
%   offsets, with the small- and large-Fbar asymptotes at p_0 = 0.75.
%
%   Gtilde' is taken by the envelope theorem from the Poisson increments, not
%   by differencing the table.
%
%   Small Fbar: e ~ (g_1 - 1) Fbar.
%   Large Fbar: W has a cusp at the peak, W(gamma_0, p) ~ G(0) - c_1 |p|, with
%   c_1 the one-sided slope at the ceiling reliance. Averaging |p| over
%   N(0, p_0^2/k) and k ~ Poisson(Fbar) gives
%       Gtilde ~ G(0) - C/sqrt(Fbar),   C = c_1 p_0 sqrt(2/pi),
%       e      ~ C/(2 G(0) sqrt(Fbar)).
%   C is computed from c_1 with no fit and checked against the tail of the
%   computed Gtilde. The curvature form C/(G(0) Fbar) needs a smooth peak,
%   which W does not have.

W  = T.W{1};  G0 = Q.G0;
Fl = logspace(log10(min(Fb)), log10(400), 240);    % to read the tail
K  = ceil(max(Fl) + 12*sqrt(max(Fl)) + 30);
np = numel(Q.p0list);
e  = nan(np, numel(Fl));  Gt = e;  gsel = e;
epk = nan(1,np);  Fpk = nan(1,np);  g1 = nan(1,np);
for i = 1:np
    rk = Q.ak(W.p, W.G, Q.p0list(i), K);
    [Gt(i,:), j, dG] = Q.mix(rk, Fl);
    gsel(i,:) = W.gam(j);
    e(i,:) = Fl.*dG./Gt(i,:);
    g1(i) = max(rk(:,1));
    [~, kk] = max(e(i,:));
    % Parabolic refinement in log Fbar.
    if kk > 1 && kk < numel(Fl)
        u = log(Fl(kk-1:kk+1));  c = polyfit(u, e(i,kk-1:kk+1), 2);
        Fpk(i) = exp(-c(2)/(2*c(1)));  epk(i) = polyval(c, log(Fpk(i)));
    else
        Fpk(i) = Fl(kk);  epk(i) = e(i,kk);
    end
end
if any(abs(gsel(:,2:end) - max(W.gam)) > 1e-12, 'all')
    fprintf('fig6c: the maximizing reliance leaves the ceiling somewhere on the sweep\n');
else
    fprintf('fig6c: the maximizing reliance sits on the ceiling for every Fbar > 0 and all three offsets\n');
end

% Cross-check against the table derivative the reduced model uses.
M = lingermodel();
edl = Fb.*M.dlogPsi(Q.psList{1}, Fb);
eme = interp1(Fl, e(1,:), Fb);
fprintf('fig6c: max |e_envelope - e_dlogPsi| at p0 = 0.75 on Fbar in [0.1, 25]: %.4f\n', ...
        max(abs(eme(Fb <= 25) - edl(Fb <= 25))));

C  = T.c1*Q.p0list(1)*sqrt(2/pi);
Ct = sqrt(Fl).*(G0 - Gt(1,:));
fprintf('fig6c: one-sided slope of W at the peak c1 = %.3f, so C = c1 p0 sqrt(2/pi) = %.3f\n', T.c1, C);
fprintf('fig6c: tail check sqrt(F)(G0 - Gtilde) at F = 25, 100, 400: %.3f %.3f %.3f\n', ...
        interp1(Fl, Ct, [25 100 400]));
fprintf('fig6c: F(G0 - Gtilde) at F = 25, 100, 400: %.2f %.2f %.2f (no 1/F limit)\n', ...
        interp1(Fl, Fl.*(G0 - Gt(1,:)), [25 100 400]));
sl = polyfit(log(Fl(Fl >= 100)), log(G0 - Gt(1, Fl >= 100)), 1);
fprintf('fig6c: log-log slope of G0 - Gtilde over F in [100, 400]: %.3f\n', sl(1));
Fpred = (C/(2*G0*(g1(1) - 1)))^(2/3);
for i = 1:np
    fprintf('fig6c: p0 = %.2f, elasticity peak e = %.3f at Fbar = %.2f\n', Q.p0list(i), epk(i), Fpk(i));
end
fprintf('fig6c: asymptote crossing predicts the p0 = 0.75 peak near Fbar = %.2f\n', Fpred);

% ---- panel --------------------------------------------------------------
h  = P.figure(P);
ax = axes('Parent',h);  hold(ax,'on');
show = Fl <= max(Fb);
Fs = Fl(show);
aC = P.grey;
% Literal sum of tents at the reference offset, with the measured peak variance.
rkv = Q.akv(W.p, W.G, Q.p0list(1), K, Q.vk.k, Q.vk.ratio, Q.vk.vinf);
[Gv, jv, dGv] = Q.mix(rkv, Fl);
ev = Fl.*dGv./Gv;
[~, kv] = max(ev);
u = log(Fl(kv-1:kv+1));  cv = polyfit(u, ev(kv-1:kv+1), 2);
Fpkv = exp(-cv(2)/(2*cv(1)));  epkv = polyval(cv, log(Fpkv));
fprintf('fig6c: summed-tent peak variance, elasticity peak e = %.3f at Fbar = %.2f\n', epkv, Fpkv);
plot(ax, Fs, (g1(1) - 1)*Fs, ':', 'Color', aC, 'LineWidth', P.ALW);
plot(ax, Fs, C./(2*G0*sqrt(Fs)), ':', 'Color', aC, 'LineWidth', P.ALW);
for i = np:-1:1
    plot(ax, Fs, e(i,show), P.p0sty{i}, 'Color', P.p0cols(i,:), 'LineWidth', P.LW);
end
vshow = show & (Fl <= 12);   % measured counts reach k = 30, so read F <~ 12
plot(ax, Fl(vshow), ev(vshow), '--', 'Color', [0.1 0.1 0.1], 'LineWidth', P.LW-1.5);
plot(ax, Fpkv, epkv, 'o', 'MarkerSize', 11, 'LineWidth', P.ALW, ...
     'MarkerFaceColor','w', 'MarkerEdgeColor', [0.1 0.1 0.1]);
for i = 1:np
    plot(ax, Fpk(i), epk(i), 'o', 'MarkerSize', 11, 'LineWidth', P.ALW, ...
         'MarkerFaceColor','w', 'MarkerEdgeColor', P.p0cols(i,:));
end
set(ax, 'XScale','log');
xlim(ax, [min(Fb) max(Fb)]);  ylim(ax, [0 0.4]);
set(ax, 'XTick', [0.1 1 10]);

fs = P.FSZ - 2;
% The peak value is convention dependent, so it is marked but not printed;
% the caption gives it.
text(ax, 0.115, 0.30, '$(g_1-1)\bar F$', 'Interpreter','latex', 'FontSize', fs, 'Color', aC);
text(ax, 5.0, 0.288, '$C/2G(0)\sqrt{\bar F}$', 'Interpreter','latex', ...
     'FontSize', fs, 'Color', aC, 'HorizontalAlignment','left');
swatches(P, ax, 6, 9, 9.8, [0.39 0.362 0.334], Q.p0list);
plot(ax, [6 9], 0.035*[1 1], '--', 'Color', [0.1 0.1 0.1], 'LineWidth', P.LW-1.5);
text(ax, 9.8, 0.035, 'summed tents', 'Interpreter','latex', 'FontSize', P.FSZ-4, ...
     'Color', [0.1 0.1 0.1], 'VerticalAlignment','middle');

xlabel(ax, 'scenting population $\bar F$', 'Interpreter','latex', 'FontSize', P.LSZ);
ylabel(ax, 'gain elasticity $e(\bar F)$', 'Interpreter','latex', 'FontSize', P.LSZ);
P.axes(P, ax);
set(ax, 'Position', [0.17 0.16 0.77 0.72]);
P.export(P, 'fig6c', outdir, h);

out = struct('epkv',epkv,'Fpkv',Fpkv,'F',Fl,'e',e,'Gt',Gt,'epk',epk,'Fpk',Fpk,'g1',g1,'C',C,'c1',T.c1, ...
             'tailslope',sl(1),'Fpred',Fpred);
end

% ------------------------------------------------------------------------
function swatches(P, ax, x0, x1, xt, ys, p0list)
%SWATCHES  Hand-drawn line-style key for the p_0 family, as axes objects.
for i = 1:numel(ys)
    plot(ax, [x0 x1], ys(i)*[1 1], P.p0sty{i}, 'Color', P.p0cols(i,:), 'LineWidth', P.LW-0.5);
    text(ax, xt, ys(i), sprintf('$p_0=%.2f$', p0list(i)), 'Interpreter','latex', ...
         'FontSize', P.FSZ-4, 'Color', P.p0cols(i,:), 'VerticalAlignment','middle');
end
end

% =========================================================================
% Panel (c), fig6d. Peak variance against scenter count.
% =========================================================================
function out = panel_c(prm, p0, outdir, klist, nmc, dx, nmed, seed, mseed)
%PANEL_C  Variance ratio k Var(peak)/p0^2 of the aggregate plume peak
%   against scenter count, against the sample-median values v_k, the model
%   value 1, and the two limits of the tent peak. Panel (c).
%
%   Peak ensemble. Scenter offsets are xi_i ~ N(0, p0^2), each casting a unit
%   tent of half-width sigma, summed on the line. The peak is the argmax of
%   the summed cue on a uniform grid of spacing DX. On a plateau of the
%   aggregate, where the maximum is attained on an interval, the peak is the
%   midpoint of the run. The grid peak is checked against the exact maximum,
%   which sits at one of the xi_i since the aggregate is piecewise linear with
%   downward kinks only at the scenters. The sample median and the sample mean
%   of the same offsets are kept for comparison. The ensemble is cached in
%   OUTDIR under a stamped key.
%
%   The Monte Carlo band is +-1.96 standard errors of the sample variance,
%   using the sample fourth moment.
%
%   For many scenters the peak is the zero of the aggregate slope, whose mean
%   is the tent-smoothed offset density. Linearizing about the origin with
%   z = sigma/p0 gives
%
%       k Var(peak)/p0^2 -> pi erf(z/sqrt 2) / (2 (1 - exp(-z^2/2))^2),
%
%   which is pi/2, the sample median, when sigma >> p0 and grows as the
%   cluster widens against the tent.
%
%   v_k for the sample median is estimated from NMED draws per k.

if nargin < 3 || isempty(outdir), outdir = 'output'; end
if nargin < 4 || isempty(klist),  klist = [1 2 3 4 5 6 8 10 13 16 20 25 30]; end
if nargin < 5 || isempty(nmc),    nmc = 40000; end
if nargin < 6 || isempty(dx),     dx = 2e-3; end
if nargin < 7 || isempty(nmed),   nmed = 400000; end
if nargin < 8 || isempty(seed),   seed = 17; end
if nargin < 9 || isempty(mseed),  mseed = 29; end

M = lingermodel();  P = linger_style();
sig = prm.sig;

% ---- peak ensemble ------------------------------------------------------
cfile = fullfile(outdir, 'fig6d_cache.mat');
key = M.stamp('fig6d-grid-argmax', sig, p0, klist, nmc, dx, seed);
D = [];
if isfile(cfile)
    C = load(cfile);
    if isfield(C, 'key') && strcmp(C.key, key), D = C.R; fprintf('fig6d: loaded stamped peak ensemble\n'); end
end
if isempty(D)
    nk = numel(klist);
    pk = zeros(nmc, nk);  pe = zeros(nmc, nk);  md = zeros(nmc, nk);  mn = zeros(nmc, nk);
    narrow = zeros(1, nk);
    for a = 1:nk
        k = klist(a);
        rng(seed + k);
        X = p0*randn(nmc, k);
        md(:,a) = median(X, 2);  mn(:,a) = mean(X, 2);
        narrow(a) = mean(max(X, [], 2) - min(X, [], 2) < sig);
        pk(:,a) = gridPeak(X, sig, dx);
        pe(:,a) = exactPeak(X, sig);
        fprintf('fig6d: k = %d done\n', k);
    end
    R = struct('k',klist,'peak',pk,'exact',pe,'median',md,'mean',mn, ...
               'narrow',narrow,'n',nmc,'dx',dx,'p0',p0,'sig',sig);
    save(cfile, 'R', 'key');
    D = R;
end

agree = mean(abs(D.peak - D.exact) <= D.dx);
atmed = mean(abs(D.peak - D.median) <= D.dx);

% ---- variance ratios ----------------------------------------------------
k = D.k;  n = D.n;
z = sig/p0;
vinf = pi*erf(z/sqrt(2))/(2*(1 - exp(-z^2/2))^2);

vr = zeros(size(k));  vse = zeros(size(k));  vmn = zeros(size(k));
for a = 1:numel(k)
    x = D.peak(:,a) - mean(D.peak(:,a));
    s2 = mean(x.^2)*n/(n-1);  m4 = mean(x.^4);
    vr(a)  = k(a)*s2/p0^2;
    vse(a) = k(a)*sqrt(max(m4 - s2^2*(n-3)/(n-1), 0)/n)/p0^2;
    vmn(a) = k(a)*var(D.mean(:,a))/p0^2;
end
% v_k from a larger median ensemble.
vk = zeros(size(k));
for a = 1:numel(k)
    rng(mseed + k(a));
    vk(a) = k(a)*var(median(randn(nmed, k(a)), 2));
end

for a = 1:numel(k)
    fprintf(['fig6d: k = %2d  tent peak %.3f +- %.3f   sample median v_k %.3f   sample mean %.3f\n' ...
             '              grid = exact peak on %.3f of draws, peak = sample median on %.3f, ' ...
             'cluster narrower than sigma on %.3f\n'], ...
            k(a), vr(a), 1.96*vse(a), vk(a), vmn(a), agree(a), atmed(a), D.narrow(a));
end
fprintf('fig6d: large-k limit of the tent peak %.3f at sigma/p0 = %.3f; median limit pi/2 = %.3f\n', ...
        vinf, z, pi/2);

% ---- panel --------------------------------------------------------------
cM = P.blues(3,:);  cV = P.greens(4,:);
h = P.figure(P);  ax = axes('Parent', h);  hold(ax, 'on');
kl = [0.8 40];
fill(ax, [k fliplr(k)], [vr + 1.96*vse, fliplr(vr - 1.96*vse)], P.blues(1,:), ...
     'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot(ax, kl, vinf*[1 1], '--', 'Color', P.grey, 'LineWidth', P.ALW);
plot(ax, kl, pi/2*[1 1], ':', 'Color', P.grey, 'LineWidth', P.ALW);
plot(ax, kl, [1 1], '-', 'Color', P.grey, 'LineWidth', P.ALW);
plot(ax, k, vk, '-', 'Color', cV, 'LineWidth', P.LW);
plot(ax, k, vk, 'o', 'Color', cV, 'MarkerFaceColor', cV, 'MarkerSize', 7);
plot(ax, k, vr, '-', 'Color', cM, 'LineWidth', P.LW);
plot(ax, k, vr, 'o', 'Color', cM, 'MarkerFaceColor', 'w', 'MarkerSize', 12, 'LineWidth', P.ALW);
set(ax, 'XScale', 'log');
xlim(ax, kl);  ylim(ax, [0 1.12*vinf]);
set(ax, 'XTick', [1 3 10 30]);
xlabel(ax, 'scenters $k$', 'Interpreter','latex', 'FontSize', P.LSZ);
ylabel(ax, 'variance ratio $k\,\mathrm{Var}(p)/p_0^2$', 'Interpreter','latex', 'FontSize', P.LSZ);
P.axes(P, ax);
set(ax, 'Position', [0.17 0.17 0.76 0.72]);
P.label(P, ax, 0.9, vinf, 'large-$k$ limit', 'Color', P.grey, 'VerticalAlignment','bottom');
P.label(P, ax, 0.9, pi/2, '$\pi/2$', 'Color', P.grey, 'VerticalAlignment','bottom');
P.label(P, ax, 38, 1, 'model, sample mean', 'Color', P.grey, 'HorizontalAlignment','right', 'VerticalAlignment','top');
P.label(P, ax, 6, interp1(k, vr, 6) + 0.25, 'tent peak', 'Color', cM, 'HorizontalAlignment','right', 'VerticalAlignment','bottom');
P.label(P, ax, 38, pi/2 + 0.05, 'sample median', 'Color', cV, 'HorizontalAlignment','right', 'VerticalAlignment','bottom');

P.export(P, 'fig6d', outdir, h);
close(h);

out = struct('k',k,'ratio',vr,'se',vse,'vk',vk,'vmean',vmn,'vinf',vinf,'z',z, ...
             'agree',agree,'atmedian',atmed,'narrow',D.narrow,'n',n,'p0',p0,'sig',sig);
end

% ------------------------------------------------------------------------
function pk = gridPeak(X, sig, dx)
%GRIDPEAK  Argmax of the summed tents on a uniform grid, plateau midpoint.
[n, k] = size(X);
pk = zeros(n, 1);  chunk = 2000;
for c0 = 1:chunk:n
    rows = c0:min(n, c0 + chunk - 1);
    Xr = X(rows,:);
    lo = floor((min(Xr(:)) - sig)/dx)*dx;  hi = ceil((max(Xr(:)) + sig)/dx)*dx;
    xg = lo:dx:hi;
    Cg = zeros(numel(rows), numel(xg));
    for i = 1:k
        Cg = Cg + max(0, 1 - abs(xg - Xr(:,i))/sig);
    end
    [cm, jm] = max(Cg, [], 2);
    tol = 1e-9*k;
    for r = 1:numel(rows)
        top = Cg(r,:) >= cm(r) - tol;
        j1 = jm(r);  j2 = jm(r);
        while j1 > 1 && top(j1-1), j1 = j1 - 1; end
        while j2 < numel(xg) && top(j2+1), j2 = j2 + 1; end
        pk(rows(r)) = 0.5*(xg(j1) + xg(j2));
    end
end
end

function pk = exactPeak(X, sig)
%EXACTPEAK  Maximum over the scenter positions, plateau midpoint. Disconnected
%   equal maxima, which occur for two scenters farther apart than sigma, are
%   resolved by a fair draw, which the grid does through its discretization.
[n, k] = size(X);
pk = zeros(n, 1);  tol = 1e-9*k;
for r = 1:n
    xi = sort(X(r,:));
    c = sum(max(0, 1 - abs(xi' - xi)/sig), 2)';
    cm = max(c);
    j = find(c >= cm - tol);
    mids = [];  j1 = j(1);  j2 = j1;
    for m = j(2:end)
        xm = 0.5*(xi(j2) + xi(m));
        if sum(max(0, 1 - abs(xm - xi)/sig)) >= cm - tol
            j2 = m;
        else
            mids(end+1) = 0.5*(xi(j1) + xi(j2)); %#ok<AGROW>
            j1 = m;  j2 = m;
        end
    end
    mids(end+1) = 0.5*(xi(j1) + xi(j2)); %#ok<AGROW>
    pk(r) = mids(randi(numel(mids)));
end
end
