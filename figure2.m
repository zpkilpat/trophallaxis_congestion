function out = figure2(varargin)
%FIGURE2 Reproduce the three panels of Figure 2 (search for one donor) and
%   assemble the lettered figure output/figs/fig2.pdf.
%
% Usage
%   figure2
%   figure2('outdir','output/figure2')
%
% Panel map
%   (a) fig2a   MFPT ratio T(x;beta)/T(x;0) over start position and turning rate
%   (b) fig2b   averaged MFPT versus cue reliance gamma at several cue offsets
%   (c) fig2c   gain surface W(gamma,p) with the optimal-reliance ridge overlaid
%
% The panels are drawn by panels_abc, below. The first-passage arc solver is afpt, the
% gain tables come from linger_tables (which calls lingermodel), and the house
% style is linger_style.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'src')));

ip = inputParser;
addParameter(ip,'outdir',fullfile('output','figure2'));
parse(ip,varargin{:});
OUTDIR = ip.Results.outdir;
if ~isfolder(OUTDIR), mkdir(OUTDIR); end

A   = afpt();
prm = A.mkprm('kbar',0.0314,'sig',0.8,'beta0',1,'sm',1,'v',1);
out = panels_abc(prm, OUTDIR);

% Lettered figure, one row of three panel canvases at the panel size.
figdir = fullfile('output','figs');
if ~isfolder(figdir), mkdir(figdir); end
T = [22 19];
tiles = [0 0 T; T(1) 0 T; 2*T(1) 0 T];
panels = fullfile(OUTDIR, strcat({'fig2a','fig2b','fig2c'}, '.fig'));
assemble_figure(fullfile(figdir,'fig2'), panels, tiles, [3*T(1) T(2)]);
end

% =========================================================================
% Panels (a)-(c), fig2a to fig2c. Search for one donor.
% =========================================================================
function out = panels_abc(prm, outdir)
%PANELS_ABC  Compute and export panels (a)-(c) of Figure 2, the single-donor search.
%
% All panels use the bounded chemotactic turning law, which admits an exact
% arc-wise first-passage calculation. Unless stated otherwise, the
% point-target panels use kbar=0.0314, sigma=0.8, beta0=s_m=v=1, matching the
% collective model in Figures 3 and 4.

if nargin < 1 || isempty(prm)
    A0  = afpt();
    prm = A0.mkprm('kbar',0.0314,'sig',0.8,'beta0',1,'sm',1,'v',1);
end
if nargin < 2 || isempty(outdir), outdir = 'output'; end
if ~isfolder(outdir), mkdir(outdir); end

A = afpt();  P = linger_style();
out = struct('prm',prm);

% Finite-width parameters for panel (a). These differ from kbar because the
% ratio plane needs a resolved interaction region. Here 2*kappa*a = 0.96.
FA = 0.4;  FK = 1.2;  v = prm.v;
out.finite = struct('a',FA,'kappa',FK,'kbar',2*FK*FA);

%% (a) fig2a --------------------------------------------------------------
xg = linspace(-pi, pi, 361);
bg = logspace(-2, 1, 121);
T0 = Tfinite(xg, 0, v, FA, FK);
R  = zeros(numel(bg), numel(xg));
for i = 1:numel(bg)
    R(i,:) = Tfinite(xg, bg(i), v, FA, FK) ./ T0;
end

% The beta = 0 profile is flat outside the target only. Inside, the MFPT
% follows cosh(lambda x) and dips toward the center, so flatness is tested
% outside the target.
out_ = abs(xg) > FA;
flat = max(abs(T0(out_)/mean(T0(out_)) - 1));
if flat > 1e-10
    warning('fig2:notFlatAtZero', ...
        ['the beta = 0 profile varies by %.1e OUTSIDE the target; a ' ...
         'ballistic searcher on a ring has no preferred start there'], flat);
end

% Color limits are symmetric about unity with half-width 1 - min(R). The
% blue side fills exactly and the red side saturates. State in the caption
% that the upper end is clipped.
h = P.figure(P, true);  ax = gca;
hs = pcolor(ax, xg, log10(bg), R);  set(hs,'EdgeColor','none');
shading(ax,'interp');  colormap(ax, P.divmap(256));
cl = 1 - min(R(:));  clim(ax, [1-cl, 1+cl]);
cb = P.cbar(P, 'time ratio to ballistic, $T(x;\beta)/T(x;0)$');
cb.Label.FontSize = P.FSZ;
hold(ax,'on');
plot(ax, [-FA -FA], log10([min(bg) max(bg)]), ':', 'Color',[.25 .25 .25], 'LineWidth',P.ALW);
plot(ax, [ FA  FA], log10([min(bg) max(bg)]), ':', 'Color',[.25 .25 .25], 'LineWidth',P.ALW);
contour(ax, xg, log10(bg), R, [1 1], '-', 'Color','k', 'LineWidth', P.ALW+0.4);
xlabel(ax,'start position $x$','Interpreter','latex','FontSize',P.LSZ);
ylabel(ax,'reversal rate magnitude $\log_{10}\beta$','Interpreter','latex','FontSize',P.LSZ);
xlim(ax,[-pi pi]);  P.axes(P);
set(ax,'Position',[0.19 0.14 0.60 0.71]);  cb.Position = [0.81 0.14 0.035 0.71];
% In-panel guides. The band where turning helps is narrow, so it is pointed
% to from the red field. The dotted lines are the donor edges.
tl = @(x,y,str,varargin) text(ax, x, y, str, 'Interpreter','latex', ...
        'FontSize', P.FSZ-2, varargin{:});
tl(-2.90, 0.55, {'turning','hurts'}, 'Color','w');
tl( 0.50, 0.30, {'$\leftarrow$ turning','\ \ \ \ helps'}, 'Color','w');
tl(-0.68, -1.30, 'donor', 'Color',[.25 .25 .25], 'Rotation',90, ...
   'HorizontalAlignment','center');
tl( 0.60, -1.55, {'rare turns,','near ballistic'}, 'Color',[.25 .25 .25]);
% Band width is measured by linear interpolation of the R = 1 crossings.
% Counting grid columns would quantize it to the column spacing.
bw = zeros(size(bg));
for i = 1:numel(bg)
    k = find(R(i,:) < 1);
    if isempty(k), bw(i) = 0; continue; end
    lo = k(1);  hi = k(end);
    xl = xg(lo);  xh = xg(hi);
    if lo > 1,          xl = interp1(R(i,lo-1:lo), xg(lo-1:lo), 1); end
    if hi < numel(xg),  xh = interp1(R(i,hi:hi+1), xg(hi:hi+1), 1); end
    bw(i) = xh - xl;
end
out.a = struct('x',xg,'beta',bg,'ratio',R,'flatdev',flat,'clim',[1-cl 1+cl], ...
               'bandwidth',bw,'Rmin',min(R(:)),'Rmax',max(R(:)));
fprintf(['fig2(a): ratio %.3f to %.3f, color clipped at %.3f;  the sub-unity ' ...
         'band runs %.4f to %.4f rad wide over the beta sweep (%.2f%% variation)\n'], ...
        min(R(:)), max(R(:)), 1+cl, min(bw), max(bw), ...
        100*(max(bw)-min(bw))/max(bw,[],'all'));
P.export(P,'fig2a',outdir,h);

%% (b) fig2b --------------------------------------------------------------
pl = [0 0.30 0.55 0.75];
gg = linspace(0, prm.gmax, 121);
gscale = prm.sm*prm.sig;  % Normalize displayed reliance; retain gamma in calculations.
u = gg/gscale;
Tc = zeros(numel(pl), numel(gg));  gstar = zeros(1,numel(pl));
for i = 1:numel(pl)
    for j = 1:numel(gg)
        Tc(i,j) = A.moments(pl(i), prm.sig, prm.beta0, gg(j), prm.sm, prm.theta, prm.v);
    end
    [~,k] = min(Tc(i,:));  gstar(i) = gg(k);
end
h = P.figure(P);  ax = gca;  hold(ax,'on');
cols = P.greens(3:6,:);
for i = 1:numel(pl)
    plot(ax, u, Tc(i,:), '-', 'Color', cols(i,:), 'LineWidth', P.LW);
    [tm,k] = min(Tc(i,:));
    plot(ax, u(k), tm, 'o', 'MarkerSize', 10, 'MarkerFaceColor','w', ...
         'MarkerEdgeColor', cols(i,:), 'LineWidth', P.ALW);
end
xlabel(ax,'reliance $\gamma/(s_m\sigma)$','Interpreter','latex','FontSize',P.LSZ);
ylabel(ax,'discovery time $\bar T(\gamma,p)$','Interpreter','latex','FontSize',P.LSZ);
xlim(ax,[0 1]);  P.axes(P);
set(ax,'XTick',0:0.25:1);
set(ax,'Position',[0.19 0.14 0.76 0.71]);
% Each curve is labeled by its offset just above it, and the markers are
% named once.
xlab = 0.775;  [~, jl] = min(abs(u - xlab));
for i = 1:numel(pl)
    text(ax, xlab, Tc(i,jl) + 12, sprintf('$p=%g$', pl(i)), 'Interpreter','latex', ...
         'FontSize', P.FSZ, 'Color', cols(i,:), 'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom');
end
text(ax, 0.0375, 62, '$\circ$ optimal reliance', 'Interpreter','latex', 'FontSize', P.FSZ, ...
     'Color', [.3 .3 .3], 'HorizontalAlignment','left', 'VerticalAlignment','bottom');
out.b = struct('p',pl,'gamma',gg,'Tbar',Tc,'gstar',gstar);
% Inset, the cue geometry. A searcher heads along the unrolled ring toward a
% tent whose peak is displaced by p from the donor.
inset_cue(h, P);
P.export(P,'fig2b',outdir,h);

%% (c) fig2c gain surface and optimal-reliance ridge -------------------
[G, W] = linger_tables(prm, linspace(0, 1.2, 61), fullfile(outdir,'figure2_gain_cache.mat'));
h = P.figure(P, true);  ax = gca;
hs = pcolor(ax, W.p, W.gam/gscale, W.G);  set(hs,'EdgeColor','none');
shading(ax,'interp');  colormap(ax, P.plasma(256));
cb = P.cbar(P, 'gain $W(\gamma,p)$');
cb.Label.FontSize = P.LSZ-4;
hold(ax,'on');
plot(ax, G.p, G.gam/gscale, '-', 'Color','w','LineWidth',P.LW+1.2);
plot(ax, G.p, G.gam/gscale, '-', 'Color',[.15 .15 .15],'LineWidth',P.LW);
contour(ax, W.p, W.gam/gscale, W.G, [1 1], '--', 'Color','w','LineWidth',P.ALW);
xlabel(ax,'cue offset $p$','Interpreter','latex','FontSize',P.LSZ);
ylabel(ax,'reliance $\gamma/(s_m\sigma)$','Interpreter','latex','FontSize',P.LSZ);
xlim(ax,[0 max(W.p)]);  ylim(ax,[0 1]);  P.axes(P);
set(ax,'YTick',0:0.25:1);
set(ax,'Position',[0.19 0.14 0.60 0.71]);  cb.Position = [0.81 0.14 0.035 0.71];
% In-panel guides for the strategy plane.
tc = @(x,y,str,varargin) text(ax, x, y/gscale, str, 'Interpreter','latex', ...
        'FontSize', P.FSZ-2, 'Color','w', varargin{:});
tc(0.08, 0.735, '$\gamma^\ast(p)/(s_m\sigma)$');
tc(0.05, 0.36, {'follow','the cue'});
tc(0.98, 0.50, {'following','the cue','hurts'}, 'HorizontalAlignment','center');
tc(0.98, 0.09, {'ignore','the cue'}, 'HorizontalAlignment','center');
tc(0.775, 0.30, '$W=1$', 'Rotation', -84, 'HorizontalAlignment','center');
out.c = struct('p',W.p,'gam',W.gam,'G',W.G,'G0',max(W.G(:)),'Gmin',min(W.G(:)));
fprintf('fig2(c): ceiling %.3f, floor %.3f\n', out.c.G0, out.c.Gmin);
P.export(P,'fig2c',outdir,h);

end

% ------------------------------------------------------------------------
function T = Tfinite(x, beta, v, a, kappa)
%TFINITE  Direction-averaged MFPT for a finite-width target, closed form.
%   The closed-form finite-target MFPT, summed over headings and halved. At beta = 0
%   it is independent of x outside the target, which the caller checks.
lam = sqrt((2*beta+kappa)*kappa)/v;
out = 2/kappa + 2*(pi-a)*lam./(kappa*tanh(lam*a));
T   = zeros(size(x));
r   = x >  a;   l = x < -a;   c = ~r & ~l;
T(r) = out + 4*pi*beta/v^2*(x(r)-a) + 2*beta/v^2*(a^2 - x(r).^2);
T(l) = out - 4*pi*beta/v^2*(x(l)+a) + 2*beta/v^2*(a^2 - x(l).^2);
T(c) = 2/kappa + 2*(pi-a)*lam/(kappa*sinh(lam*a))*cosh(lam*x(c));
T = T/2;
end

% ------------------------------------------------------------------------
function inset_cue(h, P)
%INSET_CUE  Cue geometry sketch for panel (b), in the empty upper left.
K = inset_kit();  C = K.C;
ax = K.axes(h, P, [4.5 11.9 7.0 4.3]);
y0 = 0.75;  xp = 3.85;  xd = 5.35;
K.tent(ax, xp, 1.6, y0, 2.2, 0.30, 0.93);
K.track(ax, 0.15, 6.85, y0);
plot(ax, [xd xd], [y0+0.4 3.15], ':', 'Color', C.mgrey, 'LineWidth', 1.5);
K.donor(ax, xd, y0, 0.28);
K.bee(ax, 1.15, y0, 0, 1.0, 'search');
K.bracket(ax, xp, xd, 3.3, 0.15);
K.label(ax, (xp+xd)/2, 3.42, '$p$', P, 'bottom');
end
