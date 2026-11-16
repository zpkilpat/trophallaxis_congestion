function out = figure2(varargin)
%FIGURE2 Reproduce the three panels of Figure 2 (search for one donor) and
%   assemble the lettered figure output/figs/fig2.pdf.
%
% Usage
%   figure2
%   figure2('outdir','output/figure2')
%
% Panel map
%   (a) fig2a   point-donor MFPT ratio T(x;beta)/T(x;0) against start position
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

% Point donor. With a constant reversal rate beta the mean first-passage time
% from position x, averaged over headings, is T(x;beta) = C/2 - beta*y^2/v^2
% with y = pi - |x| and C/2 = pi/(v*eta) + beta*pi^2/v^2, so that
% T(x;0) = pi/(v*eta) at every start and
%   T(x;beta)/T(x;0) = 1 + eta*beta*|x|*(2*pi - |x|)/(pi*v).
% The ring average of the ratio is 1 + 2*pi*eta*beta/(3*v).
v   = prm.v;  eta = tanh(prm.kbar/(2*v));
xg  = linspace(-pi, pi, 721);
bl  = [0.5 1 2 4];
cols = P.greens(3:6,:);
h = P.figure(P);  ax = gca;  hold(ax,'on');
Rall = zeros(numel(bl), numel(xg));  Rbar = zeros(1,numel(bl));
for i = 1:numel(bl)
    R = 1 + eta*bl(i)*abs(xg).*(2*pi - abs(xg))/(pi*v);
    Rall(i,:) = R;
    Rbar(i)   = 1 + 2*pi*eta*bl(i)/(3*v);
    plot(ax, xg, R, '-', 'Color', cols(i,:), 'LineWidth', P.LW);
    plot(ax, [-pi pi], Rbar([i i]), ':', 'Color', cols(i,:), 'LineWidth', P.ALW);
end
% Check the closed form against the ring average of the plotted curves.
for i = 1:numel(bl)
    assert(abs(trapz(xg, Rall(i,:))/(2*pi) - Rbar(i)) < 1e-4);
end
plot(ax, [0 0], [1 1.26], ':', 'Color', [.25 .25 .25], 'LineWidth', P.ALW);
xlabel(ax,'start position $x$','Interpreter','latex','FontSize',P.LSZ);
ylabel(ax,'time ratio $T(x;\beta)/T(x;0)$','Interpreter','latex','FontSize',P.LSZ);
xlim(ax,[-pi pi]);  ylim(ax,[0.995 1.26]);  P.axes(P);
set(ax,'XTick',[-pi 0 pi],'XTickLabel',{'$-\pi$','$0$','$\pi$'}, ...
    'TickLabelInterpreter','latex');
set(ax,'Position',[0.19 0.14 0.65 0.71]);
for i = 1:numel(bl)
    text(ax, pi+0.12, Rall(i,end), sprintf('$\\beta=%g$', bl(i)), ...
         'Interpreter','latex','FontSize',P.FSZ,'Color',cols(i,:), ...
         'HorizontalAlignment','left','VerticalAlignment','middle','Clipping','off');
end
text(ax, -0.12, 1.255, 'donor', 'Interpreter','latex','FontSize',P.FSZ, ...
     'Color',[.25 .25 .25], 'Rotation',90, 'VerticalAlignment','bottom','HorizontalAlignment','right');
text(ax, 0.3, 1.235, '$\cdots$ ring mean', 'Interpreter','latex', ...
     'FontSize',P.FSZ, 'Color',[.3 .3 .3], 'HorizontalAlignment','left');
out.a = struct('x',xg,'beta',bl,'ratio',Rall,'ringaverage',Rbar,'eta',eta);
fprintf('fig2(a): point donor, eta=%.5f, ring-average ratio %s\n', eta, mat2str(Rbar,5));
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
