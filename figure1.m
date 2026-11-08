function out = figure1(varargin)
%FIGURE1 Reproduce Figure 1 of
% "Chemotactic recruitment offsets congestion in honeybee trophallaxis."
%
% Usage
%   figure1
%   figure1('outdir','output/figure1')
%
% Outputs
%   fig1a.pdf/.png  search geometry and displaced recruitment cue
%   fig1b.pdf/.png  population flow and queueing schematic
%   fig1c.pdf/.png  small versus large abandonment rate
%   output/figs/fig1.pdf/.png/.fig  lettered figure
%
% The panels are drawn by draw_panels below, in the house style of
% linger_style.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'src')));

ip = inputParser;
addParameter(ip,'outdir',fullfile('output','figure1'));
parse(ip,varargin{:});
outdir = ip.Results.outdir;
if ~isfolder(outdir), mkdir(outdir); end

out = draw_panels(outdir);

% Lettered figure. Each tile is its panel canvas, so type size is kept.
figdir = fullfile('output','figs');
if ~isfolder(figdir), mkdir(figdir); end
assemble_figure(fullfile(figdir,'fig1'), out.panels, out.tiles, out.canvas, ...
    'letterxy', [0.17 0; 0.09 0; 0.16 0], 'lettersize', 30*36/44);   % prints at the letter size of Figs. 2-6 at 0.85 textwidth
end

% =========================================================================
% Panels (a)-(c), fig1a to fig1c. Model schematic.
% =========================================================================
function out = draw_panels(outdir)
%DRAW_PANELS  Model schematic, as three standalone panels.
%
%   fig1a  search geometry. Ring with searching, waiting and scenting bees,
%          the donor and the cue plume, above the unrolled line with the tent
%          cue, its offset p and its footprint sigma.
%   fig1b  population flow S -> U_j -> F_j -> D with the scenting feedback
%          and the donor-blocking self-loop.
%   fig1c  small and large abandonment rate on the ring.
%
%   All geometry is in points on one 1024 x 538 pt page, y down. Each panel
%   is one axes over its own vertical strip of that page, and the panel canvas
%   is the strip itself, so assembly at tile = canvas keeps type size and
%   position. OUT.tiles and OUT.canvas (cm) feed assemble_figure.
%
%   Gradients are vertex-colored triangle meshes, which stay vector in PDF.
%   Bees are draw_bee glyphs heading along their motion. Waiting bees queue
%   head to tail toward the donor, and scenting bees face the donor with
%   fanned wings, so the air they drive backward runs down the plume.

if nargin < 1, outdir = fullfile('output','figure1'); end
P = linger_style();

PAGE = [1024 537.93];                        % pt
XCUT = [0 405 745 1024];                     % panel strips, pt
cm   = 2.54/72;

C.blue  = [0 162 255]/255;                   % searching
C.navy  = [0 77 128]/255;                    % waiting
C.green = [29 177 0]/255;                    % scenting, fed
C.dkgrn = [1 113 0]/255;                     % cue, feedback
C.dgrey = [94 94 94]/255;                    % donor
C.mgrey = [146 146 146]/255;                 % departed
C.ring  = [213 213 213]/255;
C.alpha = 0.7;

names = {'fig1a','fig1b','fig1c'};
draw  = {@panelA, @panelB, @panelC};
tiles = zeros(3,4);
for i = 1:3
    x0 = XCUT(i);  x1 = XCUT(i+1);
    W = (x1 - x0)*cm;  H = PAGE(2)*cm;
    f = figure('Units','centimeters','Position',[1 1 W H],'Color','w', ...
               'Visible','off','InvertHardcopy','off');
    set(f,'PaperUnits','centimeters','PaperSize',[W H],'PaperPosition',[0 0 W H]);
    ax = axes(f,'Units','normalized','Position',[0 0 1 1], ...
              'XLim',[x0 x1],'YLim',[0 PAGE(2)],'YDir','reverse', ...
              'DataAspectRatio',[1 1 1],'Visible','off','Clipping','off');
    hold(ax,'on');
    % white corner marks, so the tight crop of the export keeps the full strip
    for cx = [x0 x1-0.2], for cy = [0 PAGE(2)-0.2]
        patch(ax, cx + [0 .2 .2 0], cy + [0 0 .2 .2], 'w', 'EdgeColor','none');
    end, end
    draw{i}(ax, C);
    set(ax.Children, 'Clipping', 'off');
    P.export(P, names{i}, outdir, f);
    close(f);
    tiles(i,:) = [x0*cm 0 W H];
end

out = struct('tiles',tiles,'canvas',PAGE*cm, ...
             'panels',{fullfile(outdir, strcat(names,'.fig'))});
end

% ===================================================================== (a)
function panelA(ax, C)
SA = 36;  SL = 26;                           % glyph body length, ring and line
ctr = [173.3 176.2];  R = 150;
plume(ax, C, ctr, R, [245.9 290.0], 274.1, [174.7 57.6], 216.0, ...
      [180 25], [391 58], [319 218], [0.1029 0.002384]);
ringLine(ax, C, ctr, R, 15);
plumeEdge(ax, C, ctr, R, [245.9 290.0], 274.1, [174.7 57.6], 216.0, ...
          [180 25], [391 58], [319 218]);

dn = [311 121.8];
% searching bees fly tangent to the ring, in the sense of their motion
ab = [-1.80 -2.17 -1.45 pi 2.33 0.95]';
sb = ctr + R*[cos(ab) sin(ab)];
sv = [-1 0; -1 1; 1 0; 0 -1; 1 1; -1 1];
for k = 1:size(sb,1)
    hd = tangent(ctr, sb(k,:), sv(k,:));
    glyph(ax, C, sb(k,:), hd, SA, 'search');
    if k == 2                                % the speed v, ahead of one bee
        a = sb(k,:) + 0.62*SA*hd;
        arrow(ax, a, a + 28*hd, 3, 11, 11, 'k');
        nv = (sb(k,:) - ctr)/R;
        math(ax, a(1) + 14*hd(1) + 15*nv(1), a(2) + 14*hd(2) + 15*nv(2), '$v$', 24, 'k');
    end
end
% waiting bees queue inside the ring, each facing the donor
[wb, wh] = queue(ctr, 127, dn, 0.22 + 0.25*(0:3));
for k = size(wb,1):-1:1, glyph(ax, C, wb(k,:), wh(k,:), 34, 'wait'); end
% scenting bees sit on the ring at the donor, head in, fanning down the plume
gb = [249.6 47.9; 277 65.2; 299 92.9; 326 154.2; 325.3 188.8];
for k = 1:size(gb,1)
    glyph(ax, C, gb(k,:), tangent(ctr, gb(k,:), dn - gb(k,:)), SA, 'scent');
end
bee(ax, C, dn, 14, 'donor', 2);

sans(ax, 47.0, 187.1, 'searching', 24, C.blue);
sans(ax, 189.8, 165.7, 'waiting', 24, C.navy);
sans(ax, 221.2, 118.4, 'donor', 24, C.dgrey);
sans(ax, 171.4, 85.0, 'scenting', 24, C.green);
math(ax, 121, 84, '$\beta_\pm(x)$', 24, C.blue);

% unrolled line
y0 = 464.2;
patch(ax, [198 287 240.6], [458.7 458.7 397.5], 'k', 'EdgeColor','none', ...
      'FaceVertexCData', [cue(C,0.08); cue(C,0.08); cue(C,0.965)], ...
      'FaceColor','interp');
plot(ax, [22 399], [y0 y0], '-', 'Color', C.ring, 'LineWidth', 5.6);
for xe = [22.9 399]
    patch(ax, xe + [-3 3 3 -3], [452 452 475.5 475.5], C.ring, 'EdgeColor','none');
end
bracket(ax, 200, 242, 480, 3.5, 1.8, C.mgrey);
for s = [64.6 1; 116 1; 186 -1; 369.2 -1]'
    glyph(ax, C, [s(1) y0], [s(2) 0], SL, 'search');
end
for xb = [223.5 251], glyph(ax, C, [xb y0], [1 0], SL, 'scent'); end
glyph(ax, C, [302 y0], [-1 0], SL, 'wait');
bee(ax, C, [278 y0], 9.85, 'donor', 1.8);
bracket(ax, 244, 278, 423, 4, 3, C.dgrey);
serif(ax, 266.7, 382.3, 'offset', 20, C.dgrey, 'center');
math(ax, 267.4, 404, '$p$', 24, C.dgrey);
serif(ax, 156.2, 496.7, 'footprint', 20, C.mgrey, 'left');
math(ax, 246.5, 496.5, '$\sigma$', 24, C.dgrey);
math(ax, 33.2, 503, '$-\pi$', 30, 'k');
math(ax, 394.5, 504.5, '$+\pi$', 30, 'k');
end

% ===================================================================== (b)
function panelB(ax, C)
bx = {[469.1 7.2 660.8 74.6], [473.4 112.8 656.4 184.4], ...
      [472.4 283.5 654.6 350.5], [475.6 410 658.2 477.7]};
bc = {C.blue, C.navy, C.green, C.mgrey};
for k = 1:4
    [xx, yy] = rrect(bx{k}, 12.6);
    patch(ax, xx, yy, bc{k}, 'EdgeColor','none');
end
math(ax, 557, 26, '$S$', 30, 'w');
serif(ax, 563.4, 54.4, 'hungry, searching', 24, 'w', 'center');
math(ax, 565, 134, '$U_j$', 30, 'w');
serif(ax, 564.4, 161.9, 'hungry, waiting', 24, 'w', 'center');
math(ax, 559.5, 307, '$F_j$', 30, 'k');
serif(ax, 563, 333, 'fed, lingering', 24, 'k', 'center');
math(ax, 561, 431, '$D$', 30, 'k');
serif(ax, 565.4, 459.8, 'departed', 24, 'k', 'center');

SB = 38;
glyph(ax, C, [443.7 86.4], [-0.5 -1], SB, 'search');
glyph(ax, C, [443.9 197.4], [1 0], SB, 'wait');
glyph(ax, C, [443.9 317], [1 0], SB, 'scent');
glyph(ax, C, [443.9 450], [0.35 1], SB, 'gone');

arrow(ax, [528 74.6], [528 111.7], 2.9, 12.6, 12.6, 'k');
arrow(ax, [511.6 184.4], [511.6 283.5], 2.9, 12.6, 12.6, 'k');
arrow(ax, [511.6 350.5], [511.6 408.9], 2.9, 12.6, 12.6, 'k');
serif(ax, 539.1, 88.4, 'discovery', 20, 'k', 'left');
math(ax, 571, 218, '$\displaystyle\mu_j\frac{U_j}{U_j+F_j}$', 20, 'k');
serif(ax, 520.4, 268.6, 'service (feeding)', 20, 'k', 'left');
math(ax, 495, 370.7, '$\delta$', 20, 'k');
serif(ax, 430.6, 393, 'abandon', 20, 'k', 'left');

% scenting feedback, F_j back to S
sh = spline2([653.6 300.8; 671.6 277.4; 684.2 252; 689.6 216; 690.3 108;
              690.6 70; 689.4 52; 686.5 43]);
plot(ax, sh(:,1), sh(:,2), '-', 'Color', C.dkgrn, 'LineWidth', 5.8);
patch(ax, [661.5 687.8 684.2], [39.6 30.6 54], C.dkgrn, 'EdgeColor','none');
t = serif(ax, 707.8, 160.4, 'scenting sharpens cue', 20, C.dkgrn, 'center');
t.Rotation = -90;

% donor-blocking self-loop
lp = spline2([648.5 347.4; 680 348.4; 711.2 348.7; 722 338.6; 723.8 326;
              720.2 313.4; 711.2 306.2; 700.4 304.4; 689.6 305.5; 681 308]);
dashes(ax, lp, 6, 6.2, 5, C.dkgrn);
patch(ax, [657.2 680.6 681.7], [308 296.1 318.8], C.dkgrn, 'EdgeColor','none');
serif(ax, 546.9, 368.2, 'donor-blocking', 20, C.dkgrn, 'left');
math(ax, 700, 391, '$\displaystyle\mu_j\frac{F_j}{U_j+F_j}$', 20, C.dkgrn);

math(ax, 562, 508, '$\displaystyle S=N-D-\sum_j\Big[U_j+F_j\Big]$', 19, 'k');
end

% ===================================================================== (c)
function panelC(ax, C)
% sharp cue, crowded donor
ctr = [882.7 117.7];  R = 96.8;
% the plume's inner arc is its own circle, just outside the ring
plume(ax, C, [880.9 121.5], 99.4, [932.5 171.0], 156.0, [870.1 44.2], 151.4, ...
      [891 21], [1022 43], [977 151], [-2.7229 0.003693]);
ringLine(ax, C, ctr, R, 15);
plumeEdge(ax, C, [880.9 121.5], 99.4, [932.5 171.0], 156.0, [870.1 44.2], 151.4, ...
          [891 21], [1022 43], [977 151]);
SC = 25;  SQ = 23;                           % searching, crowd glyph length
dn = [971.6 82.1];
sb = [857.8 23.9; 821.8 40.3; 877.6 214.2];
sv = [1 0; -17 20; -1 0];
for k = 1:3
    glyph(ax, C, sb(k,:), tangent(ctr, sb(k,:), sv(k,:)), SC, 'search');
end
% waiting bees queue inside the ring on both sides of the donor
[w1, h1] = queue(ctr, 80, dn, -(0.30 + 0.27*(0:2)));
[w2, h2] = queue(ctr, 80, dn, 0.30 + 0.27*(0:3));
wb = [w1; w2];  wh = [h1; h2];
for k = size(wb,1):-1:1, glyph(ax, C, wb(k,:), wh(k,:), SQ, 'wait'); end
gb = queue(ctr, R + 1, dn, [-0.56 -0.29 0.29 0.56]);
for k = 1:size(gb,1)
    glyph(ax, C, gb(k,:), tangent(ctr, gb(k,:), dn - gb(k,:)), SQ, 'scent');
end
bee(ax, C, dn, 9.1, 'donor', 1.8);
serif(ax, 871.7, 92.1, 'sharp cue;', 20, 'k', 'center');
serif(ax, 871.7, 113.0, 'crowded donor', 20, 'k', 'center');
serif(ax, 871.7, 135.0, '(\it\delta\rm\bf small)', 20, 'k', 'center', 'tex');

% no cue, slow search
ctr = [891 405.1];
ringLine(ax, C, ctr, R, 15);
sb = [866.1 311.2; 915.4 311.2; 941 321.3; 962.6 340.8; 979.6 440.2;
      962.3 470.1; 937.8 490.6; 885.7 501.6; 816 471.2; 803 446.7;
      794.8 426.9; 795.7 382.5; 829.3 327.8];
% sense of motion along the ring, +1 clockwise on the page
ss = [1 -1 1 1 -1 1 -1 1 -1 1 1 -1 -1];
for k = 1:size(sb,1)
    r = sb(k,:) - ctr;
    glyph(ax, C, sb(k,:), ss(k)*[-r(2) r(1)], SC, 'search');
end
dn = [977.9 369.9];
[wb, wh] = queue(ctr, 78, dn, 0.05 + 0.27*(0:1));
for k = size(wb,1):-1:1, glyph(ax, C, wb(k,:), wh(k,:), SC, 'wait'); end
bee(ax, C, dn, 9.1, 'donor', 1.8);
serif(ax, 877.8, 379.3, 'no cue;', 20, 'k', 'center');
serif(ax, 877.8, 400.1, 'slow search', 20, 'k', 'center');
serif(ax, 877.8, 422.0, '(\it\delta\rm\bf large)', 20, 'k', 'center', 'tex');
end

% ================================================================ helpers
function bee(ax, C, xy, rout, kind, lw)
switch kind
    case 'search', fc = C.blue;  ec = C.navy;  a = C.alpha;
    case 'wait',   fc = C.navy;  ec = C.dgrey; a = C.alpha;
    case 'scent',  fc = C.green; ec = C.navy;  a = C.alpha;
    case 'donor',  fc = C.dgrey; ec = C.mgrey; a = 1;
end
th = linspace(0, 2*pi, 73);  r = rout - lw/2;
patch(ax, xy(1) + r*cos(th), xy(2) + r*sin(th), fc, 'FaceAlpha', a, ...
      'EdgeColor', ec, 'EdgeAlpha', a, 'LineWidth', lw);
end

function glyph(ax, C, xy, hd, s, kind)
% Bee glyph of body length S whose body midpoint sits at XY, heading along
% the page vector HD. The glyph is mirror symmetric, so the reversed y axis
% needs no correction beyond taking the angle of HD in data coordinates.
switch kind
    case 'search', col = C.blue;  fan = false;
    case 'wait',   col = C.navy;  fan = false;
    case 'scent',  col = C.green; fan = true;
    case 'gone',   col = C.mgrey; fan = false;
end
d = hd/norm(hd);
o = xy + 0.07*s*d;
draw_bee(ax, o(1), o(2), atan2(d(2), d(1)), s, col, 'fan', fan, 'edge', s/32);
end

function t = tangent(ctr, xy, ref)
% unit tangent to the circle about CTR through XY, on the side of REF
r = xy - ctr;  t = [-r(2) r(1)]/norm(r);
if dot(t, ref) < 0, t = -t; end
end

function [P, H] = queue(ctr, r, dn, da)
% Positions at radius R and angles DA from the donor DN about CTR, each bee
% facing the one ahead of it and the first facing the donor.
a0 = atan2(dn(2) - ctr(2), dn(1) - ctr(1));
P = ctr + r*[cos(a0 + da(:)) sin(a0 + da(:))];
A = [dn; P(1:end-1,:)];
H = A - P;
end

function arrow(ax, p0, p1, lw, hl, hw, col)
d = (p1 - p0)/norm(p1 - p0);  n = [-d(2) d(1)];
b = p1 - hl*d;
plot(ax, [p0(1) b(1)+0.3*d(1)], [p0(2) b(2)+0.3*d(2)], '-', 'Color', col, 'LineWidth', lw);
patch(ax, [p1(1) b(1)+hw/2*n(1) b(1)-hw/2*n(1)], ...
          [p1(2) b(2)+hw/2*n(2) b(2)-hw/2*n(2)], col, 'EdgeColor','none');
end

function bracket(ax, x0, x1, y, h, lw, col)
plot(ax, [x0 x1], [y y], '-', 'Color', col, 'LineWidth', lw);
plot(ax, [x0 x0], y + [-h h], '-', 'Color', col, 'LineWidth', lw);
plot(ax, [x1 x1], y + [-h h], '-', 'Color', col, 'LineWidth', lw);
end

function ringLine(ax, C, ctr, R, lw)
th = linspace(0, 2*pi, 361);
plot(ax, ctr(1) + R*cos(th), ctr(2) + R*sin(th), '-', 'Color', C.ring, 'LineWidth', lw);
end

function c = cue(C, t)
% cue color, dark green at t = 0 fading to white at t = 1
t = min(max(t(:), 0), 1);
c = (1 - t).*C.dkgrn + t.*[1 1 1];
end

function [inner, outer] = plumePaths(ctr, R, cT, rT, cR, rR, gS, gC, gT)
% Plume bounded by the ring on the inside and two circular arcs outside.
S  = xcirc(ctr, R, cT, rT, gS);
K  = xcirc(cT, rT, cR, rR, gC);
T  = xcirc(cR, rR, ctr, R, gT);
inner = arcPts(ctr, R, S, T, 200);
outer = [arcPts(cT, rT, S, K, 120); arcPts(cR, rR, K, T, 120)];
end

function plume(ax, C, ctr, R, cT, rT, cR, rR, gS, gC, gT, tc)
[inner, outer] = plumePaths(ctr, R, cT, rT, cR, rR, gS, gC, gT);
n = 160;
A = resample(inner, n);  B = resample(outer, n);
V = [A; B];
F = zeros(2*(n-1), 3);
for k = 1:n-1
    F(2*k-1,:) = [k k+1 n+k+1];
    F(2*k,:)   = [k n+k+1 n+k];
end
t = tc(1) + tc(2)*(V(:,1) - V(:,2));
patch(ax, 'Faces', F, 'Vertices', V, 'FaceVertexCData', cue(C, t), ...
      'FaceColor','interp', 'EdgeColor','none');
end

function plumeEdge(ax, C, ctr, R, cT, rT, cR, rR, gS, gC, gT)
[inner, outer] = plumePaths(ctr, R, cT, rT, cR, rR, gS, gC, gT);
e = [outer; flipud(inner)];
plot(ax, e(:,1), e(:,2), '-', 'Color', C.dkgrn, 'LineWidth', 1.3);
end

function p = xcirc(c1, r1, c2, r2, guess)
% intersection of two circles nearest GUESS
d = norm(c2 - c1);  a = (r1^2 - r2^2 + d^2)/(2*d);
h = sqrt(max(r1^2 - a^2, 0));
m = c1 + a*(c2 - c1)/d;  n = [-(c2(2) - c1(2)) (c2(1) - c1(1))]/d;
q = [m + h*n; m - h*n];
[~, i] = min(sum((q - guess).^2, 2));
p = q(i,:);
end

function P = arcPts(c, r, p0, p1, n)
a0 = atan2(p0(2) - c(2), p0(1) - c(1));
a1 = atan2(p1(2) - c(2), p1(1) - c(1));
da = mod(a1 - a0 + pi, 2*pi) - pi;
a = a0 + da*linspace(0, 1, n)';
P = [c(1) + r*cos(a), c(2) + r*sin(a)];
end

function Q = resample(P, n)
s = [0; cumsum(sqrt(sum(diff(P).^2, 2)))];
[s, iu] = unique(s);
Q = interp1(s, P(iu,:), linspace(0, s(end), n)');
end

function Q = spline2(P)
% Catmull-Rom curve through the rows of P
P = [2*P(1,:) - P(2,:); P; 2*P(end,:) - P(end-1,:)];
t = linspace(0, 1, 20)';  Q = [];
for k = 2:size(P,1)-2
    p0 = P(k-1,:); p1 = P(k,:); p2 = P(k+1,:); p3 = P(k+2,:);
    seg = 0.5*((2*p1) + (-p0 + p2).*t + (2*p0 - 5*p1 + 4*p2 - p3).*t.^2 + ...
               (-p0 + 3*p1 - 3*p2 + p3).*t.^3);
    Q = [Q; seg(1:end-1,:)];                                   %#ok<AGROW>
end
Q = [Q; P(end-1,:)];
end

function dashes(ax, P, on, off, lw, col)
s = [0; cumsum(sqrt(sum(diff(P).^2, 2)))];
[s, iu] = unique(s);  P = P(iu,:);
for s0 = 0:(on + off):s(end)
    ss = linspace(s0, min(s0 + on, s(end)), 12)';
    Q = interp1(s, P, ss);
    plot(ax, Q(:,1), Q(:,2), '-', 'Color', col, 'LineWidth', lw);
end
end

function [x, y] = rrect(b, r)
th = linspace(0, pi/2, 16);
cx = [b(3)-r, b(1)+r, b(1)+r, b(3)-r];
cy = [b(2)+r, b(2)+r, b(4)-r, b(4)-r];
a0 = [0, -pi/2, pi, pi/2];
x = [];  y = [];
for k = 1:4
    x = [x, cx(k) + r*cos(a0(k) - th)];                        %#ok<AGROW>
    y = [y, cy(k) + r*sin(a0(k) - th)];                        %#ok<AGROW>
end
end

function t = sans(ax, x, y, str, fs, col)
t = text(ax, x, y, str, 'FontName','Arial', 'FontWeight','bold', ...
         'FontSize', fs, 'Color', col, 'Interpreter','none', ...
         'HorizontalAlignment','left', 'VerticalAlignment','middle');
end

function t = serif(ax, x, y, str, fs, col, ha, interp)
if nargin < 8, interp = 'none'; end
t = text(ax, x, y, str, 'FontName','Times New Roman', 'FontWeight','bold', ...
         'FontSize', fs, 'Color', col, 'Interpreter', interp, ...
         'HorizontalAlignment', ha, 'VerticalAlignment','middle');
end

function t = math(ax, x, y, str, fs, col)
t = text(ax, x, y, str, 'Interpreter','latex', 'FontSize', fs, 'Color', col, ...
         'HorizontalAlignment','center', 'VerticalAlignment','middle');
end
