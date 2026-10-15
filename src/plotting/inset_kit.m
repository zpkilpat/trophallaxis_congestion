function K = inset_kit()
%INSET_KIT  Primitives for the small schematic insets drawn inside data panels.
%   K = INSET_KIT() returns the schematic palette of figure1 in K.C and
%   function handles that draw on an inset axes.
%
%     ax = K.axes(h, P, box)            inset axes over BOX = [x y w h], in
%                                       centimeters of the panel canvas
%     K.track(ax, x0, x1, y)            the unrolled ring, a light grey line
%     K.donor(ax, x, y, r)              grey donor disk
%     K.tent(ax, xc, hw, y0, ht, t0, t1)  cue tent with apex at XC, half width
%                                       HW and height HT, shaded from cue
%                                       tone T0 at the base to T1 at the apex
%     K.ghost(ax, xc, hw, y0, ht)       dashed outline of an alternative tent
%     K.bracket(ax, x0, x1, y, d)       offset bracket with end ticks
%     K.bee(ax, x, y, th, s, kind)      draw_bee glyph, KIND is 'search',
%                                       'wait' or 'scent'
%     K.label(ax, x, y, str, P, va)     LaTeX label at P.FSZ-2 in donor grey
%
%   The inset axes has one data unit per centimeter and equal data aspect, so
%   glyph sizes read directly in centimeters on the 22 x 19 cm panel. It is
%   invisible with a transparent background, is created after the data axes
%   so it stacks on top through assemble_figure, and leaves the figure's
%   current axes unchanged. Every element is a patch, line, or text, so the
%   inset exports as vector.

K.C = struct('blue', [0 162 255]/255, 'navy', [0 77 128]/255, ...
             'green', [29 177 0]/255, 'dkgrn', [1 113 0]/255, ...
             'dgrey', [94 94 94]/255, 'mgrey', [146 146 146]/255, ...
             'ring', [213 213 213]/255);
C = K.C;
K.axes    = @newAxes;
K.track   = @(ax, x0, x1, y) plot(ax, [x0 x1], [y y], '-', 'Color', C.ring, 'LineWidth', 5);
K.donor   = @(ax, x, y, r) donor(ax, C, x, y, r);
K.tent    = @(ax, xc, hw, y0, ht, t0, t1) tent(ax, C, xc, hw, y0, ht, t0, t1);
K.ghost   = @(ax, xc, hw, y0, ht) plot(ax, [xc-hw xc xc+hw], [y0 y0+ht y0], '--', ...
                'Color', 0.3*C.dkgrn + 0.7, 'LineWidth', 1.2);
K.bracket = @(ax, x0, x1, y, d) plot(ax, [x0 x1 NaN x0 x0 NaN x1 x1], ...
                [y y NaN y-d y+d NaN y-d y+d], '-', 'Color', C.dgrey, 'LineWidth', 2);
K.bee     = @(ax, x, y, th, s, kind) bee(ax, C, x, y, th, s, kind);
K.label   = @label;
end

% ------------------------------------------------------------------------
function ax = newAxes(h, P, box)
ax0 = get(h, 'CurrentAxes');
ax = axes(h, 'Units','normalized', 'Position', box./[P.FIG P.FIG], ...
          'XLim',[0 box(3)], 'YLim',[0 box(4)], 'DataAspectRatio',[1 1 1], ...
          'Visible','off', 'Color','none', 'Clipping','off', 'HitTest','off');
hold(ax, 'on');
if ~isempty(ax0), set(h, 'CurrentAxes', ax0); end
end

function donor(ax, C, x, y, r)
t = linspace(0, 2*pi, 73);
patch(ax, x + r*cos(t), y + r*sin(t), C.dgrey, 'EdgeColor', C.mgrey, 'LineWidth', 1.5);
end

function tent(ax, C, xc, hw, y0, ht, t0, t1)
% The cue tone runs from dark green toward white, as in Fig. 1(a). A thin edge in
% a mid tone keeps the apex visible where the fill fades into the page.
cue = @(t) (1 - t)*C.dkgrn + t*[1 1 1];
patch(ax, [xc-hw xc+hw xc], [y0 y0 y0+ht], 'k', 'EdgeColor','none', ...
      'FaceVertexCData', [cue(t0); cue(t0); cue(t1)], 'FaceColor','interp');
plot(ax, [xc-hw xc xc+hw], [y0 y0+ht y0], '-', 'Color', cue(0.45), 'LineWidth', 1.4);
end

function bee(ax, C, x, y, th, s, kind)
% Body midpoint at (x,y), heading th, as the glyph helper of figure1 places it.
switch kind
    case 'search', col = C.blue;  fan = false;
    case 'wait',   col = C.navy;  fan = false;
    case 'scent',  col = C.green; fan = true;
end
draw_bee(ax, x + 0.07*s*cos(th), y + 0.07*s*sin(th), th, s, col, 'fan', fan, 'edge', 1.0, 'legs', false);
end

function t = label(ax, x, y, str, P, va)
if nargin < 6, va = 'middle'; end
t = text(ax, x, y, str, 'Interpreter','latex', 'FontSize', P.FSZ-2, ...
         'Color', [94 94 94]/255, 'HorizontalAlignment','center', ...
         'VerticalAlignment', va);
end
