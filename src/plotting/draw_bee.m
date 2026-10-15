function h = draw_bee(ax, x, y, theta, s, col, varargin)
%DRAW_BEE  Vector top-view honeybee glyph for schematics and figure insets.
%   H = DRAW_BEE(AX, X, Y, THETA, S, COL) draws a bee centered at (X,Y) in
%   data units of AX, heading THETA radians (0 points along +x), with body
%   length S in data units and body tint COL (RGB). Returns the graphics
%   handles. The axes should have equal data aspect ratio, or the glyph is
%   stretched.
%
%   Name-value options
%     'fan'      true draws blurred, spread wings for a fanning scenter
%     'wings'    false omits the wings, for very small glyphs
%     'legs'     false omits the legs
%     'edge'     outline width in points, default 0.8
%     'legw'     leg and antenna width in points, default 1.2*edge
%
%   The glyph is built from patches and lines only, so it exports as vector
%   through linger_style and survives assemble_figure.

ip = inputParser;
addParameter(ip, 'fan', false);
addParameter(ip, 'wings', true);
addParameter(ip, 'legs', true);
addParameter(ip, 'edge', 0.8);
addParameter(ip, 'legw', []);
parse(ip, varargin{:});
o = ip.Results;
if isempty(o.legw), o.legw = 1.2*o.edge; end

R    = [cos(theta) -sin(theta); sin(theta) cos(theta)];
tf   = @(P) (R*(P.'*s)).' + [x y];
sat  = max(col) - min(col);
warm = 0.12*min(1, sat/0.3);                 % a faint warm cast on saturated tints only
body = (1 - warm)*col + warm*[0.95 0.80 0.35];
dark = col*0.45;
lum  = [0.2126 0.7152 0.0722]*col(:);
if lum < 0.3
    band = col + 0.30*(1 - col);             % faint hair rings read on a dark body
else
    band = dark;
end
wing = [0.88 0.94 0.98];
t    = linspace(0, 2*pi, 48).';
el   = @(cx,cy,a,b,ph) [cx + a*cos(t)*cos(ph) - b*sin(t)*sin(ph), ...
                        cy + a*cos(t)*sin(ph) + b*sin(t)*cos(ph)];
h = gobjects(0);
    function add(P, varargin)
        Q = tf(P);
        h(end+1) = patch(ax, 'XData', Q(:,1), 'YData', Q(:,2), varargin{:});
    end
    function line_(P, varargin)
        Q = tf(P);
        h(end+1) = plot(ax, Q(:,1), Q(:,2), '-', varargin{:});
    end

% Six legs from the thorax, drawn first so the body covers their roots.
if o.legs
    L = {[0.14 0.05; 0.21 0.15; 0.28 0.19], ...
         [0.09 0.07; 0.08 0.20; 0.02 0.27], ...
         [0.04 0.06; -0.07 0.18; -0.20 0.24]};
    for sg = [-1 1]
        for k = 1:3
            P = L{k};  P(:,2) = sg*P(:,2);
            line_(P, 'Color', dark, 'LineWidth', o.legw);
        end
    end
end

% Fanning wings are spread and sit behind the body.
if o.wings && o.fan
    for sg = [-1 1]
        for k = -1:1
            add(el(0.02, sg*0.22, 0.32, 0.085, sg*(1.0 + 0.22*k)), ...
                'FaceColor', [0.78 0.86 0.92], 'FaceAlpha', 0.40, ...
                'EdgeColor', [0.50 0.58 0.64], 'EdgeAlpha', 0.45, 'LineWidth', 0.4*o.edge);
        end
    end
end

% Resting wings sweep back behind the body and show at its sides.
if o.wings && ~o.fan
    for sg = [-1 1]
        add(el(-0.16, sg*0.20, 0.29, 0.080, -sg*0.50), 'FaceColor', wing, ...
            'FaceAlpha', 0.65, 'EdgeColor', [0.55 0.62 0.68], 'LineWidth', 0.5*o.edge);
        add(el(-0.09, sg*0.15, 0.19, 0.055, -sg*0.40), 'FaceColor', wing, ...
            'FaceAlpha', 0.55, 'EdgeColor', [0.55 0.62 0.68], 'LineWidth', 0.5*o.edge);
    end
end

% Tapered abdomen, widest near the waist and pointed at the tip.
u   = linspace(0, 1, 40).';
xa  = @(u) 0.00 - 0.62*u;
ya  = @(u) 0.18*sin(pi*u.^0.72).^0.85;
add([xa(u) ya(u); flipud(xa(u)) -flipud(ya(u))], 'FaceColor', body, ...
    'EdgeColor', dark, 'LineWidth', o.edge);
for uc = [0.36 0.54 0.72 0.90]
    ub = linspace(uc - 0.055, min(uc + 0.055, 1), 8).';
    add([xa(ub) ya(ub); flipud(xa(ub)) -flipud(ya(ub))], 'FaceColor', band, 'EdgeColor', 'none');
end

% Thorax, head with compound eyes, elbowed antennae.
add(el(0.10, 0, 0.12, 0.115, 0), 'FaceColor', col*0.6, 'EdgeColor', dark, 'LineWidth', o.edge);
add(el(0.27, 0, 0.075, 0.095, 0), 'FaceColor', col*0.6, 'EdgeColor', dark, 'LineWidth', o.edge);
for sg = [-1 1]
    add(el(0.275, sg*0.058, 0.045, 0.030, sg*0.25), 'FaceColor', dark, 'EdgeColor', 'none');
    line_([0.33 sg*0.025; 0.39 sg*0.08; 0.45 sg*0.09], 'Color', dark, 'LineWidth', max(0.85*o.legw, 0.8*o.edge));
end
end
