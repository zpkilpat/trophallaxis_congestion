function assemble_figure(stem, panels, tiles, canvas, varargin)
%ASSEMBLE_FIGURE  Compose panel .fig files into one lettered figure.
%   ASSEMBLE_FIGURE(STEM, PANELS, TILES, CANVAS) opens each panel file in the
%   cell array PANELS, copies its axes, colorbars, and legends onto one canvas
%   of size CANVAS = [W H] in centimeters, and writes STEM.pdf, STEM.png, and
%   STEM.fig through linger_style. TILES is an n-by-4 array [x y w h] in
%   centimeters, measured from the lower-left corner of the canvas. A tile the
%   size of a panel canvas, linger_style P.FIG, keeps its type size exactly.
%
%   Name-value options
%     'letters'   cell array of panel letters, default (a), (b), (c), ... as in
%                 the captions and the APS house style
%     'letterxy'  n-by-2 offsets [dx dy] in centimeters from the upper-left
%                 corner of each tile, default [0.2 -0.2]
%     'lettersize' letter size in points, default 36, sized so the
%                 parenthesized letters clear the top tick label

P  = linger_style();
n  = numel(panels);
ip = inputParser;
addParameter(ip,'letters', cellfun(@(c) ['(' c ')'], cellstr(char('a' + (0:n-1))'), 'uni', 0));
addParameter(ip,'letterxy', repmat([0.2 -0.2], n, 1));
addParameter(ip,'lettersize', 36);
parse(ip, varargin{:});
L = ip.Results;
if size(L.letterxy,1) == 1, L.letterxy = repmat(L.letterxy, n, 1); end

W = canvas(1);  H = canvas(2);
dst = figure('Units','centimeters','Position',[1 1 W H],'Color','w', ...
             'Visible','off','InvertHardcopy','off');
% Preserve the full canvas when openfig fits the window to the screen.
setappdata(dst, 'ExportCanvasSize', canvas);
set(dst,'PaperUnits','centimeters','PaperSize',[W H], ...
        'PaperPosition',[0 0 W H]);

for i = 1:n
    src = openfig(panels{i}, 'invisible');
    set(src, 'Units', 'centimeters');
    sp = get(src, 'Position');
    kids = findobj(src.Children, 'flat', '-not', 'Type', 'uimenu');
    kids = kids(isprop(kids, 'Position'));
    set(kids, 'Units', 'normalized');
    pos0 = get(kids, 'Position');
    if ~iscell(pos0), pos0 = {pos0}; end
    typ  = get(kids, 'Type');
    if ~iscell(typ), typ = {typ}; end
    new = copyobj(kids, dst);
    % Axes first, then attached colorbars and legends, so an attached object
    % does not resize its axes after being placed.
    order = [find(strcmp(typ,'axes')); find(~strcmp(typ,'axes'))];
    t = tiles(i,:);
    sx = t(3)/sp(3);  sy = t(4)/sp(4);
    for k = order(:)'
        p = pos0{k};
        set(new(k), 'Units', 'normalized', 'Position', ...
            [(t(1) + p(1)*sp(3)*sx)/W, (t(2) + p(2)*sp(4)*sy)/H, ...
              p(3)*sp(3)*sx/W,          p(4)*sp(4)*sy/H]);
    end
    close(src);
end

% Negative tick labels get a true minus sign. Under the latex tick
% interpreter a bare '-' prints as a hyphen. Ticks are frozen first so the
% relabeled positions cannot drift.
for a = findall(dst, 'Type', 'axes')'
    if ~strcmp(get(a, 'TickLabelInterpreter'), 'latex'), continue; end
    for d = 'XYZ'
        set(a, [d 'Tick'], get(a, [d 'Tick']));
        set(a, [d 'TickLabel'], minusLabels(get(a, [d 'TickLabel'])));
    end
end
for c = findall(dst, 'Type', 'colorbar')'
    if ~strcmp(get(c, 'TickLabelInterpreter'), 'latex'), continue; end
    set(c, 'Ticks', get(c, 'Ticks'));
    set(c, 'TickLabels', minusLabels(get(c, 'TickLabels')));
end

% Letters on a transparent overlay spanning the canvas.
ov = axes(dst, 'Units','normalized', 'Position',[0 0 1 1], ...
          'Visible','off', 'XLim',[0 W], 'YLim',[0 H], 'HitTest','off');
for i = 1:n
    t = tiles(i,:);
    text(ov, t(1) + L.letterxy(i,1), t(2) + t(4) + L.letterxy(i,2), ...
         L.letters{i}, 'FontName','Arial', 'FontWeight','bold', ...
         'FontSize', L.lettersize, 'Units','data', ...
         'HorizontalAlignment','left', 'VerticalAlignment','top');
end

[outdir, name] = fileparts(stem);
if isempty(outdir), outdir = pwd; end
P.export(P, name, outdir, dst);
close(dst);
end

% ------------------------------------------------------------------------
function lab = minusLabels(lab)
if ~iscell(lab), lab = cellstr(lab); end
neg = ~cellfun(@isempty, regexp(lab, '^-[\d.]+$', 'once'));
lab(neg) = cellfun(@(s) ['$' s '$'], lab(neg), 'uni', 0);
end
