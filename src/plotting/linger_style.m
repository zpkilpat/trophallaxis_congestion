function P = linger_style()
%LINGER_STYLE  House figure style and export.
%   LaTeX interpreter throughout, no FontName set. Panels are standalone, one
%   file each, with no panel letters. Each figure driver assembles its
%   panels into one lettered figure with assemble_figure.
%
%   Every panel is exported at the same canvas size, colorbar or not. A
%   colorbar takes its room from the axes, so all panels share one outer
%   footprint and one effective type size after assembly.
%
%   Every export writes three files with a common stem, NAME.pdf (vector),
%   NAME.png (raster, P.DPI), and NAME.fig (editable MATLAB figure).
%
%   Color convention for p_0. Solid, dashed, dash-dot is p_0 = 0.75, 1.13,
%   1.52 in every panel that carries the family, darkest at the smallest
%   offset. Use P.p0cols and P.p0sty rather than indexing P.greens by hand.
%
%   CLIM_SPEEDUP is a floor for the Fig. 3(c) color limit. Fig. 3(c) takes
%   max(floor, ceil(data)) so the largest speedup is never clipped.

P.LW   = 3.5;
P.FSZ  = 26;
P.LSZ  = P.FSZ + 8;
P.ALW  = 2;
P.FIG  = [22 19];          % centimeters, identical for every panel
P.DPI  = 400;
P.CBSHRINK    = 0.82;      % axes width retained when a colorbar is added
P.CLIM_SPEEDUP = [0 60];   % percent, a floor for Fig. 3(c)

P.greens = [0.86 0.94 0.87; 0.66 0.85 0.70; 0.44 0.74 0.53;
            0.20 0.55 0.32; 0.13 0.41 0.24; 0.06 0.27 0.15];
P.blues  = [0.72 0.83 0.92; 0.36 0.58 0.80; 0.13 0.31 0.55];
P.mid    = P.greens(4,:);
P.grey   = [0.40 0.40 0.40];

% One ordering for the p_0 family, darkest first.
P.p0cols = [P.greens(6,:); P.greens(5,:); P.greens(4,:)];
P.p0sty  = {'-','--','-.'};

P.axes    = @apply_axes;
P.figure  = @new_figure;
P.export  = @export_panel;
P.reexport = @export_from_fig;
P.legend  = @legend_latex;
P.label   = @label_latex;
P.plasma  = @plasma;
P.divmap  = @divmap;
P.pivotmap = @pivotmap;
P.cbar    = @cbar_latex;
end

% ------------------------------------------------------------------------
function f = new_figure(P, ~)
% The optional second argument is accepted and ignored.
f = figure('Units','centimeters','Position',[2 2 P.FIG(1) P.FIG(2)], ...
           'Color','w');
end

function apply_axes(P, ax)
if nargin < 2, ax = gca; end
set(ax,'FontSize',P.FSZ,'LineWidth',P.ALW,'TickDir','out','Box','off', ...
       'TickLabelInterpreter','latex','Layer','top');
end

function lg = legend_latex(P, ax, hs, strs, loc)
%LEGEND_LATEX  Boxless LaTeX legend at axis font size.
if nargin < 5 || isempty(loc), loc = 'best'; end
lg = legend(ax, hs, strs, 'Interpreter','latex', 'FontSize', P.FSZ, ...
            'Location', loc, 'Box','off');
end

function t = label_latex(P, ax, x, y, str, varargin)
%LABEL_LATEX  In-axes LaTeX annotation at axis font size.
t = text(ax, x, y, str, 'Interpreter','latex', 'FontSize', P.FSZ, varargin{:});
end

function export_panel(P, name, outdir, h)
%EXPORT_PANEL  Write NAME.pdf, NAME.png and NAME.fig for figure handle H.
if nargin < 4 || isempty(h)
    h = gcf;
    warning('linger_style:exportGcf', ...
        'export of %s used gcf; pass the figure handle', name);
end
if nargin < 3 || isempty(outdir) || ~isfolder(outdir)
    fprintf('outdir missing, skipping export of %s\n', name);  return
end
stem = fullfile(outdir, name);
% The editable figure is written first. PDF and PNG are then exported from a
% fresh invisible copy of that file, which is robust in batch mode, where a
% live figure can be invalidated by its first export.
drawnow;
savefig(h, [stem '.fig']);
export_from_fig(P, stem);
fprintf('wrote %s.fig / .pdf / .png\n', name);
end

function export_from_fig(P, stem)
%EXPORT_FROM_FIG  Write STEM.pdf and STEM.png from STEM.fig.
f = openfig([stem '.fig'], 'invisible');
c = onCleanup(@() close(f));
if isappdata(f, 'ExportCanvasSize')
    canvas = getappdata(f, 'ExportCanvasSize');
    set(f, 'Units','centimeters', 'Position',[1 1 canvas]);
end
exportgraphics(f, [stem '.pdf'], 'ContentType','vector');
% Without a display the OpenGL raster drops lines, so the PNG is rendered
% with painters there, which matches the vector PDF.
try
    if ~usejava('desktop')
        print(f, [stem '.png'], '-dpng', sprintf('-r%d', P.DPI), '-painters');
    else
        exportgraphics(f, [stem '.png'], 'Resolution', P.DPI);
    end
catch err
    warning('linger_style:png', 'PNG export of %s failed: %s', stem, err.message);
end
end

function cb = cbar_latex(P, str)
%CBAR_LATEX  Colorbar with a LaTeX label. The interpreter is set before the
%   string so that LaTeX-only macros parse. The axes width is then set so a
%   colorbar panel and a plain panel share one outer footprint.
ax  = gca;
set(ax, 'Units', 'normalized');
pos = get(ax, 'Position');
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
cb.LineWidth = P.ALW;
cb.Label.Interpreter = 'latex';
cb.Label.String = str;
cb.Label.FontSize = P.LSZ;
set(ax, 'Position', [pos(1) pos(2) pos(3)*P.CBSHRINK pos(4)]);
end

% ------------------------------------------------------------------------
function C = plasma(n)
if nargin < 1, n = 256; end
A = [0.050 0.030 0.528; 0.294 0.011 0.631; 0.472 0.040 0.658;
     0.635 0.157 0.590; 0.767 0.290 0.474; 0.874 0.428 0.358;
     0.954 0.582 0.246; 0.993 0.755 0.152; 0.940 0.975 0.131];
C = interp1(linspace(0,1,size(A,1)), A, linspace(0,1,n));
end

function C = divmap(n)
%DIVMAP  Blue-white-red for data symmetric about a meaningful center.
if nargin < 1, n = 256; end
A = [0.13 0.31 0.55; 0.36 0.58 0.80; 0.78 0.86 0.94;
     1.00 1.00 1.00; 0.97 0.80 0.70; 0.87 0.45 0.32; 0.60 0.11 0.09];
C = interp1(linspace(0,1,size(A,1)), A, linspace(0,1,n));
end

function C = pivotmap(n, lo, hi, piv)
%PIVOTMAP  Diverging map whose neutral tone sits at PIV, which need not be the
%   midpoint of [LO,HI]. Warm below the pivot, cool above. Pass LO, HI, PIV in
%   the coordinate used for clim, for instance log2 of a gain so that a factor
%   of two either way gets equal visual weight.
if nargin < 1 || isempty(n), n = 256; end
f  = min(max((piv - lo)/(hi - lo), 0.02), 0.98);
n1 = max(2, round(f*n));  n2 = max(2, n - n1);
warm = [0.60 0.11 0.09; 0.87 0.45 0.32; 0.97 0.80 0.70; 1.00 1.00 1.00];
cool = [1.00 1.00 1.00; 0.78 0.86 0.94; 0.36 0.58 0.80; 0.13 0.31 0.55;
        0.03 0.16 0.32];
C = [interp1(linspace(0,1,size(warm,1)), warm, linspace(0,1,n1)); ...
     interp1(linspace(0,1,size(cool,1)), cool, linspace(0,1,n2))];
end
