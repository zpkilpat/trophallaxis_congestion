function scale_figure_fonts(K)
%SCALE_FIGURE_FONTS  Rescale the text in the lettered figures without rerunning them.
%   scale_figure_fonts(K) reopens output/figs/figN.fig (written by figureN),
%   multiplies tick, legend and in-panel text sizes by K (default 1.07), sets
%   all in-panel annotations to one common size first, leaves axis labels and
%   the bold panel letters alone, and re-exports figN.pdf and figN.png. Only
%   the saved figure files are touched, so no model code is recomputed.
%
%   Run from the repository root after figure1 ... figure9 have produced
%   their .fig files. Figure 1 is a hand-built schematic and is skipped.
if nargin < 1, K = 1.07; end
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')));
P = linger_style();
files = dir(fullfile('output','figs','fig*.fig'));
for k = 1:numel(files)
    [~, stem] = fileparts(files(k).name);
    if strcmp(stem, 'fig1'), continue; end
    path = fullfile(files(k).folder, stem);
    f = openfig([path '.fig'], 'invisible');
    keep = [];
    for a = findall(f,'Type','axes')'
        keep = [keep; a.XLabel; a.YLabel; a.Title]; %#ok<AGROW>
    end
    for c = findall(f,'Type','colorbar')'
        keep = [keep; c.Label]; %#ok<AGROW>
    end
    for o = findall(f, '-property', 'FontSize')'
        if strcmp(get(o,'Type'),'figure') || any(o == keep), continue; end
        s = get(o,'FontSize');
        if strcmp(get(o,'Type'),'text')
            if strcmp(get(o,'FontName'),'Arial') && strcmp(get(o,'FontWeight'),'bold') && s > 30
                continue                       % panel letters
            end
            if s >= 22 && s <= 31, s = 26; end % common in-panel annotation size
        end
        set(o, 'FontSize', s*K);
    end
    savefig(f, [path '.fig']);  close(f);
    P.reexport(P, path);
    fprintf('rescaled %s by %.2f\n', stem, K);
end
end
