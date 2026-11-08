%RUN_ALL Regenerate every figure panel, every lettered figure, and the
%   regression checks. Runs from a cold interpreter in the repository root.
%   The article has nine figures, Figs. 1-6 in the main text and 7-9 in the
%   appendices. figureN writes the panels of Fig. N under output/figureN/ and
%   the lettered figure to output/figs/figN.pdf. Figures 3 and 4 share the
%   single-donor gain tables, and Figs. 5 and 6 share the two-donor table,
%   each cached in the figure's own folder. selftest runs last because it
%   clears the workspace.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'src')));

tstart = tic;
figure1;  close all
figure2;  close all
figure3;  close all
figure3;  close all
figure4;  close all
figure5;  close all
figure6;  close all
figure7;  close all
figure8;  close all
figure9; close all
fprintf('\nfigures complete in %.1f minutes\n', toc(tstart)/60);
selftest; close all
fprintf('\nrun_all complete\n');
