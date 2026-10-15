function [G, W] = linger_tables(prm, pgrid, cachefile, ggrid)
%LINGER_TABLES  Build or load the single-donor gain tables.
%
%   Returns the ridge table G and the gain surface W together, so psiTable is
%   never called with a surface built on a different grid or prm from the
%   ridge table beside it.
%
%   W is a struct, not a matrix. gainSurface returns
%   struct('p',pgrid,'gam',ggrid,'G',Gs,'T0',T0,'prm',prm) with Gs sized
%   numel(ggrid) by numel(pgrid), reliance down the rows. Check the grids
%   through W.p and W.gam, since size(W,1) is always 1.
%
%   The reliance grid is part of the table. psiTable selects a reliance off
%   ggrid, so its resolution sets the resolution of the selected gamma and of
%   g1. ggrid is passed explicitly and carried in the key.
%
%   The cache stamp is stored as tkey. Give this cache a filename distinct
%   from any other cache that stores G and W under a different stamp.
%
%   numel(W.p) must equal numel(pgrid), numel(W.gam) must equal numel(ggrid),
%   and G must be recoverable from W as the argmax of W.G over its first
%   dimension.

if nargin < 2 || isempty(pgrid),     pgrid = linspace(0, 1.2, 41); end
if nargin < 3 || isempty(cachefile), cachefile = fullfile('output', 'linger_tables_cache.mat'); end
if nargin < 4 || isempty(ggrid),     ggrid = linspace(0, prm.gmax, 51); end

want = tablekey(prm, pgrid, ggrid);
G = [];  W = [];
if isfile(cachefile)
    C = load(cachefile);
    if isfield(C,'G') && isfield(C,'W') && isfield(C,'tkey') && ...
            strcmp(C.tkey, want)
        G = C.G;  W = C.W;
        fprintf('linger_tables: reusing %s\n', cachefile);
    else
        fprintf('linger_tables: %s key mismatch, rebuilding\n', cachefile);
    end
end

if isempty(G)
    M = lingermodel();
    W = M.gainSurface(prm, pgrid, ggrid);
    G = M.gainTable(prm, pgrid);
    tkey = want; %#ok<NASGU>
    save(cachefile, 'G', 'W', 'tkey', 'prm', 'pgrid', 'ggrid');
    fprintf('linger_tables: built and saved %s\n', cachefile);
end

if ~isstruct(W) || ~all(isfield(W, {'p','gam','G'}))
    error('linger_tables:notASurface', ...
        ['W is not a gain surface struct.  gainSurface returns a struct ' ...
         'with fields p, gam and G; a bare matrix here means an old cache ' ...
         'or an old call site.']);
end
if numel(W.p) ~= numel(pgrid)
    warning('linger_tables:offsetGrid', ...
        ['the surface has %d offset points against %d grid points; the ' ...
         'surface and the ridge table were not built on the same grid'], ...
        numel(W.p), numel(pgrid));
end
if numel(W.gam) ~= numel(ggrid)
    warning('linger_tables:relianceGrid', ...
        ['the surface has %d reliance points against %d requested; the key ' ...
         'did not catch a change in the reliance grid'], ...
        numel(W.gam), numel(ggrid));
end
if size(W.G,1) ~= numel(W.gam) || size(W.G,2) ~= numel(W.p)
    warning('linger_tables:surfaceOrientation', ...
        ['W.G is %dx%d against %d reliances by %d offsets; psiTable indexes ' ...
         'W.G(j,:) as a function of offset at fixed reliance and will read ' ...
         'the transpose'], size(W.G,1), size(W.G,2), numel(W.gam), numel(W.p));
end
end

% ------------------------------------------------------------------------
function k = tablekey(prm, pgrid, ggrid)
VER = 'v2-reliance-outside-average';
f = fieldnames(prm);
k = VER;
for i = 1:numel(f)
    v = prm.(f{i});
    if isnumeric(v), k = [k sprintf('|%s=%.10g', f{i}, v(:).')]; end %#ok<AGROW>
end
k = [k sprintf('|p=%.6g:%.6g:%d', pgrid(1), pgrid(end), numel(pgrid))];
k = [k sprintf('|g=%.6g:%.6g:%d', ggrid(1), ggrid(end), numel(ggrid))];
end
