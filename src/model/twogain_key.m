function k = twogain_key(prm, pgrid, kmax, ggrid)
%TWOGAIN_KEY  Stamp for the two-donor gain table.
%
%   Carries a version tag as well as the parameters, so a convention change
%   that alters no parameter still invalidates the cache. Bump VER on any
%   change to what twogain computes, not only to what it is called with.

VER = 'v2-reliance-outside-average-signed-offsets';
f = fieldnames(prm);
s = VER;
for i = 1:numel(f)
    v = prm.(f{i});
    if isnumeric(v), s = [s sprintf('|%s=%.10g', f{i}, v(:).')]; end %#ok<AGROW>
end
s = [s sprintf('|p=%.6g:%.6g:%d', pgrid(1), pgrid(end), numel(pgrid))];
s = [s sprintf('|k=%d', kmax)];
s = [s sprintf('|g=%.6g:%.6g:%d', ggrid(1), ggrid(end), numel(ggrid))];
k = s;
end
