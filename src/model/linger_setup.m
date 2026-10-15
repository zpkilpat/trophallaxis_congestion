function S = linger_setup(outdir, fast, tag)
%LINGER_SETUP  Single-donor model and gain tables shared by Figs. 4 and 5.
%   S = LINGER_SETUP(outdir, fast, tag) returns the afpt and lingermodel
%   objects, the parameter structs prm and cfg, and the gain tables G, W and
%   psList (one table per offset in P0LIST). The tables are cached in outdir as
%   [tag '_gain_cache.mat'] behind a key built from the parameters and grids.
%   The grids are part of the key, so they are chosen before the stamp, or a
%   fast run's coarse tables could be loaded by a production call.
VER    = 'v3-pcap-bisected';
P0LIST = [0.75 1.13 1.52];

A = afpt();
M = lingermodel();
prm = A.mkprm('kbar',0.0314,'sig',0.8,'beta0',1,'sm',1,'v',1);
cfg = M.mkcfg('N',120,'mu',1,'lam0',prm.lam0,'p0',0.75);

if fast
    pgrid = linspace(0,1.2,31);
    ggrid = linspace(0,prm.gmax,25);
else
    pgrid = linspace(0,1.2,61);
    ggrid = linspace(0,prm.gmax,51);
end
cachefile = fullfile(outdir,[tag '_gain_cache.mat']);
key = M.stamp(VER, prm, cfg, P0LIST, pgrid, ggrid);
haveCache = false;
if isfile(cachefile)
    C = load(cachefile);
    if isfield(C,'key') && strcmp(C.key,key)
        G = C.G; W = C.W; psList = C.psList;
        haveCache = true;
        fprintf('%s: loaded stamped gain cache (%d x %d)\n', tag, ...
                numel(pgrid), numel(ggrid));
    end
end
if ~haveCache
    fprintf('%s: building gain tables on a %d x %d grid\n', tag, ...
            numel(pgrid), numel(ggrid));
    G = M.gainTable(prm,pgrid);
    W = M.gainSurface(prm,pgrid,ggrid);
    psList = arrayfun(@(p0) M.psiTable(G,p0,W),P0LIST,'uni',0);
    save(cachefile,'G','W','psList','key');
end

fprintf('%s: lambda0 = %.5f, rho = %.3f\n', tag, prm.lam0, prm.lam0*cfg.N/cfg.mu);
fprintf('%s checks: p*/sigma = %.3f, pcap/sigma = %.3f, g1 =', tag, ...
        G.pstar/prm.sig, G.pcap/prm.sig);
fprintf(' %.3f',cellfun(@(s)s.g1,psList));
fprintf('\n');

S = struct('A',A,'M',M,'prm',prm,'cfg',cfg,'G',G,'W',W,'psList',{psList}, ...
           'ps',psList{1},'P0LIST',P0LIST);
end
