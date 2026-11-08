function S = twodonor_setup(outdir, fast, tag, needSingle)
%TWODONOR_SETUP  Two-donor arrival table and sweep shared by Figs. 5 and 6.
%   S = TWODONOR_SETUP(outdir, fast, tag, needSingle) builds or loads the exact
%   two-donor arrival table G2, caches it in outdir as [tag '_twogain_cache.mat']
%   behind a key that includes every grid and model quantity, and sweeps the
%   total recruitment benefit and the donor-specific optimum into res. With
%   needSingle it also builds the single-donor gain tables G, W and ps at the
%   same parameters.
if nargin < 4, needSingle = false; end
prm = struct('kbar',0.0314,'sig',0.8,'beta0',1,'sm',1,'v',1,'p0',0.75);
prm.gmax = prm.sm*prm.sig - 1e-6;
cfg0 = struct('N',120,'mutot',2.0,'mu1',1,'mu2',1);

if fast
    PGRID = linspace(-prm.sig,prm.sig,9);
    KMAX  = 8;
    GGRID = linspace(0,prm.gmax,7);
    rlist = [1 2 4 8 20 50];
else
    PGRID = linspace(-prm.sig,prm.sig,25);
    KMAX  = 25;
    GGRID = linspace(0,prm.gmax,13);
    rlist = [1 1.5 2 2.5 3 3.5 4 4.5 5 6 8 10 12 15 20 25 30 40 50];
end

want = twogain_key(prm,PGRID,KMAX,GGRID);
cachefile = fullfile(outdir,[tag '_twogain_cache.mat']);
G2 = [];
if isfile(cachefile)
    C = load(cachefile,'G2');
    if isfield(C,'G2') && isfield(C.G2,'key') && strcmp(C.G2.key,want)
        G2 = C.G2;
        fprintf('%s: loaded stamped two-donor gain cache\n', tag);
    end
end
if isempty(G2)
    G2 = twogain(prm,PGRID,KMAX,GGRID);
    save(cachefile,'G2');
end

M2 = twodonor();
res = M2.sweep(G2,cfg0,rlist);
save(fullfile(outdir,'two_donor_results.mat'),'res','G2');

S = struct('prm',prm,'cfg0',cfg0,'G2',G2,'res',res);
if needSingle
    A = afpt();
    M = lingermodel();
    prmS = A.mkprm('kbar',prm.kbar,'sig',prm.sig,'beta0',prm.beta0, ...
                   'sm',prm.sm,'v',prm.v);
    if fast
        pgridS = linspace(0,1.2,31);
        ggridS = linspace(0,prmS.gmax,25);
    else
        pgridS = linspace(0,1.2,61);
        ggridS = linspace(0,prmS.gmax,51);
    end
    G = M.gainTable(prmS,pgridS);
    W = M.gainSurface(prmS,pgridS,ggridS);
    S.G = G;
    S.W = W;
    S.ps = M.psiTable(G,prm.p0,W);
end
end
