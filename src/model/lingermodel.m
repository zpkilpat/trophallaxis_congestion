function M = lingermodel()
%LINGERMODEL  Recruitment and congestion at a single donor, built on afpt.
%
%   Compartments, conserved as S+U+F+D = N
%       S  hungry searchers in the bulk
%       U  hungry bees in the crowd at the donor
%       F  fed bees still lingering (the scenters)
%       D  fed bees that have abandoned, permanently out of the system
%   Random-order service. The donor engages one crowd member, drawn uniformly,
%   at rate mu, so the productive rate is mu*U/(U+F). Fed bees abandon at
%   rate delta. There is no crowd capacity. Conservation and the single server
%   supply the congestion cost, since bees held at the donor cannot search.
%
%   The cue is produced by the fed lingerers F alone. The crowd U does not
%   scent, here or in the ABM.
%
%   Reliance is chosen before the offset is drawn. A bee sets one reliance,
%   and the peak is redrawn by the scenting population each time a bee is
%   served or abandons, so the max over gamma is taken outside the average.
%   psiTable needs the full surface W = gainSurface(...) and returns the
%   selected reliance in ps.gam. The ABM uses the same ps.gam, so psiTable and
%   abm must follow one convention.
%
%   Invariants
%     Tq(delta -> inf) must equal log(1/(1-q))/lam0 to five digits.
%     S must be nonincreasing in the ABM; an abandoned bee never searches again.
%     psiTable g1 must equal max_gamma of the k=1 row of rk, not psi at Fbar=1.
%     psiTable must be nondecreasing in Fbar and bounded above by G(0).
%     optDelta must never return a Topt above the best value on its own scan.

M.gainTable   = @gainTable;
M.gainSurface = @gainSurface;
M.psiTable    = @psiTable;
M.Tq          = @Tq;
M.traj        = @traj;
M.abm         = @abm;
M.advance     = @advance;   % exposed for the SI service-clock model, abm_svc
M.mkcfg       = @mkcfg;
M.dlogPsi     = @dlogPsi;
M.psiAt       = @psiAt;
M.gamAt       = @gamAt;
M.optDelta    = @optDelta;
M.refineDelta = @refineDelta;
M.criterion   = @criterion;
M.stamp       = @stamp;
end

% ------------------------------------------------------------------------
function cfg = mkcfg(varargin)
cfg = struct('N',120,'mu',1.0,'lam0',0.005,'p0',0.75);
for k = 1:2:numel(varargin), cfg.(varargin{k}) = varargin{k+1}; end
end

function s = stamp(varargin)
%STAMP  Fingerprint of everything a cached result depends on.  Compared as a
%   string, so a cache built at other parameters is rejected. Numeric struct
%   fields are taken in sorted field order.
%
%   Callers also pass a version tag, so a convention change that alters no
%   parameter still invalidates the cache.
parts = {};
for k = 1:numel(varargin)
    v = varargin{k};
    if isstruct(v)
        f = sort(fieldnames(v));
        for i = 1:numel(f)
            x = v.(f{i});
            if isnumeric(x) || islogical(x)
                parts{end+1} = sprintf('%s=%.12g,', f{i}, double(x(:)).'); %#ok<AGROW>
            end
        end
    elseif isnumeric(v) || islogical(v)
        parts{end+1} = sprintf('%.12g,', double(v(:)).');                  %#ok<AGROW>
    elseif ischar(v)
        parts{end+1} = [v ',']; %#ok<AGROW>
    end
end
s = strjoin(parts, '|');
end

function G = gainTable(prm, pgrid)
%GAINTABLE  gamma*(p) and G(p) = Tbar(0)/Tbar(gamma*,p).  Tabulate once.
%   Graw is the unclipped ratio and the clipped G is what the model consumes.
%   atcap flags the p at which the optimum sits on the positivity bound
%   gamma = sm*sig, which is a result and not a truncation.
%
%   The ridge in Fig. 2(c) is drawn from this table. psiTable averages the full
%   surface W instead (see the header).
%
%   pstar is not read from this table. Past pstar the optimizer picks
%   gamma = 0, so Graw equals 1 to rounding and a table crossing is
%   unreliable. pstar is bisected on a predicate instead (see pstarBisect).
A = afpt();
if nargin < 2, pgrid = linspace(0, 1.2, 61); end
T0 = A.moments(max(pgrid)+2, prm.sig, prm.beta0, 0, prm.sm, prm.theta, prm.v);
np = numel(pgrid); gam = zeros(1,np); Gr = zeros(1,np); cap = false(1,np);
for i = 1:np
    obj = @(g) A.moments(pgrid(i), prm.sig, prm.beta0, g, prm.sm, prm.theta, prm.v);
    [gs, fs] = fminbnd(obj, 0, prm.gmax, optimset('TolX',1e-5));
    f0 = obj(0);
    if f0 <= fs, gs = 0; fs = f0; end
    gam(i) = gs;  Gr(i) = T0/fs;  cap(i) = gs > prm.gmax - 1e-3;
end
G = struct('p',pgrid,'gam',gam,'G',max(Gr,1),'Graw',Gr,'atcap',cap, ...
           'T0',T0,'prm',prm);
G.pstar = pstarBisect(prm, T0, min(pgrid), max(pgrid));
if ~any(cap)
    G.pcap = NaN;
else
    G.pcap = pcapBisect(prm, min(pgrid), max(pgrid));
end
j = find(pgrid < G.pstar, 1, 'last');       % how gamma* leaves: a large value
if isempty(j), G.gamBelow = NaN; else, G.gamBelow = gam(j); end   % means a jump
end

function p = pcapBisect(prm, plo, phi)
%PCAPBISECT  Offset at which gamma*(p) leaves the positivity bound.
%   Bisects the predicate "the optimum still sits on gamma = sm*sig", which is
%   monotone in p.
%
%   Read off the grid, pcap would quantize to the grid spacing, and the
%   interior-reliance window (pstar - pcap)/sig would depend on grid size.
%
%   The tolerance must match the one gainTable uses to set atcap, or the
%   bisected pcap and the flagged grid points disagree at the boundary.
A   = afpt();
tol = 1e-3;
atc = @(p) capAt(A, prm, p, tol);
if ~atc(plo) || atc(phi), p = NaN; return; end
for it = 1:60
    pm = 0.5*(plo + phi);
    if atc(pm), plo = pm; else, phi = pm; end
    if phi - plo < 1e-8, break; end
end
p = 0.5*(plo + phi);
end

function tf = capAt(A, prm, p, tol)
%CAPAT  Does the unconstrained optimum sit on the positivity bound at p.
obj = @(g) A.moments(p, prm.sig, prm.beta0, g, prm.sm, prm.theta, prm.v);
[gs, fs] = fminbnd(obj, 0, prm.gmax, optimset('TolX',1e-5));
f0 = obj(0);
if f0 <= fs, gs = 0; end
tf = gs > prm.gmax - tol;
end

function p = pstarBisect(prm, T0, plo, phi)
%PSTARBISECT  Offset at which biasing stops paying at all.
%   Bisects the predicate "some gamma > 0 beats gamma = 0", which is monotone
%   in p and converges whether gamma*(p) jumps to zero or slides to it. Past
%   pstar the gain table is flat at 1 by construction.
%
%   The predicate compares against T0, which equals obj(0) at every p, since
%   gamma = 0 removes the plume from the dynamics.
%
%   If the transition is tangential, the gain deficit near pstar goes like
%   (pstar-p)^2 and pstar is good to about 1e-4, set by the 1e-8 threshold.
%   If G.gamBelow is large, the transition is a jump and pstar is good to the
%   bisection tolerance.
A = afpt();
helps = @(p) innerMin(A, prm, p) < T0*(1 - 1e-8);
if ~helps(plo) || helps(phi), p = NaN; return; end
for it = 1:60
    pm = 0.5*(plo + phi);
    if helps(pm), plo = pm; else, phi = pm; end
    if phi - plo < 1e-8, break; end
end
p = 0.5*(plo + phi);
end

function f = innerMin(A, prm, p)
% Best strictly biased strategy. The lower bound is 1e-6 rather than 0 so
% this is the biased branch alone and never collapses onto obj(0).
% fminbnd's first output is the minimizer and the second is the minimum.
% The time is the second output.
obj = @(g) A.moments(p, prm.sig, prm.beta0, g, prm.sm, prm.theta, prm.v);
[~, f] = fminbnd(obj, 1e-6, prm.gmax, optimset('TolX',1e-7));
end

function W = gainSurface(prm, pgrid, ggrid)
%GAINSURFACE  G(p,gamma) over the whole strategy plane.
%   Used by the Fig. 2(c) heatmap and by psiTable, which needs off-ridge values
%   because the reliance it selects is generally not optimal at the realized
%   offset. The gamma grid sets the resolution of the selected reliance. 51
%   points over [0,gmax] resolve it to about 0.016, finer than the
%   interior-reliance window. Refine here, not in psiTable.
A = afpt();
if nargin < 2, pgrid = linspace(0, 1.2, 61); end
if nargin < 3, ggrid = linspace(0, prm.gmax, 51); end
T0 = A.moments(max(pgrid)+2, prm.sig, prm.beta0, 0, prm.sm, prm.theta, prm.v);
Gs = zeros(numel(ggrid), numel(pgrid));
for i = 1:numel(pgrid)
    for j = 1:numel(ggrid)
        Gs(j,i) = T0/A.moments(pgrid(i), prm.sig, prm.beta0, ggrid(j), ...
                               prm.sm, prm.theta, prm.v);
    end
end
W = struct('p',pgrid,'gam',ggrid,'G',Gs,'T0',T0,'prm',prm);
end

function ps = psiTable(G, p0, W, Fgrid)
%PSITABLE  Psi~(Fbar), the recruitment gain averaged over a stochastic
%   scenting population, with the reliance chosen once rather than per offset.
%
%       Psi~(Fbar) = max_gamma [ P(0|Fbar)*1
%                                + sum_{k>=1} P(k|Fbar) E_p[ G(gamma,|p|) ] ]
%
%   with F ~ Poisson(Fbar) and, given k >= 1 scenters each offset by
%   N(0,p0^2), aggregate peak p ~ N(0,p0^2/k). G(gamma,p) = T0/Tbar(gamma,p)
%   is a rate ratio, so the average is taken on G and not on Tbar, since
%   arrival rates from independent starts add and times do not.
%
%   The k = 0 term is 1, not G(gamma, large p). With no scenters there is no
%   plume, so the searcher is unbiased whatever its gamma.
%
%   Offsets are clamped to max(W.p) rather than to the support, so the
%   sub-unity floor beyond p = sigma is counted.
%
%   The F grid is linear where the curvature lives and logarithmic in the
%   tail, because linear interpolation on a purely log grid puts kinks in the
%   ODE right-hand side near F ~ 5.
if nargin < 4 || isempty(Fgrid)
    Fgrid = unique([0, linspace(0.02, 25, 250), logspace(log10(26), log10(150), 40)]);
end
if nargin < 3 || isempty(W)
    error('lingermodel:psiTableNeedsSurface', ...
        ['psiTable now needs the gain surface W as its third argument, so ' ...
         'that the reliance can be chosen before the offset is drawn. ' ...
         'Pass W = M.gainSurface(prm, ...); see the header.']);
end

ks = 1:400;
z  = linspace(-6, 6, 1201);
w  = exp(-z.^2/2);  w = w/trapz(z, w);
ng = numel(W.gam);
pmax = max(W.p);

% rk(j,i) = E[ G(gamma_j, |p|) ] with i scenters.  One interp1 per (gamma,k).
rk = ones(ng, numel(ks));
for i = 1:numel(ks)
    pp = min(abs(z*p0/sqrt(ks(i))), pmax);
    for j = 1:ng
        rk(j,i) = trapz(z, w .* interp1(W.p, W.G(j,:), pp, 'linear'));
    end
end

vals = ones(size(Fgrid));  gsel = zeros(size(Fgrid));
for i = 1:numel(Fgrid)
    Fb = Fgrid(i);
    if Fb <= 0, continue; end
    pk  = exp(-Fb + ks*log(Fb) - gammaln(ks+1));    % P(k|Fbar), k >= 1
    tot = exp(-Fb) + rk*pk(:);                      % ng x 1, k = 0 scores 1
    [vals(i), jb] = max(tot);
    gsel(i) = W.gam(jb);
end

% Slope at the origin.  d/dFbar of [e^-F + F e^-F r_1 + O(F^2)] at F = 0 is
% r_1 - 1, and the maximizer there maximizes r_1 alone, so g1 is a max over
% gamma of the k = 1 row, not rk(jsel,1) at some Fbar.
[g1, jg1] = max(rk(:,1));

if isfinite(G.pstar), Fsub = (p0/G.pstar)^2; else, Fsub = NaN; end
if any(diff(vals) < -1e-9)
    warning('lingermodel:psiNotMonotone', ...
        'Psi~ decreases somewhere in Fbar; check the gamma grid in W');
end
ps = struct('F',Fgrid,'psi',vals,'gam',gsel,'p0',p0, ...
            'g1',   g1, ...                         % single-scenter gain
            'gam1', W.gam(jg1), ...                 % reliance that attains it
            'psi1', interp1(Fgrid,vals,1), ...      % Poisson average at Fbar=1
            'psitop', vals(end), ...                % Psi~ at the top of the grid
            'Fsub', Fsub, ...                       % substitution threshold
            'G0',   max(W.G(:,1)));                 % true ceiling G(0)
end

function v = psiAt(ps, F)
v = interp1(ps.F, ps.psi, min(max(F,0), max(ps.F)), 'linear');
end

function g = gamAt(ps, F)
%GAMAT  Reliance a searcher sets given the current scenting population.  The
%   ABM uses it so the reduced system and the simulation share one
%   convention. It depends on how many bees are scenting, not on where their
%   plume peak landed.
g = interp1(ps.F, ps.gam, min(max(F,0), max(ps.F)), 'linear');
end

function d = dlogPsi(ps, F)
g = gradient(log(ps.psi), ps.F);
d = interp1(ps.F, g, min(max(F,0), max(ps.F)), 'linear');
end

% ------------------------------------------------------------------------
function [lhs, rhs, pays] = criterion(ps, cfg, q)
%CRITERION  Does lingering pay at all?  The large-delta expansion is
%
%     T_q(delta) = L_q/lam0 + (1/delta)*[ A_q/(mu - A_q) - (g1-1)*q*N ] + ...
%     A_q = lam0*(1-q)*N,
%
%   so lingering pays when (g1-1)*q*N*(mu - A_q) > A_q. The benefit
%   coefficient is the number of bees fed over the run, and the cost is
%   evaluated at the end of the run, where the queue must still clear.
%
%   A necessary condition is rho < 1/(1-q).
if nargin < 3, q = 0.5; end
Aq   = cfg.lam0*(1-q)*cfg.N;
lhs  = (ps.g1 - 1)*q*cfg.N*(cfg.mu - Aq);
rhs  = Aq;
pays = (cfg.mu > Aq) && (lhs > rhs);
end

% ------------------------------------------------------------------------
function [T, sol] = Tq(delta, ps, cfg, q, tmax)
%TQ  Time to feed a fraction q of the colony, from the reduced system.
%   The crowd is seeded at U = 1e-6 and the productive rate is regularized,
%   because mu*U/(U+F) is discontinuous at the origin. The limit along
%   U -> 0+ with F = 0 is mu, not 0, and an unregularized guard puts a jump
%   at t = 0 that the integrator handles differently at each delta.
if nargin < 5, tmax = 2e5; end
N = cfg.N;
opt = odeset('RelTol',1e-8,'AbsTol',1e-11,'MaxStep',tmax/1e4, ...
             'Events',@(t,y) ev(t,y,q,N));
sol = ode15s(@(t,y) rhs(t,y,delta,ps,cfg), [0 tmax], [N-1e-6; 1e-6; 0; 0], opt);
if isempty(sol.xe), T = Inf; else, T = sol.xe(1); end
end

function tr = traj(delta, ps, cfg, q, nt)
%TRAJ  The same integration as Tq, resampled, for the crowd panel.  The
%   figure and the reported T_q come from one solve.
if nargin < 5, nt = 4000; end
[T, sol] = Tq(delta, ps, cfg, q);
if ~isfinite(T), tr = struct('T',Inf); return; end
t = linspace(0, T, nt).';
Y = deval(sol, t).';
S = max(Y(:,1),0);  U = max(Y(:,2),0);  F = max(Y(:,3),0);
tr = struct('t',t,'S',S,'U',U,'F',F,'D',Y(:,4), 'crowd',U+F, ...
            'load', cfg.lam0*arrayfun(@(f) psiAt(ps,f), F).*S/cfg.mu, ...
            'T',T,'delta',delta,'q',q);
end

function dy = rhs(~, y, delta, ps, cfg)
S = max(y(1),0);  U = max(y(2),0);  F = max(y(3),0);
a   = cfg.lam0 * psiAt(ps, F) * S;
sig = cfg.mu * U / (U + F + 1e-9);
dy  = [-a; a - sig; sig - delta*F; delta*F];
end

function [val, isterm, dir] = ev(~, y, q, N)
val = (y(3) + y(4)) - q*N;  isterm = 1;  dir = 1;
end

% ------------------------------------------------------------------------
function [dopt, Topt] = refineDelta(f, dgrid, T)
%REFINEDELTA  Parabolic refinement of a grid argmin, in log10(delta), given
%   the grid values T = f(dgrid). Refining avoids staircase curves in delta*.
%
%   Topt <= min(T). A refinement worse than its own bracket means the bracket
%   is inverted, and it is rejected.
T(~isfinite(T)) = Inf;
[Tg, k] = min(T);
dopt = dgrid(k);  Topt = Tg;
if ~isfinite(Tg) || k == 1 || k == numel(dgrid), return; end
[u, Tu] = fminbnd(@(u) f(10^u), log10(dgrid(k-1)), log10(dgrid(k+1)), ...
                  optimset('TolX',1e-4));
if Tu <= Tg, dopt = 10^u;  Topt = Tu; end
end

function [dopt, Topt, Tinf] = optDelta(ps, cfg, q, dgrid)
%OPTDELTA  Optimal abandonment rate and the time it attains.  Returns dopt =
%   NaN and Topt = Tinf when no delta beats immediate departure, so callers
%   can report a speedup of exactly zero rather than a small negative one.
if nargin < 4, dgrid = logspace(-2.5, 1.0, 44); end
f = @(d) Tq(d, ps, cfg, q);
Tinf = f(1e4);
T = arrayfun(f, dgrid);
[dopt, Topt] = refineDelta(f, dgrid, T);
if ~isfinite(Topt) || Topt >= Tinf
    dopt = NaN;  Topt = Tinf;
end
end

% ------------------------------------------------------------------------
function [T, tr] = abm(delta, G, ps, cfg, q, seed)
%ABM  Full stochastic simulation.  Searchers move on the ring by exact
%   thinning with analytic trap crossings.  Each fed bee carries its own scent
%   offset ~ N(0,p0^2) and the plume peak is the empirical mean over the
%   current lingerers, so the Poisson-Gaussian averaging of psiTable emerges
%   rather than being imposed. Crowd events are drawn by thinning against
%   mu + delta*N.
%
%   A searcher sets gamma from how many bees are scenting, via ps.gam, not
%   from where their plume peak landed. This matches psiTable.
%
%   Offsets are drawn with ps.p0, not cfg.p0, so the simulation matches the
%   psiTable it is compared against.
%
%   A bee that abandons is fed and leaves the system. It does not return to
%   the bulk, so S is nonincreasing.
%
%   Crowd events resolve at the end of each thinning interval, so a bee
%   arriving mid-interval can be served at most one interval late. Intervals
%   are ~1/(mu+delta*N), so the bias is O(dt^2).
%
%   A run that hits the time cap returns NaN rather than the cap, so averages
%   over seeds are not silently truncated.
%
%   With two outputs, the compartment counts are recorded once per thinning
%   interval and returned in tr. Recording draws no random numbers and
%   changes no state, so T is identical either way.
%
%   tr.F is the integer number of bees scenting. Draw it with stairs, since
%   the count is piecewise constant between events.
%
%   The record holds one row at t = 0 and one per thinning interval, so its
%   length is roughly T*(mu+delta*N), about two thousand rows per realization
%   at these parameters. A capped run still returns its record.
if nargin < 6, seed = 0; end
rng(seed);
prm = G.prm;  N = cfg.N;  TWOPI = 2*pi;  p0 = ps.p0;
x = -pi + TWOPI*rand(N,1);
s = 2*(rand(N,1) > 0.5) - 1;
inBulk = true(N,1);
U = 0; off = zeros(1,0); nFed = 0; t = 0;
Rmax = cfg.mu + delta*N;
rec = nargout > 1;
if rec
    R = zeros(8192, 6);                 % t, S, U, F, nFed, gamma
    R(1,:) = [0, N, 0, 0, 0, 0];
    nR = 1;
end
while nFed < q*N && t < 1e7
    k = numel(off);
    if k == 0
        ppk = Inf;  gam = 0;            % no plume, so no bias, whatever gamma
    else
        ppk = mean(off);
        gam = gamAt(ps, k);
    end
    tau = -log(rand)/Rmax;
    ib = find(inBulk);
    if ~isempty(ib)
        [x(ib), s(ib), tabs] = advance(x(ib), s(ib), tau, ppk, gam, prm);
        got = ~isnan(tabs);
        if any(got), inBulk(ib(got)) = false;  U = U + sum(got); end
    end
    t = t + tau;
    r = rand*Rmax;
    if r < cfg.mu
        n = U + numel(off);
        if n > 0 && rand < U/n
            U = U - 1;  nFed = nFed + 1;
            off(end+1) = p0*randn;                 %#ok<AGROW>
        end
    elseif r < cfg.mu + delta*numel(off)
        off(randi(numel(off))) = [];               % lingerer leaves, and is done
    end
    if rec
        nR = nR + 1;
        if nR > size(R,1), R(2*size(R,1), 6) = 0; end   % double, zero filled
        R(nR,:) = [t, sum(inBulk), U, numel(off), nFed, gam];
    end
end
if nFed < q*N
    warning('lingermodel:abmCap', ...
        'ABM hit the time cap at delta = %g, seed %d, with %d of %d fed', ...
        delta, seed, nFed, ceil(q*N));
    T = NaN;
else
    T = t;
end
if rec
    R  = R(1:nR,:);
    tr = struct('t',R(:,1),'S',R(:,2),'U',R(:,3),'F',R(:,4), ...
                'nFed',R(:,5),'gam',R(:,6),'crowd',R(:,3)+R(:,4), ...
                'D',R(:,5)-R(:,4), ...              % fed and gone
                'T',T,'delta',delta,'q',q,'seed',seed,'p0',p0, ...
                'capped',~isfinite(T));
    % S+U+F+D = N at every recorded instant. A violation means a bee has been
    % created or has re-entered the bulk.
    bad = abs(tr.S + tr.U + tr.F + tr.D - N) > 0.5;
    if any(bad)
        warning('lingermodel:abmConservation', ...
            'S+U+F+D departs from N at %d of %d recorded instants', ...
            sum(bad), nR);
    end
end
end

function [x, s, tabs] = advance(x, s, tau, ppk, gam, prm)
TWOPI = 2*pi;  n = numel(x);
tabs = nan(n,1);  el = zeros(n,1);  live = true(n,1);
bmax = prm.beta0*(1 + gam/(prm.sm*prm.sig));
logq = log1p(-prm.theta);
for step = 1:10000
    idx = find(live);
    if isempty(idx), break; end
    dt = -log(rand(numel(idx),1))/bmax;
    rem = tau - el(idx);
    partial = dt >= rem;  dt(partial) = rem(partial);
    d  = prm.v*dt;  xa = x(idx);  sa = s(idx);
    u1 = zeros(size(xa));
    u1(sa > 0) = mod(-xa(sa > 0), TWOPI);
    u1(sa < 0) = mod( xa(sa < 0), TWOPI);
    u1(u1 == 0) = TWOPI;
    nc = zeros(size(d));  m = d >= u1;
    nc(m) = floor((d(m) - u1(m))/TWOPI) + 1;
    K = floor(log(rand(numel(idx),1))/logq);
    hit = K < nc;
    if any(hit)
        j = idx(hit);
        tabs(j) = el(idx(hit)) + (u1(hit) + TWOPI*K(hit))/prm.v;
        live(j) = false;
    end
    srv = ~hit;
    if any(srv)
        j  = idx(srv);
        xn = mod(xa(srv) + sa(srv).*d(srv) + pi, TWOPI) - pi;
        x(j) = xn;  el(j) = el(idx(srv)) + dt(srv);
        doflip = ~partial(srv);
        if any(doflip)
            jj = j(doflip);
            dd = mod(xn(doflip) - ppk + pi, TWOPI) - pi;
            cp = zeros(size(dd));
            cp(dd > -prm.sig & dd < 0) =  1/prm.sig;
            cp(dd >  0 & dd < prm.sig) = -1/prm.sig;
            rout = prm.beta0*(1 - s(jj).*gam.*cp/prm.sm);
            fl = rand(numel(jj),1) < rout/bmax;
            s(jj) = s(jj).*(1 - 2*fl);
        end
        live(j(partial(srv))) = false;
    end
end
if any(live)
    warning('lingermodel:advanceCap', ...
        '%d searchers unresolved after the step cap in advance', sum(live));
end
end
