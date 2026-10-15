function A = atwo()
%ATWO  Two partially absorbing point traps on the ring, each with a tent plume.
%
%   Exact Laplace-domain arc solver for the defective transforms
%       ghat_j(x,s) = E_{x}[ exp(-s T) 1{absorbed at trap j} ],
%   so ghat_j(x,0) = pi_j(x) is the splitting probability and
%   -d/ds ghat_j |_{s=0} = pi_j T_j gives the conditional mean first passage
%   time to trap j.
%
%   The rate at which a searcher joins donor j is
%       Lambda_j = pi_j / Tbar,      Tbar = sum_k pi_k T_k,
%   not pi_j/T_j. A searcher finds a donor at rate 1/Tbar, and it is donor j
%   with probability pi_j. The two forms agree only when T_1 = T_2, so
%   symmetric checks cannot distinguish them.
%
%   Traps are placed at -pi/2 and +pi/2 so neither sits on the periodic seam.
%
%   Usage
%       A = atwo();
%       [pi_,Tc,Tbar,Lam] = A.split(traps,thetas,plumes,sig,beta0,gamma);
%       [gs,Lam]          = A.optgamma(traps,thetas,plumes,sig,beta0,gmax);
%       [who,tabs]        = A.mc(n,traps,thetas,plumes,sig,beta0,gamma,seed);

A.split    = @splitTimes;
A.optgamma = @optGamma;
A.mc       = @mcTwo;
A.build    = @buildArcs;
A.assemble = @assemble;      % homogeneous arc system, for spectralgap
end

% ------------------------------------------------------------------------
function y = wrapPi(x)
y = mod(x + pi, 2*pi) - pi;
end

function c = cprimeAt(x, plumes, sig)
c = zeros(size(x));
for m = 1:numel(plumes)
    d = wrapPi(x - plumes(m));
    c = c + (d > -sig & d < 0)/sig - (d > 0 & d < sig)/sig;
end
end

function [edges, cp] = buildArcs(traps, plumes, sig)
pts = [-pi, wrapPi(traps(:)).'];
for m = 1:numel(plumes)
    pts = [pts, wrapPi(plumes(m)-sig), wrapPi(plumes(m)), wrapPi(plumes(m)+sig)]; %#ok<AGROW>
end
pts = unique(round(pts, 10));
pts = pts(pts > -pi + 1e-9);
edges = [-pi, pts, pi];
n = numel(edges) - 1;
cp = zeros(1, n);
for k = 1:n
    cp(k) = cprimeAt(0.5*(edges(k)+edges(k+1)), plumes, sig);
end
end

function [E, I] = expmInt(M, L)
tr = M(1,1) + M(2,2);
dt = M(1,1)*M(2,2) - M(1,2)*M(2,1);
mu = tr/2;  w = sqrt(mu^2 - dt);
if abs(w) < 1e-14, w = 1e-14; end
K = M - mu*eye(2);
E = exp(mu*L)*(cosh(w*L)*eye(2) + (sinh(w*L)/w)*K);
P = e1((mu+w)*L, L);  Q = e1((mu-w)*L, L);
I = 0.5*(P+Q)*eye(2) + (0.5*(P-Q)/w)*K;
end

function val = e1(z, L)
if abs(z) > 1e-6, val = L*(exp(z)-1)/z; else, val = L*(1 + z/2 + z^2/6); end
end

function [Amat, b, Is] = assemble(s, traps, thetas, plumes, sig, beta0, gamma, sm, v, jsrc)
[edges, cp] = buildArcs(traps, plumes, sig);
n  = numel(cp);
tw = wrapPi(traps);
Es = cell(1,n);  Is = cell(1,n);
for k = 1:n
    L  = edges(k+1) - edges(k);
    bp = beta0*(1 + gamma*cp(k)/sm);
    bm = beta0*(1 - gamma*cp(k)/sm);
    M  = [ (s+bm)/v, -bm/v ; bp/v, -(s+bp)/v ];
    [Es{k}, Is{k}] = expmInt(M, L);
end
Amat = complex(zeros(2*n)); b = complex(zeros(2*n,1)); r = 1;
for k = 1:n
    kn = mod(k, n) + 1;  E = Es{k};
    ck = (2*k-1):(2*k);  cn = (2*kn-1):(2*kn);
    hit = find(abs(wrapPi(edges(k+1) - tw)) < 1e-9, 1);
    if ~isempty(hit) && k < n
        th = thetas(hit);  src = th*(hit == jsrc);
        Amat(r, ck)    = E(1,:);
        Amat(r, cn(1)) = Amat(r, cn(1)) - (1-th);
        b(r) = src;  r = r + 1;
        Amat(r, cn(2)) = 1;
        Amat(r, ck)    = Amat(r, ck) - (1-th)*E(2,:);
        b(r) = src;  r = r + 1;
    else
        for j = 1:2
            Amat(r, ck)    = E(j,:);
            Amat(r, cn(j)) = Amat(r, cn(j)) - 1;
            r = r + 1;
        end
    end
end
end

function g = ghat(s, traps, thetas, plumes, sig, beta0, gamma, sm, v, jsrc)
[Amat, b, Is] = assemble(s, traps, thetas, plumes, sig, beta0, gamma, sm, v, jsrc);
y = Amat\b;  tot = 0;
for k = 1:numel(Is)
    seg = Is{k}*y((2*k-1):(2*k));
    tot = tot + seg(1) + seg(2);
end
g = tot/(4*pi);
end

function [pis, Tc, Tbar, Lam] = splitTimes(traps, thetas, plumes, sig, beta0, gamma, sm, v)
if nargin < 7, sm = 1; end
if nargin < 8, v  = 1; end
h = 1e-4;  m = numel(traps);
pis = zeros(1,m);  Tc = zeros(1,m);
for j = 1:m
    f = @(s) real(ghat(s, traps, thetas, plumes, sig, beta0, gamma, sm, v, j));
    p  = f(0);
    d1 = (8*(f(h)-f(-h)) - (f(2*h)-f(-2*h)))/(12*h);
    pis(j) = p;  Tc(j) = -d1/p;
end
Tbar = sum(pis.*Tc);           % unconditional MFPT
Lam  = pis/Tbar;               % arrival rate at each donor -- see header
end

function [gs, Lam, Tbar] = optGamma(traps, thetas, plumes, sig, beta0, gmax, sm, v)
%OPTGAMMA  Single shared reliance minimizing the unconditional MFPT.  A
%   searcher cannot know which donor it will reach, so it minimizes the time
%   to reach either.
if nargin < 7, sm = 1; end
if nargin < 8, v  = 1; end
obj = @(g) tbarOf(traps, thetas, plumes, sig, beta0, g, sm, v);
[gs, ~] = fminbnd(obj, 0, gmax, optimset('TolX',1e-4));
if obj(0) <= obj(gs), gs = 0; end
[~, ~, Tbar, Lam] = splitTimes(traps, thetas, plumes, sig, beta0, gs, sm, v);
end

function T = tbarOf(traps, thetas, plumes, sig, beta0, g, sm, v)
[~, ~, T] = splitTimes(traps, thetas, plumes, sig, beta0, g, sm, v);
end

% ------------------------------------------------------------------------
function [who, tabs] = mcTwo(n, traps, thetas, plumes, sig, beta0, gamma, sm, v, seed)
%MCTWO  Exact event-driven Monte Carlo.  Turns by Ogata thinning; trap
%   crossings enumerated analytically along each ballistic leg and processed in
%   order, so there is no timestep.
if nargin < 8,  sm = 1;   end
if nargin < 9,  v  = 1;   end
if nargin < 10, seed = 0; end
rng(seed);  TWOPI = 2*pi;
x = -pi + TWOPI*rand(n,1);
s = 2*(rand(n,1) > 0.5) - 1;
t = zeros(n,1);  who = -ones(n,1);  tabs = nan(n,1);  alive = true(n,1);
if isempty(plumes), bmax = beta0; else, bmax = beta0*(1 + gamma/(sm*sig)); end
bmax = max(bmax, 1e-12);
for step = 1:200000
    idx = find(alive);
    if isempty(idx), break; end
    xa = x(idx); sa = s(idx); ta = t(idx);
    dt = -log(rand(numel(idx),1))/bmax;
    d  = v*dt;
    hit = false(numel(idx),1);  thit = zeros(numel(idx),1);  wtr = -ones(numel(idx),1);
    for i = 1:numel(idx)
        ev = [];
        for k = 1:numel(traps)
            u = mod((traps(k) - xa(i))*sa(i), TWOPI);
            if u == 0, u = TWOPI; end
            while u <= d(i)
                ev = [ev; u, k];  %#ok<AGROW>
                u = u + TWOPI;
            end
        end
        if ~isempty(ev)
            ev = sortrows(ev, 1);
            for e = 1:size(ev,1)
                if rand < thetas(ev(e,2))
                    hit(i) = true;  thit(i) = ta(i) + ev(e,1)/v;  wtr(i) = ev(e,2);
                    break
                end
            end
        end
    end
    if any(hit)
        j = idx(hit);  tabs(j) = thit(hit);  who(j) = wtr(hit);  alive(j) = false;
    end
    srv = ~hit;
    if any(srv)
        j  = idx(srv);
        xn = wrapPi(xa(srv) + sa(srv).*d(srv));
        if isempty(plumes)
            rate = beta0*ones(numel(j),1);
        else
            rate = beta0*(1 - sa(srv).*gamma.*cprimeAt(xn, plumes, sig)/sm);
        end
        fl = rand(numel(j),1) < rate/bmax;
        x(j) = xn;  s(j) = sa(srv).*(1 - 2*fl);  t(j) = ta(srv) + dt(srv);
    end
end
end
