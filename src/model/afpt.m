function F = afpt()
%AFPT  Search layer for the trophallaxis lingering model.
%
%   Exact first-passage statistics for a velocity-jump searcher on the ring
%   [-pi,pi) with a point trap at x=0 and a tent-shaped pheromone plume whose
%   peak sits at offset p. Rates are piecewise constant on arcs, so the
%   Laplace transform of the FPT density is exact from one small linear solve.
%
%   A '-'-mover flips to '+' at rate bp(x) and a '+'-mover flips to '-' at
%   rate bm(x), with the linearized chemotactic response
%       b_{+/-}(x) = beta0 * (1 +/- gamma * c'(x) / sm).
%   Plume slope c' = +1/sig on (p-sig,p), -1/sig on (p,p+sig), 0 elsewhere.
%   Trap capture probability per crossing theta = 1 - exp(-kbar/v).
%
%   Positivity of both rates bounds the reliance at gamma <= sm*sig, mkprm's
%   gmax. At that bound the '+'-mover on the rising flank never turns, so the
%   searcher is ballistic toward the peak. The bound is physical, and
%   gainTable optima can sit on it.
%
%   Usage
%       F = afpt();
%       [T,cv] = F.moments(p,sig,beta0,gamma,sm,theta,v);
%       l1     = F.lam1(   p,sig,beta0,gamma,sm,theta,v);
%       smp    = F.mc(1e5, p,sig,beta0,gamma,sm,theta,v,seed);
%
%   prm is a struct with fields sig, beta0, sm, theta, v (see mkprm).

F.fhat    = @fhatUniform;
F.moments = @momentsFPT;
F.lam1    = @lam1v;
F.mc      = @mcFPT;
F.assemble = @assembleSys;    % homogeneous arc system, for spectralgap
F.arcs    = @buildArcs;
F.mkprm   = @mkprm;
F.expmInt = @expmInt;          % exposed for the selftest crossover sweep
end

% ------------------------------------------------------------------------
function prm = mkprm(varargin)
prm = struct('sig',0.8,'beta0',1.0,'sm',1.0,'kbar',0.0314,'v',1.0);
for k = 1:2:numel(varargin), prm.(varargin{k}) = varargin{k+1}; end
prm.theta = 1 - exp(-prm.kbar/prm.v);
prm.gmax  = prm.sm*prm.sig - 1e-6;          % ceiling for nonneg linear rates
prm.lam0  = prm.v*prm.theta/(pi*(2-prm.theta));   % unaided ballistic rate
end

function y = wrapPi(x)
y = mod(x + pi, 2*pi) - pi;
end

function [edges, cp] = buildArcs(p, sig)
% Breakpoints at the periodic seam, the trap at 0, and the three tent kinks.
pts = [-pi, 0, wrapPi(p-sig), wrapPi(p), wrapPi(p+sig)];
pts = unique(round(pts, 10));
pts = pts(pts > -pi + 1e-9);
edges = [-pi, pts(:).', pi];
n  = numel(edges) - 1;
cp = zeros(1, n);
for k = 1:n
    d = wrapPi(0.5*(edges(k) + edges(k+1)) - p);
    if d > -sig && d < 0
        cp(k) =  1/sig;      % rising flank: gradient points toward +x
    elseif d > 0 && d < sig
        cp(k) = -1/sig;      % falling flank
    end
end
end

function [E, I] = expmInt(M, L)
%EXPMINT  expm(M*L) and int_0^L expm(M*x) dx for a 2x2 (possibly complex) M.
%   The eigenvalue split mu +/- w gives E to machine precision at any w, but
%   the integral carries a factor (P-Q)/w that loses roughly one digit per
%   decade of |w*L| below 1e-2 and is wrong at the percent level at w = 0.
%   Below that crossover the augmented form is used,
%       expm([M eye(2); zeros(2,4)]*L) = [E, I; 0, eye(2)],
%   which is exact there. On a flat arc w vanishes only at s = 0 and
%   s = -2*beta0, and on a flank only at isolated s, so the branch is rare.
%
%   The two branches must agree to about 1e-15 at |w*L| = 1e-2, and both must
%   match expm of the augmented matrix. selftest section 7 checks this. If
%   that check drifts, widen the threshold rather than tightening it.
tr = M(1,1) + M(2,2);
dt = M(1,1)*M(2,2) - M(1,2)*M(2,1);
mu = tr/2;
w  = sqrt(mu^2 - dt);
if abs(w*L) < 1e-2
    X = expm([M, eye(2); zeros(2,4)]*L);
    E = X(1:2,1:2);
    I = X(1:2,3:4);
    return
end
K = M - mu*eye(2);
E = exp(mu*L) * (cosh(w*L)*eye(2) + (sinh(w*L)/w)*K);
P = e1((mu + w)*L, L);
Q = e1((mu - w)*L, L);
I = 0.5*(P + Q)*eye(2) + (0.5*(P - Q)/w)*K;
end

function val = e1(z, L)
% L * (exp(z)-1)/z, stable as z -> 0.
if abs(z) > 1e-6
    val = L*(exp(z) - 1)/z;
else
    val = L*(1 + z/2 + z^2/6);
end
end

function [A, b, edges, Is] = assembleSys(s, p, sig, beta0, gamma, sm, theta, v)
[edges, cp] = buildArcs(p, sig);
n  = numel(cp);
A  = zeros(2*n, 2*n);  A = complex(A);
b  = complex(zeros(2*n, 1));
Is = cell(1, n);
r  = 1;
Es = cell(1, n);
for k = 1:n
    L  = edges(k+1) - edges(k);
    bp = beta0*(1 + gamma*cp(k)/sm);        % rate INTO '+'
    bm = beta0*(1 - gamma*cp(k)/sm);        % rate INTO '-'
    M  = [ (s + bm)/v, -bm/v ; bp/v, -(s + bp)/v ];
    [Es{k}, Is{k}] = expmInt(M, L);
end
for k = 1:n
    kn = mod(k, n) + 1;
    E  = Es{k};
    ck = (2*k-1):(2*k);                     % unknowns at LEFT end of arc k
    cn = (2*kn-1):(2*kn);
    isTrap = abs(wrapPi(edges(k+1))) < 1e-12;
    if isTrap
        % '+' mover crosses from the left:  f+(0-) = th + (1-th) f+(0+)
        A(r, ck)    = E(1, :);
        A(r, cn(1)) = A(r, cn(1)) - (1 - theta);
        b(r) = theta;  r = r + 1;
        % '-' mover crosses from the right: f-(0+) = th + (1-th) f-(0-)
        A(r, cn(2)) = 1;
        A(r, ck)    = A(r, ck) - (1 - theta)*E(2, :);
        b(r) = theta;  r = r + 1;
    else
        for j = 1:2                          % continuity (and the seam)
            A(r, ck)    = E(j, :);
            A(r, cn(j)) = A(r, cn(j)) - 1;
            b(r) = 0;  r = r + 1;
        end
    end
end
end

function g = fhatUniform(s, p, sig, beta0, gamma, sm, theta, v)
% E[exp(-s*T)] averaged over uniform start position and both headings.
[A, b, ~, Is] = assembleSys(s, p, sig, beta0, gamma, sm, theta, v);
y = A\b;
tot = 0;
for k = 1:numel(Is)
    seg = Is{k} * y((2*k-1):(2*k));
    tot = tot + seg(1) + seg(2);
end
g = tot/(4*pi);
end

function [T, cv] = momentsFPT(p, sig, beta0, gamma, sm, theta, v)
% Mean and CV from derivatives of fhat at s=0.  Two-pass so the step size is
% scaled to the FPT itself (the FPT ranges over orders of magnitude in kbar).
f = @(s) real(fhatUniform(s, p, sig, beta0, gamma, sm, theta, v));
[T, ~] = d12(f, 1e-3);
[T, d2] = d12(f, 1e-3/max(T, 1e-12));
cv = sqrt(max(d2 - T^2, 0))/T;
end

function [T, d2] = d12(f, h)
% Central differences for -f'(0) and f''(0). The 30*f(0) term of the second
% stencil is written as 30 because absorption is certain, so f(0) = 1.
fp = f(h); fm = f(-h); fp2 = f(2*h); fm2 = f(-2*h);
T  = -(8*(fp - fm) - (fp2 - fm2))/(12*h);              % O(h^4)
d2 = (-fp2 + 16*fp - 30 + 16*fm - fm2)/(12*h^2);
end

function g = sminAt(s, p, sig, beta0, gamma, sm, theta, v)
A = assembleSys(s, p, sig, beta0, gamma, sm, theta, v);
A = A/norm(A);
sv = svd(A);
g = sv(end);
end

function l1 = lam1v(p, sig, beta0, gamma, sm, theta, v)
% Rightmost pole of fhat on the negative real axis = quasi-stationary decay
% rate.  Found as the first zero of the smallest singular value, which works
% for roots of any multiplicity (the ballistic case has a double root).
%
% A local minimum of the smallest singular value need not be a zero, so the
% refined candidate is accepted only if it is one. Otherwise NaN is returned,
% so a spurious pole does not propagate into lam1*Tbar.
T = momentsFPT(p, sig, beta0, gamma, sm, theta, v);
ss = -linspace(1e-9, 20/T, 600);
g  = arrayfun(@(s) sminAt(s, p, sig, beta0, gamma, sm, theta, v), ss);
l1 = NaN;
for i = 2:numel(ss)-1
    if g(i) < g(i-1) && g(i) < g(i+1)
        a = ss(i-1); c = ss(i+1);
        for it = 1:80                        % golden-section refine
            m1 = a + 0.382*(c - a);  m2 = a + 0.618*(c - a);
            if sminAt(m1,p,sig,beta0,gamma,sm,theta,v) < ...
               sminAt(m2,p,sig,beta0,gamma,sm,theta,v)
                c = m2;
            else
                a = m1;
            end
        end
        sroot = 0.5*(a + c);
        gmin  = sminAt(sroot,p,sig,beta0,gamma,sm,theta,v);
        if gmin > 1e-8
            warning('afpt:lam1NotRoot', ...
                ['local minimum of smin at s = %.6g is not a root ' ...
                 '(smin = %.2e); returning NaN'], sroot, gmin);
            l1 = NaN;
        else
            l1 = -sroot;
        end
        return
    end
end
end

% ------------------------------------------------------------------------
function smp = mcFPT(nwalk, p, sig, beta0, gamma, sm, theta, v, seed)
%MCFPT  Exact FPT sample.  Turns by Ogata thinning; trap crossings computed
%   analytically within each ballistic segment. No timestep, so no
%   discretization error.
%
%   Walkers alive at the step cap are dropped from the sample, which biases
%   it fast. A warning reports this.
if nargin < 9, seed = 0; end
rng(seed);
TWOPI = 2*pi;
x = -pi + TWOPI*rand(nwalk,1);
s = 2*(rand(nwalk,1) > 0.5) - 1;
t = zeros(nwalk,1);
out = nan(nwalk,1);
alive = true(nwalk,1);
bmax = beta0*(1 + gamma/(sm*sig));
logq = log1p(-theta);
nstep = 0;
for step = 1:400000
    nstep = step;
    idx = find(alive);
    if isempty(idx), break; end
    xa = x(idx); sa = s(idx); ta = t(idx);
    dt = -log(rand(numel(idx),1))/bmax;
    d  = v*dt;
    u1 = zeros(size(xa));                    % distance to first crossing of 0
    u1(sa > 0) = mod(-xa(sa > 0), TWOPI);
    u1(sa < 0) = mod( xa(sa < 0), TWOPI);
    u1(u1 == 0) = TWOPI;
    nc = zeros(size(d));
    m  = d >= u1;
    nc(m) = floor((d(m) - u1(m))/TWOPI) + 1;
    K  = floor(log(rand(numel(idx),1))/logq);    % survived crossings, Geom0
    hit = K < nc;
    if any(hit)
        j = idx(hit);
        out(j)   = ta(hit) + (u1(hit) + TWOPI*K(hit))/v;
        alive(j) = false;
    end
    srv = ~hit;
    if any(srv)
        j  = idx(srv);
        xn = wrapPi(xa(srv) + sa(srv).*d(srv));
        dd = wrapPi(xn - p);
        cp = zeros(size(xn));
        cp(dd > -sig & dd < 0) =  1/sig;
        cp(dd >  0   & dd < sig) = -1/sig;
        rateOut = beta0*(1 - sa(srv).*gamma.*cp/sm);   % rate of leaving state
        flip = rand(numel(j),1) < rateOut/bmax;
        x(j) = xn;
        s(j) = sa(srv).*(1 - 2*flip);
        t(j) = ta(srv) + dt(srv);
    end
end
if any(alive)
    warning('afpt:mcTruncated', ...
        ['%d of %d walkers unabsorbed after %d thinning steps; the returned ' ...
         'sample is biased fast'], sum(alive), nwalk, nstep);
end
smp = out(~isnan(out));
end
