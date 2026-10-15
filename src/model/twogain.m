function G2 = twogain(prm, pgrid, kmax, ggrid)
%TWOGAIN  Exact two-donor arrival rates, averaged over scenter configurations.
%
%   Builds Lambda_j(gamma, p_1, p_2) from the arc solver, averages over the
%   stochastic scenting populations at each gamma, and takes the maximization
%   over gamma after the average,
%
%     gamma*(F1,F2) = argmax_gamma  E[Lambda_1 + Lambda_2]
%     Lambdatilde_j = E[Lambda_j] at gamma*.
%
%   A searcher sets one reliance and cannot condition it on a configuration
%   redrawn by the scenting populations, so the maximization sits outside both
%   the Poisson mixture over scenter counts and the Gaussian smoothing over
%   offsets. The total rate is the objective because a searcher cannot
%   condition on which donor it will reach. Rates rather than times are
%   averaged because rates from independent starts add.
%
%   Offsets are signed. Traps sit at -pi/2 and +pi/2, so +p on donor 1
%   displaces its plume toward donor 2 and +p on donor 2 displaces its plume
%   away from donor 1. The grid must cover both directions so that
%   Lambda_2(p1,p2) = Lambda_1(-p2,-p1) holds. pgrid must stay symmetric about
%   zero.
%
%   No cue is its own state, index np+1. Beyond |p| >= sig the plume is
%   inert, so the grid covers [-sig, sig] and all Gaussian mass past the edge
%   goes to the absent state. Folding that mass onto the last grid cell would
%   inflate the gain.
%
%   Checks run by selfcheck on every build
%     no cue at either donor gives Lambda_1 = Lambda_2 = lam0 to about 2%;
%     L2(i,j,g) = L1(rev(j),rev(i),g) with rev the index reversal, at every g;
%     A1(a,b,g) = A2(b,a,g), the same statement after averaging.
%
%   Cost is numel(ggrid) times the single-gamma build. The result carries a
%   key. Recompute twogain_key and compare before reusing a cache.

if nargin < 2 || isempty(pgrid), pgrid = linspace(-prm.sig, prm.sig, 25); end
if nargin < 3 || isempty(kmax),  kmax  = 25; end
if nargin < 4 || isempty(ggrid), ggrid = linspace(0, prm.gmax, 13); end

np = numel(pgrid);  ng = numel(ggrid);  ABS = np + 1;
if mod(np,2) ~= 1 || max(abs(pgrid + fliplr(pgrid))) > 1e-12
    error('twogain:asymmetricGrid', ...
        'pgrid must be symmetric about zero with an odd number of points');
end

A      = atwo();
traps  = [-pi/2, pi/2];
thetas = [1, 1]*(1 - exp(-prm.kbar/prm.v));

% ---- Lambda_j(gamma, p1, p2) on the offset grid, gamma not optimized here ---
L1 = zeros(ABS, ABS, ng);  L2 = zeros(ABS, ABS, ng);
fprintf('twogain: building %dx%d offsets at %d reliances\n', ABS, ABS, ng);
for g = 1:ng
    for i = 1:ABS
        for j = 1:ABS
            plumes = [];
            if i <= np, plumes(end+1) = traps(1) + pgrid(i); end %#ok<AGROW>
            if j <= np, plumes(end+1) = traps(2) + pgrid(j); end %#ok<AGROW>
            [~, ~, ~, Lam] = A.split(traps, thetas, plumes, prm.sig, ...
                                     prm.beta0, ggrid(g));
            L1(i,j,g) = Lam(1);  L2(i,j,g) = Lam(2);
        end
    end
    fprintf('  gamma = %.4f done\n', ggrid(g));
end

% ---- Gaussian weights over the grid, with the tail sent to the absent state --
% Row a is k = a-1 scenters.  k = 0 is pure absence.
dp = pgrid(2) - pgrid(1);
ed = [pgrid - dp/2, pgrid(end) + dp/2];
Wf = zeros(kmax+1, ABS);
Wf(1, ABS) = 1;
for a = 1:kmax
    sd = prm.p0/sqrt(a);
    c  = 0.5*(1 + erf(ed/(sd*sqrt(2))));
    Wf(a+1, 1:np) = diff(c);
    Wf(a+1, ABS)  = 1 - sum(Wf(a+1, 1:np));
end

% ---- collapse the offset grid onto scenter counts, keeping gamma open -------
A1 = zeros(kmax+1, kmax+1, ng);  A2 = zeros(kmax+1, kmax+1, ng);
for g = 1:ng
    A1(:,:,g) = Wf*L1(:,:,g)*Wf.';
    A2(:,:,g) = Wf*L2(:,:,g)*Wf.';
end

G2 = struct('p',pgrid,'gam',ggrid,'L1',L1,'L2',L2,'A1',A1,'A2',A2, ...
            'kmax',kmax,'prm',prm,'absidx',ABS, ...
            'lam0', prm.v*thetas(1)/(pi*(2-thetas(1))));
G2.key = twogain_key(prm, pgrid, kmax, ggrid);
selfcheck(G2);
end

% ------------------------------------------------------------------------
function selfcheck(G2)
np  = numel(G2.p);  ABS = np + 1;
rev = [np:-1:1, ABS];
d1  = max(abs(G2.L2(:) - reshape(permute(G2.L1(rev,rev,:), [2 1 3]), [], 1)));
d2  = max(abs(G2.A1(:) - reshape(permute(G2.A2, [2 1 3]), [], 1)));
n0  = max(abs([G2.L1(ABS,ABS,1), G2.L2(ABS,ABS,1)] - G2.lam0))/G2.lam0;
fprintf('twogain: no-cue Lam = (%.5f, %.5f) vs lam0 = %.5f  (%.1f%%)\n', ...
        G2.L1(ABS,ABS,1), G2.L2(ABS,ABS,1), G2.lam0, 100*n0);
fprintf('twogain: swap residual  raw %.2e   averaged %.2e\n', d1, d2);
if n0 > 0.02
    warning('twogain:noCue', ...
        'no-cue rate is %.1f%% from lam0; the arc solver or theta changed', 100*n0);
end
if d1 > 1e-8*max(G2.L1(:)) || d2 > 1e-8*max(G2.A1(:))
    warning('twogain:notSymmetric', ...
        ['Lambda_2 is not the donor-swapped Lambda_1; two identical donors ' ...
         'will get unequal rates.  Check that pgrid is signed and that the ' ...
         'plume list is being attributed to the right donor.']);
end
end
