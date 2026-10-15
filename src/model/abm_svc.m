function T = abm_svc(M, delta, G, ps, cfg, q, seed, svc)
%ABM_SVC  The single-donor ABM with service drawn from an arbitrary law.
%
%   M is the lingermodel struct. The velocity-jump advance is taken from it
%   as M.advance and the reliance rule as M.gamAt, so the search layer is the
%   one abm runs.
%
%   A copy of abm with one change. Service is a scheduled interval drawn from
%   svc() rather than a propensity mu in the thinning vector. The velocity-jump
%   advance, crossing-based capture, random-order service, per-bee scent
%   offsets, and the reliance read from the number of scenters via gamAt are
%   the same as in abm.
%
%   svc is a zero-argument sampler returning one service interval. The default
%   is exponential, which by memorylessness matches abm in distribution but
%   not path by path, since the random stream differs. The exponential arm of
%   Fig. 10(c) is therefore a control.
%
%   The service clock idles. tsvc is infinite while the donor is empty and is
%   set when the aggregation next becomes occupied. This is immaterial for
%   exponential service.
%
%   A fed bee drawn for service consumes the event and produces nothing.
%   Under non-exponential service that wasted event also consumes a full
%   service interval, which is the mechanism by which shape can matter and
%   what Fig. 10(c) tests.
%
%   The search step is capped at 1/(mu + delta*N), the mean thinning interval
%   in abm, so the O(dt^2) bias from resolving crowd events at the end of an
%   interval is the same size in both.

if nargin < 7 || isempty(seed), seed = 0; end
if nargin < 8 || isempty(svc),  svc  = @() -log(rand)/cfg.mu; end
rng(seed);

prm = G.prm;  N = cfg.N;  TWOPI = 2*pi;  p0 = ps.p0;
x = -pi + TWOPI*rand(N,1);
s = 2*(rand(N,1) > 0.5) - 1;
inBulk = true(N,1);
U = 0;  off = zeros(1,0);  nFed = 0;  t = 0;

Rab  = max(delta*N, 1e-12);            % abandonment thinning ceiling
DTC  = 1/(cfg.mu + delta*N);           % search step, matching abm's spacing
tsvc = inf;

while nFed < q*N && t < 1e7
    k = numel(off);
    if k == 0
        ppk = Inf;  gam = 0;
    else
        ppk = mean(off);  gam = M.gamAt(ps, k);
    end

    if (U + k) > 0 && ~isfinite(tsvc), tsvc = t + svc(); end

    tau = -log(rand)/Rab;
    [dt, which] = min([tau, tsvc - t, DTC]);   % 1 abandon, 2 serve, 3 nothing

    ib = find(inBulk);
    if ~isempty(ib)
        [x(ib), s(ib), tabs] = M.advance(x(ib), s(ib), dt, ppk, gam, prm);
        got = ~isnan(tabs);
        if any(got), inBulk(ib(got)) = false;  U = U + sum(got); end
    end
    t = t + dt;

    switch which
        case 1                                  % thinned abandonment
            k = numel(off);
            if k > 0 && rand < k/N, off(randi(k)) = []; end
        case 2                                  % scheduled service
            n = U + numel(off);
            if n > 0 && rand < U/n
                U = U - 1;  nFed = nFed + 1;
                off(end+1) = p0*randn;          %#ok<AGROW>
            end
            if U + numel(off) > 0
                tsvc = t + svc();
            else
                tsvc = inf;
            end
    end
end

if nFed < q*N
    warning('abm_svc:cap', ...
        'hit the time cap at delta = %g, seed %d, with %d of %d fed', ...
        delta, seed, nFed, ceil(q*N));
    T = NaN;
else
    T = t;
end
end
