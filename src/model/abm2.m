function [Tq, trace] = abm2(delta, G2, cfg, q, seed, opt)
%ABM2  Event-driven two-donor individual based model.
%
%   Simulates the full velocity-jump search on the ring together with the
%   random-order service and abandonment at both donors, with no reduction.
%   The reduced system in twodonor is tested against it, so nothing is
%   borrowed from twogain except the reliance.
%
%   Each fed bee at a donor carries its own scenting mark, drawn once when it
%   is served. The plume peak is the mean of the marks present, so the
%   k^{-1/2} sharpening that twogain imposes through p ~ N(0, p0^2/k) emerges
%   here. The searcher spatial profile is whatever the dynamics produce, so
%   the coherence limit of the reduction is visible rather than imposed.
%
%   The reliance gamma is read from G2 at the current fed populations, since
%   it is a colony-level strategy. This matches the single-donor abm.
%
%   Capture at a point trap uses crossing detection, not a rate, so it is
%   exact and independent of the step size. A searcher that passes x_j in a
%   step is captured with probability q_j = 1 - exp(-kbar_j/v). Turning is
%   handled per step and is O(dt), which sets the step size.
%
%   Checks
%     with delta = Inf no bee ever scents, so Tq must match the reduced system
%       in mode 'none' to within Monte Carlo error;
%     with mu1 = mu2 the two donors must feed equal numbers on average;
%     the population must be conserved at every step, which is asserted below.
%
%   Fm shape. Each Fm{j} must stay a column. Under two subscripts, end in
%   Fm{j}(end+1,1) means size(Fm{j},1). Deleting the last element of a column
%   with Fm{j}(go) = [] returns a 1x0 row, so that append would pad and invent
%   a mark at 0, midway between the traps. Append and delete are both written
%   shape-independently.
%
%   Conservation is checked per block and per compartment. census() returns
%   the five counts and audit() reports which compartment moved, by how much,
%   against what the block should do.
%
%   A second invariant is asserted. Every bee ever fed is at a donor or has
%   departed, so
%       nfed == numel(Fm{1}) + numel(Fm{2}) + nD
%   holds at every step. This separates a spurious mark in Fm from a searcher
%   that was not removed from alive.
%
%   Returns the time at which q*N bees have been fed, or NaN if the cap is hit.

if nargin < 4 || isempty(q),    q = 0.5;   end
if nargin < 5 || isempty(seed), seed = 1;  end
if nargin < 6, opt = struct; end
if ~isfield(opt,'dt'),    opt.dt    = 0.02; end
if ~isfield(opt,'tmax'),  opt.tmax  = 4e4;  end
if ~isfield(opt,'trace'), opt.trace = false; end
% opt.audit prints the compartment vector every step. It is expensive and
% verbose. Use it with a single seed when debugging.
if ~isfield(opt,'audit'), opt.audit = false; end

rng(seed, 'twister');
prm = G2.prm;
v = prm.v;  sig = prm.sig;  beta0 = prm.beta0;  sm = prm.sm;  p0 = prm.p0;
traps = [-pi/2, pi/2];
thet  = [1 1]*(1 - exp(-prm.kbar/v));
mu    = [cfg.mu1, cfg.mu2];
N     = cfg.N;  dt = opt.dt;

% ---- searcher state -----------------------------------------------------
x  = -pi + 2*pi*rand(N,1);                % uniform position
d  = 2*(rand(N,1) > 0.5) - 1;             % uniform heading
alive = true(N,1);                        % still searching

U = [0 0];                                % hungry, waiting at each donor
Fm = {zeros(0,1), zeros(0,1)};            % scenting marks of fed bees present
nD = 0;                                   % departed, fed
nfed = 0;                                 % total fed so far
t = 0;  Tq = NaN;  step = 0;
if opt.trace, trace = zeros(0,6); else, trace = []; end

M2 = twodonor();
while t < opt.tmax
    step = step + 1;
    F = [numel(Fm{1}), numel(Fm{2})];

    % ---- reliance from the table at the current fed populations ---------
    [~, ~, gi] = M2.rates(G2, F(1), F(2), 'full');
    gam = G2.gam(gi);

    % ---- plume peaks are the mean of the marks present -----------------
    peak = nan(1,2);
    for j = 1:2
        if F(j) > 0, peak(j) = mean(Fm{j}); end
    end

    % ================= move and capture ==================================
    c0 = census(alive, U, Fm, nD);
    nseated = [0 0];
    idx = find(alive);
    if ~isempty(idx)
        xo = x(idx);
        % turning rate depends on the local cue slope and the heading
        g = cueslope(xo, peak, traps, sig);
        w = gam*g/sm;
        br = beta0*(1 - d(idx).*w);            % rate of leaving current state
        br = max(br, 0);
        flip = rand(numel(idx),1) < (1 - exp(-br*dt));
        d(idx(flip)) = -d(idx(flip));
        xn = wrapToPiLocal(xo + v*d(idx)*dt);
        % capture if the step crossed a trap. seated tracks who has already
        % been taken this step, so a bee cannot be deposited at both donors.
        % With traps pi apart and a step of v*dt this guard should never fire.
        seated = false(numel(idx),1);
        for j = 1:2
            hit = crossed(xo, xn, d(idx), traps(j), v*dt);
            hit = hit & ~seated;
            hit = hit & (rand(numel(idx),1) < thet(j));
            if any(hit)
                k = idx(hit);
                alive(k) = false;  U(j) = U(j) + numel(k);
                seated = seated | hit;
                nseated(j) = numel(k);
            end
        end
        x(idx) = xn;
    end
    audit(census(alive,U,Fm,nD), c0, ...
          [-sum(nseated), nseated(1), nseated(2), 0, 0], ...
          'move/capture', N, t, step, nfed, opt.audit);

    % ================= service, random order over U + F ===================
    c0 = census(alive, U, Fm, nD);
    nserved = [0 0];
    for j = 1:2
        if rand < 1 - exp(-mu(j)*dt)
            tot = U(j) + numel(Fm{j});
            if tot > 0 && rand < U(j)/tot         % a hungry bee was selected
                U(j) = U(j) - 1;
                % Append shape-independently (see the Fm shape note in the
                % header).
                Fm{j} = [Fm{j}(:); traps(j) + p0*randn];
                nfed = nfed + 1;
                nserved(j) = 1;
            end
        end
    end
    % U(j) falls by one and F(j) rises by one for each donor that served.
    audit(census(alive,U,Fm,nD), c0, ...
          [0, -nserved(1), -nserved(2), sum(nserved), 0], ...
          'service', N, t, step, nfed, opt.audit);

    % ================= abandonment ========================================
    c0 = census(alive, U, Fm, nD);
    nleft = 0;
    for j = 1:2
        nf = numel(Fm{j});
        if nf > 0 && isfinite(delta)
            go = rand(nf,1) < (1 - exp(-delta*dt));
            if any(go)
                % Deleting the last element of a column returns a 1x0 row, so
                % restore the column shape before the next append.
                Fm{j}(go) = [];  Fm{j} = Fm{j}(:);
                nD = nD + sum(go);  nleft = nleft + sum(go);
            end
        elseif nf > 0                              % delta = Inf, leave at once
            nD = nD + nf;  Fm{j} = zeros(0,1);  nleft = nleft + nf;
        end
    end
    audit(census(alive,U,Fm,nD), c0, ...
          [0, 0, 0, -nleft, nleft], ...
          'abandonment', N, t, step, nfed, opt.audit);

    t = t + dt;
    if opt.trace
        trace(end+1,:) = [t, sum(alive), U(1), numel(Fm{1}), U(2), numel(Fm{2})]; %#ok<AGROW>
    end

    if nfed >= q*N, Tq = t;  break;  end
end
if isnan(Tq)
    warning('abm2:timeCap', 'hit tmax = %g with %d of %d fed', ...
            opt.tmax, nfed, round(q*N));
end
end

% ------------------------------------------------------------------------
function c = census(alive, U, Fm, nD)
%CENSUS  The five compartments, in one place so every check counts the same
%   way. Order is searching, U1, U2, fed-and-present (both donors), departed.
c = [sum(alive), U(1), U(2), numel(Fm{1}) + numel(Fm{2}), nD];
end

function audit(c1, c0, expected, name, N, t, step, nfed, verbose)
%AUDIT  Compare the compartment change against what the block should do.
%   Reports the offending compartment rather than a total.
lbl = {'searching','U1','U2','F(present)','departed'};
if verbose
    fprintf('  t=%7.2f %-12s  [%s]\n', t, name, ...
            strjoin(arrayfun(@(z) sprintf('%d',z), c1, 'uni', 0), ' '));
end
if sum(c1) ~= N
    bad = find((c1 - c0) ~= expected);
    if isempty(bad)
        detail = 'no compartment moved unexpectedly, so the state was already wrong on entry';
    else
        parts = arrayfun(@(k) sprintf('%s moved %+d, expected %+d', ...
                     lbl{k}, c1(k)-c0(k), expected(k)), bad, 'uni', 0);
        detail = strjoin(parts, '; ');
    end
    error('abm2:conservation', ...
        ['%s block: %d against N = %d at t = %.2f, step %d.  %s.  ' ...
         'Compartments in [%s] were [%s] and are now [%s].'], ...
        name, sum(c1), N, t, step, detail, strjoin(lbl, ' '), ...
        strjoin(arrayfun(@(z) sprintf('%d',z), c0, 'uni', 0), ' '), ...
        strjoin(arrayfun(@(z) sprintf('%d',z), c1, 'uni', 0), ' '));
end
% Every bee ever fed is present at a donor or departed. A failure here with the
% population intact means a mark was created or destroyed without a bee moving.
if nfed ~= c1(4) + c1(5)
    error('abm2:fedCount', ...
        ['%s block: nfed = %d against F(present) + departed = %d at t = %.2f, ' ...
         'step %d.  The population is intact, so this is a scenting mark ' ...
         'added or removed without a bee changing compartment.'], ...
        name, nfed, c1(4) + c1(5), t, step);
end
end

function g = cueslope(x, peak, traps, sig)
%CUESLOPE  Slope of the summed tent cues at each position.  Slopes add where
%   the two plumes overlap. Each tent is +1/sig on its rising flank and -1/sig
%   on its falling flank, zero outside its own support.
g = zeros(size(x));
for j = 1:2
    if isnan(peak(j)), continue; end
    r = wrapToPiLocal(x - peak(j));
    g = g + (1/sig)*( (r > -sig & r <= 0) - (r > 0 & r < sig) );
end
end

function h = crossed(xo, xn, d, xt, step)
%CROSSED  Did a step from xo to xn pass the point xt, allowing for the wrap.
%   Compares signed distance travelled against distance to the target in the
%   direction of travel, which is wrap-safe for steps shorter than the ring.
gap = wrapToPiLocal((xt - xo).*d);          % forward distance to the trap
h = gap >= 0 & gap < step;
end

function y = wrapToPiLocal(x)
y = mod(x + pi, 2*pi) - pi;
end
