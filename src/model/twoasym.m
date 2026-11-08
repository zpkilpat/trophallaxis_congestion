function A = twoasym()
%TWOASYM  Simple asymptotics for the two-donor feeding time.
%
%   Two closed-form results, derived rather than fitted, compared against the
%   full individual-based model in figure6.
%
%   1. No-lingering limit. At delta = Inf no bee scents, the gain is one, and
%      every service event is productive. The colony is a pure
%      search-and-serve system with arrival rate lam0 per searcher per donor,
%      and Tinf is the exact solution of the two-donor ODE with the gain
%      removed. In the search-limited corner mu_j >> A_j it reduces to
%      L_q/(2*lam0).
%
%   2. Large-delta expansion. Repeating the single-donor calculation of
%      app:largedelta at two donors, the benefit term is unchanged because it
%      counts bees fed over the whole run, and the cost term becomes a sum over
%      donors, each evaluated at the end of the run where its queue must clear,
%
%        T_q(delta) = Tinf + (1/delta)*[ sum_j A_jq/(mu_j - A_jq)
%                                        - (g1 - 1)*q*N ] + O(delta^-2).
%
%      The asymmetry enters only through the cost sum. A_2q/(mu_2 - A_2q)
%      diverges as the slow donor approaches its capacity, which is why
%      lingering stops paying at large r.
%
%   Both results hold at large delta, outside the coherence window. The
%   reduction and the expansion agree there by construction, so they are
%   compared against the ABM. The ABM tracks the expansion at large delta and
%   departs from it as delta falls into the window.
%
%   The expansion requires A_jq < mu_j at both donors. With total capacity
%   mutot and A_jq = lam0*(1-q)*N per donor, this holds while
%   r < mutot/(lam0*(1-q)*N) - 1, which is r < 5.7 at mutot = 2,
%   lam0 = 0.005, N = 120, q = 0.5. Beyond it slope() warns. The reduction and
%   the simulation remain valid at large r.
%
%   At mu1 = mu2 the cost sum must be twice the single-donor cost at the same
%   per-donor arrival rate lam0*(1-q)*N, not half of it. Each donor absorbs
%   at lam0 per searcher, so the total discovery rate is 2*lam0.

A.Tinf   = @Tinf;
A.slope  = @slope;
A.Tlarge = @Tlarge;
A.Aend   = @Aend;
end

% ------------------------------------------------------------------------
function [T, Aq] = Tinf(G2, cfg, q)
%TINF  Exact no-lingering feeding time, from the ODE with the gain removed.
M2 = twodonor();
T  = M2.Tq(1e4, G2, cfg, 'none', q);
Aq = Aend(G2, cfg, q);
end

function Aq = Aend(G2, cfg, q)
%AEND  Arrival rate at each donor at the moment the run ends, with no cue.
%   With the gain removed the split is the no-cue split, so each donor takes
%   lam0 per searcher, and the searching pool at the threshold is (1-q)*N less
%   whatever is still queued. The leading-order pool (1-q)*N matches the order
%   of the expansion.
Aq = G2.lam0*(1 - q)*cfg.N*[1 1];
end

function s = slope(G2, cfg, q, g1)
%SLOPE  Coefficient of 1/delta in the large-delta expansion.
%   Negative means lingering pays at large delta.
mu = [cfg.mu1, cfg.mu2];
Aq = Aend(G2, cfg, q);
if any(Aq >= mu)
    warning('twoasym:overCapacity', ...
        ['donor %s is at or beyond its own capacity at the end of the run, ' ...
         'so the cost term is not defined and the expansion does not apply'], ...
        mat2str(find(Aq >= mu)));
end
cost = sum(Aq./max(mu - Aq, eps));
s    = cost - (g1 - 1)*q*cfg.N;
end

function T = Tlarge(delta, G2, cfg, q, g1, T0)
%TLARGE  The expansion itself, as a function of delta.
if nargin < 6 || isempty(T0), T0 = Tinf(G2, cfg, q); end
T = T0 + slope(G2, cfg, q, g1)./delta;
end
