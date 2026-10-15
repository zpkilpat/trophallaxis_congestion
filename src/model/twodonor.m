function M2 = twodonor()
%TWODONOR  Collective feeding at two donors sharing one abandonment rate.
%
%   State (S,U1,F1,U2,F2,D), conserved. Random-order service at each donor,
%   productive rate mu_j U_j/(U_j+F_j). Fed bees abandon at the same rate
%   delta at both donors. Nothing is tuned per donor, and there is no capacity.
%
%   Arrival rates come from the exact table built by twogain, so cross
%   suppression between the two cues is carried rather than approximated. The
%   reliance is chosen inside rates(), after the average over configurations,
%   by maximizing the total averaged rate. A searcher cannot condition its
%   reliance on which donor it will reach or on the realized offsets.
%
%   Cue modes
%     'none'  no cue.
%     'sym'   both donors receive the cue strength that the mean fed population
%             would generate. Accelerates discovery equally, cannot allocate.
%     'full'  donor-specific cue used in the manuscript.
%   The 'sym' mode and decomp entry point are retained for legacy analyses;
%   manuscript figures use sweep and do not compute the comparison.
%
%   Checks
%     at mu1 = mu2 the 'sym' and 'full' optima must coincide;
%     in 'sym' mode the two returned rates must be equal to round-off, which
%     fails if the twogain table is not donor-swap symmetric.

M2.Tq     = @Tq;
M2.best   = @bestDelta;
M2.decomp = @decompose; % legacy comparison, unused by manuscript figures
M2.sweep  = @sweep;
M2.rates  = @rates;
end

% ------------------------------------------------------------------------
function [L1, L2, gs] = rates(G2, F1, F2, mode)
%RATES  Arrival rate at each donor given the mean lingering populations.
%   Third output is the index of the selected reliance, for diagnostics.
if strcmp(mode, 'none')
    L1 = G2.lam0;  L2 = G2.lam0;  gs = 1;  return
end
if strcmp(mode, 'sym'), F1 = 0.5*(F1+F2);  F2 = F1; end
w1 = poisw(F1, G2.kmax);  w2 = poisw(F2, G2.kmax);
ng = numel(G2.gam);
t1 = zeros(1,ng);  t2 = zeros(1,ng);
for g = 1:ng
    t1(g) = w1*G2.A1(:,:,g)*w2.';
    t2(g) = w1*G2.A2(:,:,g)*w2.';
end
[~, gs] = max(t1 + t2);          % one reliance, chosen on the total rate
L1 = t1(gs);  L2 = t2(gs);
end

function w = poisw(F, kmax)
% Poisson weights on k = 0..kmax, renormalized for the truncation.
k = 0:kmax;
if F <= 0
    w = [1, zeros(1,kmax)];
else
    lw = -F + k*log(F) - gammaln(k+1);
    w  = exp(lw - max(lw));
end
w = w/sum(w);
end

function [T, sol] = Tq(delta, G2, cfg, mode, q, tmax)
if nargin < 5, q = 0.5;  end
if nargin < 6, tmax = 6e5; end
N = cfg.N;
opt = odeset('RelTol',1e-8,'AbsTol',1e-11,'MaxStep',tmax/1e4, ...
             'Events',@(t,y) ev(t,y,q,N));
sol = ode15s(@(t,y) rhs(t,y,delta,G2,cfg,mode), [0 tmax], ...
             [N-2e-6; 1e-6; 0; 1e-6; 0; 0], opt);
if isempty(sol.xe), T = Inf; else, T = sol.xe(1); end
end

function dy = rhs(~, y, delta, G2, cfg, mode)
S = max(y(1),0); U1 = max(y(2),0); F1 = max(y(3),0);
U2 = max(y(4),0); F2 = max(y(5),0);
[L1, L2] = rates(G2, F1, F2, mode);
a1 = L1*S;  a2 = L2*S;
s1 = cfg.mu1*U1/(U1+F1+1e-9);
s2 = cfg.mu2*U2/(U2+F2+1e-9);
dy = [-(a1+a2); a1-s1; s1-delta*F1; a2-s2; s2-delta*F2; delta*(F1+F2)];
end

function [val, isterm, dir] = ev(~, y, q, N)
val = (y(3)+y(5)+y(6)) - q*N;  isterm = 1;  dir = 1;
end

% ------------------------------------------------------------------------
function [ds, Ts, T] = bestDelta(G2, cfg, mode, dgrid, q)
%BESTDELTA  Optimum located off the grid.  The grid argmin is quantized at
%   about 20 percent per step, and downstream descent factors are ratios of
%   two optima, so quantization would enter them twice.
if nargin < 4, dgrid = logspace(-2.5, 1.0, 44); end
if nargin < 5, q = 0.5; end
M = lingermodel();
f = @(d) Tq(d, G2, cfg, mode, q);
T = arrayfun(f, dgrid);
T(~isfinite(T)) = Inf;
[ds, Ts] = M.refineDelta(f, dgrid, T);
end

function out = decompose(G2, cfg0, rlist, dgrid, q)
%DECOMPOSE  Split the recruitment benefit into acceleration and allocation at
%   fixed total capacity mu1 + mu2 = cfg0.mutot, sweeping the ratio r.
%
%   F1 and F2 are the fed populations at the feeding threshold, the state
%   where the event stops, not a run average or the peak. The contrast
%   condition holds at these parameters.
if nargin < 4, dgrid = logspace(-2.5, 1.0, 44); end
if nargin < 5, q = 0.5; end
n = numel(rlist);
out = struct('r',rlist,'Tnone',zeros(1,n),'dsym',zeros(1,n),'Tsym',zeros(1,n), ...
             'dfull',zeros(1,n),'Tfull',zeros(1,n),'disc',zeros(1,n), ...
             'alloc',zeros(1,n),'total',zeros(1,n),'penalty',zeros(1,n), ...
             'contrast',zeros(1,n),'F1',zeros(1,n),'F2',zeros(1,n), ...
             'q',q,'Fwhen','feeding threshold');
fprintf(['  r     Tnone   dsym    Tsym    dfull   Tfull   disc%%  alloc%%' ...
         '  tot%%   pen%%   F1     F2\n']);
for i = 1:n
    r = rlist(i);
    cfg = cfg0;  cfg.mu2 = cfg0.mutot/(1+r);  cfg.mu1 = r*cfg.mu2;
    Tn  = Tq(1e4, G2, cfg, 'none', q);
    [ds, Ts] = bestDelta(G2, cfg, 'sym',  dgrid, q);
    [df, Tf] = bestDelta(G2, cfg, 'full', dgrid, q);
    [~, sol] = Tq(df, G2, cfg, 'full', q);
    yf = sol.y(:,end);
    [L1, L2] = rates(G2, yf(3), yf(5), 'full');
    % cost of running the real cue at the acceleration-optimal abandonment rate
    Tmis = Tq(ds, G2, cfg, 'full', q);
    out.Tnone(i)=Tn; out.dsym(i)=ds; out.Tsym(i)=Ts;
    out.dfull(i)=df; out.Tfull(i)=Tf;
    out.disc(i)  = 100*(Tn-Ts)/Tn;   out.alloc(i) = 100*(Ts-Tf)/Tn;
    out.total(i) = 100*(Tn-Tf)/Tn;   out.penalty(i) = 100*(Tmis/Tf - 1);
    out.contrast(i) = L1/L2;  out.F1(i) = yf(3);  out.F2(i) = yf(5);
    fprintf('%5.1f %8.1f %7.4f %7.1f %7.4f %7.1f %6.1f %7.1f %6.1f %6.1f %6.2f %6.2f\n', ...
        r, Tn, ds, Ts, df, Tf, out.disc(i), out.alloc(i), out.total(i), ...
        out.penalty(i), out.F1(i), out.F2(i));
end
if any(out.alloc < -1e-6)
    warning('twodonor:negativeAllocation', ...
        ['allocation share is negative at r = %s; the full cue is losing to ' ...
         'the comparison cue, which means the twogain table is not ' ...
         'donor-swap symmetric'], mat2str(rlist(out.alloc < -1e-6), 3));
end
end

% ------------------------------------------------------------------------
function out = sweep(G2, cfg0, rlist, dgrid, q)
%SWEEP  Total recruitment benefit and optimal rate at fixed total capacity.
%   Only the no-cue baseline and the donor-specific cue are evaluated.
if nargin < 4, dgrid = logspace(-2.5,1.0,44); end
if nargin < 5, q = 0.5; end
n = numel(rlist);
out = struct('r',rlist,'Tnone',zeros(1,n),'dfull',zeros(1,n), ...
    'Tfull',zeros(1,n),'total',zeros(1,n),'contrast',zeros(1,n), ...
    'F1',zeros(1,n),'F2',zeros(1,n),'q',q,'Fwhen','feeding threshold');
for i = 1:n
    cfg = cfg0; cfg.mu2 = cfg0.mutot/(1+rlist(i)); cfg.mu1 = rlist(i)*cfg.mu2;
    out.Tnone(i) = Tq(1e4,G2,cfg,'none',q);
    [out.dfull(i),out.Tfull(i)] = bestDelta(G2,cfg,'full',dgrid,q);
    [~,sol] = Tq(out.dfull(i),G2,cfg,'full',q);
    yf = sol.y(:,end);
    out.F1(i) = yf(3); out.F2(i) = yf(5);
    [L1,L2] = rates(G2,yf(3),yf(5),'full');
    out.contrast(i) = L1/L2;
    out.total(i) = 100*(out.Tnone(i)-out.Tfull(i))/out.Tnone(i);
    fprintf('r=%g: Tnone=%.2f, Tfull=%.2f, delta*=%.5f, benefit=%.2f%%\n', ...
        rlist(i),out.Tnone(i),out.Tfull(i),out.dfull(i),out.total(i));
end
end
