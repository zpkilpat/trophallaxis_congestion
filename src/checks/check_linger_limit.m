function out = check_linger_limit(M, cfg, chk, rsmall)
%CHECK_LINGER_LIMIT  Indefinite lingering without recruitment, printout only.
%   Integrates the gain-free reduced dynamics at delta = 0 and reports
%   T_q(0)/T_q(inf) at the representative parameters, at the check points CHK (rows
%   [rho q]), and the small-load slope (ratio - 1)/rho at rho = RSMALL against
%   q/L_q, the slope that retains the queue lag dU/dt. The no-lingering time is
%   max(L_q/lam0, qN/mu). Quoted in the appendix.

if nargin < 3 || isempty(chk),    chk = [0.2 0.25; 0.6 0.9; 2 0.5; 4 0.1]; end
if nargin < 4 || isempty(rsmall), rsmall = 0.01; end
ps1 = struct('F', [0 1e7], 'psi', [1 1]);

rop = cfg.lam0*cfg.N/cfg.mu;  qop = 0.5;
zop = local_ratio(M, ps1, cfg, rop, qop);
fprintf('check: T_q(0)/T_q(inf) without recruitment at rho = %.3f, q = %.2f is %.4f (+%.1f%%)\n', ...
        rop, qop, zop, 100*(zop - 1));
zc = zeros(size(chk,1),1);
for j = 1:size(chk,1)
    zc(j) = local_ratio(M, ps1, cfg, chk(j,1), chk(j,2));
    fprintf('check: rho = %.2f, q = %.2f, ratio %.4f\n', chk(j,1), chk(j,2), zc(j));
end
qs = [0.1 0.25 0.5 0.75 0.9];
zs = local_ratio(M, ps1, cfg, rsmall, qs);
slope = (zs - 1)/rsmall;  pred = qs./log(1./(1 - qs));
fprintf('check: small-load slope at rho = %g, q = %s\n       reduced %s\n       q/L_q   %s\n', ...
        rsmall, mat2str(qs), mat2str(slope, 3), mat2str(pred, 3));
if any(zc <= 1) || zop <= 1
    warning('check:belowOne', 'T_q(0)/T_q(inf) without recruitment fell to or below one');
end
out = struct('rho',rop,'q',qop,'ratio',zop,'chk',chk,'chkratio',zc, ...
             'qs',qs,'slope',slope,'pred',pred);
end

% ------------------------------------------------------------------------
function z = local_ratio(M, ps1, cfg, rho, qs)
%LOCAL_RATIO  Gain-free reduced dynamics at delta = 0, one solve read at every
%   fed fraction in qs, over the exact no-lingering time max(L_q/lam0, qN/mu).
c = cfg;  c.lam0 = rho*cfg.mu/cfg.N;
qmax = max(qs);
tmax = 50*max(log(1/(1 - qmax))/c.lam0, qmax*c.N/c.mu);
[T, sol] = M.Tq(0, ps1, c, min(qmax + 0.005, 0.999), tmax);
t   = linspace(0, T, 20000);
Y   = deval(sol, t);
fed = Y(3,:) + Y(4,:);
[fed, iu] = unique(fed);
T0  = interp1(fed, t(iu), qs*c.N, 'pchip');
Ti  = max(log(1./(1 - qs))/c.lam0, qs*c.N/c.mu);
z   = T0./Ti;
end
