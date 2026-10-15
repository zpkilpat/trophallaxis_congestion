function [gap1, gap2, th1, th2] = spectralgap(prm, p, gamma)
%SPECTRALGAP  Relaxation rate of the biased search process onto a cue, app:gap.
%
%   [gap1, gap2] = spectralgap(prm)
%   [gap1, gap2, th1, th2] = spectralgap(prm, p, gamma)
%
%   The decay rates theta of the search generator are the values s = -theta
%   at which the homogeneous arc system is singular, which is the condition
%   that the monodromy matrix around the ring has eigenvalue one. One plume
%   uses the afpt system, one donor at 0 with its plume peak at p. Two plumes
%   use the atwo system, donors at -pi/2 and +pi/2, each plume displaced by p
%   in the same direction. The gap is theta_2 - theta_1.
%
%   Defaults are full reliance gamma = s_m*sigma and p = 0.3, where the paper
%   quotes 0.449 for one plume and 0.274 for two. th1 and th2 return the two
%   slowest rates, [one plume; two plumes].
%
%   Only the real axis is scanned, up to theta = 1. At beta0 = 1 the first
%   complex rates sit near Re theta = 1 (a grid search over the complex plane
%   found none below), so the two real roots found are the two slowest.
%
%   prm needs sig, beta0, sm, v, and kbar.

if nargin < 2 || isempty(p),     p = 0.3; end
if nargin < 3 || isempty(gamma), gamma = prm.sm*prm.sig; end
theta = 1 - exp(-prm.kbar/prm.v);
F = afpt();  A2 = atwo();

one = @(s) F.assemble(s, p, prm.sig, prm.beta0, gamma, prm.sm, theta, prm.v);
two = @(s) A2.assemble(s, [-pi/2 pi/2], [theta theta], [-pi/2+p pi/2+p], ...
                       prm.sig, prm.beta0, gamma, prm.sm, prm.v, 1);
r1 = slowRates(one);  r2 = slowRates(two);
gap1 = r1(2) - r1(1);  gap2 = r2(2) - r2(1);
th1 = [r1(1); r2(1)];  th2 = [r1(2); r2(2)];
end

% ------------------------------------------------------------------------
function th = slowRates(Afun)
% The two smallest theta in (0,1] at which smin(A(-theta)) vanishes. Local
% minima of smin are refined by golden section and kept only if they are roots,
% as in afpt's lam1.
f  = @(th) smin(Afun(-th));
tg = linspace(1e-6, 1, 1001);
g  = arrayfun(f, tg);
th = [];
for i = 2:numel(tg)-1
    if g(i) < g(i-1) && g(i) < g(i+1)
        a = tg(i-1);  c = tg(i+1);
        for it = 1:80
            m1 = a + 0.382*(c - a);  m2 = a + 0.618*(c - a);
            if f(m1) < f(m2), c = m2; else, a = m1; end
        end
        t = 0.5*(a + c);
        if f(t) < 1e-8, th(end+1) = t; end %#ok<AGROW>
        if numel(th) == 2, return; end
    end
end
error('spectralgap:roots', 'found %d real decay rates below 1, expected 2', numel(th));
end

function y = smin(A)
A = A/norm(A);
s = svd(A);
y = s(end);
end
