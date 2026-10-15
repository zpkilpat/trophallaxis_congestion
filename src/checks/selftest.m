%SELFTEST  Regression checks for the search and single-donor layers.
%   Analytic checks should agree to printed precision; stochastic ones to
%   within a few Monte Carlo standard errors.  Runtime ~2 min.
%
%   Expected values are stated after each block. Investigate any drift before
%   quoting a number in the paper.

clear; A = afpt(); M = lingermodel();
L = 2*pi;

%% 1. Ballistic limit, closed form for Tbar, CV and lambda1 --------------
%    Tbar = L(2-th)/(2 th);  Var = L^2/12 + L^2(1-th)/th^2;  lam1 = -log(1-th)/L
fprintf('\n--- 1. ballistic closed form (beta0 = 0, no plume) ---\n');
fprintf(' theta      Tbar        exact      CV        exact     lam1       exact\n');
for kb = [0.03 0.3 0.8 2.0]
    th = 1 - exp(-kb);
    [T, cv] = A.moments(0, 0.8, 0, 0, 1, th, 1);
    l1 = A.lam1(0, 0.8, 0, 0, 1, th, 1);
    Te  = L*(2-th)/(2*th);
    cve = sqrt(L^2/12 + L^2*(1-th)/th^2)/Te;
    fprintf('%6.4f %10.4f %10.4f %9.6f %9.6f %9.6f %9.6f\n', ...
            th, T, Te, cv, cve, l1, -log(1-th)/L);
end
% Expect exact agreement in every pair.

%% 2. Arc solver vs exact event-driven Monte Carlo -----------------------
fprintf('\n--- 2. arc solver vs Monte Carlo (2e5 walkers) ---\n');
fprintf(' kbar    p   gamma     T(MC)      +-se       T(arc)   dev(SE)  CV(MC)   CV(arc)\n');
cases = [0.8 0.0 0.0; 0.8 0.3 0.7; 0.8 0.6 0.5; 0.0314 0.2 0.79; 2.0 0.5 0.4];
for i = 1:size(cases,1)
    kb = cases(i,1); p = cases(i,2); gm = cases(i,3);
    th = 1 - exp(-kb);
    smp = A.mc(2e5, p, 0.8, 1, gm, 1, th, 1, 1);
    [Ta, cva] = A.moments(p, 0.8, 1, gm, 1, th, 1);
    se = std(smp)/sqrt(numel(smp));
    fprintf('%6.4f %4.1f %6.2f %10.4f %9.4f %10.4f %8.1f %8.4f %9.4f\n', ...
        kb, p, gm, mean(smp), se, Ta, (mean(smp)-Ta)/se, ...
        std(smp)/mean(smp), cva);
end
% Reference arc values are 14.848, 9.422, 13.252, 88.344, 9.138.
% Expect |dev| under about 2 standard errors on all five, and no
% afpt:mcTruncated warning. If one fires, the sample is biased fast.

%% 3. Search-layer structure, p*, the gain ceiling, and exponentiality ---
fprintf('\n--- 3. p* and exponentiality (kbar = 0.0314) ---\n');
prm = A.mkprm('kbar',0.0314);
pg  = linspace(0, 1.2, 41);
% G and W come from one builder so the ridge table and the surface share grids
% and prm. psiTable needs the surface to pick a reliance before the offset is
% drawn.
[G, W] = linger_tables(prm, pg);
fprintf('  G(0) = %.3f   p* = %.4f = %.3f*sigma  (sigma = %.2f)\n', ...
        G.G(1), G.pstar, G.pstar/prm.sig, prm.sig);
fprintf('  reliance sits on the positivity bound gamma = sm*sig = %.3f up to p = %.3f\n', ...
        prm.gmax, G.pcap);
fprintf(['  gamma* at the last grid point below p*: %.4f  (small => the optimum ' ...
         'slides to zero, large => it jumps)\n'], G.gamBelow);
fprintf('   p     gamma*     Gain    lam1*Tbar     CV\n');
for p = [0 0.2 0.4 0.6 0.8 1.0]
    gm = interp1(G.p, G.gam, p);
    [T, cv] = A.moments(p, prm.sig, prm.beta0, gm, prm.sm, prm.theta, prm.v);
    l1 = A.lam1(p, prm.sig, prm.beta0, gm, prm.sm, prm.theta, prm.v);
    fprintf('%4.1f %9.4f %8.3f %11.4f %8.4f\n', p, gm, interp1(G.p,G.G,p), l1*T, cv);
end
% Expect G(0) = 3.434, lam1*Tbar in [0.97,1.03], CV in [0.98,1.03] throughout.
% Expect p* = 0.7023 = 0.878*sigma, bisected rather than read off the grid,
% so it does not move when pg is refined. A change in the third digit is a
% regression. gamma* just below p* is about 0.10, so the optimum slides to
% zero, the gain approaches unity quadratically, and the 1e-8 threshold in
% pstarBisect limits p* to about 1e-4. Quote at most three digits.

%% 4. Averaged gain vs the deterministic substitution --------------------
fprintf('\n--- 4. Psitilde(Fbar) vs G(p0/sqrt(Fbar)) ---\n');
for p0 = [0.75 1.13 1.52]
    ps = M.psiTable(G, p0, W);
    fprintf(['  p0 = %.2f   g1 = %.3f (single scenter)   G(0) = %.3f (ceiling)' ...
             '   Psitilde(1) = %.3f   Fsub = %.2f\n'], ...
             p0, ps.g1, ps.G0, ps.psi1, ps.Fsub);
    fprintf('    Fbar   Psitilde   naive    excess captured\n');
    for Fb = [1 2 3 5 10 20]
        ex = interp1(ps.F, ps.psi, Fb);
        nv = max(interp1(G.p, G.G, min(p0/sqrt(Fb), max(G.p))), 1);
        fprintf('  %6.0f %10.3f %9.3f %13.1f%%\n', Fb, ex, nv, 100*(nv-1)/(ex-1));
    end
end
% Expect g1 = 1.552, 1.298, 1.159 at p0 = 0.75, 1.13, 1.52, matching the three
% curves in Fig. 4(a).  The excess these carry against the ceiling,
% (g1-1)/(G(0)-1), is 0.227, 0.123 and 0.065, so one scenting bee supplies
% between a fifteenth and a quarter of G(0) = 3.434.  Expect Fsub = 1.14, 2.59,
% 4.68 and the naive column to capture 0%, 34.7%, 71.0% and 88.4% of the excess
% at Fbar = 1, 2, 5 and 20 for p0 = 0.75. g1 is the k = 1 conditional gain,
% not Psitilde(1), which averages over a Poisson count including k = 0. This
% section uses 41 offset points and the figures use 61, which moves g1 by
% about 1e-3. The naive column captures ~0% for Fbar below Fsub = (p0/p*)^2,
% where the substituted offset lies outside the support.

%% 5. Reduced ODE, the delta -> inf limit must be exactly L_q/lam0 -------
fprintf('\n--- 5. no-lingering limit, T_q(inf) = ln(1/(1-q))/lam0 ---\n');
cfg = M.mkcfg('N',120,'mu',1,'lam0',0.005,'p0',0.75);
ps  = M.psiTable(G, cfg.p0, W);
for q = [0.25 0.5 0.9]
    T = M.Tq(1e4, ps, cfg, q);
    fprintf('  q = %.2f   ODE = %8.2f   exact = %8.2f   rel dev %8.1e\n', ...
            q, T, log(1/(1-q))/cfg.lam0, abs(T/(log(1/(1-q))/cfg.lam0)-1));
end
% Expect 57.5, 138.6, 460.5 with relative deviation under 1e-4. T_q(inf) does
% not depend on the gain, so a drift here is a different fault from a drift
% in block 4.

%% 6. Lingering criterion, evaluated ------------------------------------
fprintf('\n--- 6. does lingering pay?  (g1-1) q N (mu - A_q) > A_q ---\n');
for q = [0.25 0.5 0.9]
    [lhs, rhs, pays] = M.criterion(ps, cfg, q);
    fprintf('  q = %.2f   lhs = %8.3f   A_q = %.4f   ratio = %6.1f   -> %s\n', ...
        q, lhs, rhs, lhs/rhs, string(pays));
    % necessary condition rho < 1/(1-q), checked separately
    rho = cfg.lam0*cfg.N/cfg.mu;
    if rho >= 1/(1-q)
        fprintf('             rho = %.3f exceeds 1/(1-q) = %.3f, so no delta pays\n', ...
            rho, 1/(1-q));
    end
end
% This is eq:linger_criterion, taken from M.criterion so the test and the
% model state the same inequality.
% Expect true at every q, with ratio = 20.2, 77.3 and 934 at q = 0.25, 0.5, 0.9.

%% 7. expmInt across the confluent crossover ----------------------------
%   The arc propagator splits along the eigenvalues mu +/- w.  The integral
%   carries (P-Q)/w, which loses a digit per decade of |w*L| below 1e-2, so
%   expmInt switches to the augmented form there. This checks both branches
%   and checks that they agree at the crossover. A gap there would put a
%   discontinuity into Tbar as a function of s.
fprintf('\n--- 7. expmInt across the confluent crossover ---\n');
fprintf('   |w*L|      err(E)     err(I)\n');
%   The references are independent of both branches, expm for E and adaptive
%   quadrature for I. The augmented matrix is not used as a reference, since
%   below the crossover it would compare the confluent branch with itself.
Lseg = 1.7;  mu = 0.4;  k = 0.3;  b = 1.0;
for wl = [1e0 1e-1 1e-2 1e-3 1e-6 0]
    w  = wl/Lseg;
    Mm = mu*eye(2) + [k b; (w^2-k^2)/b -k];      % w is the eigenvalue split
    [E, I] = A.expmInt(Mm, Lseg);
    Eref = expm(Mm*Lseg);
    Iref = integral(@(x) expm(Mm*x), 0, Lseg, 'ArrayValued', true, ...
                    'AbsTol', 1e-14, 'RelTol', 1e-12);
    fprintf(' %9.1e  %10.2e %10.2e\n', wl, ...
        norm(E - Eref)/norm(Eref), norm(I - Iref)/norm(Iref));
end
% Expect both columns at 1e-14 or below at every |w*L|, including 0. The
% crossover row is the last on the closed form and runs about 1e-14. The rows
% below it use the confluent branch.

%% 8. Spectral gaps of the biased search generator, app:gap ---------------
%   The gap is computed, not quoted. At full reliance and p = 0.3 it should
%   print 0.449 for one plume and 0.274 for two. At zero reliance the slowest
%   rate of one plume should be 1/Tbar(0) of eq:gamma0.
fprintf('\n--- 8. spectral gaps (full reliance, p = 0.3) ---\n');
prm8 = struct('sig',0.8,'beta0',1,'sm',1,'v',1,'kbar',0.0314);
[g1, g2, t1, t2] = spectralgap(prm8);
fprintf('  one plume   theta1 = %.5f  theta2 = %.5f  gap = %.3f\n', t1(1), t2(1), g1);
fprintf('  two plumes  theta1 = %.5f  theta2 = %.5f  gap = %.3f\n', t1(2), t2(2), g2);
[~, ~, t1z] = spectralgap(prm8, 0.3, 0);
q0 = 1 - exp(-prm8.kbar/prm8.v);
T0 = pi*(2-q0)/(prm8.v*q0) + 2*prm8.beta0*pi^2/(3*prm8.v^2);
fprintf('  zero reliance  1/theta1 = %.2f  against eq:gamma0 Tbar(0) = %.2f\n', 1/t1z(1), T0);
% Expect gaps 0.449 and 0.274, and 1/theta1 within a fraction of a percent of Tbar(0).
