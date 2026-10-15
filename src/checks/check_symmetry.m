function out = check_symmetry(M, ps, G2, Fgrid, hlog)
%CHECK_SYMMETRY  Stability margin of the symmetric two-donor state, printout only.
%   At identical donors and F1 = F2 = Fbar, the self log-elasticity
%   d log L1/d log F1 and the cross term d log L1/d log F2 from twodonor rates
%   on the cached table, by central differences in log F. The antisymmetric
%   mode grows only when self - cross exceeds one (Eq. symmetry_criterion).
%   The single-donor elasticity peak is printed for comparison. Quoted in
%   the main text and the symmetry appendix.

if nargin < 4 || isempty(Fgrid), Fgrid = logspace(-2, log10(12), 90); end
if nargin < 5 || isempty(hlog),  hlog = 0.01; end
M2 = twodonor();
nF = numel(Fgrid);  eh = exp(hlog);
self = zeros(1,nF);  cross = zeros(1,nF);
for i = 1:nF
    F = Fgrid(i);
    self(i)  = (log(M2.rates(G2, F*eh, F, 'full')) - log(M2.rates(G2, F/eh, F, 'full')))/(2*hlog);
    cross(i) = (log(M2.rates(G2, F, F*eh, 'full')) - log(M2.rates(G2, F, F/eh, 'full')))/(2*hlog);
end
diffc = self - cross;
e1 = Fgrid.*M.dlogPsi(ps, Fgrid);
[smax, is] = max(self);  [cmin, ic] = min(cross);  [dmax, id] = max(diffc);  [emax, ie] = max(e1);
ipos = find(cross < 0, 1, 'last');
fprintf('check: self elasticity peaks at %.3f at Fbar = %.2f\n', smax, Fgrid(is));
if ipos < nF
    fprintf('check: cross term minimum %.3f at Fbar = %.2f, positive for Fbar > %.2f (max %.3f)\n', ...
            cmin, Fgrid(ic), Fgrid(ipos + 1), max(cross));
else
    fprintf('check: cross term minimum %.3f at Fbar = %.2f, negative throughout\n', cmin, Fgrid(ic));
end
fprintf('check: self - cross peaks at %.3f at Fbar = %.2f, short of the threshold by a factor %.2f\n', ...
        dmax, Fgrid(id), 1/dmax);
fprintf('check: single-donor elasticity peaks at %.3f at Fbar = %.2f\n', emax, Fgrid(ie));
out = struct('F',Fgrid,'self',self,'cross',cross,'diff',diffc,'e1',e1, ...
             'smax',smax,'cmin',cmin,'dmax',dmax,'Fdmax',Fgrid(id),'emax',emax);
end
