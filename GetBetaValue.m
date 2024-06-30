function [beta,betap] = GetBetaValue(quantities)
pq = quantities.pq;
pT = quantities.pT;
pr = quantities.pr;
pp = quantities.pp;
ps = quantities.ps;
pjphi = quantities.pjphi;
ms = quantities.ms;
mchi = quantities.mchi;
msstep = ms(2:end)-ms(1:end-1);
mchistep = mchi(2:end)-mchi(1:end-1);
parea = mchistep*msstep';
% jacobian
pJ = pq.*pr.^2./pT;
% volum integration of pressure
Intp = sum(sum(pp.*pJ.*ps.*parea));
% volum integration of B2
IntB2 = sum(sum(pJ./pr.^2.*ps.*parea));
% volum integration of p/r
IntR = sum(sum(pp.*pJ./pr.*ps.*parea));
IntE = sum(sum(pjphi.*pJ./pr.*ps.*parea));

beta = 2*Intp/IntB2;
betap = 8*pi*IntR/(IntE^2);
disp(betap);