figure(5);
subplot(2,3,1);mesh(ps,pchi,pjphi.^2./ppsiGradNorm.^2);
subplot(2,3,2);mesh(ps,pchi,pjphi./pr.*plnpsiGradNormDpsi);
subplot(2,3,3);mesh(ps,pchi,ppDpsi.*plogrDpsi);
subplot(2,3,4);mesh(ps,pchi,pjphi.^2./ppsiGradNorm.^2-...
pjphi./pr.*plnpsiGradNormDpsi);
subplot(2,3,5);mesh(ps,pchi,pjphi.^2./ppsiGradNorm.^2-...
pjphi./pr.*plnpsiGradNormDpsi-ppDpsi.*plogrDpsi);

