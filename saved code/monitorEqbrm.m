% figure(3);
% subplot(3,3,1);mesh(cr,cz,fpsi_rz(cr,cz));
% subplot(3,3,2);mesh(cr,cz,fpsiDr_rz(cr,cz));
% subplot(3,3,3);mesh(cr,cz,fpsiDz_rz(cr,cz));
% subplot(3,3,4);mesh(cr,cz,fpsiGradNorm_rz(cr,cz));
% subplot(3,3,5);mesh(cr,cz,fpsiGrad2Dr_rz(cr,cz));
% subplot(3,3,6);mesh(cr,cz,fpsiGrad2Dz_rz(cr,cz));
% subplot(3,3,7);mesh(cr,cz,fpsiGrad2Dpsi_rz(cr,cz));
% subplot(3,3,8);mesh(cr,cz,flogpsiGradNormDpsi_rz(cr,cz));
% subplot(3,3,9);mesh(cr,cz,flogrDpsi_rz(cr,cz));
figure(3);
subplot(3,3,1);mesh(ps,pchi,ppsi);
subplot(3,3,2);mesh(ps,pchi,ppsiGradNorm);
subplot(3,3,3);mesh(ps,pchi,ppDpsi);
subplot(3,3,4);mesh(ps,pchi,pjphi);
subplot(3,3,5);mesh(ps,pchi,pqDs);
subplot(3,3,6);mesh(ps,pchi,pbetachi);
subplot(3,3,7);mesh(ps,pchi,pH);
subplot(3,3,8);mesh(ps,pchi,plogrDpsi);
subplot(3,3,9);mesh(ps,pchi,e);

figure(4);
subplot(2,3,1);mesh(ps,pchi,pjphi.^2./ppsiGradNorm.^2.*ppsi);
subplot(2,3,2);mesh(ps,pchi,pjphi./pr.*plogpsiGradNormDpsi.*ppsi);
subplot(2,3,3);mesh(ps,pchi,-ppDpsi.*plogrDpsi.*ppsi);
subplot(2,3,4);mesh(ps,pchi,pjphi.^2./ppsiGradNorm.^2.*ppsi+...
pjphi./pr.*plogpsiGradNormDpsi.*ppsi);
subplot(2,3,5);mesh(ps,pchi,pjphi.^2./ppsiGradNorm.^2.*ppsi+...
pjphi./pr.*plogpsiGradNormDpsi.*ppsi-ppDpsi.*plogrDpsi.*ppsi);

