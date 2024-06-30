pksi1 = pT./(pq.*ppsiGradNorm).*pXR;
pksi2 = pr.*ppsiGradNorm./(2*ps*Psi_s).*(pVR+pYR-pbetachi.*pXR);
pksi3 = pr.*pT./(2*ps*Psi_s).*pYR;
% the compontant of ksi in r,z direction
ppsiDr = double(subs(fpsiDr_rz,{r,z},{pr,pz}));
ppsiDz = double(subs(fpsiDz_rz,{r,z},{pr,pz}));
pksir = (pksi1.*ppsiDr-pksi2.*ppsiDz)./ppsiGradNorm;
pksiz = (pksi1.*ppsiDz+pksi2.*ppsiDr)./ppsiGradNorm;
figure(2);
quiver3(pr,zeros(size(pr)),pz,pksir,pksi3,pksiz);