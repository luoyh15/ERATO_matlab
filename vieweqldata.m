d=load('eqldata.mat');
d=d.totaldata;
figure(5);
subplot(2,4,1);mesh(d.pr(:,2:end-1),d.pz(:,2:end-1),d.pbetachi(:,2:end-1));
subplot(2,4,2);mesh(d.pr(:,2:end-1),d.pz(:,2:end-1),d.plnr2Dchi(:,2:end-1));
subplot(2,4,3);mesh(d.pr(:,2:end-1),d.pz(:,2:end-1),d.plnr2Ds(:,2:end-1));
subplot(2,4,4);mesh(d.pr(:,2:end-1),d.pz(:,2:end-1),d.plnrDpsi(:,2:end-1));
subplot(2,4,5);mesh(d.pr(:,2:end-1),d.pz(:,2:end-1),d.pjphi(:,2:end-1));
subplot(2,4,6);mesh(d.pr(:,2:end-1),d.pz(:,2:end-1),d.ppsiGrad2(:,2:end-1));
subplot(2,4,7);mesh(d.pr(:,2:end-1),d.pz(:,2:end-1),d.plnpsiGradNormDpsi(:,2:end-1));
subplot(2,4,8);mesh(d.ps(:,2:end-1),d.pchi(:,2:end-1),d.pbetachi(:,2:end-1));

figure(6);
plot(d.pr(:,2:end-1),d.pz(:,2:end-1));
hold on;
plot(d.pr(:,2:end-1)',d.pz(:,2:end-1)');