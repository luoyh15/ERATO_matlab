range_s = floor(linspace(1,n_s,10));
range_chi = floor(linspace(1,n_s,15));
figure(1);
plot(r_sc(:,range_s),z_sc(:,range_s),'k');
hold on;
plot(r_sc(range_chi,:)',z_sc(range_chi,:)','k');
hold off;
xlim([0.5,1.4]);
text(1.22,0.45,'$$\psi=$$ constant','Interpreter','latex');
text(1.22,0.55,'$$\chi=$$ constant','Interpreter','latex');
% annotation('arrow','X',[0.32,0.5],'Y',[0.6,0.4])

ylabel('$$z$$','Interpreter','latex');
xlabel('$$r$$','Interpreter','latex');
title(['$$\psi-\chi$$','  coordinate system'],'Interpreter','latex');