figure(2);
contour(cr,cz,psi_rz,linspace(0,Psi_s,8),'k');
axis('equal');
xlim([0.5,1.3]);
ylim([0,0.8]);

ylabel('$$z$$','Interpreter','latex');
xlabel('$$r$$','Interpreter','latex');
title('Solovev equilibrium,$$E=2,a=1/3,q_0=0.3$$','Interpreter','latex');