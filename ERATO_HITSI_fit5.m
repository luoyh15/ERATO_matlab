function eqlstate=ERATO_HITSI_fit5(filename)
% load the equilibrium state data from DCON output
data = load(filename);

psi_p = max(max(data.psi));
data.psi = flip(psi_p-data.psi);
data.R = flip(data.R(:,1:data.ntheta/2+1)',2);
data.Z = flip(data.Z(:,1:data.ntheta/2+1)',2);
data.P = flip(data.P);
data.T = flip(data.T);
data.q = flip(data.q);
R = data.R(end,1);
rmin = min(min(data.R));
rmax = max(max(data.R));
zmax = max(max(data.Z));
% (s,chi) coordinate system
n_s = 2^7;
n_chi = n_s/2+1;
% ps = 1/(2*n_s):1/n_s:1-1/(2*n_s);
pchi = 0:pi/(n_chi-1):pi;

% psi
% ppsi = psi_p*ps.^2;
ppsi = data.psi;
ps = sqrt(ppsi/psi_p);
% pressure and its derivative as a function of psi
fp_psi = fit(data.psi',data.P','spline');
pp = fp_psi(ppsi)';
ppDpsi = differentiate(fp_psi,ppsi)';
% T and its derivative as a function of psi
fT_psi = fit(data.psi',data.T','spline');
pT = fT_psi(ppsi)';
pTDpsi = differentiate(fT_psi,ppsi)';
% safety factor and its derivative
fq_psi = fit(data.psi',data.q','spline');
pq = fq_psi(ppsi)';
pqDpsi = differentiate(fq_psi,ppsi)';
%----initionalize
%r z
pr = zeros(length(pchi),length(ps));
pz = zeros(length(pchi),length(ps));
% the non-orthogonality betachi
pbetachi = zeros(length(pchi),length(ps));
% lnr2Dchi
plnr2Dchi = zeros(length(pchi),length(ps));
% lnr2Ds
plnr2Ds = zeros(length(pchi),length(ps));
% jphi
pjphi = zeros(length(pchi),length(ps));
% psiGrad2
ppsiGrad2 = zeros(length(pchi),length(ps));
% lnpsiGradNormDpsi
plnpsiGradNormDpsi = zeros(length(pchi),length(ps));
% lnrDpsi
plnrDpsi = zeros(length(pchi),length(ps));



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % the analytic area near magnetic axis
n_anl = floor(sqrt(n_s));
Msk_anl = data.psi<ppsi(n_anl);
psi_anl = data.psi(Msk_anl);
r_anl = data.R(:,Msk_anl);
z_anl = data.Z(:,Msk_anl);
rmin_anl = min(r_anl(:));
rmax_anl = max(r_anl(:));
zmax_anl = max(z_anl(:));
zmin_anl = min(z_anl(:));
Mpsi_anl = repmat(psi_anl,data.ntheta/2+1,1);
fpsi_r2z2 = fit([(r_anl(:)-R).^2,(z_anl(:)).^2],Mpsi_anl(:),'poly33');
syms r z spsi rho theta positive
% the analytic form of psi and its relative quantities
fpsi_rz =  sum(coeffvalues(fpsi_r2z2).*[1,r,z,r^2,r*z,z^2,r^3,r^2*z,r*z^2,z^3].^2);
fpsiDr_rz =  diff(fpsi_rz,r);
fpsiDz_rz =  diff(fpsi_rz,z);
fpsiGradNorm_rz = sqrt(fpsiDr_rz.^2+fpsiDz_rz.^2);
fpsiGrad2Dr_rz = diff(fpsiGradNorm_rz.^2,r);
fpsiGrad2Dz_rz = diff(fpsiGradNorm_rz.^2,z);
fpsiGrad2Dpsi_rz = (fpsiGrad2Dr_rz.*fpsiDr_rz+fpsiGrad2Dz_rz.*fpsiDz_rz)...
    ./fpsiGradNorm_rz.^2;
flnpsiGradNormDpsi_rz = fpsiGrad2Dpsi_rz./(2*fpsiGradNorm_rz.^2);
flnrDpsi_rz = 1./r.*fpsiDr_rz./fpsiGradNorm_rz.^2;
% rho as a function of theta and psi
fpsi_rhotheta = subs(fpsi_rz,{r,z},{rho*cos(theta),rho*sin(theta)});
fpsiGradNorm_rhotheta = subs(fpsiGradNorm_rz,{r,z},{rho*cos(theta),rho*sin(theta)});
fpsiGrad2Dpsi_rhotheta = subs(fpsiGrad2Dpsi_rz,{r,z},{rho*cos(theta),rho*sin(theta)});
flnpsiGradNormDpsi_rhotheta = subs(flnpsiGradNormDpsi_rz,{r,z},{rho*cos(theta),rho*sin(theta)});
flnrDpsi_rhotheta = subs(flnrDpsi_rz,{r,z},{rho*cos(theta),rho*sin(theta)});

[frho_thetapsi,~,~] = solve(fpsi_rhotheta-spsi,rho,'ReturnConditions',true,'MaxDegree',3);



for i = 2:n_anl
   
    % z function on the psi=const surface
    frho = subs(frho_thetapsi,spsi,ppsi(i));
    frhomsk = real(double(subs(frho,theta,pi)))>0&...
        real(double(subs(frho,theta,pi)))<=rmax_anl-R&...
        abs(imag(double(subs(frho,theta,pi))))<1e-20;
    if frhomsk==0
        pause(1);
    end
    rho_path = frho(frhomsk);
    rhoDtheta_path = diff(rho_path,theta);
    % plasma pressure and it's derivative
    pDpsi_path = ppDpsi(i);
    % T and it's derivative
    T_path = pT(i);
    TDpsi_path = pTDpsi(i);
    % norm of psi gradient
    psiGradNorm_path = subs(fpsiGradNorm_rhotheta,rho,rho_path);
    % the toroidal current density
    jphi_path = r.*pDpsi_path+TDpsi_path.*T_path./r;
    % the derivative of square of psi gradient
    psiGrad2Dpsi_path = subs(fpsiGrad2Dpsi_rhotheta,rho,rho_path);
    % art length
    dl = sqrt(rho_path.^2+rhoDtheta_path.^2);
    % the integral kernal of q
    kernel_q = subs(T_path./(r.*psiGradNorm_path),r,R+rho_path*cos(theta));
    %     pq(i) = double(1/pi*int(-fint_q,r_right,r_left));
    fint_q = matlabFunction(kernel_q*dl);
    pq(i) = 1/pi*integral(@(t) fint_q(t),0,pi);
    % the integral kernal of chi
    kernel_chi = subs(T_path./(pq(i)*r.*psiGradNorm_path),r,R+rho_path*cos(theta));
    %     chi_r = int(-fint_chi,r_right,r);
    fint_chi = matlabFunction(kernel_chi*dl);
    chi_t = @(t) integral(@(tt) fint_chi(tt),0,t);
    % dqdpsi
    kernel_qDpsi = subs(((-jphi_path.*r-psiGrad2Dpsi_path)./psiGradNorm_path.^2+...
        TDpsi_path./T_path)*kernel_chi,r,R+rho_path*cos(theta));
    %     pqDpsi(i) = double(pq(i)/pi*int(-fint_qDpsi,r_right,r_left));
    %     pqDpsi(i) = real(pqDpsi(i));
    fint_qDpsi = matlabFunction(kernel_qDpsi*dl);
    pqDpsi(i) = pq(i)/pi*integral(@(rr) fint_qDpsi(rr),0,pi);
%     pqDpsi(i) = real(pqDpsi(i));
    % the crosponding r of each mesh
%     ptheta(1,i) = 0;
%     ptheta(end,i) = pi;
%     % betachi
%     kernel_betachi = subs(((-jphi_path.*r-psiGrad2Dpsi_path)./psiGradNorm_path.^2+...
%         TDpsi_path./T_path-pqDpsi(i)/pq(i))*kernel_chi,r,R+rho_path*cos(theta));
%     fint_betachi = matlabFunction(kernel_betachi);
%     for j = 2:length(pchi)-1
%         % r,z
%         pr(j,i) = fzero(@(r) chi_r(r)-pchi(j),[r_left,R+(r_right-R)*cos(pi/(4*n_chi))]);
%         pz(j,i) = double(subs(z_path,r,pr(j,i)));
%         % betachi
%         pbetachi(j,i) = 2*ps(i)*psi_p*integral(@(rr) fint_betachi(rr),r_right,pr(j,i));
%         % lnr2Dchi
%         plnr2Dchi(j,i) = double(subs(2./r.*1./fint_chi,r,pr(j,i)));
%     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the aera not near the magnetic axis
% for i = n_anl+1:n_s
% figure(4);
% hold on;
% ax1 = gca;
% chi_st = zeros(size(data.R));
% betachi_st = zeros(size(data.R));
% for i = n_anl+1:n_s
%     %     psi_range = [(ppsi(i)+ppsi(i-1))/2,(ppsi(i)+ppsi(i+1))/2];
%     %     Msk_i = data.psi>psi_range(1)&data.psi<psi_range(2);
%     %     Msk_i = [false,false,data.psi(4:end)<ppsi(i)&data.psi(1:end-3)>ppsi(i)];
%     psi_i = data.psi(i+(-1:1));
%     r_i = data.R(:,i+(-1:1));
%     z_i = data.Z(:,i+(-1:1));  
%     % plasma pressure and it's derivative
%     pDpsi_path = ppDpsi(i);
%     % T and it's derivative
%     T_path = pT(i);
%     TDpsi_path = pTDpsi(i);
%     % the toroidal current density
%     jphi_path = r.*pDpsi_path+TDpsi_path.*T_path./r;
%     
%     ichi = 2;
%     for j = 2:size(r_i,1)-1
%         Mpsi_ij = repmat(psi_i,3,1);
%         r_ij = r_i(j+(-1:1),:);
%         z_ij = z_i(j+(-1:1),:);
%         fpsi_r2z2 = fit([r_ij(:).^2,z_ij(:).^2],Mpsi_ij(:),'poly22');
%         % the analytic form of psi and its relative quantities
%         fpsi_rz =  sum(coeffvalues(fpsi_r2z2).*[1,r,z,r^2,r*z,z^2].^2);
%         fpsiDr_rz =  diff(fpsi_rz,r);
%         fpsiDz_rz =  diff(fpsi_rz,z);
%         fpsiGradNorm_rz = sqrt(fpsiDr_rz.^2+fpsiDz_rz.^2);
%         fpsiGrad2Dr_rz = diff(fpsiGradNorm_rz.^2,r);
%         fpsiGrad2Dz_rz = diff(fpsiGradNorm_rz.^2,z);
%         fpsiGrad2Dpsi_rz = (fpsiGrad2Dr_rz.*fpsiDr_rz+fpsiGrad2Dz_rz.*fpsiDz_rz)...
%             ./fpsiGradNorm_rz.^2;
%         flnpsiGradNormDpsi_rz = fpsiGrad2Dpsi_rz./(2*fpsiGradNorm_rz.^2);
%         flnrDpsi_rz = 1./r.*fpsiDr_rz./fpsiGradNorm_rz.^2;
%        
%         % the psi=constant curve
%         fbound_up = polyfit(r_ij(:,2),z_ij(:,2),1);
%         
%         if abs(fbound_up(1))<1
%             % z as a function of r
%             [fz,~,~] = solve(fpsi_rz-ppsi(i),z,'ReturnConditions',true);
%             fzmsk = double(subs(fz,r,r_ij(2,2)))>=(1-1e-3)*z_ij(2,2)...
%                 &double(subs(fz,r,r_ij(2,2)))<=(1+1e-3)*z_ij(2,2);
%             %             display(fzmsk);
%             if sum(fzmsk)~=1
%                 pause(1);
%             end
%             fz = fz(fzmsk);
%             % right boundary of the fit area
%             r_right = r_ij(1,2);
%             % left boundary of the fit area
%             if j==size(r_i,1)-1
%                 r_left = r_ij(3,2);
%             else
%                 r_left = r_ij(2,2);
%             end
%             
%             fr = r;
%             dr = 1;
%             dz = diff(fz,r);
%             dl = -sqrt(dz.^2+dr.^2);
%             
%             lright = r_right;
%             lleft = r_left;
%             var = r;
%             % plot the constant psi surface
%             rr = linspace(r_left,r_right,10);
%             plot(ax1,rr,double(subs(fz,r,rr)));
%         else
%             % r as a function of z
%             [fr,~,~] = solve(fpsi_rz-ppsi(i),r,'ReturnConditions',true);
%             frmsk = double(subs(fr,z,z_ij(2,2)))>=(1-1e-3)*r_ij(2,2)...
%                 &double(subs(fr,z,z_ij(2,2)))<=(1+1e-3)*r_ij(2,2);
%             %             display(frmsk);
%             if sum(frmsk)~=1
%                 pause(1);
%             end
%             fr = fr(frmsk);
%             % right boundary of the fit area
%             z_right = z_ij(1,2);
%             % left boundary of the fit area
%             if j==size(r_i,1)-1
%                 z_left = z_ij(3,2);
%             else
%                 z_left = z_ij(2,2);
%             end
%             
%             fz = z;
%             dr = diff(fr,z);
%             dz = 1;
%             dl = -sqrt(dz.^2+dr.^2)*sign(fbound_up(1));
%             
%             lright = z_right;
%             lleft = z_left;
%             var = z;
%             % plot the constant psi surface
%             zz = linspace(z_left,z_right,10);
%             
%             plot(ax1,double(subs(fr,z,zz)),zz);
%         end
%         
% %         % the integral kernal of q
% %         fint_q = T_path./(r.*fpsiGradNorm_rz).*dl;
% %         %     pq(i) = double(1/pi*int(-fint_q,r_right,r_left));
% %         kernel_q = matlabFunction(subs(fint_q,{r,z},{fr,fz}));
% %         pq(i) = pq(i)+1/pi*integral(@(lx) kernel_q(lx),lright,lleft);
% %         %dqdpsi
% %         fint_qDpsi = 1/pi*((-jphi_path.*r-fpsiGrad2Dpsi_rz)./fpsiGradNorm_rz.^2+...
% %             TDpsi_path./T_path).*T_path./(r.*fpsiGradNorm_rz).*dl;
% %         kernel_qDpsi = matlabFunction(subs(fint_qDpsi,{r,z},{fr,fz}));
% %         pqDpsi(i) = pqDpsi(i)+integral(@(lx) kernel_qDpsi(lx),lright,lleft);
%        
%         % the integral kernal of chi
%         fint_chi = T_path./(data.q(i)*r.*fpsiGradNorm_rz).*dl;
%         %     chi_r = int(-fint_chi,r_right,r);
%         kernel_chi = matlabFunction(subs(fint_chi,{r,z},{fr,fz}));
%         dchi = integral(@(lx) kernel_chi(lx),lright,lleft);
%         chi_st(j,i) = chi_st(j-1,i)+dchi;
%         % betachi
%         fint_betachi = ((-jphi_path.*r-fpsiGrad2Dpsi_rz)./fpsiGradNorm_rz.^2+...
%             TDpsi_path./T_path-pqDpsi(i)/pq(i)).*fint_chi;
%         kernel_betachi = matlabFunction(subs(fint_betachi,{r,z},{fr,fz}));
%         dbetachi = 2*ps(i)*psi_p*integral(@(lx) kernel_betachi(lx),lright,lleft);
%         betachi_st(j,i) = betachi_st(j-1,i)+dbetachi;
%         while ichi<n_chi&&pchi(ichi)<=chi_st(j,i)
%             chi_r = @(x) integral(@(lx) kernel_chi(lx),lright,x)+chi_st(j-1,i);
%             x = fzero(@(x) chi_r(x)-pchi(ichi),[lleft,lright]);
%             pr(ichi,i) = double(subs(fr,var,x));
%             pz(ichi,i) = double(subs(fz,var,x));
%             %betachi_chi
%             pbetachi(ichi,i) = betachi_st(j-1,i)+...
%                 2*ps(i)*psi_p*integral(@(lx) kernel_betachi(lx),lright,x);          
%             % lnr2Dchi
%             plnr2Dchi(ichi,i) = double(subs(2*pq(i)*fpsiGradNorm_rz/T_path...
%                 *dr/dl,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             % lnr2Ds
%             plnr2Ds(ichi,i) = -pbetachi(ichi,i)*plnr2Dchi(ichi,i)+(4*ps(i)*psi_p)*...
%                 double(subs(flnrDpsi_rz,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             % jphi
%             pjphi(ichi,i) = double(subs(jphi_path,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             % psiGrad2
%             ppsiGrad2(ichi,i) = double(subs(fpsiGradNorm_rz.^2,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             % lnpsiGradNormDpsi
%             plnpsiGradNormDpsi(ichi,i) = double(subs(flnpsiGradNormDpsi_rz,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             % lnrDpsi
%             plnrDpsi(ichi,i) = double(subs(flnrDpsi_rz,{r,z},{pr(ichi,i),pz(ichi,i)}));
% 
%             ichi = ichi+1;
%         end
%         if j==2
%             ichi = 1;
%             pr(ichi,i) = double(subs(fr,var,lright));
%             pz(ichi,i) = double(subs(fz,var,lright));
%             %betachi_chi
%             pbetachi(ichi,i) = 0;         
%             % lnr2Dchi
%             plnr2Dchi(ichi,i) = 0;
%             % lnr2Ds
%             plnr2Ds(ichi,i) = -pbetachi(ichi,i)*plnr2Dchi(ichi,i)+(4*ps(i)*psi_p)*...
%                 double(subs(flnrDpsi_rz,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             % jphi
%             pjphi(ichi,i) = double(subs(jphi_path,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             % psiGrad2
%             ppsiGrad2(ichi,i) = double(subs(fpsiGradNorm_rz.^2,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             % lnpsiGradNormDpsi
%             plnpsiGradNormDpsi(ichi,i) = double(subs(flnpsiGradNormDpsi_rz,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             % lnrDpsi
%             plnrDpsi(ichi,i) = double(subs(flnrDpsi_rz,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             
%             ichi = ichi+1;
%         elseif j==size(r_i,1)-1
%             ichi = size(r_i,1);
%             pr(ichi,i) = double(subs(fr,var,lleft));
%             pz(ichi,i) = double(subs(fz,var,lleft));
%             %betachi_chi
%             pbetachi(ichi,i) = 0;         
%             % lnr2Dchi
%             plnr2Dchi(ichi,i) = 0;
%             % lnr2Ds
%             plnr2Ds(ichi,i) = -pbetachi(ichi,i)*plnr2Dchi(ichi,i)+(4*ps(i)*psi_p)*...
%                 double(subs(flnrDpsi_rz,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             % jphi
%             pjphi(ichi,i) = double(subs(jphi_path,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             % psiGrad2
%             ppsiGrad2(ichi,i) = double(subs(fpsiGradNorm_rz.^2,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             % lnpsiGradNormDpsi
%             plnpsiGradNormDpsi(ichi,i) = double(subs(flnpsiGradNormDpsi_rz,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             % lnrDpsi
%             plnrDpsi(ichi,i) = double(subs(flnrDpsi_rz,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             
%         end
%     end
% end
[ps,pchi] = meshgrid(ps,pchi);
% pressure
pp = repmat(pp,n_chi,1);
% psi derivative of p
ppDpsi = repmat(ppDpsi,n_chi,1);
% flux function
pT = repmat(pT,n_chi,1);
% s derivative of T
pTDs = pTDpsi.*(2*ps*psi_p);
% safety factor
pq = repmat(pq,n_chi,1);
% s derivative of q
pqDs = pqDpsi.*(2*ps*psi_p);
% H defined at the ERATO paper
pH = 2*pjphi.*ppsi.*pr./(ps.*ppsiGrad2)+pTDs./pT-pqDs./pq;
% mass density
prho = 1;

eqlstate.ps = ps;
eqlstate.pchi = pchi;
eqlstate.pp = pp;
eqlstate.ppDpsi = ppDpsi;
eqlstate.pT = pT;
eqlstate.pTDs = pTDs;
eqlstate.pq = pq;
eqlstate.pqDs = pqDs;
eqlstate.pH = pH;
eqlstate.prho = prho;
eqlstate.pr = pr;
eqlstate.pz = pz;
eqlstate.pbetachi = pbetachi;
eqlstate.plnr2Dchi = plnr2Dchi;
eqlstate.plnr2Ds = plnr2Ds;
eqlstate.plnrDpsi = plnrDpsi;
eqlstate.pjphi = pjphi;
eqlstate.ppsiGrad2 =ppsiGrad2;
eqlstate.plnpsiGradNormDpsi = plnpsiGradNormDpsi;
eqlstate.psi_p = psi_p;
