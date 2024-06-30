function eqlstate=GetQuantitiesFromPSITri_old()
%% load the equilibrium state data from PSI-Tri output
data_gs = GetGSprofiles();
% total plasma toroidal magnetic flux 
psi_p = data_gs.psi_p;
% pressure and its derivative as a function of psi
fp_psi = fit(data_gs.psi,data_gs.P,'spline');
fpDpsi_psi = fit(data_gs.psi,data_gs.dPDpsi,'spline');
% toroidal current density and its derivative as a function of psi
fT_psi  = fit(data_gs.psi,data_gs.T,'spline');
fTDpsi_psi = fit(data_gs.psi,data_gs.dTDpsi,'spline');
% safety factor as a function of psi
data_q  = Getqprofiles();
fq_psi = fit(data_q.psi,data_q.q,'spline');


%% (s,chi) coordinate system
n_s = 24;
n_chi = n_s+1;
ps = (1/(2*n_s):1/n_s:1-1/(2*n_s))';
pchi = (0:pi/(n_chi-1):pi)';

% ppsi = psi_p*ps.^2;
ppsi = psi_p*ps.^2;
% pressure and its derivative as a function of psi
pp = fp_psi(ppsi);
ppDpsi = fpDpsi_psi(ppsi);
% T and its derivative as a function of psi
pT = fT_psi(ppsi);
pTDpsi = fTDpsi_psi(ppsi);
% safety factor and its derivative
pq2 = fq_psi(ppsi);
pqDpsi = differentiate(fq_psi,ppsi);

%----initionalize
pq = zeros(size(ps));
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

ppsiDr = zeros(length(pchi),length(ps));
ppsiDz = zeros(length(pchi),length(ps));

% get path functions
PF = GetPathFunctions(ppsi);
% PF = load('PF.mat');
% PF = PF.PF;
syms r z positive

figure(4);
hold on;
ax1 = gca;
for i = 2:n_s
    % plasma pressure and it's derivative
    p_path = pp(i);
    pDpsi_path = ppDpsi(i);
    % T and it's derivative
    T_path = pT(i);
    TDpsi_path = pTDpsi(i);
    % the toroidal current density
    jphi_path = r.*pDpsi_path+TDpsi_path.*T_path./r;
    
    ichi = 2;
    for j = 1:size(PF{i}.functions)
 
        % the analytic form of psi and its relative quantities
        fpsi_rz =  sum(PF{i}.functions(j,:).*[1,r,z,r^2,r*z,z^2].^2);
        fpsiDr_rz =  diff(fpsi_rz,r);
        fpsiDz_rz =  diff(fpsi_rz,z);
        fpsiGradNorm_rz = sqrt(fpsiDr_rz.^2+fpsiDz_rz.^2);
        fpsiGrad2Dr_rz = diff(fpsiGradNorm_rz.^2,r);
        fpsiGrad2Dz_rz = diff(fpsiGradNorm_rz.^2,z);
        fpsiGrad2Dpsi_rz = (fpsiGrad2Dr_rz.*fpsiDr_rz+fpsiGrad2Dz_rz.*fpsiDz_rz)...
            ./fpsiGradNorm_rz.^2;
        flnpsiGradNormDpsi_rz = fpsiGrad2Dpsi_rz./(2*fpsiGradNorm_rz.^2);
        flnrDpsi_rz = 1./r.*fpsiDr_rz./fpsiGradNorm_rz.^2;
        
        if abs(PF{i}.lines(j,1))<abs(PF{i}.lines(j,2))
            % z as a function of r
            fz = -(PF{i}.lines(j,1)*r+PF{i}.lines(j,3))/PF{i}.lines(j,2);
            % left boundary of the fit area
            r_left = PF{i}.points(j,1);
            % right boundary of the fit area
            r_right = PF{i}.points(j+1,1);
            
            fr = r;
            dr = 1;
            dz = diff(fz,r);
            dl = -sqrt(dz.^2+dr.^2);
            
            lright = r_right;
            lleft = r_left;
            var = r;
            % plot the constant psi surface
            rr = linspace(r_left,r_right,10);
            plot(ax1,rr,double(subs(fz,r,rr)));
        else
            % r as a function of z
            fr = -(PF{i}.lines(j,2)*z+PF{i}.lines(j,3))/PF{i}.lines(j,1);
            % left boundary of the fit area
            z_left = PF{i}.points(j,2);
            % right boundary of the fit area
            z_right = PF{i}.points(j+1,2);
            
            fz = z;
            dr = diff(fr,z);
            dz = 1;
            dl = sqrt(dz.^2+dr.^2)*sign(PF{i}.lines(j,2)/PF{i}.lines(j,1));
            
            lright = z_right;
            lleft = z_left;
            var = z;
            % plot the constant psi surface
            zz = linspace(z_left,z_right,10);
            
            plot(ax1,double(subs(fr,z,zz)),zz);
        end
        
        % the integral kernal of q
        fint_q = T_path./(r.*fpsiGradNorm_rz).*dl;
        %     pq(i) = double(1/pi*int(-fint_q,r_right,r_left));
        kernel_q = matlabFunction(subs(fint_q,{r,z},{fr,fz}));
        pq(i) = pq(i)+1/pi*integral(@(lx) kernel_q(lx),lleft,lright);
        %dqdpsi
%         fint_qDpsi = 1/pi*((-jphi_path.*r-fpsiGrad2Dpsi_rz)./fpsiGradNorm_rz.^2+...
%             TDpsi_path./T_path).*T_path./(r.*fpsiGradNorm_rz).*dl;
%         kernel_qDpsi = matlabFunction(subs(fint_qDpsi,{r,z},{fr,fz}));
%         pqDpsi(i) = pqDpsi(i)+integral(@(lx) kernel_qDpsi(lx),lright,lleft);
       
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
%             ppsiDr(ichi,i) = double(subs(fpsiDr_rz,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             ppsiDz(ichi,i) = double(subs(fpsiDz_rz,{r,z},{pr(ichi,i),pz(ichi,i)}));
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
%             ppsiDr(ichi,i) = double(subs(fpsiDr_rz,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             ppsiDz(ichi,i) = double(subs(fpsiDz_rz,{r,z},{pr(ichi,i),pz(ichi,i)}));
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
%             ppsiDr(ichi,i) = double(subs(fpsiDr_rz,{r,z},{pr(ichi,i),pz(ichi,i)}));
%             ppsiDz(ichi,i) = double(subs(fpsiDz_rz,{r,z},{pr(ichi,i),pz(ichi,i)}));
%         end
    end
end
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
eqlstate.ppsiDr = ppsiDr;
eqlstate.ppsiDz = ppsiDz;
