function eqlstate=GetQuantitiesFromPSITri(psi_p,fT_psi,fp_psi,fq_psi,ps,pchi,ms,mchi,PF)

%% (s,chi) coordinate system
n_s = length(ps);
n_chi = length(pchi);
% ps = (1/(2*n_s):1/n_s:1-1/(2*n_s))';
% pchi = (0:pi/(n_chi-1):pi)';
% ms  = (0:1/n_s:1);
% mchi = (-pi/(2*(n_chi-1)):pi/(n_chi-1):pi+pi/(2*(n_chi-1)));
% [ps,pchi,ms,mchi] = NonEqualMesh(psi_p,fq_psi,n_s,n_chi,n);
% ppsi = psi_p*ps.^2;
ppsi = psi_p*ps.^2;
% pressure and its derivative as a function of psi
pp = fp_psi(ppsi);
ppDpsi = differentiate(fp_psi,ppsi);
% T and its derivative as a function of psi
pT = fT_psi(ppsi);
pTDpsi = differentiate(fT_psi,ppsi);
% safety factor and its derivative
pq = fq_psi(ppsi);
pqDpsi = differentiate(fq_psi,ppsi);

%----initionalize
% pq = zeros(size(ps));
% pqDpsi = zeros(size(ps));
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
% PF = GetPathFunctions_squre(ppsi);
% PF = load('PF-14.mat');
% PF = PF.PF;
% figure(4);
% hold on;
% ax1 = gca;

for i = 1:n_s
    % plasma pressure and it's derivative
    
    pDpsi_path = ppDpsi(i);
    % T and it's derivative
    T_path = pT(i);
    TDpsi_path = pTDpsi(i);
    % the toroidal current density
    jphi_path = @(r,z) r.*pDpsi_path+TDpsi_path.*T_path./r;
    % the number of segment functions
    n_f = length(PF{i}.fpsis);
    % calculate safety factor q
    %     for j = 1:n_f
    %         % the function of psi and its relative quantities
    %         fpsiGradNorm_rz = PF{i}.fpsiGradNorms{j};
    %         % information of integration path
    %         fr = PF{i}.frs{j};
    %         fz = PF{i}.fzs{j};
    %         dl = PF{i}.dls{j};
    %         Interval = PF{i}.Intervals(j,:);
    %         % plot the constant psi surface
    %         xx = linspace(Interval(1),Interval(2),10);
    %         rr = fr(xx);
    %         zz = fz(xx);
    %         plot3(ax1,rr,zz,jphi_path(rr,zz));
    %         % the integral kernal of q
    %         fint_q = @(r,z) T_path./(r.*fpsiGradNorm_rz(r,z));
    %         pq(i) = pq(i)+1/(pi)*integral(@(x) fint_q(fr(x),fz(x)).*dl(x),Interval(1),Interval(2));
    %         %dqdpsi
    %         %         fint_qDpsi = @(r,z) 1/pi*((-jphi_path(r,z).*r-fpsiGrad2Dpsi_rz(r,z))./fpsiGradNorm_rz(r,z).^2+...
    %         %             TDpsi_path./T_path).*T_path./(r.*fpsiGradNorm_rz(r,z));
    %         %         pqDpsi(i) = pqDpsi(i)+integral(@(x) fint_qDpsi(fr(x),fz(x)).*dl(x),Interval(1),Interval(2));
    %
    %     end
    % q derivative of psi
    %     fq_psi = fit(ppsi,pq,'spline');
    %     pqDpsi = differentiate(fq_psi,ppsi);
    
    % some temprate quantities
    chi_temp = zeros(n_f+1,1);
    betachi_temp = zeros(n_f+1,1);
    i_chi = 1;
    
    for j = 1:n_f
        fpsiGradNorm_rz = PF{i}.fpsiGradNorms{j};
        fr = PF{i}.frs{j};
        fz = PF{i}.fzs{j};
        dl = PF{i}.dls{j};
        Interval = PF{i}.Intervals(j,:);
        % the integral kernal of chi
        fint_chi = @(r,z) T_path./(pq(i)*r.*fpsiGradNorm_rz(r,z));
        dchi = integral(@(x) fint_chi(fr(x),fz(x)).*dl(x),Interval(1),Interval(2));
        chi_temp(j+1) = chi_temp(j)+real(dchi);% the integration will have a very very small image part
    end
    % normalizaion of the chi_temp make sure [0,pi]
    chi_norm_coff = pi*(1+1e-8)/chi_temp(end);%make sure chi_temp(end)>pi
    chi_temp = chi_norm_coff*chi_temp;
    for j = 1:n_f
        % the function of psi and its relative quantities
        fpsiDr_rz = PF{i}.fpsiDrs{j};
        fpsiDz_rz = PF{i}.fpsiDzs{j};
        fpsiGradNorm_rz = PF{i}.fpsiGradNorms{j};
        fpsiGrad2Dpsi_rz = PF{i}.fpsiGrad2Dpsis{j};
        flnpsiGradNormDpsi_rz = PF{i}.flnpsiGradNormDpsis{j};
        flnrDpsi_rz = PF{i}.flnrDpsis{j};
        
        fr = PF{i}.frs{j};
        fz = PF{i}.fzs{j};
        dl = PF{i}.dls{j};
        dr = PF{i}.drs{j};
        Interval = PF{i}.Intervals(j,:);
        
        % dchi: the chi's element
        fint_chi = @(r,z) chi_norm_coff*T_path./(pq(i)*r.*fpsiGradNorm_rz(r,z));
        % betachi
        fint_betachi = @(r,z) 2*ps(i)*psi_p*((-jphi_path(r,z).*r-fpsiGrad2Dpsi_rz(r,z))./fpsiGradNorm_rz(r,z).^2+...
            TDpsi_path./T_path-pqDpsi(i)/pq(i)).*fint_chi(r,z);
        dbetachi = integral(@(x) fint_betachi(fr(x),fz(x)).*dl(x),Interval(1),Interval(2));
        betachi_temp(j+1) = betachi_temp(j)+real(dbetachi);% the integration will have a very very small image part
        while i_chi<n_chi+1&&pchi(i_chi)<=chi_temp(j+1)
            if i_chi == 1
                xtemp = Interval(1);
                pr(i_chi,i) = real(fr(xtemp));
                pz(i_chi,i) = 0;
                % betachi_chi
                pbetachi(i_chi,i) = 0;
                % lnr2Dchi
                plnr2Dchi(i_chi,i) = 0;
                
            elseif i_chi == n_chi
                xtemp = Interval(2);
                pr(i_chi,i) = real(fr(xtemp));
                pz(i_chi,i) = 0;
                % betachi_chi
                pbetachi(i_chi,i) = 0;
                % lnr2Dchi
                plnr2Dchi(i_chi,i) = 0;
            else
                chi_r = @(x) real(integral(@(xx) fint_chi(fr(xx),fz(xx)).*dl(xx),Interval(1),x)+chi_temp(j));
                xtemp = fzero(@(x) chi_r(x)-pchi(i_chi),[Interval(1),Interval(2)]);
                
                pr(i_chi,i) = real(fr(xtemp));
                pz(i_chi,i) = real(fz(xtemp));
                %betachi_chi
                pbetachi(i_chi,i) = betachi_temp(j)+...
                    real(integral(@(x) fint_betachi(fr(x),fz(x)).*dl(x),Interval(1),xtemp));
                % lnr2Dchi
                plnr2Dchi(i_chi,i) = 2*pq(i)*fpsiGradNorm_rz(pr(i_chi,i),pz(i_chi,i))/T_path...
                    *real(dr(xtemp)/dl(xtemp));
            end
            % lnr2Ds
            plnr2Ds(i_chi,i) = -pbetachi(i_chi,i)*plnr2Dchi(i_chi,i)+(4*ps(i)*psi_p)*...
                flnrDpsi_rz(pr(i_chi,i),pz(i_chi,i));
            % jphi
            pjphi(i_chi,i) = jphi_path(pr(i_chi,i),pz(i_chi,i));
            % psiGrad2
            ppsiGrad2(i_chi,i) = fpsiGradNorm_rz(pr(i_chi,i),pz(i_chi,i)).^2;
            % lnpsiGradNormDpsi
            plnpsiGradNormDpsi(i_chi,i) = flnpsiGradNormDpsi_rz(pr(i_chi,i),pz(i_chi,i));
            % lnrDpsi
            plnrDpsi(i_chi,i) = flnrDpsi_rz(pr(i_chi,i),pz(i_chi,i));
            
            ppsiDr(i_chi,i) = fpsiDr_rz(pr(i_chi,i),pz(i_chi,i));
            ppsiDz(i_chi,i) = fpsiDz_rz(pr(i_chi,i),pz(i_chi,i));
            
            i_chi = i_chi+1;
        end
        
    end
end
[ps,pchi] = meshgrid(ps,pchi);
% pressure
pp = repmat(pp',n_chi,1);
% psi derivative of p
ppDpsi = repmat(ppDpsi',n_chi,1);
% flux function
pT = repmat(pT',n_chi,1);
% s derivative of T
pTDs = pTDpsi'.*(2*ps*psi_p);
% safety factor
pq = repmat(pq',n_chi,1);
% s derivative of q
pqDs = pqDpsi'.*(2*ps*psi_p);
% H defined at the ERATO paper
pH = 2*pjphi.*ppsi'.*pr./(ps.*ppsiGrad2)+pTDs./pT-pqDs./pq;
% mass density
prho = 1;


eqlstate.ms = ms;
eqlstate.mchi = mchi;
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
