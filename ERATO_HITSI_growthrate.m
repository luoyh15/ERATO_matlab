% function gamma2 = ERATO_HITSI_growthrate()
% data = load('eqldata3.mat');
% data = data.eqldata3;
% m=2^3;
% spart = 1+m:2*m:129;
% chipart = 1:m:65;
% quantities.ms = data.ps(1,1:2*m:129);
% quantities.mchi = data.pchi(1+m/2:m:65,1)';
% quantities.mchi = [-quantities.mchi(1),quantities.mchi,2*pi-quantities.mchi(end)];
% quantities.psi_p = data.psi_p;
% quantities.ps = data.ps(chipart,spart);
% quantities.pchi = data.pchi(chipart,spart);
% quantities.pr = data.pr(chipart,spart);
% quantities.pz = data.pz(chipart,spart);
% quantities.pbetachi = data.pbetachi(chipart,spart);
% quantities.plnr2Dchi = data.plnr2Dchi(chipart,spart);
% quantities.plnr2Ds = data.plnr2Ds(chipart,spart);
% quantities.pjphi = data.pjphi(chipart,spart);
% quantities.pp = data.pp(chipart,spart);
% quantities.ppDpsi = data.ppDpsi(chipart,spart);
% quantities.pT = data.pT(chipart,spart);
% quantities.pq = data.pq(chipart,spart);
% quantities.ppsiGrad2 = data.ppsiGrad2(chipart,spart);
% quantities.pH = data.pH(chipart,spart);
% quantities.plnpsiGradNormDpsi = data.plnpsiGradNormDpsi(chipart,spart);
% quantities.plnrDpsi = data.plnrDpsi(chipart,spart);
% quantities.ppsiDr = data.ppsiDr(chipart,spart);
% quantities.ppsiDz = data.ppsiDz(chipart,spart);

% % preprocess of the psi and mesh data
% [TR,psidata] = PreprocessPsiData();
% % get all the information about fit function
% psiF = PsiFitFunctionOrder2(TR,psidata);
% % load the equilibrium state data from PSI-Tri output
% [psi_p,fT_psi,fp_psi,fq_psi] = GetGSprofiles();
        
n_series = 2;
% n_s_series = 14:10:54;
n_s_series = 14;
n_fig = 1;
gamma2_series = zeros(length(n_series),length(n_s_series));
for j = 1:length(n_s_series)
    for i = 1:length(n_series)
        % toroidal mode number
        n=n_series(i);
        % mesh number
        n_s = n_s_series(j);
        n_chi = n_s+1;
        % non-equal mesh for s coordinate
        [ps,pchi,ms,mchi] = NonEqualMesh(psi_p,fq_psi,n_s,n_chi,n);
        ppsi  = psi_p*ps.^2;
        % get path functions
        PF = GetPathFunctions(ppsi,TR,psidata,psiF);
        
        % get all the quantities
        quantities = GetQuantitiesFromPSITri(psi_p,fT_psi,fp_psi,fq_psi,ps,pchi,ms,mchi,PF);
        figure(n_fig);
        gamma2 = -GrowthRateOfModeN(quantities,n);
        gamma2_series(i,j) = gamma2;
        if gamma2~=0
            title(['n=',int2str(n),',  mesh=',int2str(n_s)])
            n_fig = n_fig+1;
        end
    end
end



