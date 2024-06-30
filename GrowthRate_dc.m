
% % preprocess of the psi and mesh data
% psidata = PreprocessPsiData_dc();
% % get all the information about fit function
% psiF = PsiFitFunctionOrder2_dc(psidata);
% % load the equilibrium state data from PSI-Tri output
% [psi_p,fT_psi,fp_psi,fq_psi] = GetGSprofiles();
%         
n_series = 2;
n_s_series = 24:10:74;
% n_s_series = 84;
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
        s_max = sqrt(max(psidata.psi(:))/psi_p);
        [ps,pchi,ms,mchi] = NonEqualMesh_bc(psi_p,s_max,fq_psi,n_s,n_chi,n);
        ppsi  = psi_p*ps.^2;
        % get path functions
        PF = GetPathFunctions_dc(ppsi,psidata,psiF);
        
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




