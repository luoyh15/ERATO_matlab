function psiF = PsiFitFunctionOrder2_dc(psidata)
n_psi = psidata.npsi;
n_theta = psidata.ntheta/2;
% initialization
psiF.sectorsF = cell(n_psi,n_theta);

%% treat the region near magnetic axis as analytic
% the analysis region near magnetic axis
n_ana = floor(sqrt(n_psi)/2);
psi_ana = psidata.psi(n_ana+1,1);% threshold of analytic region
index_ana = psidata.psi<psi_ana;
% index_ana(end,2:end) = false;
r_ana = psidata.R(index_ana);
z_ana = psidata.Z(index_ana);
psivalue_ana = psidata.psi(index_ana);
[r_ana,index,~] = unique(r_ana,'rows');
z_ana = z_ana(index);
psivalue_ana = psivalue_ana(index);
% using poly33 to fit the analytic area
fpsi = fit([r_ana(:).^2,z_ana(:).^2],psivalue_ana(:),'poly33');
% the max distance from magnetic axis using for calculating rho coordinate
rmax_ana = max(r_ana(:));
rmin_ana = min(r_ana(:));
zmax_ana = max(z_ana(:));
% find the coordinate of magnetic axis
[xmin, ~] = fminsearch(@(x)fpsi(x(1).^2,x(2).^2), [(rmax_ana+rmin_ana)/2,0]);
r_magaxis = xmin(1);
% the max rho used to calculate the intesection of flux surface path
rhomax_ana = sqrt(max(rmax_ana-r_magaxis,r_magaxis-rmin_ana).^2+zmax_ana^2);

psiF.psi_ana = psi_ana;
psiF.analyticF = fpsi;
psiF.r_magaxis = r_magaxis;
psiF.rhomax_ana = rhomax_ana;
%% for the region not analytic
% number of points needed to fit is 6
for i = n_ana+1:n_psi
    for j = 1:n_theta
        [p_in,psi_in] = get4vertex(i,j);
        [p_add,psi_add] = getNeighbors(i,j);
        
        condmin = 10^10;
        for m = 1:size(p_add,1)-1
            for n = m+1:size(p_add,1)
                
                p_total = [p_in;p_add([m,n],:)];
                
                condA = condVdMatrix(p_total);
                if condA<condmin
                    condmin = condA;
                    p_fit_coordinates = p_total;
                    psi_fit_value = [psi_in;psi_add([m,n],:)];
                end
                
            end
        end
        fpsi = fit([p_fit_coordinates(:,1),p_fit_coordinates(:,2)],psi_fit_value,'poly22');
        coeff.p00 = fpsi.p00;coeff.p10 = fpsi.p10; coeff.p01 = fpsi.p01;
        coeff.p20 = fpsi.p20;coeff.p11 = fpsi.p11; coeff.p02 = fpsi.p02;
        psiF.sectorsF{i,j} = coeff;
    end
end
    function [p_neighbors,psi_neighbors] = getNeighbors(i_psi,i_theta)
        psi_index_down = i_psi-1;
        theta_index_down = i_theta+(0:1);
        rtemp_down = psidata.R(psi_index_down,theta_index_down);
        ztemp_down = psidata.Z(psi_index_down,theta_index_down);
        psi_down = psidata.psi(psi_index_down,theta_index_down);
        
        psi_index_side = i_psi+(0:1);
        theta_index_side = i_theta+[-1,2];
        theta_index_side(theta_index_side==0) = psidata.ntheta;
            
        rtemp_side = psidata.R(psi_index_side,theta_index_side);
        ztemp_side = psidata.Z(psi_index_side,theta_index_side);
        psi_side = psidata.psi(psi_index_side,theta_index_side);
        
        p_neighbors = [rtemp_down(:),ztemp_down(:);rtemp_side(:),ztemp_side(:)];
        psi_neighbors = [psi_down(:);psi_side(:)];
    end

    function [p_in,psi_in] = get4vertex(i_psi,i_theta)
        psi_index = i_psi+(0:1);
        theta_index = i_theta+(0:1);
        rtemp = psidata.R(psi_index,theta_index);
        ztemp = psidata.Z(psi_index,theta_index);
        psitemp = psidata.psi(psi_index,theta_index);
        p_in = [rtemp(:),ztemp(:)];
        psi_in = psitemp(:);
        
    end
    function condA = condVdMatrix(p_coordinates)
        % condition value of the vd matrix which used to decide the
        % accurecy of fit
        r_M = p_coordinates(:,1);
        z_M = p_coordinates(:,2);
        A = [ones(6,1),r_M,z_M,r_M.^2,r_M.*z_M,z_M.^2];
        condA = cond(A);
    end
end