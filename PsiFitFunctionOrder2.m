function psiF = PsiFitFunctionOrder2(TR,psidata)
% number of all the triangulars
n_triangulars = size(TR.ConnectivityList,1);
% initialization
psiF.triangleF = cell(n_triangulars,1);
%% treat the region near magnetic axis as analytic
% the analysis region near magnetic axis
psi_ana = 0.1;% threshold of analytic region
points_ana = TR.Points(psidata<psi_ana,:);
psidata_ana = psidata(psidata<psi_ana);
% using poly33 to fit the analytic area
fpsi = fit([points_ana(:,1).^2,points_ana(:,2).^2],psidata_ana(:),'poly33');
% find the coordinate of magnetic axis
[xmin, ~] = fminsearch(@(x)fpsi(x(1).^2,x(2).^2), [0.3,0]);
r_magaxis = xmin(1);
% z_magaxis = xmin(2);
% the max distance from magnetic axis using for calculating rho coordinate
rmax_ana = max(points_ana(:,1));
rmin_ana = min(points_ana(:,1));
zmax_ana = max(points_ana(:,2));
% the max rho used to calculate the intesection of flux surface path
rhomax_ana = sqrt(max(rmax_ana-r_magaxis,r_magaxis-rmin_ana).^2+zmax_ana^2);

psiF.psi_ana = psi_ana;
psiF.analyticF = fpsi;
psiF.r_magaxis = r_magaxis;
psiF.rhomax_ana = rhomax_ana;
%% find all the triangulars that z=0 line passes
% (outside of the magnetic axis),
% the 'polodial' coordinate chi will start from 0 inside these triangulars

tri_midplane = [];%list contains all the according trianglars
for i = 1:n_triangulars
    p_index =  TR.ConnectivityList(i,:);
    p_location = TR.Points(p_index,:);
    zmax = max(p_location(:,2));
    zmin = min(p_location(:,2));
    rmin = min(p_location(:,1));
    if zmax>0&&zmin<0&&rmin>=r_magaxis
        tri_midplane = [tri_midplane;i];
        %         patch(p_location(:,1),p_location(:,2),'b');
    end
end
psiF.tri_midplane = tri_midplane;
%% for the region not analytic
np = 6; % number of points needed to fit
for i = 1:n_triangulars
    tri_index = i;
    p_total = TR.ConnectivityList(tri_index,:);
    tritemp = tri_index;
    i_p = np;
    while 1
        tritemp = TR.neighbors(tritemp(:));
        tritemp = tritemp(tritemp<n_triangulars);
        p_add = TR.ConnectivityList(tritemp(:),:);
        p_total = unique([p_total(:);p_add(:)],'rows','stable');
        if length(p_total)>=np
            while i_p<=length(p_total)
                p_fit_index = p_total([1:np-1,i_p]);
                p_fit_coordinates = TR.Points(p_fit_index,:);
                condA = condVdMatrix(p_fit_coordinates);
                if condA<10^7
                    break;
                end
                i_p = i_p+1;
            end
            if i_p<=length(p_total)
                break;
            end
        end
    end
    psi_fit_value = psidata(p_fit_index);
    fpsi = fit([p_fit_coordinates(:,1),p_fit_coordinates(:,2)],psi_fit_value,'poly22');
    psiF.triangleF{tri_index} = fpsi;
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