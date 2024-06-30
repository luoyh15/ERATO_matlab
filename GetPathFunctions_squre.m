function pF = GetPathFunctions_squre(psi)
% this function gets the path functions of the flux surfaces of value psi
% input: psi
%        the psi value of every flux surfaces
%        make sure psi is a column vector
% output: pF
%         cell array same size with psi, each cell containing the
%         infomation of path function of the corresponding psi value
%         pF{i}.functions ([m,3] double):containing all the path functions'
%         coefficients which are the form psi = c1*r + c2*z + c3
%         pF{i}.points ([m+1,2] double):containing the start and end points
%         of corresponding path function
%         m is the number of segment functions of that path

%% input psi must be a column vector
if size(psi,2)~=1
    error('make sure ''psi'' is a column vector');
end

%% read data from Psi-tri output files
% all the points or vertexes coordinates
points = h5read('mesh.h5','/R');
points = points';
% the three vertexes numbering of all the triangulars
triangulars = h5read('mesh.h5','/LC');
triangulars = triangulars'+1;% note that matlab's numbering start from 1
% the psi value of all the points/vertexes
psidata = h5read('scalar_dump.h5','/Psi0018');


%% preprocess of the raw data
% construct the triangulation mesh using matlab's triangulation structure
TR = triangulation(double(triangulars),points);
% number of all the triangulars
n_triangulars = size(triangulars,1);
% renormalize the psi data in order to consist with ERATO framework
psi_p = max(psidata);
psidata = (psi_p-psidata)/psi_p;

% number of flux surface that need calculated
n_psi = size(psi,1);
% initialization
pF  = cell(size(psi));

% figure(1);
% hold on;
%% treat the region near magnetic axis as analytic
% the analysis region near magnetic axis
n_ana = floor(sqrt(n_psi));
points_ana = points(psidata<psi(n_ana),:);
psidata_ana = psidata(psidata<psi(n_ana));
% % using poly33 to fit the analytic area
% fpsi = fit([points_ana(:,1),points_ana(:,2)],psidata_ana(:),'poly33');
% using poly22 to fit the analytic area
fpsi = fit([points_ana(:,1).^2,points_ana(:,2).^2],psidata_ana(:),'poly33');


% find the coordinate of magnetic axis
[xmin, ~] = fminsearch(@(x)fpsi(x(1).^2,x(2).^2), [0.3,0]);
r_magaxis = xmin(1);
% z_magaxis = xmin(2);
% the max distance from magnetic axis using for calculating rho coordinate
rmax_ana = max(points_ana(:,1));
rmin_ana = min(points_ana(:,1));
zmax_ana = max(points_ana(:,2));
rhomax_ana = sqrt(max(rmax_ana-r_magaxis,r_magaxis-rmin_ana).^2+zmax_ana^2);


%% find all the triangulars that z=0 line passes
% (outside of the magnetic axis),
% the 'polodial' coordinate chi will start inside these triangulars

tri_midplane = [];%list contains all the according trianglars
for i = 1:n_triangulars
    p_index =  TR.ConnectivityList(i,:);
    p_location = TR.Points(p_index,:);
    zmax = max(p_location(:,2));
    zmin = min(p_location(:,2));
    rmin = min(p_location(:,1));
    if zmax>0&&zmin<=0&&rmin>=r_magaxis
        tri_midplane = [tri_midplane;i];
        %         patch(p_location(:,1),p_location(:,2),'b');
    end
end
%% get the path functions of one specific psi value
% r = sym('r','positive');
% z = sym('z','positive');
for j = 1:n_psi
    psii = psi(j);
    % find the triangular that containing psii value
    tri_path = [];% list containing all the triangular index pass through
    edges_path = [];% list containing all the edges that pass through
    
    Interval_path = [];
    var_path = [];% list containing all the integrate variables
    point_path = [];% list containing the start and end points of the corresponding line
    for i = 1:size(tri_midplane,1)
        p_index =  TR.ConnectivityList(tri_midplane(i),:);
        psi_value = psidata(p_index);
        psimax = max(psi_value);
        psimin = min(psi_value);
        if psimax>=psii&&psimin<psii% the triangular contains value 'psii'
            tri_neighbors = TR.neighbors(tri_midplane(i));
            p_nb_index = TR.ConnectivityList(tri_neighbors,:);
            p_nb_index = unique(p_nb_index(:),'rows');
            psi_nb_value = psidata(p_nb_index);
            p_nb_coordinates = TR.Points(p_nb_index,:);
            fpsi = fit([p_nb_coordinates(:,1).^2,p_nb_coordinates(:,2).^2],psi_nb_value,'poly22');
            fpath = @(r,z) fpsi(r.^2,z.^2)-psii;
            
            z_start = 0;
            r_start = fzero(@(r) fpath(r,0),[min(p_nb_coordinates(:,1)),max(p_nb_coordinates(:,1))]);
            tri_temp = TR.pointLocation(r_start,z_start);
            if tri_temp==tri_midplane(i)
                tri_path = tri_temp;
%                 fpsi_path = [fpsi.p00,fpsi.p10,fpsi.p01,fpsi.p20,fpsi.p11,fpsi.p02];
                point_path = [r_start,z_start];% store first point
                
                p1 = p_index(psi_value==psimax);
                p2 = p_index(psi_value==psimin);
                p3 = p_index(~((p_index==p1)|(p_index==p2)));
                
                r1 = TR.Points(p1,1);z1 = TR.Points(p1,2);
                r2 = TR.Points(p2,1);z2 = TR.Points(p2,2);
                r3 = TR.Points(p3,1);z3 = TR.Points(p3,2);
                
%                 patch([r1,r2,r3],[z1,z2,z3],'r');
                
                [r_end,z_end] = FindIntesectionPoint(fpath,[r1,z1],[r2,z2]);
                % find the edge 'psi=psii' surface line intesect
                if r_end>min(r1,r2)&&r_end<=max(r1,r2)&&z_end>0
                    edges_path = [p1,p2];
                elseif psidata(p3)>=psii
                    [r_end,z_end] = FindIntesectionPoint(fpath,[r3,z3],[r2,z2]);
                    edges_path = [p3,p2];
                else
                    [r_end,z_end] = FindIntesectionPoint(fpath,[r1,z1],[r3,z3]);
                    edges_path = [p1,p3];
                end
                % determine which coordinate as integrat variable
                [fr,fz,dl,Interval] = SolvePathFunction(fpsi,psii,point_path(end,:),[r_end,z_end]);
                % get function fpsi's derivatives using in the calculating relevent quantities
                [fpsiDr,fpsiDz,fpsiGradNorm,fpsiGrad2Dpsi,flnpsiGradNormDpsi,flnrDpsi] = GetPsiDerivative(fpsi);
                
                i_count = 1;% to count the number of triangulars in the path
                fpsi_path{i_count,1} = @(r,z) fpsi(r.^2,z.^2); 
                fpsiDr_path{i_count,1} = fpsiDr;
                fpsiDz_path{i_count,1} = fpsiDz;
                fpsiGradNorm_path{i_count,1} = fpsiGradNorm;
                fpsiGrad2Dpsi_path{i_count,1} = fpsiGrad2Dpsi;
                flnpsiGradNormDpsi_path{i_count,1} = flnpsiGradNormDpsi;
                flnrDpsi_path{i_count,1} = flnrDpsi;
                fr_path{i_count,1} = fr;
                fz_path{i_count,1} = fz;
                dl_path{i_count,1} = dl;
                Interval_path = [Interval_path;Interval];
                
%                 var_path = [var_path;var_temp];
                point_path = [point_path;r_end,z_end];% store next point
%                 plot(point_path([end-1,end],1),point_path([end-1,end],2),'k');
                break;
            end
        end
    end
    
    while 1
        
        % find next triangular that attached to the last edge
        tri_attached = TR.edgeAttachments(edges_path(end,:));
        tri_attached = tri_attached{1};
        if tri_attached(1)==tri_path(end)
            tri_next = tri_attached(2);
        else
            tri_next = tri_attached(1);
        end
        p_index =  TR.ConnectivityList(tri_next,:);
        
        tri_neighbors = TR.neighbors(tri_next);
        p_nb_index = TR.ConnectivityList(tri_neighbors,:);
        p_nb_index = unique(p_nb_index(:),'rows');
        psi_nb_value = psidata(p_nb_index);
        p_nb_coordinates = TR.Points(p_nb_index,:);
        fpsi = fit([p_nb_coordinates(:,1).^2,p_nb_coordinates(:,2).^2],psi_nb_value,'poly22');
        fpath = @(r,z) fpsi(r.^2,z.^2)-psii;
        %         patch(p_coordinates(:,1),p_coordinates(:,2),'r');
        p1 = edges_path(end,1);
        p2 = edges_path(end,2);
        % make sure psidata(p1)>psidata(p2)
        if psidata(p1)<psidata(p2)
            p = p2;
            p2 = p1;
            p1 = p;
        end
        p3 = p_index(~((p_index==p1)|(p_index==p2)));
        
        r1 = TR.Points(p1,1);z1 = TR.Points(p1,2);
        r2 = TR.Points(p2,1);z2 = TR.Points(p2,2);
        r3 = TR.Points(p3,1);z3 = TR.Points(p3,2);
%         patch([r1,r2,r3],[z1,z2,z3],'r');
        % determine if the line cross z=0
        if ~(z1>0&&z2>0&&z3>0)
            z_temp = 0;
            r_temp = fzero(@(r) fpath(r,0),[min(p_nb_coordinates(:,1)),max(p_nb_coordinates(:,1))]);
            tri_temp = TR.pointLocation(r_temp,z_temp);
            if tri_temp==tri_next
                % r as a function of z
                [fr,fz,dl,Interval] = SolvePathFunction(fpsi,psii,point_path(end,:),[r_temp,z_temp]);
                % get function fpsi's derivatives using in the calculating relevent quantities
                [fpsiDr,fpsiDz,fpsiGradNorm,fpsiGrad2Dpsi,flnpsiGradNormDpsi,flnrDpsi] = GetPsiDerivative(fpsi);
                
                i_count = i_count+1;
                fpsi_path{i_count,1} = @(r,z) fpsi(r.^2,z.^2);  
                fpsiDr_path{i_count,1} = fpsiDr;
                fpsiDz_path{i_count,1} = fpsiDz;
                fpsiGradNorm_path{i_count,1} = fpsiGradNorm;
                fpsiGrad2Dpsi_path{i_count,1} = fpsiGrad2Dpsi;
                flnpsiGradNormDpsi_path{i_count,1} = flnpsiGradNormDpsi;
                flnrDpsi_path{i_count,1} = flnrDpsi;
                
                fr_path{i_count,1} = fr;
                fz_path{i_count,1} = fz;
                dl_path{i_count,1} = dl;
                Interval_path = [Interval_path;Interval];
                
                point_path = [point_path;r_temp,z_temp];
%                 plot(point_path([end-1,end],1),point_path([end-1,end],2),'k');
                break;
            end
        end
%         patch([r1,r2,r3],[z1,z2,z3],'r');
        % find next edge that intesect 'psi=psii' surface line
        if psidata(p3)>=psii
            [r_end,z_end] = FindIntesectionPoint(fpath,[r3,z3],[r2,z2]);
            edges_next = [p3,p2];
        else
            [r_end,z_end] = FindIntesectionPoint(fpath,[r1,z1],[r3,z3]);
            edges_next = [p1,p3];
        end
        % determine which coordinate as integrat variable
        [fr,fz,dl,Interval] = SolvePathFunction(fpsi,psii,point_path(end,:),[r_end,z_end]);
        % get function fpsi's derivatives using in the calculating relevent quantities
        [fpsiDr,fpsiDz,fpsiGradNorm,fpsiGrad2Dpsi,flnpsiGradNormDpsi,flnrDpsi] = GetPsiDerivative(fpsi);
        
        i_count = i_count+1;
        fpsi_path{i_count,1} = @(r,z) fpsi(r.^2,z.^2);
        fpsiDr_path{i_count,1} = fpsiDr;
        fpsiDz_path{i_count,1} = fpsiDz;
        fpsiGradNorm_path{i_count,1} = fpsiGradNorm;
        fpsiGrad2Dpsi_path{i_count,1} = fpsiGrad2Dpsi;
        flnpsiGradNormDpsi_path{i_count,1} = flnpsiGradNormDpsi;
        flnrDpsi_path{i_count,1} = flnrDpsi;
        
        fr_path{i_count,1} = fr;
        fz_path{i_count,1} = fz;
        dl_path{i_count,1} = dl;
        Interval_path = [Interval_path;Interval];
        
%         var_path = [var_path;var_temp];
        point_path = [point_path;r_end,z_end];
        edges_path = [edges_path;edges_next];
        tri_path = [tri_path;tri_next];
        
%         plot(point_path([end-1,end],1),point_path([end-1,end],2),'k');
    end
    pF{j}.fpsis = fpsi_path;
    pF{j}.fpsiDrs = fpsiDr_path;
    pF{j}.fpsiDzs = fpsiDz_path;
    pF{j}.fpsiGradNorms = fpsiGradNorm_path;
    pF{j}.fpsiGrad2Dpsis = fpsiGrad2Dpsi_path;
    pF{j}.flnpsiGradNormDpsis = flnpsiGradNormDpsi_path;
    pF{j}.flnrDpsis = flnrDpsi_path;
    
    pF{j}.frs = fr_path;
    pF{j}.fzs = fz_path;
    pF{j}.dls = dl_path;
    pF{j}.Intervals = Interval_path;
    
    pF{j}.vars = var_path;
    pF{j}.points = point_path;
end

    function [r_end,z_end] = FindIntesectionPoint(fpath,p1,p2)
        a = p2(2)-p1(2);
        b = -(p2(1)-p1(1));
        c = -a*p1(1)-b*p1(2);
        if b~=0
            fz_temp = @(r) -(a*r+c)/b;
            r_end = fzero(@(r) fpath(r,fz_temp(r)),[p1(1),p2(1)]);
            z_end = fz_temp(r_end);
        else
            fr_temp = @(z) -(b*z+c)/a;
            z_end = fzero(@(z) fpath(fr_temp(z),z),[p1(2),p2(2)]);
            r_end = fr_temp(z_end);
        end
    end

    function [fpsiDr_rz,fpsiDz_rz,fpsiGradNorm_rz,fpsiGrad2Dpsi_rz,flnpsiGradNormDpsi_rz,flnrDpsi_rz] = GetPsiDerivative(fpsi)
        c2 = fpsi.p10;c3 = fpsi.p01;
        c4 = fpsi.p20;c5 = fpsi.p11;c6 = fpsi.p02;
        fpsiDr_r2z2 = @(r,z) c2 + 2*c4*r + c5*z;
        fpsiDz_r2z2 = @(r,z) c3 + c5*r + 2*c6*z;
        fpsiDr_rz = @(r,z) 2*r.*fpsiDr_r2z2(r.^2,z.^2);
        fpsiDz_rz = @(r,z) 2*z.*fpsiDz_r2z2(r.^2,z.^2);
        fpsiGradNorm_rz = @(r,z) sqrt(fpsiDr_rz(r,z).^2+fpsiDz_rz(r,z).^2);
%         fpsiGrad2Dr_rz = @(r,z) (4*c4*(c2 + 2*c4*r + c5*z) + 2*c5*(c3 + c5*r + 2*c6*z))./(2*((c2 + 2*c4*r + c5*z).^2 + (c3 + c5*r + 2*c6*z).^2).^(1/2));
%         fpsiGrad2Dz_rz = @(r,z) (2*c5*(c2 + 2*c4*r + c5*z) + 4*c6*(c3 + c5*r + 2*c6*z))./(2*((c2 + 2*c4*r + c5*z).^2 + (c3 + c5*r + 2*c6*z).^2).^(1/2));
        fpsiGrad2Dr_rz = @(r,z) 8*r.*(2*c4*r.^2 + c5*z.^2 + c2).^2 + 32*c4*r.^3.*(2*c4*r.^2 + c5*z.^2 + c2) + 16*c5*r.*z.^2.*(c5*r.^2 + 2*c6*z.^2 + c3);
        fpsiGrad2Dz_rz = @(r,z) 8*z.*(c5*r.^2 + 2*c6*z.^2 + c3).^2 + 32*c6*z.^3.*(c5*r.^2 + 2*c6*z.^2 + c3) + 16*c5*r.^2.*z.*(2*c4*r.^2 + c5*z.^2 + c2);

        fpsiGrad2Dpsi_rz = @(r,z) (fpsiGrad2Dr_rz(r,z).*fpsiDr_rz(r,z)+fpsiGrad2Dz_rz(r,z).*fpsiDz_rz(r,z))./fpsiGradNorm_rz(r,z).^2;
        flnpsiGradNormDpsi_rz = @(r,z) fpsiGrad2Dpsi_rz(r,z)./(2*fpsiGradNorm_rz(r,z).^2);
        flnrDpsi_rz = @(r,z) 1./r.*fpsiDr_rz(r,z)./fpsiGradNorm_rz(r,z).^2;
        
        
        
    end

    function [fr,fz,dl,Interval] = SolvePathFunction(fpsi,psi,p1,p2)
        if abs(p2(1)-p1(1))<abs(p2(2)-p1(2))
            c1 = fpsi.p00-psi;c2 = fpsi.p10;c3 = fpsi.p01;
            c4 = fpsi.p20;c5 = fpsi.p11;c6 = fpsi.p02;
            fr2_1 = @(z) -(c2 + c5*z - ((c2 + c5*z).^2 - c4*(4*c6*z.^2 + 4*c3*z + 4*c1)).^(1/2))/(2*c4);
            fr2_2 = @(z) -(c2 + c5*z + ((c2 + c5*z).^2 - c4*(4*c6*z.^2 + 4*c3*z + 4*c1)).^(1/2))/(2*c4);
            fr1 = @(z) sqrt(fr2_1(z.^2));
            fr2 = @(z) sqrt(fr2_2(z.^2));
 
            if abs(fr1(p2(2))-p2(1))<(p2(1)*1e-6)
                fr = fr1;
                fr2dz2 = @(z) -(c5 + (c4*(4*c3 + 8*c6*z) - 2*c5*(c2 + c5*z))./(2*((c2 + c5*z).^2 - c4*(4*c6*z.^2 + 4*c3*z + 4*c1)).^(1/2)))/(2*c4);
            elseif abs(fr2(p2(2))-p2(1))<(p2(1)*1e-6)
                fr = fr2;
                fr2dz2 = @(z) -(c5 - (c4*(4*c3 + 8*c6*z) - 2*c5*(c2 + c5*z))./(2*((c2 + c5*z).^2 - c4*(4*c6*z.^2 + 4*c3*z + 4*c1)).^(1/2)))/(2*c4);
            else
                error('path function not exist');
            end
            fz = @(z) z;
            frdz = @(z) z./fr(z).*fr2dz2(z.^2);
            Interval = [p1(2),p2(2)];
            
            if p1(2)>p2(2)
                dl = @(z) -sqrt(1+frdz(z).^2);
            else
                dl = @(z) sqrt(1+frdz(z).^2);
            end
        else
            c1 = fpsi.p00-psi;c2 = fpsi.p10;c3 = fpsi.p01;
            c4 = fpsi.p20;c5 = fpsi.p11;c6 = fpsi.p02;
            fz2_1 = @(r) -(c3 + c5*r - ((c3 + c5*r).^2 - c6*(4*c4*r.^2 + 4*c2*r + 4*c1)).^(1/2))/(2*c6);
            fz2_2 = @(r) -(c3 + c5*r + ((c3 + c5*r).^2 - c6*(4*c4*r.^2 + 4*c2*r + 4*c1)).^(1/2))/(2*c6);
            fz1 = @(r) sqrt(fz2_1(r.^2));
            fz2 = @(r) sqrt(fz2_2(r.^2));
 
            if abs(fz1(p2(1))-p2(2))<(p2(2)*1e-6)
                fz = fz1;
                fz2dr2 = @(r) -(c5 + (c6*(4*c2 + 8*c4*r) - 2*c5*(c3 + c5*r))./(2*((c3 + c5*r).^2 - c6*(4*c4*r.^2 + 4*c2*r + 4*c1)).^(1/2)))/(2*c6);
            elseif abs(fz2(p2(1))-p2(2))<(p2(2)*1e-6)
                fz = fz2;
                fz2dr2 = @(r) -(c5 - (c6*(4*c2 + 8*c4*r) - 2*c5*(c3 + c5*r))./(2*((c3 + c5*r).^2 - c6*(4*c4*r.^2 + 4*c2*r + 4*c1)).^(1/2)))/(2*c6);
            else
                error('path function not exist');
            end
            fr = @(r) r;
            fzdr = @(r) r./fz(r).*fz2dr2(r.^2);
            Interval = [p1(1),p2(1)];
            if p1(1)>p2(1)
                dl = @(r) -sqrt(fzdr(r).^2+1);
            else
                dl = @(r) sqrt(fzdr(r).^2+1);
            end
        end
    end
        
end