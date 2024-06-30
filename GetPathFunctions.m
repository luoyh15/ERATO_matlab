function pF = GetPathFunctions(psi,TR,psidata,psiF)
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

% number of flux surface that need calculated
n_psi = size(psi,1);
% initialization
pF  = cell(size(psi));
% treat the region near magnetic axis as analytic
n_ana = length(psi(psi<psiF.psi_ana));
fpsi = psiF.analyticF;
rhomax_ana = psiF.rhomax_ana;
r_magaxis = psiF.r_magaxis;
% number of lines that for the segment path function
n_lines = 4;
for j = 1:n_ana
    psii = psi(j);
    %     fpath = @(rho,theta) fpsi((rho*cos(theta)+r_magaxis),(rho*sin(theta)))-psii;
    theta_total = linspace(0,pi,n_lines);
    % the first point of this path surface
    p_start = FindIntesectionPointAna(fpsi,psii,0,r_magaxis,rhomax_ana);
    point_path = p_start;
    % initialization
    fpsi_path = [];% list containing all the segment psi functions
    fpsiDr_path = [];
    fpsiDz_path = [];
    fpsiGradNorm_path = [];
    fpsiGrad2Dpsi_path = [];
    flnpsiGradNormDpsi_path = [];
    flnrDpsi_path = [];
    fr_path = [];
    fz_path = [];
    dl_path = [];
    dr_path = [];
    Interval_path = [];
    
    i_count = 0;% to count the number of segments
    for i= 2:n_lines
        theta_temp = theta_total(i);
        p_next = FindIntesectionPointAna(fpsi,psii,theta_temp,r_magaxis,rhomax_ana);
        % get the information of the current part of path function
        [fr,fz,dl,dr,Interval] = SolvePathFunctionAna(fpsi,psii,point_path(end,:),p_next);
        % get function fpsi's derivatives using in the calculating relevent quantities
        [fpsiDr,fpsiDz,fpsiGradNorm,fpsiGrad2Dpsi,flnpsiGradNormDpsi,flnrDpsi] = GetPsiDerivativeAna(fpsi);
        
        i_count = i_count+1;
        fpsi_path{i_count,1} = fpsi;
        fpsiDr_path{i_count,1} = fpsiDr;
        fpsiDz_path{i_count,1} = fpsiDz;
        fpsiGradNorm_path{i_count,1} = fpsiGradNorm;
        fpsiGrad2Dpsi_path{i_count,1} = fpsiGrad2Dpsi;
        flnpsiGradNormDpsi_path{i_count,1} = flnpsiGradNormDpsi;
        flnrDpsi_path{i_count,1} = flnrDpsi;
        fr_path{i_count,1} = fr;
        fz_path{i_count,1} = fz;
        dl_path{i_count,1} = dl;
        dr_path{i_count,1} = dr;
        Interval_path(i_count,:) = Interval;
        
        point_path = [point_path;p_next];% store next point
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
    pF{j}.drs = dr_path;
    pF{j}.Intervals = Interval_path;
    
    pF{j}.points = point_path;
end

%% get the path functions of one specific psi value not analytic
for j = n_ana+1:n_psi
% for j = 78:n_psi
    psii = psi(j);
    % find the triangular that containing psii value
    tri_next = [];% list containing all the triangular index pass through
    edge_next = [];% list containing all the edges that pass through
    % initialization
    fpsi_path = [];% list containing all the segment psi functions
    fpsiDr_path = [];
    fpsiDz_path = [];
    fpsiGradNorm_path = [];
    fpsiGrad2Dpsi_path = [];
    flnpsiGradNormDpsi_path = [];
    flnrDpsi_path = [];
    fr_path = [];
    fz_path = [];
    dl_path = [];
    dr_path = [];
    Interval_path = [];
    point_path = [];% list containing the start and end points of the corresponding line
    tri_midplane = psiF.tri_midplane;
    
    [fpsi,p_start,p_next,tri_next,edge_next] = FindFirstTriangle(psiF,TR,psidata,psii,tri_midplane);
    point_path = p_start;
    % get the information of the path function
    [fr,fz,dl,dr,Interval] = SolvePathFunction(fpsi,psii,point_path(end,:),p_next);
    % get function fpsi's derivatives using in the calculating relevent quantities
    [fpsiDr,fpsiDz,fpsiGradNorm,fpsiGrad2Dpsi,flnpsiGradNormDpsi,flnrDpsi] = GetPsiDerivative(fpsi);
    
    i_count = 1;% to count the number of triangulars in the path
    fpsi_path{i_count,1} = fpsi;
    fpsiDr_path{i_count,1} = fpsiDr;
    fpsiDz_path{i_count,1} = fpsiDz;
    fpsiGradNorm_path{i_count,1} = fpsiGradNorm;
    fpsiGrad2Dpsi_path{i_count,1} = fpsiGrad2Dpsi;
    flnpsiGradNormDpsi_path{i_count,1} = flnpsiGradNormDpsi;
    flnrDpsi_path{i_count,1} = flnrDpsi;
    fr_path{i_count,1} = fr;
    fz_path{i_count,1} = fz;
    dl_path{i_count,1} = dl;
    dr_path{i_count,1} = dr;
    Interval_path(i_count,:) = Interval;
    
    %                 var_path = [var_path;var_temp];
    point_path = [point_path;p_next];% store next point
    
    while 1
        %         if isempty(tri_next)
        %             error('can not find next triangle');
        %         end
        % get psi function by fit to 'poly22'
        fpsi = psiF.triangleF{tri_next};
        % get  next trianglar infomation
        [p_next,tri_next,edge_next] = FindNextTriangular(fpsi,psii,TR,psidata,tri_next,edge_next);
        % determine if is the last triangular
        if isempty(tri_next)
            % r as a function of z
            [fr,fz,dl,dr,Interval] = SolvePathFunction(fpsi,psii,point_path(end,:),p_next);
            % get function fpsi's derivatives using in the calculating relevent quantities
            [fpsiDr,fpsiDz,fpsiGradNorm,fpsiGrad2Dpsi,flnpsiGradNormDpsi,flnrDpsi] = GetPsiDerivative(fpsi);
            
            i_count = i_count+1;
            fpsi_path{i_count,1} = fpsi;
            fpsiDr_path{i_count,1} = fpsiDr;
            fpsiDz_path{i_count,1} = fpsiDz;
            fpsiGradNorm_path{i_count,1} = fpsiGradNorm;
            fpsiGrad2Dpsi_path{i_count,1} = fpsiGrad2Dpsi;
            flnpsiGradNormDpsi_path{i_count,1} = flnpsiGradNormDpsi;
            flnrDpsi_path{i_count,1} = flnrDpsi;
            
            fr_path{i_count,1} = fr;
            fz_path{i_count,1} = fz;
            dl_path{i_count,1} = dl;
            dr_path{i_count,1} = dr;
            Interval_path(i_count,:) = Interval;
            
            point_path = [point_path;p_next];
            
            break;
        end
        
        % determine which coordinate as integrat variable
        [fr,fz,dl,dr,Interval] = SolvePathFunction(fpsi,psii,point_path(end,:),p_next);
        % get function fpsi's derivatives using in the calculating relevent quantities
        [fpsiDr,fpsiDz,fpsiGradNorm,fpsiGrad2Dpsi,flnpsiGradNormDpsi,flnrDpsi] = GetPsiDerivative(fpsi);
        
        i_count = i_count+1;
        fpsi_path{i_count,1} = fpsi;
        fpsiDr_path{i_count,1} = fpsiDr;
        fpsiDz_path{i_count,1} = fpsiDz;
        fpsiGradNorm_path{i_count,1} = fpsiGradNorm;
        fpsiGrad2Dpsi_path{i_count,1} = fpsiGrad2Dpsi;
        flnpsiGradNormDpsi_path{i_count,1} = flnpsiGradNormDpsi;
        flnrDpsi_path{i_count,1} = flnrDpsi;
        
        fr_path{i_count,1} = fr;
        fz_path{i_count,1} = fz;
        dl_path{i_count,1} = dl;
        dr_path{i_count,1} = dr;
        Interval_path(i_count,:) = Interval;
        
        point_path = [point_path;p_next];
        
        
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
    pF{j}.drs = dr_path;
    pF{j}.Intervals = Interval_path;
    
    pF{j}.points = point_path;
end
end

function [fpsi,p_start,p_next,tri_next,edge_next] = FindFirstTriangle(psiF,TR,psidata,psii,tri_midplane)
tri_store = [];
line_store = [];
for i = 1:size(tri_midplane,1)
    p_indexes =  TR.ConnectivityList(tri_midplane(i),:);
    psi_values = psidata(p_indexes);
    if max(psi_values)>=psii&&min(psi_values)<psii% the triangular contains value 'psii'
        % fit psi function in the triangular
        %             fpsi = PsiFitFunctionOrder2(TR,psidata,tri_midplane(i));
        fpsi = psiF.triangleF{tri_midplane(i)};
        tri_index = tri_midplane(i);
        % three edge's intesect point with z=0, one is empty
        point1 = FindIntesectionPoint(@(r,z) z,TR.Points(p_indexes(1),:),TR.Points(p_indexes(2),:));
        point2 = FindIntesectionPoint(@(r,z) z,TR.Points(p_indexes(1),:),TR.Points(p_indexes(3),:));
        point3 = FindIntesectionPoint(@(r,z) z,TR.Points(p_indexes(2),:),TR.Points(p_indexes(3),:));
        line = unique([point1;point2;point3],'rows');
        if size(line,1)~=2% there is just one point
            continue;
        end
        tri_store = [tri_store;tri_midplane(i)];
        line_store = [line_store;{line}];
        % start point
        p_start = FindIntesectionPoint(@(r,z)fpsi(r,z)-psii,line(1,:),line(2,:));
        if ~isempty(p_start)% start point in the triangular
            % get fist edge that intersect with the line 'psi = psii'
            p1 = p_indexes(psi_values==max(psi_values));
            p2 = p_indexes(psi_values==min(psi_values));
            p3 = p_indexes(~((p_indexes==p1)|(p_indexes==p2)));
            if psidata(p3)>=psii
                if TR.Points(p1,2)>TR.Points(p3,2)
                    p_next = FindIntesectionPoint(@(r,z)fpsi(r,z)-psii,TR.Points(p1,:),TR.Points(p2,:));
                    edge_next = [p2,p1];% make sure psi value of the two points is in the order
                else
                    p_next = FindIntesectionPoint(@(r,z)fpsi(r,z)-psii,TR.Points(p3,:),TR.Points(p2,:));
                    edge_next = [p2,p3];
                end
            else
                if TR.Points(p2,2)>TR.Points(p3,2)
                    p_next = FindIntesectionPoint(@(r,z)fpsi(r,z)-psii,TR.Points(p1,:),TR.Points(p2,:));
                    edge_next = [p2,p1];% make sure psi value of the two points is in the order
                else
                    p_next = FindIntesectionPoint(@(r,z)fpsi(r,z)-psii,TR.Points(p1,:),TR.Points(p3,:));
                    edge_next = [p3,p1];
                end
            end
            %                 plot(point_path([end-1,end],1),point_path([end-1,end],2),'k');
            break;
        end
    end
end
if isempty(p_start)% need to use linear fit to calculate the start point
    for i = 1:length(tri_store)
        fpsi = psiF.triangleF{tri_store(i)};
        tri_index = tri_store(i);
        p_indexes = TR.ConnectivityList(tri_store(i),:);
        psi_values = psidata(p_indexes);
        fpsi_linear = PsiFitFunctionOrder1(TR,psidata,tri_store(i));
        line = line_store{i};
        % start point
        p_start = FindIntesectionPoint(@(r,z)fpsi_linear(r,z)-psii,line(1,:),line(2,:));
        if ~isempty(p_start)% start point in the triangular
            % get fist edge that intersect with the line 'psi = psii'
            p1 = p_indexes(psi_values==max(psi_values));
            p2 = p_indexes(psi_values==min(psi_values));
            p3 = p_indexes(~((p_indexes==p1)|(p_indexes==p2)));
            if psidata(p3)>=psii
                if TR.Points(p1,2)>TR.Points(p3,2)
                    p_next = FindIntesectionPoint(@(r,z)fpsi_linear(r,z)-psii,TR.Points(p1,:),TR.Points(p2,:));
                    edge_next = [p2,p1];% make sure psi value of the two points is in the order
                else
                    p_next = FindIntesectionPoint(@(r,z)fpsi_linear(r,z)-psii,TR.Points(p3,:),TR.Points(p2,:));
                    edge_next = [p2,p3];
                end
            else
                if TR.Points(p2,2)>TR.Points(p3,2)
                    p_next = FindIntesectionPoint(@(r,z)fpsi_linear(r,z)-psii,TR.Points(p1,:),TR.Points(p2,:));
                    edge_next = [p2,p1];% make sure psi value of the two points is in the order
                else
                    p_next = FindIntesectionPoint(@(r,z)fpsi_linear(r,z)-psii,TR.Points(p1,:),TR.Points(p3,:));
                    edge_next = [p3,p1];
                end
            end
            %                 plot(point_path([end-1,end],1),point_path([end-1,end],2),'k');
            break;
        end
    end
end
if isempty(p_start)
    error('no start point found.')
end
% find next triangular that attached to the last edge
tri_attached = TR.edgeAttachments(edge_next(end,:));
tri_attached = tri_attached{1};
if tri_attached(1)==tri_index
    tri_next = tri_attached(2);
else
    tri_next = tri_attached(1);
end
end

function fpsi = PsiFitFunctionOrder1(TR,psidata,tri_index)
p_total = TR.ConnectivityList(tri_index,:);
p_fit_index = p_total;
psi_fit_value = psidata(p_fit_index);
p_fit_coordinates = TR.Points(p_fit_index,:);
fpsi = fit([p_fit_coordinates(:,1),p_fit_coordinates(:,2)],psi_fit_value,'poly11');
end

function p = FindIntesectionPoint(fpath,p1,p2)
a = p2(2)-p1(2);
b = -(p2(1)-p1(1));
c = -a*p1(1)-b*p1(2);
if a==0&&b==0
    error('please input two different points of a line');
elseif abs(b)>abs(a)
    fz_temp = @(r) -(a*r+c)/b;
    % if there is no intersection between fpath and line[p1,p2]
    if fpath(p1(1),fz_temp(p1(1)))*fpath(p2(1),fz_temp(p2(1)))>0
        p = [];
    else
        r = fzero(@(r) fpath(r,fz_temp(r)),[p1(1),p2(1)]);
        z = fz_temp(r);
        p = [r,z];
    end
else
    fr_temp = @(z) -(b*z+c)/a;
    if fpath(fr_temp(p1(2)),p1(2))*fpath(fr_temp(p2(2)),p2(2))>0
        p = [];
    else
        z = fzero(@(z) fpath(fr_temp(z),z),[p1(2),p2(2)]);
        r = fr_temp(z);
        p = [r,z];
    end
end
end

function  p = FindIntesectionPointAna(fpsi,psii,theta,r_magaxis,rhomax_ana)
fpath = @(rho,theta) fpsi((rho*cos(theta)+r_magaxis).^2,(rho*sin(theta)).^2)-psii;
rho = fzero(@(rho) fpath(rho,theta),[0,rhomax_ana]);
r = rho*cos(theta)+r_magaxis;
z = rho*sin(theta);
p = [r,z];
end
function [p_next,tri_next,edge_next] = FindNextTriangular(fpsi,psii,TR,psidata,tri_index,edge)
p_indexes = TR.ConnectivityList(tri_index,:);
p1 = edge(2);% the bigger psi value point
p2 = edge(1);% the smaller psi value point
p3 = p_indexes(~((p_indexes==p1)|(p_indexes==p2)));
if ~(TR.Points(p1,2)>0&&TR.Points(p2,2)>0&&TR.Points(p3,2)>0)%
    % three edge's intesect point with z=0, one is empty
    Ip1 = FindIntesectionPoint(@(r,z) z,TR.Points(p1,:),TR.Points(p2,:));
    Ip2 = FindIntesectionPoint(@(r,z) z,TR.Points(p1,:),TR.Points(p3,:));
    Ip3 = FindIntesectionPoint(@(r,z) z,TR.Points(p2,:),TR.Points(p3,:));
    Iline = unique([Ip1;Ip2;Ip3],'rows');
    if size(Iline,1)==2
        % end point
        p_end = FindIntesectionPoint(@(r,z)fpsi(r,z)-psii,Iline(1,:),Iline(2,:));
        if ~isempty(p_end)% the last triangular encounted
            tri_next = [];
            edge_next = [];
            p_next = p_end;
            return;
        end
    end
end
if psidata(p3)>=psii
    p_next = FindIntesectionPoint(@(r,z)fpsi(r,z)-psii,TR.Points(p3,:),TR.Points(p2,:));
    %                 if isempty(p_next)% need to use first order to calculate the intersection
    %                     p_coordinates = TR.Points(p_indexes,:);
    %                     fpsi_linear = fit([p_coordinates(:,1),p_coordinates(:,2)],psi_values,'poly11');
    %                     p_next = FindIntesectionPoint(@(r,z)fpsi_linear(r,z)-psii,TR.Points(p3,:),TR.Points(p2,:));
    %                 end
    edge_next = [p2,p3];
else
    p_next = FindIntesectionPoint(@(r,z)fpsi(r,z)-psii,TR.Points(p1,:),TR.Points(p3,:));
    %                 if isempty(p_next)% need to use first order to calculate the intersection
    %                     p_coordinates = TR.Points(p_indexes,:);
    %                     fpsi_linear = fit([p_coordinates(:,1),p_coordinates(:,2)],psi_values,'poly11');
    %                     p_next = FindIntesectionPoint(@(r,z)fpsi_linear(r,z)-psii,TR.Points(p1,:),TR.Points(p3,:));
    %                 end
    edge_next = [p3,p1];
end
% find next triangular that attached to the last edge
tri_attached = TR.edgeAttachments(edge_next(end,:));
tri_attached = tri_attached{1};
if tri_attached(1)==tri_index
    tri_next = tri_attached(2);
else
    tri_next = tri_attached(1);
end
end

function [fpsiDr_rz,fpsiDz_rz,fpsiGradNorm_rz,fpsiGrad2Dpsi_rz,flnpsiGradNormDpsi_rz,flnrDpsi_rz] = GetPsiDerivativeAna(fpsi)
c2 = fpsi.p10;c3 = fpsi.p01;
c4 = fpsi.p20;c5 = fpsi.p11;c6 = fpsi.p02;
c7 = fpsi.p30;c8 = fpsi.p21;c9 = fpsi.p12;c10 = fpsi.p03;
fpsiDr_rz = @(r,z)c2.*r.*2.0+c4.*r.^3.*4.0+c7.*r.^5.*6.0+c5.*r.*z.^2.*2.0+c9.*r.*z.^4.*2.0+c8.*r.^3.*z.^2.*4.0;
fpsiDz_rz = @(r,z)c3.*z.*2.0+c6.*z.^3.*4.0+c10.*z.^5.*6.0+c5.*r.^2.*z.*2.0+c8.*r.^4.*z.*2.0+c9.*r.^2.*z.^3.*4.0;
fpsiGradNorm_rz = @(r,z) sqrt(fpsiDr_rz(r,z).^2+fpsiDz_rz(r,z).^2);
fpsiGrad2Dr_rz = @(r,z)(c2.*r.*2.0+c4.*r.^3.*4.0+c7.*r.^5.*6.0+c5.*r.*z.^2.*2.0+c9.*r.*z.^4.*2.0+c8.*r.^3.*z.^2.*4.0).*(c2.*2.0+c4.*r.^2.*1.2e1+c7.*r.^4.*3.0e1+c5.*z.^2.*2.0+c9.*z.^4.*2.0+c8.*r.^2.*z.^2.*1.2e1).*2.0+(c8.*r.^3.*z.*8.0+c9.*r.*z.^3.*8.0+c5.*r.*z.*4.0).*(c3.*z.*2.0+c6.*z.^3.*4.0+c10.*z.^5.*6.0+c5.*r.^2.*z.*2.0+c8.*r.^4.*z.*2.0+c9.*r.^2.*z.^3.*4.0).*2.0;
fpsiGrad2Dz_rz = @(r,z)(c3.*z.*2.0+c6.*z.^3.*4.0+c10.*z.^5.*6.0+c5.*r.^2.*z.*2.0+c8.*r.^4.*z.*2.0+c9.*r.^2.*z.^3.*4.0).*(c3.*2.0+c5.*r.^2.*2.0+c8.*r.^4.*2.0+c6.*z.^2.*1.2e1+c10.*z.^4.*3.0e1+c9.*r.^2.*z.^2.*1.2e1).*2.0+(c8.*r.^3.*z.*8.0+c9.*r.*z.^3.*8.0+c5.*r.*z.*4.0).*(c2.*r.*2.0+c4.*r.^3.*4.0+c7.*r.^5.*6.0+c5.*r.*z.^2.*2.0+c9.*r.*z.^4.*2.0+c8.*r.^3.*z.^2.*4.0).*2.0;

fpsiGrad2Dpsi_rz = @(r,z) (fpsiGrad2Dr_rz(r,z).*fpsiDr_rz(r,z)+fpsiGrad2Dz_rz(r,z).*fpsiDz_rz(r,z))./fpsiGradNorm_rz(r,z).^2;
flnpsiGradNormDpsi_rz = @(r,z) fpsiGrad2Dpsi_rz(r,z)./(2*fpsiGradNorm_rz(r,z).^2);
flnrDpsi_rz = @(r,z) 1./r.*fpsiDr_rz(r,z)./fpsiGradNorm_rz(r,z).^2;
end
function [fpsiDr_rz,fpsiDz_rz,fpsiGradNorm_rz,fpsiGrad2Dpsi_rz,flnpsiGradNormDpsi_rz,flnrDpsi_rz] = GetPsiDerivative(fpsi)
c2 = fpsi.p10;c3 = fpsi.p01;
c4 = fpsi.p20;c5 = fpsi.p11;c6 = fpsi.p02;
fpsiDr_rz = @(r,z) c2 + 2*c4*r + c5*z;
fpsiDz_rz = @(r,z) c3 + c5*r + 2*c6*z;
fpsiGradNorm_rz = @(r,z) sqrt(fpsiDr_rz(r,z).^2+fpsiDz_rz(r,z).^2);
fpsiGrad2Dr_rz = @(r,z) 4*c4*(c2 + 2*c4*r + c5*z) + 2*c5*(c3 + c5*r + 2*c6*z);
fpsiGrad2Dz_rz = @(r,z) 2*c5*(c2 + 2*c4*r + c5*z) + 4*c6*(c3 + c5*r + 2*c6*z);

fpsiGrad2Dpsi_rz = @(r,z) (fpsiGrad2Dr_rz(r,z).*fpsiDr_rz(r,z)+fpsiGrad2Dz_rz(r,z).*fpsiDz_rz(r,z))./fpsiGradNorm_rz(r,z).^2;
flnpsiGradNormDpsi_rz = @(r,z) fpsiGrad2Dpsi_rz(r,z)./(2*fpsiGradNorm_rz(r,z).^2);
flnrDpsi_rz = @(r,z) 1./r.*fpsiDr_rz(r,z)./fpsiGradNorm_rz(r,z).^2;
end

function [fr,fz,dl,dr,Interval] = SolvePathFunctionAna(fpsi,psi,p1,p2)
% make sure the fit function is in the form of 'poly33'
c1 = fpsi.p00-psi;c2 = fpsi.p10;c3 = fpsi.p01;
c4 = fpsi.p20;c5 = fpsi.p11;c6 = fpsi.p02;
c7 = fpsi.p30;c8 = fpsi.p21;c9 = fpsi.p12;c10 = fpsi.p03;
if abs(p2(1)-p1(1))<abs(p2(2)-p1(2))
    % the three roots of the 'poly33' function
    fr2_1 = @(z) (c4.*(-1.0./3.0)-c8.*z.*(1.0./3.0))./c7+(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).*1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(1.0./3.0)+((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(1.0./3.0);
    fr2_2 = @(z) (c4.*(-1.0./3.0)-c8.*z.*(1.0./3.0))./c7-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).*1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(1.0./3.0).*(1.0./2.0)-((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(1.0./3.0).*(1.0./2.0)-sqrt(3.0).*((1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).*1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(1.0./3.0)-((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(1.0./3.0)).*5.0e-1i;
    fr2_3 = @(z) (c4.*(-1.0./3.0)-c8.*z.*(1.0./3.0))./c7-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).*1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(1.0./3.0).*(1.0./2.0)-((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(1.0./3.0).*(1.0./2.0)+sqrt(3.0).*((1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).*1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(1.0./3.0)-((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(1.0./3.0)).*5.0e-1i;
    fr1 = @(z) sqrt(fr2_1(z.^2));
    fr2 = @(z) sqrt(fr2_2(z.^2));
    fr3 = @(z) sqrt(fr2_3(z.^2));
    if abs(fr1(p2(2))-p2(1))<(p2(1)*1e-6)
        fr = fr1;
        fr2dz2 = @(z) -((c5.*(1.0./3.0)+c9.*z.*(2.0./3.0))./c7-1.0./c7.^2.*c8.*(c4+c8.*z).*(2.0./9.0)).*1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(1.0./3.0)+1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(2.0./3.0).*(1.0./sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3).*(((c5.*(1.0./3.0)+c9.*z.*(2.0./3.0))./c7-1.0./c7.^2.*c8.*(c4+c8.*z).*(2.0./9.0)).*(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^2.*3.0+((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).*((c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7-1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)-1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)+1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*2.0).*(1.0./2.0)-(c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7+1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)+1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)-1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*(1.0./3.0)-(c8.*(1.0./3.0))./c7-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).*1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(4.0./3.0).*(1.0./sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3).*(((c5.*(1.0./3.0)+c9.*z.*(2.0./3.0))./c7-1.0./c7.^2.*c8.*(c4+c8.*z).*(2.0./9.0)).*(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^2.*3.0+((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).*((c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7-1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)-1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)+1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*2.0).*(1.0./2.0)-(c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7+1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)+1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)-1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*(1.0./3.0);
    elseif abs(fr2(p2(2))-p2(1))<(p2(1)*1e-6)
        fr = fr2;
        fr2dz2 = @(z) ((c5.*(1.0./3.0)+c9.*z.*(2.0./3.0))./c7-1.0./c7.^2.*c8.*(c4+c8.*z).*(2.0./9.0)).*1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(1.0./3.0).*(1.0./2.0)-1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(2.0./3.0).*(1.0./sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3).*(((c5.*(1.0./3.0)+c9.*z.*(2.0./3.0))./c7-1.0./c7.^2.*c8.*(c4+c8.*z).*(2.0./9.0)).*(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^2.*3.0+((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).*((c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7-1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)-1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)+1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*2.0).*(1.0./2.0)-(c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7+1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)+1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)-1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*(1.0./6.0)-(c8.*(1.0./3.0))./c7+sqrt(3.0).*(((c5.*(1.0./3.0)+c9.*z.*(2.0./3.0))./c7-1.0./c7.^2.*c8.*(c4+c8.*z).*(2.0./9.0)).*1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(1.0./3.0)+1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(2.0./3.0).*(1.0./sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3).*(((c5.*(1.0./3.0)+c9.*z.*(2.0./3.0))./c7-1.0./c7.^2.*c8.*(c4+c8.*z).*(2.0./9.0)).*(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^2.*3.0+((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).*((c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7-1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)-1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)+1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*2.0).*(1.0./2.0)-(c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7+1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)+1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)-1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*(1.0./3.0)+(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).*1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(4.0./3.0).*(1.0./sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3).*(((c5.*(1.0./3.0)+c9.*z.*(2.0./3.0))./c7-1.0./c7.^2.*c8.*(c4+c8.*z).*(2.0./9.0)).*(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^2.*3.0+((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).*((c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7-1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)-1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)+1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*2.0).*(1.0./2.0)-(c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7+1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)+1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)-1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*(1.0./3.0)).*5.0e-1i+(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).*1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(4.0./3.0).*(1.0./sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3).*(((c5.*(1.0./3.0)+c9.*z.*(2.0./3.0))./c7-1.0./c7.^2.*c8.*(c4+c8.*z).*(2.0./9.0)).*(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^2.*3.0+((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).*((c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7-1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)-1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)+1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*2.0).*(1.0./2.0)-(c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7+1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)+1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)-1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*(1.0./6.0);
    elseif abs(fr3(p2(2))-p2(1))<(p2(1)*1e-6)
        fr = fr3;
        fr2dz2 = @(z) ((c5.*(1.0./3.0)+c9.*z.*(2.0./3.0))./c7-1.0./c7.^2.*c8.*(c4+c8.*z).*(2.0./9.0)).*1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(1.0./3.0).*(1.0./2.0)-1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(2.0./3.0).*(1.0./sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3).*(((c5.*(1.0./3.0)+c9.*z.*(2.0./3.0))./c7-1.0./c7.^2.*c8.*(c4+c8.*z).*(2.0./9.0)).*(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^2.*3.0+((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).*((c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7-1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)-1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)+1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*2.0).*(1.0./2.0)-(c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7+1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)+1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)-1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*(1.0./6.0)-(c8.*(1.0./3.0))./c7-sqrt(3.0).*(((c5.*(1.0./3.0)+c9.*z.*(2.0./3.0))./c7-1.0./c7.^2.*c8.*(c4+c8.*z).*(2.0./9.0)).*1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(1.0./3.0)+1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(2.0./3.0).*(1.0./sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3).*(((c5.*(1.0./3.0)+c9.*z.*(2.0./3.0))./c7-1.0./c7.^2.*c8.*(c4+c8.*z).*(2.0./9.0)).*(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^2.*3.0+((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).*((c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7-1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)-1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)+1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*2.0).*(1.0./2.0)-(c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7+1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)+1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)-1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*(1.0./3.0)+(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).*1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(4.0./3.0).*(1.0./sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3).*(((c5.*(1.0./3.0)+c9.*z.*(2.0./3.0))./c7-1.0./c7.^2.*c8.*(c4+c8.*z).*(2.0./9.0)).*(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^2.*3.0+((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).*((c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7-1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)-1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)+1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*2.0).*(1.0./2.0)-(c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7+1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)+1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)-1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*(1.0./3.0)).*5.0e-1i+(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).*1.0./((c1.*(-1.0./2.0)-c3.*z.*(1.0./2.0)-c6.*z.^2.*(1.0./2.0)-c10.*z.^3.*(1.0./2.0))./c7+sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3)-1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)+1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^(4.0./3.0).*(1.0./sqrt(((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).^2-(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^3).*(((c5.*(1.0./3.0)+c9.*z.*(2.0./3.0))./c7-1.0./c7.^2.*c8.*(c4+c8.*z).*(2.0./9.0)).*(1.0./c7.^2.*(c4+c8.*z).^2.*(1.0./9.0)-(c2.*(1.0./3.0)+c5.*z.*(1.0./3.0)+c9.*z.^2.*(1.0./3.0))./c7).^2.*3.0+((c1.*(1.0./2.0)+c3.*z.*(1.0./2.0)+c6.*z.^2.*(1.0./2.0)+c10.*z.^3.*(1.0./2.0))./c7+1.0./c7.^3.*(c4+c8.*z).^3.*(1.0./2.7e1)-1.0./c7.^2.*(c4+c8.*z).*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)).*((c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7-1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)-1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)+1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*2.0).*(1.0./2.0)-(c3.*(1.0./2.0)+c6.*z+c10.*z.^2.*(3.0./2.0))./c7+1.0./c7.^2.*(c4+c8.*z).*(c5+c9.*z.*2.0).*(1.0./6.0)+1.0./c7.^2.*c8.*(c2+c5.*z+c9.*z.^2).*(1.0./6.0)-1.0./c7.^3.*c8.*(c4+c8.*z).^2.*(1.0./9.0)).*(1.0./6.0);
    else
        error('path function not exist');
    end
    frdz = @(z) z./fr(z).*fr2dz2(z.^2);
    fz = @(z) z;
    dr = frdz;
    Interval = [p1(2),p2(2)];
    
    if p1(2)>p2(2)
        dl = @(z) -sqrt(1+frdz(z).^2);
    else
        dl = @(z) sqrt(1+frdz(z).^2);
    end
else
    fz2_1 = @(r)(c6.*(-1.0./3.0)-c9.*r.*(1.0./3.0))./c10+(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).*1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(1.0./3.0)+(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(1.0./3.0);
    fz2_2 = @(r)(c6.*(-1.0./3.0)-c9.*r.*(1.0./3.0))./c10-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).*1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(1.0./3.0).*(1.0./2.0)-(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(1.0./3.0).*(1.0./2.0)-sqrt(3.0).*((1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).*1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(1.0./3.0)-(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(1.0./3.0)).*5.0e-1i;
    fz2_3 = @(r)(c6.*(-1.0./3.0)-c9.*r.*(1.0./3.0))./c10-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).*1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(1.0./3.0).*(1.0./2.0)-(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(1.0./3.0).*(1.0./2.0)+sqrt(3.0).*((1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).*1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(1.0./3.0)-(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(1.0./3.0)).*5.0e-1i;
    fz1 = @(r) sqrt(fz2_1(r.^2));
    fz2 = @(r) sqrt(fz2_2(r.^2));
    fz3 = @(r) sqrt(fz2_3(r.^2));
    if abs(fz1(p2(1))-p2(2))<(p2(2)*1e-6)
        fz = fz1;
        fz2dr2 = @(r)(c9.*(-1.0./3.0))./c10-((c5.*(1.0./3.0)+c8.*r.*(2.0./3.0))./c10-c9.*1.0./c10.^2.*(c6+c9.*r).*(2.0./9.0)).*1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(1.0./3.0)+1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(2.0./3.0).*((((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).*((c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10-1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)-c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)+c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*2.0+(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^2.*((c5.*(1.0./3.0)+c8.*r.*(2.0./3.0))./c10-c9.*1.0./c10.^2.*(c6+c9.*r).*(2.0./9.0)).*3.0).*1.0./sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3).*(1.0./2.0)-(c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10+1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)+c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)-c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*(1.0./3.0)-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).*1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(4.0./3.0).*((((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).*((c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10-1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)-c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)+c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*2.0+(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^2.*((c5.*(1.0./3.0)+c8.*r.*(2.0./3.0))./c10-c9.*1.0./c10.^2.*(c6+c9.*r).*(2.0./9.0)).*3.0).*1.0./sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3).*(1.0./2.0)-(c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10+1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)+c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)-c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*(1.0./3.0);
    elseif abs(fz2(p2(1))-p2(2))<(p2(2)*1e-6)
        fz = fz2;
        fz2dr2 = @(r)(c9.*(-1.0./3.0))./c10+((c5.*(1.0./3.0)+c8.*r.*(2.0./3.0))./c10-c9.*1.0./c10.^2.*(c6+c9.*r).*(2.0./9.0)).*1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(1.0./3.0).*(1.0./2.0)+sqrt(3.0).*(((c5.*(1.0./3.0)+c8.*r.*(2.0./3.0))./c10-c9.*1.0./c10.^2.*(c6+c9.*r).*(2.0./9.0)).*1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(1.0./3.0)+1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(2.0./3.0).*((((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).*((c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10-1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)-c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)+c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*2.0+(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^2.*((c5.*(1.0./3.0)+c8.*r.*(2.0./3.0))./c10-c9.*1.0./c10.^2.*(c6+c9.*r).*(2.0./9.0)).*3.0).*1.0./sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3).*(1.0./2.0)-(c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10+1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)+c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)-c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*(1.0./3.0)+(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).*1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(4.0./3.0).*((((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).*((c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10-1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)-c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)+c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*2.0+(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^2.*((c5.*(1.0./3.0)+c8.*r.*(2.0./3.0))./c10-c9.*1.0./c10.^2.*(c6+c9.*r).*(2.0./9.0)).*3.0).*1.0./sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3).*(1.0./2.0)-(c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10+1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)+c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)-c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*(1.0./3.0)).*5.0e-1i-1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(2.0./3.0).*((((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).*((c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10-1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)-c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)+c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*2.0+(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^2.*((c5.*(1.0./3.0)+c8.*r.*(2.0./3.0))./c10-c9.*1.0./c10.^2.*(c6+c9.*r).*(2.0./9.0)).*3.0).*1.0./sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3).*(1.0./2.0)-(c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10+1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)+c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)-c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*(1.0./6.0)+(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).*1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(4.0./3.0).*((((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).*((c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10-1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)-c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)+c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*2.0+(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^2.*((c5.*(1.0./3.0)+c8.*r.*(2.0./3.0))./c10-c9.*1.0./c10.^2.*(c6+c9.*r).*(2.0./9.0)).*3.0).*1.0./sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3).*(1.0./2.0)-(c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10+1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)+c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)-c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*(1.0./6.0);
    elseif abs(fz3(p2(1))-p2(2))<(p2(2)*1e-6)
        fz = fz3;
        fz2dr2 = @(r)(c9.*(-1.0./3.0))./c10+((c5.*(1.0./3.0)+c8.*r.*(2.0./3.0))./c10-c9.*1.0./c10.^2.*(c6+c9.*r).*(2.0./9.0)).*1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(1.0./3.0).*(1.0./2.0)-sqrt(3.0).*(((c5.*(1.0./3.0)+c8.*r.*(2.0./3.0))./c10-c9.*1.0./c10.^2.*(c6+c9.*r).*(2.0./9.0)).*1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(1.0./3.0)+1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(2.0./3.0).*((((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).*((c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10-1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)-c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)+c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*2.0+(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^2.*((c5.*(1.0./3.0)+c8.*r.*(2.0./3.0))./c10-c9.*1.0./c10.^2.*(c6+c9.*r).*(2.0./9.0)).*3.0).*1.0./sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3).*(1.0./2.0)-(c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10+1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)+c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)-c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*(1.0./3.0)+(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).*1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(4.0./3.0).*((((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).*((c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10-1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)-c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)+c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*2.0+(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^2.*((c5.*(1.0./3.0)+c8.*r.*(2.0./3.0))./c10-c9.*1.0./c10.^2.*(c6+c9.*r).*(2.0./9.0)).*3.0).*1.0./sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3).*(1.0./2.0)-(c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10+1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)+c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)-c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*(1.0./3.0)).*5.0e-1i-1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(2.0./3.0).*((((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).*((c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10-1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)-c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)+c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*2.0+(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^2.*((c5.*(1.0./3.0)+c8.*r.*(2.0./3.0))./c10-c9.*1.0./c10.^2.*(c6+c9.*r).*(2.0./9.0)).*3.0).*1.0./sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3).*(1.0./2.0)-(c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10+1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)+c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)-c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*(1.0./6.0)+(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).*1.0./(sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3)-(c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10-1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)+1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^(4.0./3.0).*((((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).*((c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10-1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)-c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)+c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*2.0+(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^2.*((c5.*(1.0./3.0)+c8.*r.*(2.0./3.0))./c10-c9.*1.0./c10.^2.*(c6+c9.*r).*(2.0./9.0)).*3.0).*1.0./sqrt(((c1.*(1.0./2.0)+c2.*r.*(1.0./2.0)+c4.*r.^2.*(1.0./2.0)+c7.*r.^3.*(1.0./2.0))./c10+1.0./c10.^3.*(c6+c9.*r).^3.*(1.0./2.7e1)-1.0./c10.^2.*(c6+c9.*r).*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)).^2-(1.0./c10.^2.*(c6+c9.*r).^2.*(1.0./9.0)-(c3.*(1.0./3.0)+c5.*r.*(1.0./3.0)+c8.*r.^2.*(1.0./3.0))./c10).^3).*(1.0./2.0)-(c2.*(1.0./2.0)+c4.*r+c7.*r.^2.*(3.0./2.0))./c10+1.0./c10.^2.*(c5+c8.*r.*2.0).*(c6+c9.*r).*(1.0./6.0)+c9.*1.0./c10.^2.*(c3+c5.*r+c8.*r.^2).*(1.0./6.0)-c9.*1.0./c10.^3.*(c6+c9.*r).^2.*(1.0./9.0)).*(1.0./6.0);
    else
        error('path function not exist');
    end
    fzdr = @(r) r./fz(r).*fz2dr2(r.^2);
    fr = @(r) r;
    dr = @(r) 1;
    Interval = [p1(1),p2(1)];
    
    if p1(1)>p2(1)
        dl = @(r) -sqrt(1+fzdr(r).^2);
    else
        dl = @(r) sqrt(1+fzdr(r).^2);
    end
end
end

function [fr,fz,dl,dr,Interval] = SolvePathFunction(fpsi,psi,p1,p2)
% make sure the fit function is in the form of 'poly22'
c1 = fpsi.p00-psi;c2 = fpsi.p10;c3 = fpsi.p01;
c4 = fpsi.p20;c5 = fpsi.p11;c6 = fpsi.p02;
if abs(p2(1)-p1(1))<abs(p2(2)-p1(2))
    if abs(p2(2)-p1(2))<1e-4% two points are too close
        fr = @(z) p2(1);
        frdz = @(z) 0;
    elseif c4==0
        fr = @(z)-(c1+c3.*z+c6.*z.^2)./(c2+c5.*z);
        frdz = @(z)-(c3+c6.*z.*2.0)./(c2+c5.*z)+c5.*1.0./(c2+c5.*z).^2.*(c1+c3.*z+c6.*z.^2);
    else
        fr1 = @(z) -(c2 + c5*z - ((c2 + c5*z).^2 - c4*(4*c6*z.^2 + 4*c3*z + 4*c1)).^(1/2))/(2*c4);
        fr2 = @(z) -(c2 + c5*z + ((c2 + c5*z).^2 - c4*(4*c6*z.^2 + 4*c3*z + 4*c1)).^(1/2))/(2*c4);
        
        if abs(fr1(p2(2))-p2(1))<(p2(1)*1e-3)
            fr = fr1;
            frdz = @(z) -(c5 + (c4*(4*c3 + 8*c6*z) - 2*c5*(c2 + c5*z))./(2*((c2 + c5*z).^2 - c4*(4*c6*z.^2 + 4*c3*z + 4*c1)).^(1/2)))/(2*c4);
        elseif abs(fr2(p2(2))-p2(1))<(p2(1)*1e-3)
            fr = fr2;
            frdz = @(z) -(c5 - (c4*(4*c3 + 8*c6*z) - 2*c5*(c2 + c5*z))./(2*((c2 + c5*z).^2 - c4*(4*c6*z.^2 + 4*c3*z + 4*c1)).^(1/2)))/(2*c4);
        else
            error('path function not exist');
        end
    end
    fz = @(z) z;
    dr = frdz;
    Interval = [p1(2),p2(2)];
    
    if p1(2)>p2(2)
        dl = @(z) -sqrt(1+frdz(z).^2);
    else
        dl = @(z) sqrt(1+frdz(z).^2);
    end
else
    if abs(p2(1)-p1(1))<1e-4% two points are too close
        fz = @(r) p2(2);
        fzdr = @(r) 0;
    elseif c6==0
        fz = @(r)-(c1+c2.*r+c4.*r.^2)./(c3+c5.*r);
        fzdr = @(r)-(c2+c4.*r.*2.0)./(c3+c5.*r)+c5.*1.0./(c3+c5.*r).^2.*(c1+c2.*r+c4.*r.^2);
    else
        fz1 = @(r) -(c3 + c5*r - ((c3 + c5*r).^2 - c6*(4*c4*r.^2 + 4*c2*r + 4*c1)).^(1/2))/(2*c6);
        fz2 = @(r) -(c3 + c5*r + ((c3 + c5*r).^2 - c6*(4*c4*r.^2 + 4*c2*r + 4*c1)).^(1/2))/(2*c6);
        
        if abs(fz1(p2(1))-p2(2))<(p2(2)*1e-3)
            fz = fz1;
            fzdr = @(r) -(c5 + (c6*(4*c2 + 8*c4*r) - 2*c5*(c3 + c5*r))./(2*((c3 + c5*r).^2 - c6*(4*c4*r.^2 + 4*c2*r + 4*c1)).^(1/2)))/(2*c6);
        elseif abs(fz2(p2(1))-p2(2))<(p2(2)*1e-3)
            fz = fz2;
            fzdr = @(r) -(c5 - (c6*(4*c2 + 8*c4*r) - 2*c5*(c3 + c5*r))./(2*((c3 + c5*r).^2 - c6*(4*c4*r.^2 + 4*c2*r + 4*c1)).^(1/2)))/(2*c6);
        else
            error('path function not exist');
        end
    end
    fr = @(r) r;
    dr = @(r) 1;
    Interval = [p1(1),p2(1)];
    if p1(1)>p2(1)
        dl = @(r) -sqrt(fzdr(r).^2+1);
    else
        dl = @(r) sqrt(fzdr(r).^2+1);
    end
end
end