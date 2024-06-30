function GrowthRate = ERATO5(Es,as,q0s,ns)
%% set global variable
global R E a Psi_s
%% the parameters determine the equilibrium
R = 1;  B = 1; T0 = B*R; gamma = 5/3;
%%%% parameters determine the equilibrium
% E = 2; a = 1/3; q0 = 0.8; n = 1;
E = Es; a = as; q0 = q0s; n = ns;

Psi_s = E*a^2*B/(q0*2); %the psi value on the plasma surface
% the boundary of plasma
r_min = sqrt(R^2-2*a*R); r_max = sqrt(R^2+2*a*R);
r_temp = R*(1-4*a^2/R^2)^(1/4);
z_max = sqrt((4*R^2*a^2-(r_temp^2-R^2)^2)*E^2/(4*r_temp^2));
rho_max = sqrt(max(R-r_min,r_max-R)^2+z_max^2);
%% r z coordinate system
%the dimensions of matrixes to fit
n_r = 5000;
n_z = 2500;
r = linspace(r_min,r_max,n_r);
z = linspace(0,z_max,n_z);
[cr,cz] = meshgrid(r,z);
% psi matrix for fit
psi_rz = Psi_s/(a^2*R^2)*(cr.^2.*cz.^2/E^2+(cr.^2-R^2).^2/4);
% the derivatives of psi
[psiDr_rz,psiDz_rz] = gradient(psi_rz,r,z);
psiDz_rz(1,:) = 0;
% the norm of psi gradient
psiGradNorm_rz = sqrt(psiDr_rz.^2+psiDz_rz.^2);
% the derivatives of the square of psi gradient
[psiGrad2Dr_rz,psiGrad2Dz_rz] = gradient(psiGradNorm_rz.^2,r,z);
% the gradient psi direction's derivative of psi gradient square
psiGrad2Dpsi_rz = (psiGrad2Dr_rz.*psiDr_rz+psiGrad2Dz_rz.*psiDz_rz)./psiGradNorm_rz;
% the gradient psi direction's derivative of log psi gradient norm
logpsiGradNormDpsi_rz = psiGrad2Dpsi_rz./(2*psiGradNorm_rz.^2);
% the gradient psi direction's derivative of log r
logrDpsi_rz = 1./cr.*psiGrad2Dr_rz./psiGradNorm_rz;
%% s as an independent variable for flux surface quantities
n_s = floor(n_r/10);
s = linspace(0,1,n_s);
psi = s.^2*Psi_s;
% T
T_s = T0*ones(size(psi));
% derivative of T
TDpsi_s = gradient(T_s,psi);

% pressure
p_s = Psi_s*(1+E^2)*(Psi_s-psi)/(a^2*E^2*R^2);
% derivative of p
pDpsi_s = gradient(p_s,psi);

% toroidal current density
% psi_in = psi_rz;
% psi_in(psi_in>1) = 1;
% jphi_rz = cr.*interp1(psi,pDpsi_s,psi_in)+...
%     1./cr.*interp1(psi,T_s,psi_in).*interp1(psi,TDpsi_s,psi_in);
jphi_rz = -Psi_s*(E^2+1)*cr/(a^2*R^2*E^2);
%% s theta coordinate system

% theta as an independent coordinate
n_theta = n_s;
theta = linspace(0,pi,n_theta);

[cs,ctheta] = meshgrid(psi,theta);

% the corresponding rho matrix
rho_st = zeros(size(ctheta));
for i = 2:n_s% exclude the magnetic axis situation
    psi_i = psi(i);
    
    % the start and end point of the loop integration
    path_left = fzero(@(r) fpsi_rz(r,0)-psi_i,[r_min*0.9,R]);
    path_right = fzero(@(r) fpsi_rz(r,0)-psi_i,[R,r_max*1.1]);
    
    % the initialize vector of rho
    rho_st(1,i) = path_right-R;
    rho_st(end,i) = R-path_left;
    for i_path = 2:1:(n_theta-1)% exclude the start and end point at which z = 0
        rho_st(i_path,i) = fzero(@(rho) ...
            fpsi_rz(R+rho*cos(ctheta(i_path)),...
            rho*sin(ctheta(i_path)))-psi_i,...
            [0,rho_max]);
    end
end


% safety factor
% integrate over the constant psi path
% the corresponding r and z
r_st = R+rho_st.*cos(ctheta);
z_st = rho_st.*sin(ctheta);
r_st(r_st>max(max(cr))) = max(max(cr));
r_st(r_st<min(min(cr))) = min(min(cr));
z_st(z_st>max(max(cz))) = max(max(cz));
z_st(z_st<min(min(cz))) = min(min(cz));
% T
T_st = T0*ones(size(cs));
% the norm of gradient of psi
psiGradNorm_st = interp2(cr,cz,psiGradNorm_rz,r_st,z_st);
% length of path element
dl_st = zeros(size(cs));
dl_st(2:end,:) = sqrt((rho_st(2:end,:)-rho_st(1:end-1,:)).^2+...
    (rho_st(2:end,:).*(ctheta(2:end,:)-ctheta(1:end-1,:))).^2);
% calculate q from loop integration
q_s = 1/pi*sum(T_st./(r_st.*psiGradNorm_st).*dl_st);
q_s(1) = q_s(2);

figure(1);
plot(s,q_s);
% derivative of q
qDpsi_s = gradient(q_s,psi);


% the corresponding chi
q_st = ones(n_theta,1)*q_s;
% calculate integral kernel
kernel = T_st./(q_st.*r_st.*psiGradNorm_st).*dl_st;
% initialize of chi_st
chi_st = zeros(size(ctheta));
% calculate the chi value of this point
for k = 2:(n_theta)
    chi_st(k,:) = sum(kernel(1:k,:));
end
chi_st(:,1) = ctheta(:,1);
% make sure the biggest chi at fixed s big than pi
chi_st(end,chi_st(end,:)<pi) = pi;
%% s chi coordinate system
% chi as a independent variable
n_chi = n_s;
chi = linspace(0,pi,n_chi);

[cs,cchi] = meshgrid(s,chi);

% the corresponding rho theta
rho_sc = zeros(size(cs));
theta_sc = zeros(size(cs));
theta_sc(:,1) = cchi(:,1);
for i = 2:n_s
    theta_sc(:,i) = interp1(chi_st(:,i),ctheta(:,i),cchi(:,i));
    rho_sc(:,i) = interp1(ctheta(:,i),rho_st(:,i),theta_sc(:,i));
end

% the corresponding r z
r_sc = R+rho_sc.*cos(theta_sc);
z_sc = rho_sc.*sin(theta_sc);

figure(2);
plot(r_sc,z_sc);
hold on;
plot(r_sc',z_sc');
hold off;
% the gradient psi at every point
psiDr_sc = interp2(cr,cz,psiDr_rz,r_sc,z_sc);
psiDz_sc = interp2(cr,cz,psiDz_rz,r_sc,z_sc);
% the norm of gradient of psi
psiGradNorm_sc = interp2(cr,cz,psiGradNorm_rz,r_sc,z_sc);

% non-orthogonality
betachi_sc = zeros(size(cs));
dpsi = psi(:,3:end)-psi(:,2:end-1);
dr = dpsi./psiGradNorm_sc(:,2:end-1).^2.*psiDr_sc(:,2:end-1);
dz = dpsi./psiGradNorm_sc(:,2:end-1).^2.*psiDz_sc(:,2:end-1);
rtemp = r_sc(:,2:end-1)+dr;
ztemp = z_sc(:,2:end-1)+dz;
thetatemp = acos((rtemp-R)./sqrt((rtemp-R).^2+ztemp.^2));
chitemp = interp2(cs,ctheta,chi_st,cs(:,3:end),thetatemp);
betachi_sc(:,2:end-1) = (chitemp-cchi(:,2:end-1))./(cs(:,3:end)-cs(:,2:end-1));
betachi_sc(1,:) = 0;
betachi_sc(end,:) = 0;
betachi_sc(:,end) = betachi_sc(:,end-1);

% s chi derivatives of ln(r^2)
logr2Ds_sc = zeros(size(cs));
logr2Dchi_sc = zeros(size(cs));
logr2Ds_sc(:,1:end-1) = (log(r_sc(:,2:end).^2)-log(r_sc(:,1:end-1).^2))...
    ./(cs(:,2:end)-cs(:,1:end-1));
logr2Ds_sc(:,end) = logr2Ds_sc(:,end-1);
logr2Dchi_sc(1:end-1,:) = (log(r_sc(2:end,:).^2)-log(r_sc(1:end-1,:).^2))...
    ./(cchi(2:end,:)-cchi(1:end-1,:));
logr2Dchi_sc(end,:) = logr2Dchi_sc(end-1,:);

% s derivative of T
TDs_s = zeros(size(s));
TDs_s(1:end-1) = (T_s(2:end)-T_s(1:end-1))./(s(2:end)-s(1:end-1));
TDs_s(end) = TDs_s(end-1);

% s derivative of q
qDs_s = zeros(size(s));
qDs_s(1:end-1) = (q_s(2:end)-q_s(1:end-1))./(s(2:end)-s(1:end-1));
qDs_s(end) = qDs_s(end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct of matrixes A and B
%% s and chi for the mesh
n_s = 24;
n_chi = 25;
ms = 0:1/n_s:1;
mchi = -pi/(2*(n_chi-1)):pi/(n_chi-1):pi+pi/(2*(n_chi-1));
%% initialize matrix A
n_x = (3*n_s+1)*(2*n_chi+2);% the totle number of unknow variables
%A = zeros(n_x,n_x);
%% matrix form of the seven dependent components
%initialize
n_cell = n_s*(2*n_chi-2);
XR = spalloc(n_cell,n_x,4*n_cell);
XI = spalloc(n_cell,n_x,4*n_cell);
DXDchiR = spalloc(n_cell,n_x,4*n_cell);
DXDchiI = spalloc(n_cell,n_x,4*n_cell);
DXDsR = spalloc(n_cell,n_x,4*n_cell);
DXDsI = spalloc(n_cell,n_x,4*n_cell);
YR = spalloc(n_cell,n_x,2*n_cell);
YI = spalloc(n_cell,n_x,2*n_cell);
DYDchiR = spalloc(n_cell,n_x,2*n_cell);
DYDchiI = spalloc(n_cell,n_x,2*n_cell);
VR = spalloc(n_cell,n_x,2*n_cell);
VI = spalloc(n_cell,n_x,2*n_cell);
DVDchiR = spalloc(n_cell,n_x,2*n_cell);
DVDchiI = spalloc(n_cell,n_x,2*n_cell);
% the real part of the matrixes
for i = 1:n_s
    for j = 1:n_chi
        % the coefficients of the discretization form
        i_cell = (i-1)*(2*n_chi-2)+j;
        i_xij = (2*n_chi+2)*(3*i-3)+2*(j-1)+1;
        % XR
        XR_coef = [1, 1, 1, 1]/4;
        XR(i_cell,i_xij) = XR_coef(1);
        XR(i_cell,i_xij+2) = XR_coef(2);
        XR(i_cell,i_xij+3*(2*n_chi+2)) = XR_coef(3);
        XR(i_cell,i_xij+3*(2*n_chi+2)+2) = XR_coef(4);
        % XI
        XI_coef = [1, 1, 1, 1]/4;
        XI(i_cell,i_xij+1) = XI_coef(1);
        XI(i_cell,i_xij+3) = XI_coef(2);
        XI(i_cell,i_xij+3*(2*n_chi+2)+1) = XI_coef(3);
        XI(i_cell,i_xij+3*(2*n_chi+2)+3) = XI_coef(4);
        % DXDchiR
        DXDchiR_coef = [-1, 1, -1, 1]/(2*(mchi(j+1)-mchi(j)));
        DXDchiR(i_cell,i_xij) = DXDchiR_coef(1);
        DXDchiR(i_cell,i_xij+2) = DXDchiR_coef(2);
        DXDchiR(i_cell,i_xij+3*(2*n_chi+2)) = DXDchiR_coef(3);
        DXDchiR(i_cell,i_xij+3*(2*n_chi+2)+2) = DXDchiR_coef(4);
        % DXDchiI
        DXDchiI_coef = [-1, 1, -1, 1]/(2*(mchi(j+1)-mchi(j)));
        DXDchiI(i_cell,i_xij+1) = DXDchiI_coef(1);
        DXDchiI(i_cell,i_xij+3) = DXDchiI_coef(2);
        DXDchiI(i_cell,i_xij+3*(2*n_chi+2)+1) = DXDchiI_coef(3);
        DXDchiI(i_cell,i_xij+3*(2*n_chi+2)+3) = DXDchiI_coef(4);
        % DXDsR
        DXDsR_coef = [-1, -1, 1, 1]/(2*(ms(i+1)-ms(i)));
        DXDsR(i_cell,i_xij) = DXDsR_coef(1);
        DXDsR(i_cell,i_xij+2) = DXDsR_coef(2);
        DXDsR(i_cell,i_xij+3*(2*n_chi+2)) = DXDsR_coef(3);
        DXDsR(i_cell,i_xij+3*(2*n_chi+2)+2) = DXDsR_coef(4);
        % DXDsI
        DXDsI_coef = [-1, -1, 1, 1]/(2*(ms(i+1)-ms(i)));
        DXDsI(i_cell,i_xij+1) = DXDsI_coef(1);
        DXDsI(i_cell,i_xij+3) = DXDsI_coef(2);
        DXDsI(i_cell,i_xij+3*(2*n_chi+2)+1) = DXDsI_coef(3);
        DXDsI(i_cell,i_xij+3*(2*n_chi+2)+3) = DXDsI_coef(4);
        
        i_xij = i_xij + (2*n_chi+2);
        % YR
        YR_coef = [1, 1]/2;
        YR(i_cell,i_xij) = YR_coef(1);
        YR(i_cell,i_xij+2) = YR_coef(2);
        % YI
        YI_coef = [1, 1]/2;
        YI(i_cell,i_xij+1) = YI_coef(1);
        YI(i_cell,i_xij+3) = YI_coef(2);
        % DYDchiR
        DYDchiR_coef = [-1, 1]/(mchi(j+1)-mchi(j));
        DYDchiR(i_cell,i_xij) = DYDchiR_coef(1);
        DYDchiR(i_cell,i_xij+2) = DYDchiR_coef(2);
        % DYDchiI
        DYDchiI_coef = [-1, 1]/(mchi(j+1)-mchi(j));
        DYDchiI(i_cell,i_xij+1) = DYDchiI_coef(1);
        DYDchiI(i_cell,i_xij+3) = DYDchiI_coef(2);
        
        i_xij = i_xij + (2*n_chi+2);
        % VR
        VR_coef = [1, 1]/2;
        VR(i_cell,i_xij) = VR_coef(1);
        VR(i_cell,i_xij+2) = VR_coef(2);
        % VI
        VI_coef = [1, 1]/2;
        VI(i_cell,i_xij+1) = VI_coef(1);
        VI(i_cell,i_xij+3) = VI_coef(2);
        % DVDchiR
        DVDchiR_coef = [-1, 1]/(mchi(j+1)-mchi(j));
        DVDchiR(i_cell,i_xij) = DVDchiR_coef(1);
        DVDchiR(i_cell,i_xij+2) = DVDchiR_coef(2);
        % DVDchiI
        DVDchiI_coef = [-1, 1]/(mchi(j+1)-mchi(j));
        DVDchiI(i_cell,i_xij+1) = DVDchiI_coef(1);
        DVDchiR(i_cell,i_xij+3) = DVDchiI_coef(2);
        
        if j~=1&&j~=n_chi
            i_cell = (i-1)*(2*n_chi-2)+2*n_chi-j;
            i_xij = (2*n_chi+2)*(3*i-3)+2*(j-1)+1;
            % the symmetric part of XR
            XRS_coef = -[1, 1, 1, 1]/4;
            XR(i_cell,i_xij) = XRS_coef(1);
            XR(i_cell,i_xij+2) = XRS_coef(2);
            XR(i_cell,i_xij+3*(2*n_chi+2)) = XRS_coef(3);
            XR(i_cell,i_xij+3*(2*n_chi+2)+2) = XRS_coef(4);
            % the symmetric part of XI
            XIS_coef = [1, 1, 1, 1]/4;
            XI(i_cell,i_xij+1) = XIS_coef(1);
            XI(i_cell,i_xij+3) = XIS_coef(2);
            XI(i_cell,i_xij+3*(2*n_chi+2)+1) = XIS_coef(3);
            XI(i_cell,i_xij+3*(2*n_chi+2)+3) = XIS_coef(4);
            % the symmetric part of DXDchiR
            DXDchiRS_coef = -[-1, 1, -1, 1]/(2*(mchi(j)-mchi(j+1)));
            DXDchiR(i_cell,i_xij) = DXDchiRS_coef(1);
            DXDchiR(i_cell,i_xij+2) = DXDchiRS_coef(2);
            DXDchiR(i_cell,i_xij+3*(2*n_chi+2)) = DXDchiRS_coef(3);
            DXDchiR(i_cell,i_xij+3*(2*n_chi+2)+2) = DXDchiRS_coef(4);
            % DXDchiI
            DXDchiIS_coef = [-1, 1, -1, 1]/(2*(mchi(j)-mchi(j+1)));
            DXDchiI(i_cell,i_xij+1) = DXDchiIS_coef(1);
            DXDchiI(i_cell,i_xij+3) = DXDchiIS_coef(2);
            DXDchiI(i_cell,i_xij+3*(2*n_chi+2)+1) = DXDchiIS_coef(3);
            DXDchiI(i_cell,i_xij+3*(2*n_chi+2)+3) = DXDchiIS_coef(4);
            % DXDsR
            DXDsRS_coef = -[-1, -1, 1, 1]/(2*(ms(i+1)-ms(i)));
            DXDsR(i_cell,i_xij) = DXDsRS_coef(1);
            DXDsR(i_cell,i_xij+2) = DXDsRS_coef(2);
            DXDsR(i_cell,i_xij+3*(2*n_chi+2)) = DXDsRS_coef(3);
            DXDsR(i_cell,i_xij+3*(2*n_chi+2)+2) = DXDsRS_coef(4);
            % DXDsI
            DXDsIS_coef = [-1, -1, 1, 1]/(2*(ms(i+1)-ms(i)));
            DXDsI(i_cell,i_xij+1) = DXDsIS_coef(1);
            DXDsI(i_cell,i_xij+3) = DXDsIS_coef(2);
            DXDsI(i_cell,i_xij+3*(2*n_chi+2)+1) = DXDsIS_coef(3);
            DXDsI(i_cell,i_xij+3*(2*n_chi+2)+3) = DXDsIS_coef(4);
            
            i_xij = i_xij + (2*n_chi+2);
            % YR
            YRS_coef = [1, 1]/2;
            YR(i_cell,i_xij) = YRS_coef(1);
            YR(i_cell,i_xij+2) = YRS_coef(2);
            % YI
            YIS_coef = -[1, 1]/2;
            YI(i_cell,i_xij+1) = YIS_coef(1);
            YI(i_cell,i_xij+3) = YIS_coef(2);
            % DYDchiR
            DYDchiRS_coef = [-1, 1]/(mchi(j)-mchi(j+1));
            DYDchiR(i_cell,i_xij) = DYDchiRS_coef(1);
            DYDchiR(i_cell,i_xij+2) = DYDchiRS_coef(2);
            % DYDchiI
            DYDchiIS_coef = -[-1, 1]/(mchi(j)-mchi(j+1));
            DYDchiI(i_cell,i_xij+1) = DYDchiIS_coef(1);
            DYDchiI(i_cell,i_xij+3) = DYDchiIS_coef(2);
            
            i_xij = i_xij + (2*n_chi+2);
            % VR
            VRS_coef = [1, 1]/2;
            VR(i_cell,i_xij) = VRS_coef(1);
            VR(i_cell,i_xij+2) = VRS_coef(2);
            % VI
            VIS_coef = -[1, 1]/2;
            VI(i_cell,i_xij+1) = VIS_coef(1);
            VI(i_cell,i_xij+3) = VIS_coef(2);
            % DVDchiR
            DVDchiRS_coef = [-1, 1]/(mchi(j)-mchi(j+1));
            DVDchiR(i_cell,i_xij) = DVDchiRS_coef(1);
            DVDchiR(i_cell,i_xij+2) = DVDchiRS_coef(2);
            % DVDchiI
            DVDchiIS_coef = -[-1, 1]/(mchi(j)-mchi(j+1));
            DVDchiI(i_cell,i_xij+1) = DVDchiIS_coef(1);
            DVDchiR(i_cell,i_xij+3) = DVDchiIS_coef(2);
        end
    end
end

%% the normalized equilibrium quantities in the middle of each cell
% s and chi the physical components
ps = 1/(2*n_s):1/n_s:1-1/(2*n_s);
pchi = [0:pi/(n_chi-1):pi,pi-pi/(n_chi-1):-pi/(n_chi-1):pi/(n_chi-1)];
[ps,pchi] = meshgrid(ps,pchi);
% psi
ppsi = ps.^2*Psi_s;
% r z
pr = interp2(cs,cchi,r_sc,ps,pchi);
pz = interp2(cs,cchi,z_sc,ps,pchi);
% the norm of gradient of psi
ppsiGradNorm = interp2(cr,cz,psiGradNorm_rz,pr,pz);
% mass density
prho = 1;
% plasma pressure
pp = interp1(s,p_s,ps);
% psi derivative of p
ppDpsi = interp1(s,pDpsi_s,ps);
% flux function
pT = interp1(s,T_s,ps);
% s derivative of T
pTDs = interp1(s,TDs_s,ps);
% toroidal current density
pjphi = interp2(cr,cz,jphi_rz,pr,pz);
% safety factor
pq = interp1(s,q_s,ps);
% s derivative of q
pqDs = interp1(s,qDs_s,ps);
% poloidal magnetic field
pBp = ppsi.*pr.^2./(q0*ppsiGradNorm.^2);
% non-orthogonality
pbetachi = [interp2(cs,cchi,betachi_sc,ps(1:n_chi,:),pchi(1:n_chi,:));...
    -interp2(cs,cchi,betachi_sc,ps(n_chi+1:end,:),pchi(n_chi+1:end,:))];
% dlog(r^2)/ds and dlog(r^2)/dchi
plogr2Ds = interp2(cs,cchi,logr2Ds_sc,ps,pchi);
plogr2Dchi = [interp2(cs,cchi,logr2Dchi_sc,ps(1:n_chi,:),pchi(1:n_chi,:));...
    -interp2(cs,cchi,logr2Dchi_sc,ps(n_chi+1:end,:),pchi(n_chi+1:end,:))];
% H defined at the ERATO paper
pH = 2*pjphi.*ppsi.*pr./(ps.*ppsiGradNorm.^2)+pTDs./pT-pqDs./pq;

% the gradient psi direction's derivative of ln psi gradient norm
plogpsiGradNormDpsi = interp2(cr,cz,logpsiGradNormDpsi_rz,pr,pz);
% the gradient psi direction's derivative of log r
plogrDpsi = interp2(cr,cz,logrDpsi_rz,pr,pz);
% K defined at the ERATO paper
pK = 2*ppsi/q0.*(pjphi.^2./ppsiGradNorm.^2+pjphi.^2./pr.*plogpsiGradNormDpsi...
    -ppDpsi.*plogrDpsi);


%% the matrix form of Is
I1R = 1./pq(:).*DXDchiR-n*XI;
I1I = 1./pq(:).*DXDchiI+n*XR;
I2R = DXDsR+DVDchiR;
I2I = DXDsI+DVDchiI;
I3R = pH(:).*XR-pbetachi(:).*n.*pq(:).*XI...
    -n*pq(:).*VI+pbetachi(:).*DXDchiR+DXDsR;
I3I = pH(:).*XI+pbetachi(:).*n.*pq(:).*XR...
    +n*pq(:).*VR+pbetachi(:).*DXDchiI+DXDsI;
I4R = plogr2Ds(:).*XR+plogr2Dchi(:).*(VR+pq(:).*YR)...
    -n*pq(:).^2.*YI+DXDsR+DVDchiR+pq(:).*DYDchiR;
I4I = plogr2Ds(:).*XI+plogr2Dchi(:).*(VI+pq(:).*YI)...
    +n*pq(:).^2.*YR+DXDsI+DVDchiI+pq(:).*DYDchiI;
I5R = XR;
I5I = XI;
%% the coefficient a b c d e f g h
J = pq.*pr.^2./pT;
a = 2*pq.^2.*ppsi.*pr.^4./(J.^3.*ppsiGradNorm.^2);
b = pT.^2.*pr.^2./(2*Psi_s*J);
c = ppsiGradNorm.^2.*pr.^2./(2*Psi_s*J);
d = pr.^4*gamma.*pp./(2*Psi_s*J);
e = pK.*2.*pr.^4*q0./J;
f = 2*prho.*ppsi.*pT.*pr.^2./(pq.*ppsiGradNorm.^2);
g = prho.*ppsiGradNorm.^2.*pq.*pr.^4./(2*pT*Psi_s);
h = prho.*pr.^4.*pT.*pq/(2*Psi_s);



%% the matrix A
A =   I1R'*(a(:)./ps(:).*I1R) + I1I'*(a(:)./ps(:).*I1I)...
    + I2R'*(b(:)./ps(:).*I2R) + I2I'*(b(:)./ps(:).*I2I)...
    + I3R'*(c(:)./ps(:).*I3R) + I3I'*(c(:)./ps(:).*I3I)...
    + I4R'*(d(:)./ps(:).*I4R) + I4I'*(d(:)./ps(:).*I4I)...
    - I5R'*(e(:)./ps(:).*I5R) - I5I'*(e(:)./ps(:).*I5I);
%% the matrix B
B =  XR'*(f(:)./ps(:).*XR) + XI'*(f(:)./ps(:).*XI) +...
    (VR+YR-pbetachi(:).*XR)'*(g(:)./ps(:).*(VR+YR-pbetachi(:).*XR))+...
    (VI+YI-pbetachi(:).*XI)'*(g(:)./ps(:).*(VI+YI-pbetachi(:).*XI))+...
    YR'*(h(:)./ps(:).*YR) + YI'*(h(:)./ps(:).*YI);
% figure(1);
% subplot(2,4,1);spy(A);
% subplot(2,4,5);spy(B);
% A = full(A);B = full(B);
% disp([rank(A),rank(B)]);
% A = sparse(A);
% B = sparse(B);

%% symmetry conditions
U = diag(ones(n_x,1));
for i = 1:n_s+1
    % symmetry conditions of XR
    U((2*n_chi+2)*(3*i-3)+1,(2*n_chi+2)*(3*i-3)+1) = 1;
    U((2*n_chi+2)*(3*i-3)+1,(2*n_chi+2)*(3*i-3)+3) = 1;
    U((2*n_chi+2)*(3*i-3)+2*n_chi+1,(2*n_chi+2)*(3*i-3)+2*n_chi-1) = 1;
    U((2*n_chi+2)*(3*i-3)+2*n_chi+1,(2*n_chi+2)*(3*i-3)+2*n_chi+1) = 1;
    % symmetry conditions of XI
    U((2*n_chi+2)*(3*i-3)+2,(2*n_chi+2)*(3*i-3)+2) = -1;
    U((2*n_chi+2)*(3*i-3)+2,(2*n_chi+2)*(3*i-3)+4) = 1;
    U((2*n_chi+2)*(3*i-3)+2*n_chi+2,(2*n_chi+2)*(3*i-3)+2*n_chi) = -1;
    U((2*n_chi+2)*(3*i-3)+2*n_chi+2,(2*n_chi+2)*(3*i-3)+2*n_chi+2) = 1;
end
for i = 1:n_s
    % symmetry conditions of VR
    U((2*n_chi+2)*(3*i-2)+1,(2*n_chi+2)*(3*i-2)+1) = -1;
    U((2*n_chi+2)*(3*i-2)+1,(2*n_chi+2)*(3*i-2)+3) = 1;
    U((2*n_chi+2)*(3*i-2)+2*n_chi+1,(2*n_chi+2)*(3*i-2)+2*n_chi-1) = -1;
    U((2*n_chi+2)*(3*i-2)+2*n_chi+1,(2*n_chi+2)*(3*i-2)+2*n_chi+1) = 1;
    % symmetry conditions of VI
    U((2*n_chi+2)*(3*i-2)+2,(2*n_chi+2)*(3*i-2)+2) = 1;
    U((2*n_chi+2)*(3*i-2)+2,(2*n_chi+2)*(3*i-2)+4) = 1;
    U((2*n_chi+2)*(3*i-2)+2*n_chi+2,(2*n_chi+2)*(3*i-2)+2*n_chi) = 1;
    U((2*n_chi+2)*(3*i-2)+2*n_chi+2,(2*n_chi+2)*(3*i-2)+2*n_chi+2) = 1;
    % symmetry conditions of YR
    U((2*n_chi+2)*(3*i-1)+1,(2*n_chi+2)*(3*i-1)+1) = -1;
    U((2*n_chi+2)*(3*i-1)+1,(2*n_chi+2)*(3*i-1)+3) = 1;
    U((2*n_chi+2)*(3*i-1)+2*n_chi+1,(2*n_chi+2)*(3*i-1)+2*n_chi-1) = -1;
    U((2*n_chi+2)*(3*i-1)+2*n_chi+1,(2*n_chi+2)*(3*i-1)+2*n_chi+1) = 1;
    % symmetry conditions of YI
    U((2*n_chi+2)*(3*i-1)+2,(2*n_chi+2)*(3*i-1)+2) = 1;
    U((2*n_chi+2)*(3*i-1)+2,(2*n_chi+2)*(3*i-1)+4) = 1;
    U((2*n_chi+2)*(3*i-1)+2*n_chi+2,(2*n_chi+2)*(3*i-1)+2*n_chi) = 1;
    U((2*n_chi+2)*(3*i-1)+2*n_chi+2,(2*n_chi+2)*(3*i-1)+2*n_chi+2) = 1;
end
% transformation using U
Uinv = inv(U);
At = Uinv'*A*Uinv;
Bt = Uinv'*B*Uinv;

% subplot(2,4,2);spy(At);
% subplot(2,4,6);spy(Bt);
% disp([rank(At),rank(Bt)]);
% force the new tilda variables to zero
w2inv = 1e-20; % which make the condition fullfilled when w2~=10^20
for i = 1:n_s+1
    %tilda A
    % symmetry conditions of XR
    At((2*n_chi+2)*(3*i-3)+1,:) = 0;
    At(:,(2*n_chi+2)*(3*i-3)+1) = 0;
    At((2*n_chi+2)*(3*i-3)+1,(2*n_chi+2)*(3*i-3)+1) = 1;
    At((2*n_chi+2)*(3*i-3)+2*n_chi+1,:) = 0;
    At(:,(2*n_chi+2)*(3*i-3)+2*n_chi+1) = 0;
    At((2*n_chi+2)*(3*i-3)+2*n_chi+1,(2*n_chi+2)*(3*i-3)+2*n_chi+1) = 1;
    
    At((2*n_chi+2)*(3*i-3)+3,:) = 0;
    At(:,(2*n_chi+2)*(3*i-3)+3) = 0;
    At((2*n_chi+2)*(3*i-3)+3,(2*n_chi+2)*(3*i-3)+3) = 1;
   
    % symmetry conditions of XI
    At((2*n_chi+2)*(3*i-3)+2,:) = 0;
    At(:,(2*n_chi+2)*(3*i-3)+2) = 0;
    At((2*n_chi+2)*(3*i-3)+2,(2*n_chi+2)*(3*i-3)+2) = 1;
    At((2*n_chi+2)*(3*i-3)+2*n_chi+2,:) = 0;
    At(:,(2*n_chi+2)*(3*i-3)+2*n_chi+2) = 0;
    At((2*n_chi+2)*(3*i-3)+2*n_chi+2,(2*n_chi+2)*(3*i-3)+2*n_chi+2) = 1;
    
    % tilda B
    % symmetry conditions of XR
    Bt((2*n_chi+2)*(3*i-3)+1,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-3)+1) = 0;
    Bt((2*n_chi+2)*(3*i-3)+1,(2*n_chi+2)*(3*i-3)+1) = w2inv;
    Bt((2*n_chi+2)*(3*i-3)+2*n_chi+1,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-3)+2*n_chi+1) = 0;
    Bt((2*n_chi+2)*(3*i-3)+2*n_chi+1,(2*n_chi+2)*(3*i-3)+2*n_chi+1) = w2inv;
    
    Bt((2*n_chi+2)*(3*i-3)+3,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-3)+3) = 0;
    Bt((2*n_chi+2)*(3*i-3)+3,(2*n_chi+2)*(3*i-3)+3) = w2inv;
    % symmetry conditions of XI
    Bt((2*n_chi+2)*(3*i-3)+2,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-3)+2) = 0;
    Bt((2*n_chi+2)*(3*i-3)+2,(2*n_chi+2)*(3*i-3)+2) = w2inv;
    Bt((2*n_chi+2)*(3*i-3)+2*n_chi+2,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-3)+2*n_chi+2) = 0;
    Bt((2*n_chi+2)*(3*i-3)+2*n_chi+2,(2*n_chi+2)*(3*i-3)+2*n_chi+2) = w2inv;
end


for i = 1:n_s
    % symmetry conditions of VR
    At((2*n_chi+2)*(3*i-2)+1,:) = 0;
    At(:,(2*n_chi+2)*(3*i-2)+1) = 0;
    At((2*n_chi+2)*(3*i-2)+1,(2*n_chi+2)*(3*i-2)+1) = 1;
    At((2*n_chi+2)*(3*i-2)+2*n_chi+1,:) = 0;
    At(:,(2*n_chi+2)*(3*i-2)+2*n_chi+1) = 0;
    At((2*n_chi+2)*(3*i-2)+2*n_chi+1,(2*n_chi+2)*(3*i-2)+2*n_chi+1) = 1;
    % symmetry conditions of VI
    At((2*n_chi+2)*(3*i-2)+2,:) = 0;
    At(:,(2*n_chi+2)*(3*i-2)+2) = 0;
    At((2*n_chi+2)*(3*i-2)+2,(2*n_chi+2)*(3*i-2)+2) = 1;
    At((2*n_chi+2)*(3*i-2)+2*n_chi+2,:) = 0;
    At(:,(2*n_chi+2)*(3*i-2)+2*n_chi+2) = 0;
    At((2*n_chi+2)*(3*i-2)+2*n_chi+2,(2*n_chi+2)*(3*i-2)+2*n_chi+2) = 1;
    
    At((2*n_chi+2)*(3*i-2)+4,:) = 0;
    At(:,(2*n_chi+2)*(3*i-2)+4) = 0;
    At((2*n_chi+2)*(3*i-2)+4,(2*n_chi+2)*(3*i-2)+4) = 1;
    % symmetry conditions of YR
    At((2*n_chi+2)*(3*i-1)+1,:) = 0;
    At(:,(2*n_chi+2)*(3*i-1)+1) = 0;
    At((2*n_chi+2)*(3*i-1)+1,(2*n_chi+2)*(3*i-1)+1) = 1;
    At((2*n_chi+2)*(3*i-1)+2*n_chi+1,:) = 0;
    At(:,(2*n_chi+2)*(3*i-1)+2*n_chi+1) = 0;
    At((2*n_chi+2)*(3*i-1)+2*n_chi+1,(2*n_chi+2)*(3*i-1)+2*n_chi+1) = 1;
    % symmetry conditions of YI
    At((2*n_chi+2)*(3*i-1)+2,:) = 0;
    At(:,(2*n_chi+2)*(3*i-1)+2) = 0;
    At((2*n_chi+2)*(3*i-1)+2,(2*n_chi+2)*(3*i-1)+2) = 1;
    At((2*n_chi+2)*(3*i-1)+2*n_chi+2,:) = 0;
    At(:,(2*n_chi+2)*(3*i-1)+2*n_chi+2) = 0;
    At((2*n_chi+2)*(3*i-1)+2*n_chi+2,(2*n_chi+2)*(3*i-1)+2*n_chi+2) = 1;
    
    At((2*n_chi+2)*(3*i-1)+4,:) = 0;
    At(:,(2*n_chi+2)*(3*i-1)+4) = 0;
    At((2*n_chi+2)*(3*i-1)+4,(2*n_chi+2)*(3*i-1)+4) = 1;
    % tilda B
    % symmetry conditions of VR
    Bt((2*n_chi+2)*(3*i-2)+1,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-2)+1) = 0;
    Bt((2*n_chi+2)*(3*i-2)+1,(2*n_chi+2)*(3*i-2)+1) = w2inv;
    Bt((2*n_chi+2)*(3*i-2)+2*n_chi+1,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-2)+2*n_chi+1) = 0;
    Bt((2*n_chi+2)*(3*i-2)+2*n_chi+1,(2*n_chi+2)*(3*i-2)+2*n_chi+1) = w2inv;
    % symmetry conditions of VI
    Bt((2*n_chi+2)*(3*i-2)+2,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-2)+2) = 0;
    Bt((2*n_chi+2)*(3*i-2)+2,(2*n_chi+2)*(3*i-2)+2) = w2inv;
    Bt((2*n_chi+2)*(3*i-2)+2*n_chi+2,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-2)+2*n_chi+2) = 0;
    Bt((2*n_chi+2)*(3*i-2)+2*n_chi+2,(2*n_chi+2)*(3*i-2)+2*n_chi+2) = w2inv;
    
    Bt((2*n_chi+2)*(3*i-2)+4,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-2)+4) = 0;
    Bt((2*n_chi+2)*(3*i-2)+4,(2*n_chi+2)*(3*i-2)+4) = w2inv;
    % symmetry conditions of YR
    Bt((2*n_chi+2)*(3*i-1)+1,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-1)+1) = 0;
    Bt((2*n_chi+2)*(3*i-1)+1,(2*n_chi+2)*(3*i-1)+1) = w2inv;
    Bt((2*n_chi+2)*(3*i-1)+2*n_chi+1,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-1)+2*n_chi+1) = 0;
    Bt((2*n_chi+2)*(3*i-1)+2*n_chi+1,(2*n_chi+2)*(3*i-1)+2*n_chi+1) = w2inv;
    % symmetry conditions of YI
    Bt((2*n_chi+2)*(3*i-1)+2,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-1)+2) = 0;
    Bt((2*n_chi+2)*(3*i-1)+2,(2*n_chi+2)*(3*i-1)+2) = w2inv;
    Bt((2*n_chi+2)*(3*i-1)+2*n_chi+2,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-1)+2*n_chi+2) = 0;
    Bt((2*n_chi+2)*(3*i-1)+2*n_chi+2,(2*n_chi+2)*(3*i-1)+2*n_chi+2) = w2inv;
    
    Bt((2*n_chi+2)*(3*i-1)+4,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-1)+4) = 0;
    Bt((2*n_chi+2)*(3*i-1)+4,(2*n_chi+2)*(3*i-1)+4) = w2inv;
end

% subplot(2,4,3);spy(At);
% subplot(2,4,7);spy(Bt);
% disp([rank(At),rank(Bt)]);
%% regularity conditions at s = 0
% the procejure is same as symmetry conditions
for j = 1:n_chi+1
    % tilda A
    % XR
    At(2*(j-1)+1,:) = 0;
    At(:,2*(j-1)+1) = 0;
    At(2*(j-1)+1,2*(j-1)+1) = 1;
    % XI
    At(2*(j-1)+2,:) = 0;
    At(:,2*(j-1)+2) = 0;
    At(2*(j-1)+2,2*(j-1)+2) = 1;
    % tilda B
    % XR
    Bt(2*(j-1)+1,:) = 0;
    Bt(:,2*(j-1)+1) = 0;
    Bt(2*(j-1)+1,2*(j-1)+1) = w2inv;
    % XI
    Bt(2*(j-1)+2,:) = 0;
    Bt(:,2*(j-1)+2) = 0;
    Bt(2*(j-1)+2,2*(j-1)+2) = w2inv;
end

%% boundary conditions at plasma surface s = 1
% a wall straight on the plasma surface
%the procejure is same as symmetry conditions
for j = 1:n_chi+1
    % tilda A
    % XR
    At(3*n_s*(2*n_chi+2)+2*(j-1)+1,:) = 0;
    At(:,3*n_s*(2*n_chi+2)+2*(j-1)+1) = 0;
    At(3*n_s*(2*n_chi+2)+2*(j-1)+1,3*n_s*(2*n_chi+2)+2*(j-1)+1) = 1;
    % XI
    At(3*n_s*(2*n_chi+2)+2*(j-1)+2,:) = 0;
    At(:,3*n_s*(2*n_chi+2)+2*(j-1)+2) = 0;
    At(3*n_s*(2*n_chi+2)+2*(j-1)+2,3*n_s*(2*n_chi+2)+2*(j-1)+2) = 1;
    % tilda B
    % XR
    Bt(3*n_s*(2*n_chi+2)+2*(j-1)+1,:) = 0;
    Bt(:,3*n_s*(2*n_chi+2)+2*(j-1)+1) = 0;
    Bt(3*n_s*(2*n_chi+2)+2*(j-1)+1,3*n_s*(2*n_chi+2)+2*(j-1)+1) = w2inv;
    % XI
    Bt(3*n_s*(2*n_chi+2)+2*(j-1)+2,:) = 0;
    Bt(:,3*n_s*(2*n_chi+2)+2*(j-1)+2) = 0;
    Bt(3*n_s*(2*n_chi+2)+2*(j-1)+2,3*n_s*(2*n_chi+2)+2*(j-1)+2) = w2inv;
end
% subplot(2,4,4);spy(At);
% subplot(2,4,8);spy(Bt);
% disp([rank(At),rank(Bt)]);
% At = sparse(At);
% Bt = sparse(Bt);
w2 = eig(Bt\At);
disp(min(w2));
GrowthRate = w2;


% figure(1);
% plot(r_sc,z_sc);
% hold on
% plot(r_sc',z_sc')
% figure(2);
% mesh(cs,cchi,betachi_sc);
