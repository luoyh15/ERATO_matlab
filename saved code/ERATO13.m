function GrowthRate = ERATO13(E,a,q0,n,n_r)
% the calculation of non-orthogonality
%% the parameters determine the equilibrium
R = 1;  B = 1; T0 = B*R; gamma = 5/3;
% parameters determine the equilibrium
% E = 2; a = 1/3; q0 = 0.8; n = 1;
% E = Es; a = as; q0 = q0s; n = ns;

Psi_s = E*a^2*B/(q0*2); %the psi value on the plasma surface
syms r z ;
fpsi_rz = Psi_s/(a^2*R^2)*(z.^2.*r.^2/E^2+(r.^2-R^2)^2/4);
fpsiDr_rz =  diff(fpsi_rz,r);
fpsiDz_rz =  diff(fpsi_rz,z);
fpsiGradNorm_rz = sqrt(fpsiDr_rz.^2+fpsiDz_rz.^2);
fpsiGrad2Dr_rz = diff(fpsiGradNorm_rz.^2,r);
fpsiGrad2Dz_rz = diff(fpsiGradNorm_rz.^2,z);
fpsiGrad2Dpsi_rz = (fpsiGrad2Dr_rz.*fpsiDr_rz+fpsiGrad2Dz_rz.*fpsiDz_rz)...
    ./fpsiGradNorm_rz.^2;
flogpsiGradNormDpsi_rz = fpsiGrad2Dpsi_rz./(2*fpsiGradNorm_rz.^2);
flogrDpsi_rz = 1./r.*fpsiDr_rz./fpsiGradNorm_rz.^2;

fp_psi = @(psi) 2*Psi_s*(1+E^2)*(Psi_s-psi)/(a^2*E^2*R^4);
fpDpsi_psi = @(psi) -2*Psi_s*(1+E^2)/(a^2*E^2*R^4)*ones(size(psi));

% fp_psi = matlabFunction(fp_psi);
% fpDpsi_psi = matlabFunction(fpDpsi_psi);

fpsi_rz = matlabFunction(fpsi_rz);
fpsiDr_rz = matlabFunction(fpsiDr_rz);
fpsiDz_rz = matlabFunction(fpsiDz_rz);
fpsiGradNorm_rz = matlabFunction(fpsiGradNorm_rz);
fpsiGrad2Dr_rz = matlabFunction(fpsiGrad2Dr_rz);
fpsiGrad2Dz_rz = matlabFunction(fpsiGrad2Dz_rz);
fpsiGrad2Dpsi_rz = matlabFunction(fpsiGrad2Dpsi_rz);
flogpsiGradNormDpsi_rz = matlabFunction(flogpsiGradNormDpsi_rz);
flogrDpsi_rz = matlabFunction(flogrDpsi_rz);

% the boundary of plasma
r_min = sqrt(R^2-2*a*R); r_max = sqrt(R^2+2*a*R);
r_temp = R*(1-4*a^2/R^2)^(1/4);
z_max = sqrt((4*R^2*a^2-(r_temp^2-R^2)^2)*E^2/(4*r_temp^2));
rho_max = sqrt(max(R-r_min,r_max-R)^2+z_max^2);
%% r z coordinate system
%the dimensions of matrixes to fit
% n_r = 200;
n_z = floor(n_r/2);
r = linspace(r_min,r_max,n_r);
z = linspace(0,z_max,n_z);
[cr,cz] = meshgrid(r,z);% psi matrix for fit
psi_rz = fpsi_rz(cr,cz);
% the derivatives of psi
psiDr_rz = fpsiDr_rz(cr,cz);
psiDz_rz = fpsiDz_rz(cr,cz);
% the norm of psi gradient
psiGradNorm_rz = fpsiGradNorm_rz(cr,cz);
% the derivatives of the square of psi gradient
psiGrad2Dr_rz = fpsiGrad2Dr_rz(cr,cz);
psiGrad2Dz_rz = fpsiGrad2Dz_rz(cr,cz);
% the gradient psi direction's derivative of psi gradient square
psiGrad2Dpsi_rz = fpsiGrad2Dpsi_rz(cr,cz);
% the gradient psi direction's derivative of log psi gradient norm
lnpsiGradNormDpsi_rz = flogpsiGradNormDpsi_rz(cr,cz);
% the gradient psi direction's derivative of log r
lnrDpsi_rz = flogrDpsi_rz(cr,cz);
% toroidal current density
jphi_rz = -cr.*fpDpsi_psi(psi_rz);
%% s as an independent variable for flux surface quantities
n_s = floor(n_r/6);
s = linspace(0,1,n_s);
psi = s.^2*Psi_s;
% T
T_s = T0*ones(size(psi));
% derivative of T
TDpsi_s = gradient(T_s,psi);

% pressure
p_s = fp_psi(psi);
% derivative of p
pDpsi_s = fpDpsi_psi(psi);

%% s theta coordinate system
interpmethod = 'linear';
% theta as an independent coordinate
n_theta = n_s*3;
theta = linspace(0,pi,n_theta);

[cs,ctheta] = meshgrid(s,theta);

% the corresponding rho matrix
rho_st = zeros(size(ctheta));
for i = 2:n_s% exclude the magnetic axis situation
    psi_i = psi(i);
    
    % the start and end point of the loop integration
    path_left = fzero(@(r) fpsi_rz(r,0)-psi_i,[r_min*0.9,R]);
    path_right = fzero(@(r) fpsi_rz(r,0)-psi_i,[R,r_max*1.1]);
    
    % the initialize vector of rho
    rho_st(1,i) = R-path_left;
    rho_st(end,i) = path_right-R;
    for i_path = 2:1:(n_theta-1)% exclude the start and end point at which z = 0
        rho_st(i_path,i) = fzero(@(rho) ...
            fpsi_rz(R-rho*cos(ctheta(i_path)),...
            rho*sin(ctheta(i_path)))-psi_i,...
            [0,rho_max]);
    end
end


% safety factor
% integrate over the constant psi path
% the corresponding r and z
r_st = R-rho_st.*cos(ctheta);
z_st = rho_st.*sin(ctheta);
r_st(r_st>max(max(cr))) = max(max(cr));
r_st(r_st<min(min(cr))) = min(min(cr));
z_st(z_st>max(max(cz))) = max(max(cz));
z_st(z_st<min(min(cz))) = min(min(cz));
% T
T_st = T0*ones(size(cs));
% the norm of gradient of psi
psiGradNorm_st = interp2(cr,cz,psiGradNorm_rz,r_st,z_st,interpmethod);
% length of path element
dl_st = zeros(size(cs));
dl_st(2:end,:) = sqrt((rho_st(2:end,:)-rho_st(1:end-1,:)).^2+...
    ((rho_st(2:end,:)+rho_st(1:end-1,:))/2.*...
    (ctheta(2:end,:)-ctheta(1:end-1,:))).^2);
% calculate q from loop integration
q_s = 1/pi*sum(T_st./(r_st.*psiGradNorm_st).*dl_st);
q_s(1) = q_s(2);
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
% make sure the biggest chi at fixed s not smaller than pi
chi_st(end,chi_st(end,:)<pi) = pi;
%% s chi coordinate system
% chi as a independent variable
n_chi = floor(n_theta/3);
chi = linspace(0,pi,n_chi);

[cs,cchi] = meshgrid(s,chi);

% the corresponding rho theta
rho_sc = zeros(size(cs));
theta_sc = zeros(size(cs));
theta_sc(:,1) = cchi(:,1);
for i = 2:n_s
    theta_sc(:,i) = interp1(chi_st(:,i),ctheta(:,i),cchi(:,i),interpmethod);
    rho_sc(:,i) = interp1(ctheta(:,i),rho_st(:,i),theta_sc(:,i),interpmethod);
end

% the corresponding r z
r_sc = R-rho_sc.*cos(theta_sc);
z_sc = rho_sc.*sin(theta_sc);

% the gradient psi at every point
psiDr_sc = interp2(cr,cz,psiDr_rz,r_sc,z_sc,interpmethod);
psiDz_sc = interp2(cr,cz,psiDz_rz,r_sc,z_sc,interpmethod);
% the norm of gradient of psi
psiGradNorm_sc = interp2(cr,cz,psiGradNorm_rz,r_sc,z_sc,interpmethod);
% jphi
jphi_sc = interp2(cr,cz,jphi_rz,r_sc,z_sc,interpmethod);
% T and dTdpsi
T_sc = T0*ones(size(cs));
TDpsi_sc = 0;
% q and dqdpsi
q_sc = ones(n_chi,1)*q_s;
qDpsi_sc = ones(n_chi,1)*qDpsi_s;
logpsiGradNormDpsi_sc = interp2(cr,cz,lnpsiGradNormDpsi_rz,r_sc,z_sc,interpmethod);
psiGrad2Dpsi_sc = interp2(cr,cz,psiGrad2Dpsi_rz,r_sc,z_sc,interpmethod);

psiGradNormDs_sc = zeros(size(cs));
psiGradNormDs_sc(:,1:end-1) = (psiGradNorm_sc(:,2:end)-psiGradNorm_sc(:,1:end-1))...
    ./(cs(:,2:end)-cs(:,1:end-1));
psiGradNormDs_sc(:,end) = psiGradNormDs_sc(:,end-1);

% non-orthogonality
betachi_sc = zeros(size(cs));
dpsi = psi(:,3:end)-psi(:,2:end-1);
dr = dpsi./psiGradNorm_sc(:,2:end-1).^2.*psiDr_sc(:,2:end-1);
dz = dpsi./psiGradNorm_sc(:,2:end-1).^2.*psiDz_sc(:,2:end-1);
rtemp = r_sc(:,2:end-1)+dr;
ztemp = z_sc(:,2:end-1)+dz;
thetatemp = acos((R-rtemp)./sqrt((rtemp-R).^2+ztemp.^2));
chitemp = zeros(size(thetatemp));
for i = 1:n_s-2
    chitemp(:,i) = interp1(theta,chi_st(:,i+1),thetatemp(:,i),interpmethod);
end
% chitemp = interp2(cs,ctheta,chi_st,cs(:,3:end),thetatemp,interpmethod);
betachi_sc(:,2:end-1) = (chitemp-cchi(:,2:end-1))./(cs(:,3:end)-cs(:,2:end-1));
betachi_sc(1,:) = 0;
betachi_sc(end,:) = 0;
betachi_sc(:,end) = betachi_sc(:,end-1);

% psiGradNormtemp = interp2(cr,cz,psiGradNorm_rz,rtemp,ztemp,interpmethod);
% npsiGrad2Dpsi_sc = zeros(size(cs));
% npsiGrad2Dpsi_sc(:,2:end-1) = (psiGradNormtemp.^2-psiGradNorm_sc(:,2:end-1).^2)./dpsi;

% figure(1);
% mesh(cs,cchi,betachi_sc);
% 
% betachi_sc2 = zeros(size(cs));
% kernel = (r_sc.*jphi_sc-psiGrad2Dpsi_sc...
%     +(TDpsi_sc./T_sc-qDpsi_sc./q_sc).*psiGradNorm_sc.^2)*pi/(n_chi-1);
% for i = 2:n_chi
%     betachi_sc2(i,:) = sum(kernel(1:i-1,:),1);
% end
% betachi_sc2 = 2*cs*Psi_s.*betachi_sc2./psiGradNorm_sc.^2;
% figure(2);
% mesh(cs,cchi,betachi_sc2);

% s chi derivatives of ln(r^2)
lnr2Ds_sc = zeros(size(cs));
lnr2Dchi_sc = zeros(size(cs));
lnr2Ds_sc(:,1:end-1) = (log(r_sc(:,2:end).^2)-log(r_sc(:,1:end-1).^2))...
    ./(cs(:,2:end)-cs(:,1:end-1));
lnr2Ds_sc(:,end) = lnr2Ds_sc(:,end-1);
lnr2Dchi_sc(1:end-1,:) = (log(r_sc(2:end,:).^2)-log(r_sc(1:end-1,:).^2))...
    ./(cchi(2:end,:)-cchi(1:end-1,:));
lnr2Dchi_sc(end,:) = lnr2Dchi_sc(end-1,:);

% s chi derivatives of ln(psiGradNorm)
lnpsiGradNormDs_sc = zeros(size(cs));
lnpsiGradNormDchi_sc = zeros(size(cs));
lnpsiGradNormDs_sc(:,1:end-1) = (log(psiGradNorm_sc(:,2:end).^2)...
    -log(psiGradNorm_sc(:,1:end-1).^2))./(cs(:,2:end)-cs(:,1:end-1));
lnpsiGradNormDs_sc(:,end) = lnpsiGradNormDs_sc(:,end-1);
lnpsiGradNormDchi_sc(1:end-1,:) = (log(psiGradNorm_sc(2:end,:).^2)...
    -log(psiGradNorm_sc(1:end-1,:).^2))./(cchi(2:end,:)-cchi(1:end-1,:));
lnpsiGradNormDchi_sc(end,:) = lnpsiGradNormDchi_sc(end-1,:);


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
n_s = 14;
n_chi = 15;
% n_s = floor(n_s/2);
% n_chi = ceil(n_chi/2);
ms = 0:1/n_s:1;
mchi = -pi/(2*(n_chi-1)):pi/(n_chi-1):pi+pi/(2*(n_chi-1));
% the totle number of unknow variables
n_x = (3*n_s+1)*(2*n_chi+2);
%total number of cells
n_cell = n_s*n_chi;


%% matrix form of the seven dependent components
i_cell = 1:n_cell;
xij = (2*n_chi+2)*(3*(1:n_s)-3)+2*((1:n_chi)'-1)+1;
i_xij = xij(:);
% XR
XR_coef = [1, 1, 1, 1]/4;
XR = sparse(i_cell,i_xij,XR_coef(1),n_cell,n_x)+...
    sparse(i_cell,i_xij+2,XR_coef(2),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2),XR_coef(3),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+2,XR_coef(4),n_cell,n_x);
% XI
XI_coef = [1, 1, 1, 1]/4;
XI = sparse(i_cell,i_xij+1,XI_coef(1),n_cell,n_x)+...
    sparse(i_cell,i_xij+3,XI_coef(2),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+1,XI_coef(3),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+3,XI_coef(4),n_cell,n_x);
% DXDchiR
mchistep = (mchi(2:end)-mchi(1:end-1))'*ones(n_s,1)';
DXDchiR_coef = [-1, 1, -1, 1]'./(2*mchistep(:)');
DXDchiR = sparse(i_cell,i_xij,DXDchiR_coef(1,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+2,DXDchiR_coef(2,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2),DXDchiR_coef(3,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+2,DXDchiR_coef(4,:),n_cell,n_x);
% DXDchiI
DXDchiI_coef = [-1, 1, -1, 1]'./(2*mchistep(:)');
DXDchiI = sparse(i_cell,i_xij+1,DXDchiI_coef(1,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3,DXDchiI_coef(2,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+1,DXDchiI_coef(3,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+3,DXDchiI_coef(4,:),n_cell,n_x);
% DXDsR
msstep = ones(n_chi,1)*(ms(2:end)-ms(1:end-1));
DXDsR_coef = [-1, -1, 1, 1]'./(2*msstep(:)');
DXDsR = sparse(i_cell,i_xij,DXDsR_coef(1,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+2,DXDsR_coef(2,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2),DXDsR_coef(3,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+2,DXDsR_coef(4,:),n_cell,n_x);
% DXDsI
DXDsI_coef = [-1, -1, 1, 1]'./(2*msstep(:)');
DXDsI = sparse(i_cell,i_xij+1,DXDsI_coef(1,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3,DXDsI_coef(2,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+1,DXDsI_coef(3,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+3,DXDsI_coef(4,:),n_cell,n_x);

i_xij = i_xij + (2*n_chi+2);
% VR
VR_coef = [1, 1]/2;
VR = sparse(i_cell,i_xij,VR_coef(1),n_cell,n_x)+...
    sparse(i_cell,i_xij+2,VR_coef(2),n_cell,n_x);
% VI
VI_coef = [1, 1]/2;
VI = sparse(i_cell,i_xij+1,VI_coef(1),n_cell,n_x)+...
    sparse(i_cell,i_xij+3,VI_coef(2),n_cell,n_x);
% DVDchiR
DVDchiR_coef = [-1, 1]'./mchistep(:)';
DVDchiR = sparse(i_cell,i_xij,DVDchiR_coef(1,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+2,DVDchiR_coef(2,:),n_cell,n_x);
% DVDchiI
DVDchiI_coef = [-1, 1]'./mchistep(:)';
DVDchiI = sparse(i_cell,i_xij,DVDchiI_coef(1,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+2,DVDchiI_coef(2,:),n_cell,n_x);

i_xij = i_xij + (2*n_chi+2);
% YR
YR_coef = [1, 1]/2;
YR = sparse(i_cell,i_xij,YR_coef(1),n_cell,n_x)+...
    sparse(i_cell,i_xij+2,YR_coef(2),n_cell,n_x);
% YI
YI_coef = [1, 1]/2;
YI = sparse(i_cell,i_xij+1,YI_coef(1),n_cell,n_x)+...
    sparse(i_cell,i_xij+3,YI_coef(2),n_cell,n_x);
% DYDchiR
DYDchiR_coef = [-1, 1]'./mchistep(:)';
DYDchiR = sparse(i_cell,i_xij,DYDchiR_coef(1,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+2,DYDchiR_coef(2,:),n_cell,n_x);
% DYDchiI
DYDchiI_coef = [-1, 1]'./mchistep(:)';
DYDchiI = sparse(i_cell,i_xij+1,DYDchiI_coef(1,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3,DYDchiI_coef(2,:),n_cell,n_x);

%% the normalized equilibrium quantities in the middle of each cell
% s and chi the physical components
ps = 1/(2*n_s):1/n_s:1-1/(2*n_s);
pchi = 0:pi/(n_chi-1):pi;
[ps,pchi] = meshgrid(ps,pchi);
% psi
ppsi = ps.^2*Psi_s;
% r z
pr = interp2(cs,cchi,r_sc,ps,pchi);
pz = interp2(cs,cchi,z_sc,ps,pchi);
% the norm of gradient of psi
ppsiGradNorm = interp2(cr,cz,psiGradNorm_rz,pr,pz,interpmethod);
% mass density
prho = 1;
% plasma pressure
pp = interp1(s,p_s,ps,interpmethod);
% psi derivative of p
ppDpsi = interp1(s,pDpsi_s,ps,interpmethod);
% flux function
pT = interp1(s,T_s,ps,interpmethod);
% s derivative of T
pTDs = interp1(s,TDs_s,ps,interpmethod);
% toroidal current density
pjphi = interp2(cr,cz,jphi_rz,pr,pz,interpmethod);
pjphi(:,end)=0;
% safety factor
pq = interp1(s,q_s,ps,interpmethod);
% s derivative of q
pqDs = interp1(s,qDs_s,ps,interpmethod);
% poloidal magnetic field
pBp = ppsi.*pr.^2./(q0*ppsiGradNorm.^2);
% non-orthogonality
pbetachi = interp2(cs,cchi,betachi_sc,ps,pchi,interpmethod);

% dlog(r^2)/ds and dlog(r^2)/dchi
plnr2Ds = interp2(cs,cchi,lnr2Ds_sc,ps,pchi,interpmethod);
plnr2Dchi = interp2(cs,cchi,lnr2Dchi_sc,ps,pchi,interpmethod);
% H defined at the ERATO paper
pH = -2*pjphi.*ppsi.*pr./(ps.*ppsiGradNorm.^2)+pTDs./pT-pqDs./pq;

% the gradient psi direction's derivative of ln psi gradient norm
plnpsiGradNormDs = interp2(cs,cchi,lnpsiGradNormDs_sc,ps,pchi,interpmethod);
plnpsiGradNormDchi = interp2(cs,cchi,lnpsiGradNormDchi_sc,ps,pchi,interpmethod);
plnpsiGradNormDpsi = (plnpsiGradNormDs+pbetachi.*plnpsiGradNormDchi)...
    ./(2*2*ps*Psi_s);
% plnpsiGradNormDpsi = interp2(cr,cz,logpsiGradNormDpsi_rz,pr,pz,interpmethod);
% the gradient psi direction's derivative of log r
plogrDpsi = interp2(cr,cz,lnrDpsi_rz,pr,pz,interpmethod);
plogrDpsi=1/2*(plnr2Ds+pbetachi.*plnr2Dchi)./(2*ps*Psi_s);
% K defined at the ERATO paper
pK = 2*ppsi/q0.*(pjphi.^2./ppsiGradNorm.^2-pjphi./pr.*plnpsiGradNormDpsi...
    -ppDpsi.*plogrDpsi);
% figure(3);mesh(pr,pz,pK);hold on;
% figure(2);mesh(ps,pchi,plogpsiGradNormDpsi);hold on;
%% the matrix form of Is
I1R = 1./pq(:).*DXDchiR-n*XI;
I1I = 1./pq(:).*DXDchiI+n*XR;
I2R = DXDsR+DVDchiR;
I2I = DXDsI+DVDchiI;
I3R = pH(:).*XR-pbetachi(:).*n.*pq(:).*XI...
    +n*pq(:).*VI+pbetachi(:).*DXDchiR+DXDsR;
I3I = pH(:).*XI+pbetachi(:).*n.*pq(:).*XR...
    -n*pq(:).*VR+pbetachi(:).*DXDchiI+DXDsI;
I4R = plnr2Ds(:).*XR+plnr2Dchi(:).*(VR+YR)...
    -n*pq(:).*YI+DXDsR+DVDchiR+DYDchiR;
I4I = plnr2Ds(:).*XI+plnr2Dchi(:).*(VI+YI)...
    +n*pq(:).*YR+DXDsI+DVDchiI+DYDchiI;
I5R = XR;
I5I = XI;
%% the coefficient a b c d e f g h
J = pq.*pr.^2./pT;
acoeff = ones(n_chi,n_s);
acoeff([1,end],:) = 1/2;
a = 2*pq.^2.*ppsi.*pr.^4./(J.^3.*ppsiGradNorm.^2).*acoeff;
b = pT.^2.*pr.^2./(2*Psi_s*J).*acoeff;
c = ppsiGradNorm.^2.*pr.^2.*J./(2*Psi_s).*acoeff;
d = pr.^4*gamma.*pp./(2*Psi_s*J).*acoeff;
e = 2*pK.*pr.^4*q0./J.*acoeff;
f = 2*prho.*ppsi.*pT.*pr.^2./(pq.*ppsiGradNorm.^2).*acoeff;
g = prho.*ppsiGradNorm.^2.*pq.*pr.^4./(2*pT*Psi_s).*acoeff;
h = prho.*pr.^4.*pT.*pq/(2*Psi_s).*acoeff;

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
% subplot(2,4,1);spy(A);
% subplot(2,4,5);spy(B);
% disp([rank(A),rank(B)]);
% A = sparse(A);
% B = sparse(B);
%% symmetry conditions

i_X = (2*n_chi+2)*(3*(1:n_s+1)-3)+1;
i_V = (2*n_chi+2)*(3*(1:n_s)-2)+1;
i_Y = (2*n_chi+2)*(3*(1:n_s)-1)+1;
U = sparse(1:n_x,1:n_x,1,n_x,n_x)+...
    sparse(i_X,i_X+2,1,n_x,n_x)+...%XR
    sparse(i_X+2*n_chi,i_X+2*n_chi-2,1,n_x,n_x)+...
    sparse(i_X+1,i_X+3,-1,n_x,n_x)+...%XI
    sparse(i_X+2*n_chi+1,i_X+2*n_chi-1,-1,n_x,n_x)+...
    sparse(i_V,i_V+2,-1,n_x,n_x)+...%VR
    sparse(i_V+2*n_chi,i_V+2*n_chi-2,-1,n_x,n_x)+...
    sparse(i_V+1,i_V+3,1,n_x,n_x)+...%VI
    sparse(i_V+2*n_chi+1,i_V+2*n_chi-1,1,n_x,n_x)+...
    sparse(i_Y,i_Y+2,-1,n_x,n_x)+...%YR
    sparse(i_Y+2*n_chi,i_Y+2*n_chi-2,-1,n_x,n_x)+...
    sparse(i_Y+1,i_Y+3,-1,n_x,n_x)+...%YI
    sparse(i_Y+2*n_chi+1,i_Y+2*n_chi-1,-1,n_x,n_x);

% transformation using U
Uinv = inv(U);
At = Uinv'*A*Uinv;
Bt = Uinv'*B*Uinv;

% subplot(2,4,2);spy(At);
% subplot(2,4,6);spy(Bt);
% disp([rank(At),rank(Bt)]);
% force the new tilda variables to zero
w2inv = 1e-20; % which make the condition fullfilled when w2~=10^20

i_X = (2*n_chi+2)*(3*(1:n_s+1)-3)+1;
i_V = (2*n_chi+2)*(3*(1:n_s)-2)+1;
i_Y = (2*n_chi+2)*(3*(1:n_s)-1)+1;
i_x = [i_X,i_V,i_Y];
i_x = [i_x,i_x+1,i_x+2*n_chi,i_x+2*n_chi+1];
% symmetry conditions
%tilda A
At(i_x,:) = 0;
At(:,i_x) = 0;
At = At + sparse(i_x,i_x,1,n_x,n_x);
% tilda B
Bt(i_x,:) = 0;
Bt(:,i_x) = 0;
Bt = Bt + sparse(i_x,i_x,w2inv,n_x,n_x);

% the addition condition to make matrix not ill conditioned
i_x = [i_X+2,i_V+3,i_Y+3];
%tilda A
At(i_x,:) = 0;
At(:,i_x) = 0;
At = At + sparse(i_x,i_x,1,n_x,n_x);
% tilda B
Bt(i_x,:) = 0;
Bt(:,i_x) = 0;
Bt = Bt + sparse(i_x,i_x,w2inv,n_x,n_x);

% subplot(2,4,3);spy(At);
% subplot(2,4,7);spy(Bt);
% disp([rank(At),rank(Bt)]);
%% regularity conditions at s = 0
% the procejure is same as symmetry conditions
i_x = 2*(1:n_chi+1)-1 ;
i_x = [i_x,i_x+1];
% tilda A
At(i_x,:) = 0;
At(:,i_x) = 0;
At = At + sparse(i_x,i_x,1,n_x,n_x);
% tilda B
Bt(i_x,:) = 0;
Bt(:,i_x) = 0;
Bt = Bt + sparse(i_x,i_x,w2inv,n_x,n_x);

%% boundary conditions at plasma surface s = 1
% a wall straight on the plasma surface
%the procejure is same as symmetry conditions
i_x =3*n_s*(2*n_chi+2)+2*(1:n_chi+1)-1 ;
i_x = [i_x,i_x+1];
% tilda A
At(i_x,:) = 0;
At(:,i_x) = 0;
At = At + sparse(i_x,i_x,1,n_x,n_x);
% tilda B
Bt(i_x,:) = 0;
Bt(:,i_x) = 0;
Bt = Bt + sparse(i_x,i_x,w2inv,n_x,n_x);

% subplot(2,4,4);spy(At);
% subplot(2,4,8);spy(Bt);
% disp([rank(At),rank(Bt)]);
w2 = eig(full(Bt\At));
% mu0 = -2;
% [x,w2] = eigSolver(At,Bt,1e-8,mu0,rand(n_x,1));
% w2 = w2+mu0;
disp(min(w2));
% figure(1);
% plot(w2);
GrowthRate =w2;


