


% parameters determine the equilibrium
E = 2; a = 1/3; q0 = 0.3; n = 2; n_r = 500;

%% equilibrium state
R = 1;  B = 1; T0 = B*R; gamma = 5/3;

Psi_s = E*a^2*B/(q0*2); %the psi value on the plasma surface
syms r z
fpsi_rz = Psi_s/(a^2*R^4)*(z^2*r^2/E^2+(r^2-R^2)^2/4);
fpsiDr_rz =  diff(fpsi_rz,r);
fpsiDz_rz =  diff(fpsi_rz,z);
fpsiGradNorm_rz = sqrt(fpsiDr_rz.^2+fpsiDz_rz.^2);
fpsiGrad2Dr_rz = diff(fpsiGradNorm_rz.^2,r);
fpsiGrad2Dz_rz = diff(fpsiGradNorm_rz.^2,z);
fpsiGrad2Dpsi_rz = (fpsiGrad2Dr_rz.*fpsiDr_rz+fpsiGrad2Dz_rz.*fpsiDz_rz)...
    ./fpsiGradNorm_rz.^2;
flnpsiGradNormDpsi_rz = fpsiGrad2Dpsi_rz./(2*fpsiGradNorm_rz.^2);
flnrDpsi_rz = 1./r.*fpsiDr_rz./fpsiGradNorm_rz.^2;

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
flnpsiGradNormDpsi_rz = matlabFunction(flnpsiGradNormDpsi_rz);
flnrDpsi_rz = matlabFunction(flnrDpsi_rz);

% the boundary of plasma
r_min = R*sqrt(1-2*a); r_max = R*sqrt(1+2*a);
r_temp = R*(1-4*a^2)^(1/4);
z_max = sqrt((4*R^4*a^2-(r_temp^2-R^2)^2)*E^2/(4*r_temp^2));
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
lnpsiGradNormDpsi_rz = flnpsiGradNormDpsi_rz(cr,cz);
% the gradient psi direction's derivative of log r
lnrDpsi_rz = flnrDpsi_rz(cr,cz);
% toroidal current density
jphi_rz = cr.*fpDpsi_psi(psi_rz);
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
interpmethod = 'spline';
% theta as an independent coordinate
n_theta = n_s*5;
theta = linspace(0,pi,n_theta);

[cs,ctheta] = meshgrid(s,theta);

% the corresponding rho matrix
% the initialize vector of rho
rho_st = zeros(size(ctheta));
for i = 2:n_s% exclude the magnetic axis situation
    psi_i = psi(i);
    
    % the start and end point of the loop integration
    path_left = fzero(@(r) fpsi_rz(r,0)-psi_i,[r_min*0.9,R]);
    path_right = fzero(@(r) fpsi_rz(r,0)-psi_i,[R,r_max*1.1]);
       
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
psiGradNorm_st = interp2(cr,cz,psiGradNorm_rz,r_st,z_st,interpmethod);
% length of path element
dl_st = zeros(size(cs));
dl_st(2:end,:) = sqrt((rho_st(2:end,:)-rho_st(1:end-1,:)).^2+...
    ((rho_st(2:end,:)+rho_st(1:end-1,:))/2.*...
    (ctheta(2:end,:)-ctheta(1:end-1,:))).^2);
% calculate q from loop integration
q_s = 1/pi*sum(T_st./(r_st.*psiGradNorm_st).*dl_st);
q_s(1) = q_s(2);

fq_psi = polyfit(psi,q_s,5);
q_s = polyval(fq_psi,psi);
% derivative of q
fqDpsi_psi = polyder(fq_psi);
% qDs_s = zeros(size(q_s));
% qDs_s(1:end-1) = (q_s(2:end)-q_s(1:end-1))./(s(2:end)-s(1:end-1));
% qDs_s(end) = qDs_s(end-1);
qDpsi_s = polyval(fqDpsi_psi,psi);
qDs_s = qDpsi_s*2.*s*Psi_s;

% the corresponding chi
q_st = ones(n_theta,1)*q_s;
% calculate integral kernel
kernel = T_st./(q_st.*r_st.*psiGradNorm_st).*dl_st;
% initialize of chi_st
chi_st = zeros(size(ctheta));
% calculate the chi value of this point
for k = 2:n_theta
    chi_st(k,:) = sum(kernel(1:k,:));
end
chi_st(:,1) = chi_st(:,2);
% make sure the biggest chi at fixed s is not smaller than pi
chi_st(end,chi_st(end,:)<pi) = pi;
%% s chi coordinate system
% chi as a independent variable
n_chi = floor(n_theta/5);
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
r_sc = R+rho_sc.*cos(theta_sc);
z_sc = rho_sc.*sin(theta_sc);

% the gradient psi at every point
psiDr_sc = interp2(cr,cz,psiDr_rz,r_sc,z_sc,interpmethod);
psiDz_sc = interp2(cr,cz,psiDz_rz,r_sc,z_sc,interpmethod);
% the norm of gradient of psi
psiGradNorm_sc = interp2(cr,cz,psiGradNorm_rz,r_sc,z_sc,interpmethod);

% non-orthogonality
betachi_sc = zeros(size(cs));
dpsi = psi(:,3:end)-psi(:,2:end-1);
dr = dpsi./psiGradNorm_sc(:,2:end-1).^2.*psiDr_sc(:,2:end-1);
dz = dpsi./psiGradNorm_sc(:,2:end-1).^2.*psiDz_sc(:,2:end-1);
rtemp = r_sc(:,2:end-1)+dr;
ztemp = z_sc(:,2:end-1)+dz;
thetatemp = acos((rtemp-R)./sqrt((rtemp-R).^2+ztemp.^2));
chitemp = zeros(size(thetatemp));
for i = 1:n_s-2
    chitemp(:,i) = interp1(theta,chi_st(:,i+1),thetatemp(:,i),interpmethod);
end
% chitemp = interp2(cs,ctheta,chi_st,cs(:,3:end),thetatemp,interpmethod);
betachi_sc(:,2:end-1) = (chitemp-cchi(:,2:end-1))./(cs(:,3:end)-cs(:,2:end-1));
betachi_sc(1,:) = 0;
betachi_sc(end,:) = 0;
betachi_sc(:,end) = betachi_sc(:,end-1);
% figure(1);mesh(cs,cchi,betachi_sc);

% jphi
jphi_sc = interp2(cr,cz,jphi_rz,r_sc,z_sc,interpmethod);
% T and dTdpsi
T_sc = T0*ones(size(cs));
TDpsi_sc = 0;
% logpsiGradNormDpsi_sc = interp2(cr,cz,lnpsiGradNormDpsi_rz,r_sc,z_sc,interpmethod);
psiGrad2Dpsi_sc = interp2(cr,cz,psiGrad2Dpsi_rz,r_sc,z_sc,interpmethod);
% q and dqdpsi
q_sc = ones(n_chi,1)*q_s;
qDpsi_sc = ones(n_chi,1)*qDpsi_s;

kernel = ((-r_sc.*jphi_sc-psiGrad2Dpsi_sc)./psiGradNorm_sc.^2+...
     TDpsi_sc./T_sc)*pi/(n_chi-1);
qDpsi_sc = q_sc/pi.*sum(kernel,1);

betachi_sc = zeros(size(cs));
kernel = 2*cs*Psi_s.*((-r_sc.*jphi_sc-psiGrad2Dpsi_sc)./psiGradNorm_sc.^2+...
     TDpsi_sc./T_sc-qDpsi_sc./q_sc)*pi/(n_chi-1);
for i = 2:n_chi
    betachi_sc(i,:) = sum(kernel(1:i-1,:),1);
end
% figure(2);
% mesh(cs,cchi,betachi_sc);

% s chi derivatives of ln(r^2)
lnr2Ds_sc = zeros(size(cs));
lnr2Dchi_sc = zeros(size(cs));
lnr2Ds_sc(:,2:end-1) = (log(r_sc(:,3:end).^2)-log(r_sc(:,1:end-2).^2))...
    ./(cs(:,3:end)-cs(:,1:end-2));
lnr2Ds_sc(:,end) = lnr2Ds_sc(:,end-1);
lnr2Ds_sc(:,1) = lnr2Ds_sc(:,2);
lnr2Dchi_sc(2:end-1,:) = (log(r_sc(3:end,:).^2)-log(r_sc(1:end-2,:).^2))...
    ./(cchi(3:end,:)-cchi(1:end-2,:));
lnr2Dchi_sc(end,:) = lnr2Dchi_sc(end-1,:);
lnr2Dchi_sc(1,:) = lnr2Dchi_sc(2,:);

% s derivative of T
TDs_s = zeros(size(s));
TDs_s(1:end-1) = (T_s(2:end)-T_s(1:end-1))./(s(2:end)-s(1:end-1));
TDs_s(end) = TDs_s(end-1);

% s derivative of q
% qDs_s = zeros(size(s));
% qDs_s(2:end-1) = (q_s(3:end)-q_s(1:end-2))./(s(3:end)-s(1:end-2));
% qDs_s(end) = qDs_s(end-1);
% qDs_s = qDpsi_s.*2.*s*Psi_s;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct of matrixes A and B


