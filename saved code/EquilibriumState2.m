%% set global variable
global R E a Psi_s
%% the parameters determine the equilibrium
R = 1;  B = 1; T0 = B*R;
lamda = 1; gamma = 5/3;
%%%% parameters determine the equilibrium 
E = 2; a = 1/3; q0 = 0.8; n = 1;

Psi_s = E*a^2*B/(q0*2); %the psi value on the plasma surface
% the boundary of plasma
r_min = sqrt(R^2-2*a*R); r_max = sqrt(R^2+2*a*R);
r_temp = R*(1-4*a^2/R^2)^(1/4);
z_max = sqrt((4*R^2*a^2-(r_temp^2-R^2)^2)*E^2/(4*r_temp^2));
rho_max = sqrt(max(R-r_min,r_max-R)^2+z_max^2);
%% r z coordinate system
%the dimensions of matrixes to fit
n_r = 500;
n_z = 250;
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



% figure(1);
% plot(r_sc,z_sc);
% hold on
% plot(r_sc',z_sc')
% figure(2);
% mesh(cs,cchi,betachi_sc);
