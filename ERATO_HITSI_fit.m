function ERATO_HITSI_fit(filename)
% load the equilibrium state data from DCON output
data = load(filename);

psi_p = max(max(data.psi));
data.psi = psi_p-data.psi;
data.R = data.R(:,1:data.ntheta/2+1);
data.Z = data.Z(:,1:data.ntheta/2+1);
R = data.R(end,1);
zmax = max(max(data.Z));
% (s,chi) coordinate system
n_s = 14;
n_chi = n_s+1;
ps = 1/(2*n_s):1/n_s:1-1/(2*n_s);
pchi = 0:pi/(n_chi-1):pi;

% psi
ppsi = psi_p*ps.^2;
% pressure and its derivative as a function of psi
fp_psi = fit(data.psi',data.P','spline');
pp = fp_psi(ppsi);
ppDpsi = differentiate(fp_psi,ppsi);
% T and its derivative as a function of psi
fT_psi = fit(data.psi',data.T','spline');
pT = fT_psi(ppsi);
pTDpsi = differentiate(fT_psi,ppsi);
%----initionalize 
% safety factor
pq = zeros(size(ppsi));
%psi derivative of q
pqDpsi = zeros(size(ppsi));
%r z
pr = zeros(length(pchi),length(ps));
pz = zeros(length(pchi),length(ps));
% the non-orthogonality betachi
pbetachi = zeros(length(pchi),length(ps));
% lnr2Dchi
plnr2Dchi = zeros(length(pchi),length(ps));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the analytic area near magnetic axis
n_anl = floor(sqrt(n_s));
Msk_anl = data.psi<ppsi(n_anl);
psi_anl = data.psi(Msk_anl);
r_anl = data.R(Msk_anl,:);
z_anl = data.Z(Msk_anl,:);
rmin_anl = min(min(r_anl));
rmax_anl = max(max(r_anl));
zmax_anl = max(max(z_anl));
Mpsi_anl = repmat(psi_anl',1,data.ntheta/2+1);
fpsi_r2z2 = fit([r_anl(:).^2,z_anl(:).^2],Mpsi_anl(:),'poly33');
syms r z spsi positive
% the analytic form of psi and its relative quantities
fpsi_rz =  sum(coeffvalues(fpsi_r2z2).*[1,r,z,r^2,r*z,z^2,r^3,r^2*z,r*z^2,z^3].^2);
fpsiDr_rz =  diff(fpsi_rz,r);
fpsiDz_rz =  diff(fpsi_rz,z);
fpsiGradNorm_rz = sqrt(fpsiDr_rz.^2+fpsiDz_rz.^2);
fpsiGrad2Dr_rz = diff(fpsiGradNorm_rz.^2,r);
fpsiGrad2Dz_rz = diff(fpsiGradNorm_rz.^2,z);
fpsiGrad2Dpsi_rz = (fpsiGrad2Dr_rz.*fpsiDr_rz+fpsiGrad2Dz_rz.*fpsiDz_rz)...
    ./fpsiGradNorm_rz.^2;
flnpsiGradNormDpsi_rz = fpsiGrad2Dpsi_rz./(2*fpsiGradNorm_rz.^2);
flnrDpsi_rz = 1./r.*fpsiDr_rz./fpsiGradNorm_rz.^2;
% z as a function of r and psi
[fz_rpsi,~,~] = solve(fpsi_rz-spsi,z,'ReturnConditions',true,'MaxDegree',3);

% the left and right point of the 'psi = constant' path
[fr0_psi,~,~] = solve(subs(fpsi_rz-spsi,z,0),r,'ReturnConditions',true);

for i = 1:n_anl 
    % the left and right point of the 'psi = constant' path
    fr0 = double(subs(fr0_psi,spsi,ppsi(i)));
    r_right = fr0(fr0<1.1*rmax_anl&fr0>R);
    r_left = fr0(fr0<R&fr0>0.9*rmin_anl);
    % z function on the psi=const surface
    fz = subs(fz_rpsi,spsi,ppsi(i));
    fzmsk = double(subs(fz,r,R))>0&double(subs(fz,r,R))<=1.1*zmax_anl;
    z_path = fz(fzmsk);
    zDr_path = diff(z_path,r);
    % plasma pressure and it's derivative
    pDpsi_path = ppDpsi(i);
    % T and it's derivative
    T_path = pT(i);
    TDpsi_path = pTDpsi(i);
    % norm of psi gradient
    psiGradNorm_path = subs(fpsiGradNorm_rz,z,z_path);
    % the toroidal current density
    jphi_path = r.*pDpsi_path+TDpsi_path.*T_path./r;
    % the derivative of square of psi gradient
    psiGrad2Dpsi_path = subs(fpsiGrad2Dpsi_rz,z,z_path);
    % the integral kernal of q
    fint_q = -T_path./(r.*psiGradNorm_path).*sqrt(1+zDr_path.^2);
%     pq(i) = double(1/pi*int(-fint_q,r_right,r_left));
    kernel_q = matlabFunction(fint_q);
    pq(i) = 1/pi*integral(@(rr) kernel_q(rr),r_right,r_left);  
    % the integral kernal of chi
    fint_chi = -T_path./(pq(i)*r.*psiGradNorm_path).*sqrt(1+zDr_path.^2);
%     chi_r = int(-fint_chi,r_right,r);
    kernel_chi = matlabFunction(fint_chi);
    chi_r = @(r) integral(@(rr) kernel_chi(rr),r_right,r);
    % dqdpsi
    fint_qDpsi = ((-jphi_path.*r-psiGrad2Dpsi_path)./psiGradNorm_path.^2+...
        TDpsi_path./T_path).*fint_chi;
%     pqDpsi(i) = double(pq(i)/pi*int(-fint_qDpsi,r_right,r_left));
%     pqDpsi(i) = real(pqDpsi(i));
    kernel_qDpsi = matlabFunction(fint_qDpsi);
    pqDpsi(i) = pq(i)/pi*integral(@(rr) kernel_qDpsi(rr),r_right,r_left);
    pqDpsi(i) = real(pqDpsi(i));
    % the crosponding r of each mesh
    pr(1,i) = r_right;
    pr(end,i) = r_left;
    % betachi
    fint_betachi = ((-jphi_path.*r-psiGrad2Dpsi_path)./psiGradNorm_path.^2+...
        TDpsi_path./T_path-pqDpsi(i)/pq(i)).*fint_chi;
    kernel_betachi = matlabFunction(fint_betachi);
    for j = 2:length(pchi)-1
        % r,z
        pr(j,i) = fzero(@(r) chi_r(r)-pchi(j),[r_left,R+(r_right-R)*cos(pi/(4*n_chi))]);
        pz(j,i) = double(subs(z_path,r,pr(j,i)));
        % betachi
        pbetachi(j,i) = 2*ps(i)*psi_p*integral(@(rr) kernel_betachi(rr),r_right,pr(j,i));
        % lnr2Dchi
        plnr2Dchi(j,i) = double(subs(2./r.*1./fint_chi,r,pr(j,i)));      
    end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the aera not near the magnetic axis
for i = n_anl+1:n_s
%     psi_range = [(ppsi(i)+ppsi(i-1))/2,(ppsi(i)+ppsi(i+1))/2];
%     Msk_i = data.psi>psi_range(1)&data.psi<psi_range(2);
    Msk_i = [false,false,data.psi(4:end)<ppsi(i)&data.psi(1:end-3)>ppsi(i)];  
    psi_i = data.psi(Msk_i);
    r_i = data.R(Msk_i,:);
    z_i = data.Z(Msk_i,:);
    for j = 1:size(r_i,2)-1
        Mpsi_ij = repmat(psi_i',1,2);
        r_ij = r_i(:,j+(0:1));
        z_ij = z_i(:,j+(0:1));
        fpsi_r2z2 = fit([r_ij(:).^2,z_ij(:).^2],Mpsi_ij(:),'poly22');
        % the analytic form of psi and its relative quantities
        fpsi_rz =  sum(coeffvalues(fpsi_r2z2).*[1,r,z,r^2,r*z,z^2].^2);
        fpsiDr_rz =  diff(fpsi_rz,r);
        fpsiDz_rz =  diff(fpsi_rz,z);
        fpsiGradNorm_rz = sqrt(fpsiDr_rz.^2+fpsiDz_rz.^2);
        fpsiGrad2Dr_rz = diff(fpsiGradNorm_rz.^2,r);
        fpsiGrad2Dz_rz = diff(fpsiGradNorm_rz.^2,z);
        fpsiGrad2Dpsi_rz = (fpsiGrad2Dr_rz.*fpsiDr_rz+fpsiGrad2Dz_rz.*fpsiDz_rz)...
            ./fpsiGradNorm_rz.^2;
        flnpsiGradNormDpsi_rz = fpsiGrad2Dpsi_rz./(2*fpsiGradNorm_rz.^2);
        flnrDpsi_rz = 1./r.*fpsiDr_rz./fpsiGradNorm_rz.^2;
        
        
        % z as a function of r and psi
        [fz,~,~] = solve(fpsi_rz-ppsi(i),z,'ReturnConditions',true);
        fzmsk = double(subs(fz,r,r_ij(2,1)))>0&double(subs(fz,r,r_ij(2,1)))<=zmax;
        fz = fz(fzmsk);
        % right boundary of the fit area
        fbound_right = polyfit(r_ij(1:2,1),z_ij(1:2,1),1);
        fr_right = matlabFunction(fz-fbound_right*[r,1]');
        r_right = fzero(fr_right,[r_ij(2,1),r_ij(1,1)]);
        % left boundary of the fit area
        fbound_left = polyfit(r_ij(1:2,end),z_ij(1:2,end),1);
        fr_left = matlabFunction(fz-fbound_left*[r,1]');
        r_left = fzero(fr_left,[r_ij(2,end),r_ij(1,end)]);
        
    end
end
    