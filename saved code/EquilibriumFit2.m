% parameters determine the equilibrium
E = 2; a = 1/3; q0 = 0.3; n = 2; n_r = 100; n_s = 14;

%% equilibrium state
R = 1;  B0 = 1; T0 = B0*R; gamma = 5/3;

Psi_s = E*a^2*R^2*B0/(2*q0); %the psi value on the plasma surface
syms r z psis positive
fpsi_rz = Psi_s/(a^2*R^4)*(z^2*r^2/E^2+(r^2-R^2)^2/4);

fpsiDr_rz = diff(fpsi_rz,r);
fpsiDz_rz = diff(fpsi_rz,z);
fpsiGradNorm_rz = sqrt(fpsiDr_rz.^2+fpsiDz_rz.^2);
% pressure functin
fp_rz = 2*Psi_s*(1+E^2)*(Psi_s-fpsi_rz)/(a^2*E^2*R^4);
fpDr_rz = diff(fp_rz,r);
fpDz_rz = diff(fp_rz,z);
fpDpsi_rz = (fpDr_rz.*fpsiDr_rz+fpDz_rz.*fpsiDz_rz)./fpsiGradNorm_rz.^2;
% T function
fT_rz = T0;
fTDr_rz = diff(fT_rz,r);
fTDz_rz = diff(fT_rz,z);
fTDpsi_rz = (fTDr_rz.*fpsiDr_rz+fTDz_rz.*fpsiDz_rz)./fpsiGradNorm_rz.^2;
% the toroidal current density
fjphi_rz = r.*fpDpsi_rz+fTDpsi_rz.*fT_rz./r;

% the (r,z) coordinate
n_z = ceil(n_r/2);
cr = linspace(0.5,1.5,n_r);
cz = linspace(0,1,n_z);
[cr,cz] = meshgrid(cr,cz);
Mpsi_rz = double(subs(fpsi_rz,{r,z},{cr,cz}));
r_max = max(max(cr));
r_min = min(min(cr));
z_max = max(max(cz));
% s coordinate
n_chi = n_s+1;
ps = 1/(2*n_s):1/n_s:1-1/(2*n_s);
pchi = 0:pi/(n_chi-1):pi;
% psi
ppsi = ps.^2*Psi_s;

% initialization
%safety factor
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

% fit near the magnetic axis
n_in = ceil(sqrt(n_s));
psi_in = ps(n_in).^2*Psi_s;
Msk_in = Mpsi_rz<psi_in;
fit_in = fit([cr(Msk_in).^2,cz(Msk_in).^2],Mpsi_rz(Msk_in),'poly33');
fpsi_rz = fit_in(r.^2,z.^2);
fpsiDr_rz =  diff(fpsi_rz,r);
fpsiDz_rz =  diff(fpsi_rz,z);
fpsiGradNorm_rz = sqrt(fpsiDr_rz.^2+fpsiDz_rz.^2);
fpsiGrad2Dr_rz = diff(fpsiGradNorm_rz.^2,r);
fpsiGrad2Dz_rz = diff(fpsiGradNorm_rz.^2,z);
fpsiGrad2Dpsi_rz = (fpsiGrad2Dr_rz.*fpsiDr_rz+fpsiGrad2Dz_rz.*fpsiDz_rz)...
    ./fpsiGradNorm_rz.^2;
flnpsiGradNormDpsi_rz = fpsiGrad2Dpsi_rz./(2*fpsiGradNorm_rz.^2);
flnrDpsi_rz = 1./r.*fpsiDr_rz./fpsiGradNorm_rz.^2;

for i = 1:n_in-1
    fpath = matlabFunction(fpsi_rz-ppsi(i));
    r_left = fzero(@(r) fpath(r,0),[r_min,R]);
    r_right = fzero(@(r) fpath(r,0),[R,r_max]);
    
    [fz_r,~,~] = solve(fpsi_rz-ppsi(i),z,'ReturnConditions',true,'MaxDegree',3);
    % the corresponding z
    zMsk = double(subs(fz_r,r,R))>0&double(subs(fz_r,r,R))<z_max&...
        abs(imag(double(subs(fz_r,r,R))))<1e-100;
    z_path = fz_r(zMsk);
    zDr_path = diff(z_path,r);   
    % the integral kernal of q    
    % the quantites along the 'psi = constant' path
    T_path = subs(fT_rz,z,z_path);
    TDpsi_path = subs(fTDpsi_rz,z,z_path);
    psiGradNorm_path = subs(fpsiGradNorm_rz,z,z_path);
    jphi_path = subs(fjphi_rz,z,z_path);
    psiGrad2Dpsi_path = subs(fpsiGrad2Dpsi_rz,z,z_path);
    % the integral kernal of q
    fint_q = -T_path./(r.*psiGradNorm_path).*sqrt(1+zDr_path.^2);
%     pq(i) = double(1/pi*int(-fint_q,r_right,r_left));
    kernel_q = matlabFunction(fint_q);
    pq(i) = real(1/pi*integral(@(rr) kernel_q(rr),r_right,r_left,'AbsTol',1e-6,'RelTol',1e-3));     
    % the integral kernal of chi
%     fint_chi = fint_q/pq(i);
%     chi_r = int(-fint_chi,r_right,r);
    kernel_chi = matlabFunction(fint_q/pq(i));
    chi_r = @(r) integral(@(rr) kernel_chi(rr),r_right,r);
    % dqdpsi
    fint_qDpsi = ((-jphi_path.*r-psiGrad2Dpsi_path)./psiGradNorm_path.^2+...
        TDpsi_path./T_path).*fint_chi;
%     pqDpsi(i) = double(pq(i)/pi*int(-fint_qDpsi,r_right,r_left));
%     pqDpsi(i) = real(pqDpsi(i));
    kernel_qDpsi = matlabFunction(fint_qDpsi);
    pqDpsi(i) = pq(i)/pi*integral(@(rr) kernel_qDpsi(rr),r_right,r_left,'AbsTol',1e-6,'RelTol',1e-3);
    pqDpsi(i) = real(pqDpsi(i));
    
    % the crosponding r of each mesh
    pr(1,i) = r_right;
    pr(end,i) = r_left;
    for j = 2:length(pchi)-1
        % r,z
        pr(j,i) = fzero(@(r) chi_r(r)-pchi(j),[r_left,R+(r_right-R)*cos(pi/(4*n_chi))]);
        pz(j,i) = double(subs(fz_rpsi,{psis,r},{ppsi(i),pr(j,i)}));
        % betachi
        fint_betachi = ((-jphi_path.*r-psiGrad2Dpsi_path)./psiGradNorm_path.^2+...
        TDpsi_path./T_path-pqDpsi(i)/pq(i)).*fint_chi;
        kernel_betachi = matlabFunction(simplify(fint_betachi));
        pbetachi(j,i) = 2*ps(i)*Psi_s*integral(@(rr) kernel_betachi(rr),r_right,pr(j,i));
        % lnr2Dchi
        plnr2Dchi(j,i) = double(subs(2./r.*1./fint_chi,r,pr(j,i)));      
    end  
    
    
    
end

