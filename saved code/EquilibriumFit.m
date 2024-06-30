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
% presure functin
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
fpsi_in = fit_in(r.^2,z.^2);

for i = 1:n_in-1
    fpsi_rz = fpsi_in;
    fpsiDr_rz =  diff(fpsi_rz,r);
    fpsiDz_rz =  diff(fpsi_rz,z);
    fpsiGradNorm_rz = sqrt(fpsiDr_rz.^2+fpsiDz_rz.^2);
    fpsiGrad2Dr_rz = diff(fpsiGradNorm_rz.^2,r);
    fpsiGrad2Dz_rz = diff(fpsiGradNorm_rz.^2,z);
    fpsiGrad2Dpsi_rz = (fpsiGrad2Dr_rz.*fpsiDr_rz+fpsiGrad2Dz_rz.*fpsiDz_rz)...
        ./fpsiGradNorm_rz.^2;
    flnpsiGradNormDpsi_rz = fpsiGrad2Dpsi_rz./(2*fpsiGradNorm_rz.^2);
    flnrDpsi_rz = 1./r.*fpsiDr_rz./fpsiGradNorm_rz.^2;
    
    fpath = matlabFunction(fpsi_rz-ppsi(i));
    r_left = fzero(@(r) fpath(r,0),[r_min,R]);
    r_right = fzero(@(r) fpath(r,0),[R,r_max]);
    % the integral kernal of q
    fint_q = matlabFunction(1/pi*fT_rz./(r.*fpsiGradNorm_rz));
    pq(i) = IntegralPath(fpath,fint_q,r_left,0,r_right,0,z_max,1e-6);
    
    % the integral kernal of chi
    fint_chi = -T_path./(pq(i)*r.*psiGradNorm_path).*sqrt(1+zDr_path.^2);
%     chi_r = int(-fint_chi,r_right,r);
    kernel_chi = matlabFunction(simplify(fint_chi));
    chi_r = @(r) integral(@(rr) kernel_chi(rr),r_right,r);
   
    

end

