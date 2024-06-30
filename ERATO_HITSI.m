function ERATO_HITSI(filename)
% load the equilibrium state data from DCON output
data = load(filename);

psi_p = max(max(data.psi));
ppsi = psi_p-data.psi;
R = data.R(end,1);
% (s,chi) coordinate system
% n_s = 14;
% n_chi = n_s+1;
% ps = 1/(2*n_s):1/n_s:1-1/(2*n_s);
n_chi = 15;
ps = sqrt(ppsi/psi_p);
pchi = 0:pi/(n_chi-1):pi;

syms r z positive
for i = 1:length(ps)
    r_path = data.R(i,:)';
    z_path = data.Z(i,:)';
    fz = fit(r_path,z_path,'spline');
    
    
% the analytic form of psi and its relative quantities
fpsi_rz =  sum(coeffvalues(fpsi_r2z2).*[1,r,z,r^2,r*z,z^2,r^3,r^2*z,r*z^2,z^3]);
fpsiDr_rz =  diff(fpsi_rz,r);
fpsiDz_rz =  diff(fpsi_rz,z);
fpsiGradNorm_rz = sqrt(fpsiDr_rz.^2+fpsiDz_rz.^2);
fpsiGrad2Dr_rz = diff(fpsiGradNorm_rz.^2,r);
fpsiGrad2Dz_rz = diff(fpsiGradNorm_rz.^2,z);
fpsiGrad2Dpsi_rz = (fpsiGrad2Dr_rz.*fpsiDr_rz+fpsiGrad2Dz_rz.*fpsiDz_rz)...
    ./fpsiGradNorm_rz.^2;
flnpsiGradNormDpsi_rz = fpsiGrad2Dpsi_rz./(2*fpsiGradNorm_rz.^2);
flnrDpsi_rz = 1./r.*fpsiDr_rz./fpsiGradNorm_rz.^2;
% pressure and its derivative as a function of psi
fp_psi = fit(data.psi',data.P','spline');
pp = fp_psi(ppsi);
ppDpsi = differentiate(fp_psi,ppsi);
% T and its derivative as a function of psi
fT_psi = fit(data.psi',data.T','spline');
pT = fT_psi(ppsi);
pTDpsi = differentiate(fT_psi,ppsi);
% the toroidal current density
fjphi_rz = r.*fpDpsi_psi(fpsi_rz)...
    +fTDpsi_psi(fpsi_rz).*fT_psi(fpsi_rz)./r;

for i = 1:n_anl   
    
    
end
    