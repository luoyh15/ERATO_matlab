function 
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