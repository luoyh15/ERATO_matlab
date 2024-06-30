function psi = fpsi_rz(r,z)

% global cr cz M_psi
% %% magnetic flux profile or the psi coordinate
% psi = interp2(cr,cz,M_psi,r,z,'spline',2);
global Psi_s a R E
%% psi
psi = Psi_s/(a^2*R^2)*(r.^2.*z.^2/E^2+(r.^2-R^2).^2/4);