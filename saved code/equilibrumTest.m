syms Psi_s a R r z E
psi=Psi_s/(a^2*R^4)*((r^2)*z^2/E^2+(r^2-R^2)^2/4);
psiDr=diff(psi,r);
jphi=-(r*diff(psiDr/r,r)+diff(psi,z,2))/r;
dpdpsi=simplify(jphi/r);