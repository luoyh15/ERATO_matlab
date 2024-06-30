% parameters determine the equilibrium
E = 2; a = 1/3; q0 = 0.3; n = 2; n_r = 500;

%% equilibrium state
R = 1;  B = 1; T0 = B*R; gamma = 5/3;

Psi_s = E*a^2*B/(q0*2); %the psi value on the plasma surface
syms r z psi positive
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

% the function defines the 'psi = constant' path
fz_rpsi = solve(fpsi_rz-psi,z);
fzDr_rpsi = diff(fz_rpsi,r);

% the left and right point of the 'psi = constant' path
fr0_psi = solve(subs(fpsi_rz-psi,z,0),r);
frleft_psi = fr0_psi(2);
frright_psi = fr0_psi(1);

%% s chi coordinate system
n_s = 14;
n_chi = 15;
cs = linspace(0,1,n_s);
cpsi = cs.^2*Psi_s;
cchi = linspace(0,pi,n_chi);
% betachi
betachi = zeros(n_chi,n_s);

%% equilibrium quantites
%the coordinate
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
for i = 1:length(ps)
    % the start and end point of the loop integration
    r_left = double(subs(frleft_psi,psi,ppsi(i)));
    r_right = double(subs(frright_psi,psi,ppsi(i)));
    % the corresponding z
    z_path = subs(fz_rpsi,psi,ppsi(i));
    % the quantites along the 'psi = constant' path
    T_path = subs(fT_rz,z,z_path);
    TDpsi_path = subs(fTDpsi_rz,z,z_path);
    psiGradNorm_path = subs(fpsiGradNorm_rz,z,z_path);
    zDr_path = subs(fzDr_rpsi,psi,ppsi(i));
    jphi_path = subs(fjphi_rz,z,z_path);
    psiGrad2Dpsi_path = subs(fpsiGrad2Dpsi_rz,z,z_path);
    % the integral kernal of q
    fint_q = T_path./(r.*psiGradNorm_path).*sqrt(1+zDr_path.^2);
%     pq(i) = double(1/pi*int(-fint_q,r_right,r_left));
    kernel_q = matlabFunction(fint_q);
    pq(i) = 1/pi*integral(@(rr) -kernel_q(rr),r_right,r_left);     
    % the integral kernal of chi
    fint_chi = T_path./(pq(i)*r.*psiGradNorm_path).*sqrt(1+zDr_path.^2);
%     chi_r = int(-fint_chi,r_right,r);
    kernel_chi = matlabFunction(simplify(fint_chi));
    chi_r = @(r) integral(@(rr) -kernel_chi(rr),r_right,r);
    % dqdpsi
    fint_qDpsi = ((-jphi_path.*r-psiGrad2Dpsi_path)./psiGradNorm_path.^2+...
        TDpsi_path./T_path).*fint_chi;
%     pqDpsi(i) = double(pq(i)/pi*int(-fint_qDpsi,r_right,r_left));
%     pqDpsi(i) = real(pqDpsi(i));
    kernel_qDpsi = matlabFunction(simplify(fint_qDpsi));
    pqDpsi(i) = pq(i)/pi*integral(@(rr)-kernel_qDpsi(rr),r_right,r_left);
    pqDpsi(i) = real(pqDpsi(i));
    
    % the crosponding r of each mesh
    pr(1,i) = r_right;
    pr(end,i) = r_left;
    for j = 2:length(pchi)-1
        % r,z
        pr(j,i) = fzero(@(r) chi_r(r)-pchi(j),[r_left,R+(r_right-R)*cos(pi/n_chi)]);
        pz(j,i) = double(subs(fz_rpsi,{psi,r},{ppsi(i),pr(j,i)}));
        % betachi
        fint_betachi = ((-jphi_path.*r-psiGrad2Dpsi_path)./psiGradNorm_path.^2+...
        TDpsi_path./T_path-pqDpsi(i)/pq(i)).*fint_chi;
        kernel_betachi = matlabFunction(simplify(fint_betachi));
        pbetachi(j,i) = integral(@(rr) -kernel_betachi(rr),r_right,pr(j,i));
        % lnr2Dchi
        plnr2Dchi(j,i) = double(subs(2./r.*1./fint_chi,r,pr(j,i)));
       
    end
        
    
end
[ps,pchi] = meshgrid(ps,pchi);
% lnr2Ds
plnr2Ds = -pbetachi.*plnr2Dchi+...
    double(subs(flnrDpsi_rz,{r,z},{pr,pz}).*(4*ps*Psi_s));
% jphi
pjphi = double(subs(fjphi_rz,{r,z},{pr,pz}));
% pressure
pp = double(subs(fp_rz,{r,z},{pr,pz}));
% psi derivative of p
ppDpsi = double(subs(fpDpsi_rz,{r,z},{pr,pz}));
% flux function
pT = double(subs(fT_rz,{r,z},{pr,pz}));
% s derivative of T
pTDs = double(subs(fTDpsi_rz,{r,z},{pr,pz})).*(2*ps*Psi_s);
% safety factor
pq = ones(n_chi,1)*pq;
% s derivative of q
pqDs = pqDpsi.*(2*ps*Psi_s);
% psiGrad2
ppsiGrad2 = double(subs(fpsiGradNorm_rz.^2,{r,z},{pr,pz}));
% H defined at the ERATO paper
pH = 2*pjphi.*ppsi.*pr./(ps.*ppsiGrad2)+pTDs./pT-pqDs./pq;

% lnpsiGradNormDpsi
plnpsiGradNormDpsi = double(subs(flnpsiGradNormDpsi_rz,{r,z},{pr,pz}));
% lnrDpsi
plnrDpsi = double(subs(flnrDpsi_rz,{r,z},{pr,pz}));
% mass density
prho = 1;
mesh(pr,pz,pbetachi);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct of matrixes A and B


