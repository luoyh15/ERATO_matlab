function GrowthRate = ERATO17(E,a,q0,n,n_s)


% parameters determine the equilibrium
% E = 2; a = 1/3; q0 = 0.3; n = 2; n_r = 500;

%% equilibrium state
R = 1;  B = 1; T0 = B*R; gamma = 5/3;

Psi_s = E*a^2*B/(q0*2); %the psi value on the plasma surface
syms r z psis positive
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
fz_rpsi = solve(fpsi_rz-psis,z);
fzDr_rpsi = diff(fz_rpsi,r);

% the left and right point of the 'psi = constant' path
fr0_psi = solve(subs(fpsi_rz-psis,z,0),r);
frleft_psi = fr0_psi(2);
frright_psi = fr0_psi(1);

%% s chi coordinate system
% n_s = 14;
n_chi = n_s+1;

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
    r_left = double(subs(frleft_psi,psis,ppsi(i)));
    r_right = double(subs(frright_psi,psis,ppsi(i)));
    % the corresponding z
    z_path = subs(fz_rpsi,psis,ppsi(i));
    % the quantites along the 'psi = constant' path
    T_path = subs(fT_rz,z,z_path);
    TDpsi_path = subs(fTDpsi_rz,z,z_path);
    psiGradNorm_path = subs(fpsiGradNorm_rz,z,z_path);
    zDr_path = subs(fzDr_rpsi,psis,ppsi(i));
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
        pr(j,i) = fzero(@(r) chi_r(r)-pchi(j),[r_left,R+(r_right-R)*cos(pi/(4*n_chi))]);
        pz(j,i) = double(subs(fz_rpsi,{psis,r},{ppsi(i),pr(j,i)}));
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
% mesh(pr,pz,pbetachi);
% modified n
alpha = 0.822;
beta = 0.142;
pn = n./(1-alpha*(n*pq/n_chi).^2-beta*(n*pq/n_chi).^4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct of matrixes A and B
% function GrowthRate = GetGrowthrate(n_s,n_chi)
% load('equilibrium.mat');

%% s and chi for the mesh
% n_s = 24;
% n_chi = 25;
% n_s = floor(n_s/2);
% n_chi = ceil(n_chi/2);
ms = 0:1/n_s:1;
mchi = -pi/(2*(n_chi-1)):pi/(n_chi-1):pi+pi/(2*(n_chi-1));
% the totle number of unknow variables
n_x = (3*n_s+1)*(2*n_chi+2);
%total number of cells
n_cell = n_s*n_chi;


%% matrix form of the seven dependent components
i_cell = 1:n_cell;
xij = (2*n_chi+2)*(3*(1:n_s)-3)+2*((1:n_chi)'-1)+1;
i_xij = xij(:);
% XR
X_coef = [1, 1, 1, 1]/4;
XR = sparse(i_cell,i_xij  ,X_coef(1),n_cell,n_x)+...
    sparse(i_cell,i_xij+2,X_coef(2),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)  ,X_coef(3),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+2,X_coef(4),n_cell,n_x);
% XI
XI = sparse(i_cell,i_xij+1,X_coef(1),n_cell,n_x)+...
    sparse(i_cell,i_xij+3,X_coef(2),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+1,X_coef(3),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+3,X_coef(4),n_cell,n_x);
% DXDchiR
mchistep = (mchi(2:end)-mchi(1:end-1))'*ones(n_s,1)';
DXDchi_coef = [-1, 1, -1, 1]'./(2*mchistep(:)');
DXDchiR = sparse(i_cell,i_xij  ,DXDchi_coef(1,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+2,DXDchi_coef(2,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)  ,DXDchi_coef(3,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+2,DXDchi_coef(4,:),n_cell,n_x);
% DXDchiI
DXDchiI = sparse(i_cell,i_xij+1,DXDchi_coef(1,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3,DXDchi_coef(2,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+1,DXDchi_coef(3,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+3,DXDchi_coef(4,:),n_cell,n_x);
% DXDsR
msstep = ones(n_chi,1)*(ms(2:end)-ms(1:end-1));
DXDs_coef = [-1, -1, 1, 1]'./(2*msstep(:)');
DXDsR = sparse(i_cell,i_xij  ,DXDs_coef(1,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+2,DXDs_coef(2,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)  ,DXDs_coef(3,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+2,DXDs_coef(4,:),n_cell,n_x);
% DXDsI
DXDsI = sparse(i_cell,i_xij+1,DXDs_coef(1,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3,DXDs_coef(2,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+1,DXDs_coef(3,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3*(2*n_chi+2)+3,DXDs_coef(4,:),n_cell,n_x);

i_xij = i_xij + (2*n_chi+2);
% VR
V_coef = [1, 1]/2;
VR = sparse(i_cell,i_xij  ,V_coef(1),n_cell,n_x)+...
    sparse(i_cell,i_xij+2,V_coef(2),n_cell,n_x);
% VI
VI = sparse(i_cell,i_xij+1,V_coef(1),n_cell,n_x)+...
    sparse(i_cell,i_xij+3,V_coef(2),n_cell,n_x);
% DVDchiR
DVDchi_coef = [-1, 1]'./mchistep(:)';
DVDchiR = sparse(i_cell,i_xij  ,DVDchi_coef(1,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+2,DVDchi_coef(2,:),n_cell,n_x);
% DVDchiI
DVDchiI = sparse(i_cell,i_xij+1,DVDchi_coef(1,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3,DVDchi_coef(2,:),n_cell,n_x);

i_xij = i_xij + (2*n_chi+2);
% YR
Y_coef = [1, 1]/2;
YR = sparse(i_cell,i_xij  ,Y_coef(1),n_cell,n_x)+...
    sparse(i_cell,i_xij+2,Y_coef(2),n_cell,n_x);
% YI
YI = sparse(i_cell,i_xij+1,Y_coef(1),n_cell,n_x)+...
    sparse(i_cell,i_xij+3,Y_coef(2),n_cell,n_x);
% DYDchiR
DYDchi_coef = [-1, 1]'./mchistep(:)';
DYDchiR = sparse(i_cell,i_xij  ,DYDchi_coef(1,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+2,DYDchi_coef(2,:),n_cell,n_x);
% DYDchiI
DYDchiI = sparse(i_cell,i_xij+1,DYDchi_coef(1,:),n_cell,n_x)+...
    sparse(i_cell,i_xij+3,DYDchi_coef(2,:),n_cell,n_x);

%% the coefficient a b c d e f g h
J = pq.*pr.^2./pT;
acoeff = ones(n_chi,n_s);
% acoeff([1,end],:) = 1/2;
a = 2*pq.^2.*ppsi.*pr.^4./(J.^3.*ppsiGrad2).*acoeff;
b = pT.^2.*pr.^2./(2*Psi_s*J).*acoeff;
c = ppsiGrad2.*pr.^2./(2*Psi_s*J).*acoeff;
d = pr.^4*gamma.*pp./(2*Psi_s*J).*acoeff;
e = 4*pr.^4./J.*(pjphi.^2./ppsiGrad2.*ppsi+...
    pjphi./pr.*plnpsiGradNormDpsi.*ppsi...
    -ppDpsi.*plnrDpsi.*ppsi).*acoeff;
f = 2*prho.*ppsi.*pT.*pr.^2./(pq.*ppsiGrad2).*acoeff;
g = prho.*ppsiGrad2.*pq.*pr.^4./(2*pT*Psi_s).*acoeff;
h = prho.*pr.^4.*pT.*pq/(2*Psi_s).*acoeff;
%% the matrix form of Is
I1R = 1./pq(:).*DXDchiR-pn(:).*XI;
I1I = 1./pq(:).*DXDchiI+pn(:).*XR;
I2R = DXDsR+DVDchiR;
I2I = DXDsI+DVDchiI;
I3R = pH(:).*XR-pbetachi(:).*pn(:).*pq(:).*XI...
    +pn(:).*pq(:).*VI+pbetachi(:).*DXDchiR+DXDsR;
I3I = pH(:).*XI+pbetachi(:).*pn(:).*pq(:).*XR...
    -pn(:).*pq(:).*VR+pbetachi(:).*DXDchiI+DXDsI;
I4R = plnr2Ds(:).*XR+plnr2Dchi(:).*(VR+YR)...
    -pn(:).*pq(:).*YI+DXDsR+DVDchiR+DYDchiR;
I4I = plnr2Ds(:).*XI+plnr2Dchi(:).*(VI+YI)...
    +pn(:).*pq(:).*YR+DXDsI+DVDchiI+DYDchiI;
I5R = XR;
I5I = XI;

%% the matrix A
A =   I1R'*(a(:)./ps(:).*I1R) + I1I'*(a(:)./ps(:).*I1I)...
    + I2R'*(b(:)./ps(:).*I2R) + I2I'*(b(:)./ps(:).*I2I)...
    + I3R'*(c(:)./ps(:).*I3R) + I3I'*(c(:)./ps(:).*I3I)...
    + I4R'*(d(:)./ps(:).*I4R) + I4I'*(d(:)./ps(:).*I4I)...
    - I5R'*(e(:)./ps(:).*I5R) - I5I'*(e(:)./ps(:).*I5I);
%% the matrix B
B =  XR'*(f(:)./ps(:).*XR) + XI'*(f(:)./ps(:).*XI) +...
    (VR+YR-pbetachi(:).*XR)'*(g(:)./ps(:).*(VR+YR-pbetachi(:).*XR))+...
    (VI+YI-pbetachi(:).*XI)'*(g(:)./ps(:).*(VI+YI-pbetachi(:).*XI))+...
    YR'*(h(:)./ps(:).*YR) + YI'*(h(:)./ps(:).*YI);
% subplot(2,4,1);spy(A);
% subplot(2,4,5);spy(B);
% disp([rank(A),rank(B)]);
% A = sparse(A);
% B = sparse(B);
%% symmetry conditions

i_X = (2*n_chi+2)*(3*(1:n_s+1)-3)+1;
i_V = (2*n_chi+2)*(3*(1:n_s)-2)+1;
i_Y = (2*n_chi+2)*(3*(1:n_s)-1)+1;
U = sparse(1:n_x,1:n_x,1,n_x,n_x)+...
    sparse(i_X,i_X+2,1,n_x,n_x)+...%XR
    sparse(i_X+2*n_chi,i_X+2*n_chi-2,1,n_x,n_x)+...
    sparse(i_X+1,i_X+3,-1,n_x,n_x)+...%XI
    sparse(i_X+2*n_chi+1,i_X+2*n_chi-1,-1,n_x,n_x)+...
    sparse(i_V,i_V+2,-1,n_x,n_x)+...%VR
    sparse(i_V+2*n_chi,i_V+2*n_chi-2,-1,n_x,n_x)+...
    sparse(i_V+1,i_V+3,1,n_x,n_x)+...%VI
    sparse(i_V+2*n_chi+1,i_V+2*n_chi-1,1,n_x,n_x)+...
    sparse(i_Y,i_Y+2,-1,n_x,n_x)+...%YR
    sparse(i_Y+2*n_chi,i_Y+2*n_chi-2,-1,n_x,n_x)+...
    sparse(i_Y+1,i_Y+3,1,n_x,n_x)+...%YI
    sparse(i_Y+2*n_chi+1,i_Y+2*n_chi-1,1,n_x,n_x);
%     sparse(i_X+2,i_X,-1,n_x,n_x)+...
%     sparse(i_X+2*n_chi-2,i_X+2*n_chi,-1,n_x,n_x)+...
%     sparse(i_X+3,i_X+1,1,n_x,n_x)+...
%     sparse(i_X+2*n_chi-1,i_X+2*n_chi+1,1,n_x,n_x)+...
%     sparse(i_V+2,i_V,1,n_x,n_x)+...
%     sparse(i_V+2*n_chi-2,i_V+2*n_chi,1,n_x,n_x)+...
%     sparse(i_V+3,i_V+1,-1,n_x,n_x)+...
%     sparse(i_V+2*n_chi-1,i_V+2*n_chi+1,-1,n_x,n_x)+...
%     sparse(i_Y+2,i_Y,1,n_x,n_x)+...
%     sparse(i_Y+2*n_chi-2,i_Y+2*n_chi,1,n_x,n_x)+...
%     sparse(i_Y+3,i_Y+1,-1,n_x,n_x)+...
%     sparse(i_Y+2*n_chi-1,i_Y+2*n_chi+1,-1,n_x,n_x);


% transformation using U
Uinv = inv(U);
At = Uinv'*A*Uinv;
Bt = Uinv'*B*Uinv;

% subplot(2,4,2);spy(At);
% subplot(2,4,6);spy(Bt);
% disp([rank(At),rank(Bt)]);
% force the new tilda variables to zero
w2inv = 1; % which make the condition fullfilled when w2~=10^20
rAold = rank(full(At));
rBold = rank(full(Bt));
figure(3);
clf();
subplot(1,2,1);
xlim([1,n_s+1]);
ylim([1,n_chi+1]);
hold on;
subplot(1,2,2);
xlim([1,n_s+1]);
ylim([1,n_chi+1]);
hold on;
dis=1/10;
% symmetry conditions
for j=1:n_s+1
    for i=[1,n_chi+1]
        i_X = (2*i+2)*(3*(j)-3)+1;
        i_V = (2*i+2)*(3*(j)-2)+1;
        i_Y = (2*i+2)*(3*(j)-1)+1;
        %XR
        i_x=i_X;
        At(i_x,:) = 0; At(:,i_x) = 0;At(i_x,i_x) = 1;
        rAnew = rank(full(At));
        if rAnew~=rAold
            subplot(1,2,1);text(j,i+dis,'XR');rAold = rAnew;
        end
        Bt(i_x,:) = 0;Bt(:,i_x) = 0;Bt(i_x,i_x) = 1;rBnew = rank(full(Bt));
        if rBnew~=rBold
            subplot(1,2,2);text(j,i+dis,'XR');rBold = rBnew;
        end
        %XI
        i_x=i_X+1;
        At(i_x,:) = 0; At(:,i_x) = 0;At(i_x,i_x) = 1;
        rAnew = rank(full(At));
        if rAnew~=rAold
            subplot(1,2,1);text(j,i-dis,'XI');rAold = rAnew;
        end
        Bt(i_x,:) = 0;Bt(:,i_x) = 0;Bt(i_x,i_x) = 1;rBnew = rank(full(Bt));
        if rBnew~=rBold
            subplot(1,2,2);text(j,i-dis,'XI');rBold = rBnew;
        end
        if j==n_s+1
            continue;
        end
        %VR
        i_x=i_V;
        At(i_x,:) = 0; At(:,i_x) = 0;At(i_x,i_x) = 1;
        rAnew = rank(full(At));
        if rAnew~=rAold
            subplot(1,2,1);text(j+1/2-dis,i+dis,'VR');rAold = rAnew;
        end
        Bt(i_x,:) = 0;Bt(:,i_x) = 0;Bt(i_x,i_x) = 1;rBnew = rank(full(Bt));
        if rBnew~=rBold
            subplot(1,2,2);text(j+1/2-dis,i+dis,'VR');rBold = rBnew;
        end
        %VI
        i_x=i_V+1;
        At(i_x,:) = 0; At(:,i_x) = 0;At(i_x,i_x) = 1;
        rAnew = rank(full(At));
        if rAnew~=rAold
            subplot(1,2,1);text(j+1/2-dis,i-dis,'VI');rAold = rAnew;
        end
        Bt(i_x,:) = 0;Bt(:,i_x) = 0;Bt(i_x,i_x) = 1;rBnew = rank(full(Bt));
        if rBnew~=rBold
            subplot(1,2,2);text(j+1/2-dis,i-dis,'VI');rBold = rBnew;
        end
        %YR
        i_x=i_Y;
        At(i_x,:) = 0; At(:,i_x) = 0;At(i_x,i_x) = 1;
        rAnew = rank(full(At));
        if rAnew~=rAold
            subplot(1,2,1);text(j+1/2+dis,i+dis,'YR');rAold = rAnew;
        end
        Bt(i_x,:) = 0;Bt(:,i_x) = 0;Bt(i_x,i_x) = 1;rBnew = rank(full(Bt));
        if rBnew~=rBold
            subplot(1,2,2);text(j+1/2+dis,i+dis,'YR');rBold = rBnew;
        end
        %YI
        i_x=i_Y+1;
        At(i_x,:) = 0; At(:,i_x) = 0;At(i_x,i_x) = 1;
        rAnew = rank(full(At));
        if rAnew~=rAold
            subplot(1,2,1);text(j+1/2+dis,i-dis,'YI');rAold = rAnew;
        end
        Bt(i_x,:) = 0;Bt(:,i_x) = 0;Bt(i_x,i_x) = 1;rBnew = rank(full(Bt));
        if rBnew~=rBold
            subplot(1,2,2);text(j+1/2+dis,i-dis,'YI');rBold = rBnew;
        end
    end
end
% the addition condition to make matrix B not ill conditioned
% i_x = [i_X+2,i_V+3,i_Y+3];
% %tilda A
% At(i_x,:) = 0;
% At(:,i_x) = 0;
% At = At + sparse(i_x,i_x,1,n_x,n_x);
% % tilda B
% Bt(i_x,:) = 0;
% Bt(:,i_x) = 0;
% Bt = Bt + sparse(i_x,i_x,w2inv,n_x,n_x);

% subplot(2,4,3);spy(At);
% subplot(2,4,7);spy(Bt);
% disp([rank(At),rank(Bt)]);
%% regularity conditions at s = 0
% the procejure is same as symmetry conditions
i_x = 2*(1:n_chi+1)-1 ;
i_x = [i_x,i_x+1];
% tilda A
At(i_x,:) = 0;
At(:,i_x) = 0;
At = At + sparse(i_x,i_x,1,n_x,n_x);
% tilda B
Bt(i_x,:) = 0;
Bt(:,i_x) = 0;
Bt = Bt + sparse(i_x,i_x,w2inv,n_x,n_x);

%% boundary conditions at plasma surface s = 1
% a wall straight on the plasma surface
%the procejure is same as symmetry conditions
i_x =3*n_s*(2*n_chi+2)+2*(1:n_chi+1)-1 ;
i_x = [i_x,i_x+1];
% tilda A
At(i_x,:) = 0;
At(:,i_x) = 0;
At = At + sparse(i_x,i_x,1,n_x,n_x);
% tilda B
Bt(i_x,:) = 0;
Bt(:,i_x) = 0;
Bt = Bt + sparse(i_x,i_x,w2inv,n_x,n_x);

% subplot(2,4,4);spy(At);
% subplot(2,4,8);spy(Bt);
% disp([rank(At),rank(Bt)]);
w2 = eig(full(Bt\At));
disp(min(w2));
% mu0 = -8;
% [x,w2] = eigSolver(At,Bt,1e-8,mu0,rand(n_x,1));
% w2 = w2+mu0;
% disp(w2);
% figure(1);
% plot(w2);
GrowthRate =w2;


