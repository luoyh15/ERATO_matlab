%% s and chi for the mesh
n_s = 4;
n_chi = 5;
ms = 0:1/n_s:1;
mchi = -pi/(2*(n_chi-1)):pi/(n_chi-1):pi+pi/(2*(n_chi-1));
%% initialize matrix A
n_x = (3*n_s+1)*(2*n_chi+2);% the totle number of unknow variables
%A = zeros(n_x,n_x);
%% matrix form of the seven dependent components
%initialize
n_cell = n_s*(n_chi);
XR = zeros(n_cell,n_x);
XI = zeros(n_cell,n_x);
DXDchiR = zeros(n_cell,n_x);
DXDchiI = zeros(n_cell,n_x);
DXDsR = zeros(n_cell,n_x);
DXDsI = zeros(n_cell,n_x);
YR = zeros(n_cell,n_x);
YI = zeros(n_cell,n_x);
DYDchiR = zeros(n_cell,n_x);
DYDchiI = zeros(n_cell,n_x);
VR = zeros(n_cell,n_x);
VI = zeros(n_cell,n_x);
DVDchiR = zeros(n_cell,n_x);
DVDchiI = zeros(n_cell,n_x);
% the real part of the matrixes
for i = 1:n_s
    for j = 1:n_chi
        % the coefficients of the discretization form
        XR_coef = [1, 1, 1, 1]/4;
        XR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*(j-1)+1) = XR_coef(1);
        XR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*j+1) = XR_coef(2);
        XR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*(j-1)+1) = XR_coef(3);
        XR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*j+1) = XR_coef(4);
       
        DXDchiR_coef = [-1, 1, -1, 1]/(2*(mchi(j+1)-mchi(j)));
        DXDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*(j-1)+1) = DXDchiR_coef(1);
        DXDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*j+1) = DXDchiR_coef(2);
        DXDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*(j-1)+1) = DXDchiR_coef(3);
        DXDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*j+1) = DXDchiR_coef(4);
        
        DXDsR_coef = [-1, -1, 1, 1]/(2*(ms(i+1)-ms(i)));
        DXDsR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*(j-1)+1) = DXDsR_coef(1);
        DXDsR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*j+1) = DXDsR_coef(2);
        DXDsR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*(j-1)+1) = DXDsR_coef(3);
        DXDsR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*j+1) = DXDsR_coef(4);
        
        YR_coef = [1, 1]/2;
        YR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-2)+2*(j-1)+1) = YR_coef(1);
        YR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-2)+2*j+1) = YR_coef(2);
        
        DYDchiR_coef = [-1, 1]/(mchi(j+1)-mchi(j));
        DYDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-2)+2*(j-1)+1) = DYDchiR_coef(1);
        DYDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-2)+2*j+1) = DYDchiR_coef(2);
        
        VR_coef = [1, 1]/2;
        VR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-1)+2*(j-1)+1) = VR_coef(1);
        VR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-1)+2*j+1) = VR_coef(2);
        
        DVDchiR_coef = [-1, 1]/(mchi(j+1)-mchi(j));
        DVDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-1)+2*(j-1)+1) = DVDchiR_coef(1);
        DVDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-1)+2*j+1) = DVDchiR_coef(2);
       
    end
end
% the image part of the seven dependent variables
XI(:,2:2:end) = XR(:,1:2:end);
DXDchiI(:,2:2:end) = DXDchiR(:,1:2:end);
DXDsI(:,2:2:end) = DXDsR(:,1:2:end);
YI(:,2:2:end) = YR(:,1:2:end);
DYDchiI(:,2:2:end) = DYDchiR(:,1:2:end);
VI(:,2:2:end) = VR(:,1:2:end);
DVDchiI(:,2:2:end) = DVDchiR(:,1:2:end);

%% the normalized equilibrium quantities in the middle of each cell
% s and chi the physical components
ps = 1/(2*n_s):1/n_s:1-1/(2*n_s);
pchi = 0:pi/(n_chi-1):pi;
[ps,pchi] = meshgrid(ps,pchi);
% psi
ppsi = ps.^2*Psi_s;
% r z
pr = interp2(cs,cchi,r_sc,ps,pchi);
pz = interp2(cs,cchi,z_sc,ps,pchi);
% the norm of gradient of psi
ppsiGradNorm = interp2(cr,cz,psiGradNorm_rz,pr,pz);
% mass density
prho = 1;
% plasma pressure 
pp = interp1(s,p_s,ps);
% psi derivative of p
ppDpsi = interp1(s,pDpsi_s,ps);
% flux function
pT = interp1(s,T_s,ps);
% s derivative of T
pTDs = interp1(s,TDs_s,ps);
% toroidal current density
pjphi = interp2(cr,cz,jphi_rz,pr,pz);
% safety factor
pq = interp1(s,q_s,ps);
% s derivative of q
pqDs = interp1(s,qDs_s,ps);
% poloidal magnetic field
pBp = ppsi.*pr.^2./(q0*ppsiGradNorm.^2);
% non-orthogonality
pbetachi = interp2(cs,cchi,betachi_sc,ps,pchi);
% dlog(r^2)/ds and dlog(r^2)/dchi
plogr2Ds = interp2(cs,cchi,logr2Ds_sc,ps,pchi);
plogr2Dchi = interp2(cs,cchi,logr2Dchi_sc,ps,pchi);
% H defined at the ERATO paper
pH = 2*pjphi.*ppsi.*pr./(ps.*ppsiGradNorm.^2)+pTDs./pT-pqDs./pq;

% the gradient psi direction's derivative of ln psi gradient norm
plogpsiGradNormDpsi = interp2(cr,cz,logpsiGradNormDpsi_rz,pr,pz);
% the gradient psi direction's derivative of log r
plogrDpsi = interp2(cr,cz,logrDpsi_rz,pr,pz);
% K defined at the ERATO paper
pK = 2*ppsi/q0.*(pjphi.^2./ppsiGradNorm.^2+pjphi.^2./pr.*plogpsiGradNormDpsi...
    -ppDpsi.*plogrDpsi);


%% the matrix form of Is
I1R = 1./pq(:).*DXDchiR-n*XI;
I1I = 1./pq(:).*DXDchiI+n*XR;
I2R = DXDsR+DVDchiR;
I2I = DXDsI+DVDchiI;
I3R = pH(:).*XR-pbetachi(:).*n.*pq(:).*XI...
    -n*pq(:).*VI+pbetachi(:).*DXDchiR+DXDsR;
I3I = pH(:).*XI+pbetachi(:).*n.*pq(:).*XR...
    +n*pq(:).*VR+pbetachi(:).*DXDchiI+DXDsI;
I4R = plogr2Ds(:).*XR+plogr2Dchi(:).*(VR+pq(:).*YR)...
    -n*pq(:).^2.*YI+DXDsR+DVDchiR+pq(:).*DYDchiR;
I4I = plogr2Ds(:).*XI+plogr2Dchi(:).*(VI+pq(:).*YI)...
    +n*pq(:).^2.*YR+DXDsI+DVDchiI+pq(:).*DYDchiI;
I5R = XR;
I5I = XI;
%% the coefficient a b c d e f g h
J = pq.*pr.^2./pT;
a = 2*pq.^2.*ppsi.*pr.^4./(J.^3.*ppsiGradNorm.^2);
b = pT.^2.*pr.^2./(2*Psi_s*J);
c = ppsiGradNorm.^2.*pr.^2./(2*Psi_s*J);
d = pr.^4*gamma.*pp./(2*Psi_s*J);
e = pK.*2.*pr.^4*q0./J;
f = 2*prho.*ppsi.*pT.*pr.^2./(pq.*ppsiGradNorm.^2);
g = prho.*ppsiGradNorm.^2.*pq.*pr.^4./(2*pT*Psi_s);
h = prho.*pr.^4.*pT.*pq/(2*Psi_s);


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
subplot(2,4,1);spy(A);
subplot(2,4,5);spy(B);
disp([rank(A),rank(B)]);
% A = sparse(A);
% B = sparse(B);
A = A+diag(diag(A));
B = B+diag(diag(B));
%% symmetry conditions
U = diag(ones(n_x,1));
for i = 1:n_s+1
    % symmetry conditions of XR
    U((2*n_chi+2)*(3*i-3)+1,(2*n_chi+2)*(3*i-3)+1) = 1; 
    U((2*n_chi+2)*(3*i-3)+1,(2*n_chi+2)*(3*i-3)+3) = 1;
    U((2*n_chi+2)*(3*i-3)+2*n_chi+1,(2*n_chi+2)*(3*i-3)+2*n_chi-1) = 1; 
    U((2*n_chi+2)*(3*i-3)+2*n_chi+1,(2*n_chi+2)*(3*i-3)+2*n_chi+1) = 1;
    % symmetry conditions of XI
    U((2*n_chi+2)*(3*i-3)+2,(2*n_chi+2)*(3*i-3)+2) = -1; 
    U((2*n_chi+2)*(3*i-3)+2,(2*n_chi+2)*(3*i-3)+4) = 1;
    U((2*n_chi+2)*(3*i-3)+2*n_chi+2,(2*n_chi+2)*(3*i-3)+2*n_chi) = -1; 
    U((2*n_chi+2)*(3*i-3)+2*n_chi+2,(2*n_chi+2)*(3*i-3)+2*n_chi+2) = 1;
end
for i = 1:n_s
    % symmetry conditions of VR
    U((2*n_chi+2)*(3*i-2)+1,(2*n_chi+2)*(3*i-2)+1) = -1; 
    U((2*n_chi+2)*(3*i-2)+1,(2*n_chi+2)*(3*i-2)+3) = 1;
    U((2*n_chi+2)*(3*i-2)+2*n_chi+1,(2*n_chi+2)*(3*i-2)+2*n_chi-1) = -1; 
    U((2*n_chi+2)*(3*i-2)+2*n_chi+1,(2*n_chi+2)*(3*i-2)+2*n_chi+1) = 1;
    % symmetry conditions of VI
    U((2*n_chi+2)*(3*i-2)+2,(2*n_chi+2)*(3*i-2)+2) = 1; 
    U((2*n_chi+2)*(3*i-2)+2,(2*n_chi+2)*(3*i-2)+4) = 1;
    U((2*n_chi+2)*(3*i-2)+2*n_chi+2,(2*n_chi+2)*(3*i-2)+2*n_chi) = 1; 
    U((2*n_chi+2)*(3*i-2)+2*n_chi+2,(2*n_chi+2)*(3*i-2)+2*n_chi+2) = 1;
    % symmetry conditions of YR
    U((2*n_chi+2)*(3*i-1)+1,(2*n_chi+2)*(3*i-1)+1) = -1; 
    U((2*n_chi+2)*(3*i-1)+1,(2*n_chi+2)*(3*i-1)+3) = 1;
    U((2*n_chi+2)*(3*i-1)+2*n_chi+1,(2*n_chi+2)*(3*i-1)+2*n_chi-1) = -1; 
    U((2*n_chi+2)*(3*i-1)+2*n_chi+1,(2*n_chi+2)*(3*i-1)+2*n_chi+1) = 1;
    % symmetry conditions of YI
    U((2*n_chi+2)*(3*i-1)+2,(2*n_chi+2)*(3*i-1)+2) = 1; 
    U((2*n_chi+2)*(3*i-1)+2,(2*n_chi+2)*(3*i-1)+4) = 1;
    U((2*n_chi+2)*(3*i-1)+2*n_chi+2,(2*n_chi+2)*(3*i-1)+2*n_chi) = 1; 
    U((2*n_chi+2)*(3*i-1)+2*n_chi+2,(2*n_chi+2)*(3*i-1)+2*n_chi+2) = 1;
end
% transformation using U
Uinv = inv(U);
At = Uinv'*A*Uinv;
Bt = Uinv'*B*Uinv;

subplot(2,4,2);spy(At);
subplot(2,4,6);spy(Bt);
disp([rank(At),rank(Bt)]);
% force the new tilda variables to zero
w2inv = 1; % which make the condition fullfilled when w2~=10^20
for i = 1:n_s+1
    %tilda A
    % symmetry conditions of XR
    At((2*n_chi+2)*(3*i-3)+1,:) = 0;
    At(:,(2*n_chi+2)*(3*i-3)+1) = 0;
    At((2*n_chi+2)*(3*i-3)+1,(2*n_chi+2)*(3*i-3)+1) = 1;
    At((2*n_chi+2)*(3*i-3)+2*n_chi+1,:) = 0;
    At(:,(2*n_chi+2)*(3*i-3)+2*n_chi+1) = 0;
    At((2*n_chi+2)*(3*i-3)+2*n_chi+1,(2*n_chi+2)*(3*i-3)+2*n_chi+1) = 1;
    % symmetry conditions of XI
    At((2*n_chi+2)*(3*i-3)+2,:) = 0;
    At(:,(2*n_chi+2)*(3*i-3)+2) = 0;
    At((2*n_chi+2)*(3*i-3)+2,(2*n_chi+2)*(3*i-3)+2) = 1;
    At((2*n_chi+2)*(3*i-3)+2*n_chi+2,:) = 0;
    At(:,(2*n_chi+2)*(3*i-3)+2*n_chi+2) = 0;
    At((2*n_chi+2)*(3*i-3)+2*n_chi+2,(2*n_chi+2)*(3*i-3)+2*n_chi+2) = 1;
    
    % tilda B
    % symmetry conditions of XR
    Bt((2*n_chi+2)*(3*i-3)+1,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-3)+1) = 0;
    Bt((2*n_chi+2)*(3*i-3)+1,(2*n_chi+2)*(3*i-3)+1) = w2inv;
    Bt((2*n_chi+2)*(3*i-3)+2*n_chi+1,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-3)+2*n_chi+1) = 0;
    Bt((2*n_chi+2)*(3*i-3)+2*n_chi+1,(2*n_chi+2)*(3*i-3)+2*n_chi+1) = w2inv;
    % symmetry conditions of XI
    Bt((2*n_chi+2)*(3*i-3)+2,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-3)+2) = 0;
    Bt((2*n_chi+2)*(3*i-3)+2,(2*n_chi+2)*(3*i-3)+2) = w2inv;
    Bt((2*n_chi+2)*(3*i-3)+2*n_chi+2,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-3)+2*n_chi+2) = 0;
    Bt((2*n_chi+2)*(3*i-3)+2*n_chi+2,(2*n_chi+2)*(3*i-3)+2*n_chi+2) = w2inv;
end


for i = 1:n_s
    % symmetry conditions of VR
    At((2*n_chi+2)*(3*i-2)+1,:) = 0;
    At(:,(2*n_chi+2)*(3*i-2)+1) = 0;
    At((2*n_chi+2)*(3*i-2)+1,(2*n_chi+2)*(3*i-2)+1) = 1;
    At((2*n_chi+2)*(3*i-2)+2*n_chi+1,:) = 0;
    At(:,(2*n_chi+2)*(3*i-2)+2*n_chi+1) = 0;
    At((2*n_chi+2)*(3*i-2)+2*n_chi+1,(2*n_chi+2)*(3*i-2)+2*n_chi+1) = 1;
    % symmetry conditions of VI
    At((2*n_chi+2)*(3*i-2)+2,:) = 0;
    At(:,(2*n_chi+2)*(3*i-2)+2) = 0;
    At((2*n_chi+2)*(3*i-2)+2,(2*n_chi+2)*(3*i-2)+2) = 1;
    At((2*n_chi+2)*(3*i-2)+2*n_chi+2,:) = 0;
    At(:,(2*n_chi+2)*(3*i-2)+2*n_chi+2) = 0;
    At((2*n_chi+2)*(3*i-2)+2*n_chi+2,(2*n_chi+2)*(3*i-2)+2*n_chi+2) = 1;
    % symmetry conditions of YR
    At((2*n_chi+2)*(3*i-1)+1,:) = 0;
    At(:,(2*n_chi+2)*(3*i-1)+1) = 0;
    At((2*n_chi+2)*(3*i-1)+1,(2*n_chi+2)*(3*i-1)+1) = 1;
    At((2*n_chi+2)*(3*i-1)+2*n_chi+1,:) = 0;
    At(:,(2*n_chi+2)*(3*i-1)+2*n_chi+1) = 0;
    At((2*n_chi+2)*(3*i-1)+2*n_chi+1,(2*n_chi+2)*(3*i-1)+2*n_chi+1) = 1;
    % symmetry conditions of YI
    At((2*n_chi+2)*(3*i-1)+2,:) = 0;
    At(:,(2*n_chi+2)*(3*i-1)+2) = 0;
    At((2*n_chi+2)*(3*i-1)+2,(2*n_chi+2)*(3*i-1)+2) = 1;
    At((2*n_chi+2)*(3*i-1)+2*n_chi+2,:) = 0;
    At(:,(2*n_chi+2)*(3*i-1)+2*n_chi+2) = 0;
    At((2*n_chi+2)*(3*i-1)+2*n_chi+2,(2*n_chi+2)*(3*i-1)+2*n_chi+2) = 1;
    % tilda B
    % symmetry conditions of VR
    Bt((2*n_chi+2)*(3*i-2)+1,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-2)+1) = 0;
    Bt((2*n_chi+2)*(3*i-2)+1,(2*n_chi+2)*(3*i-2)+1) = w2inv;
    Bt((2*n_chi+2)*(3*i-2)+2*n_chi+1,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-2)+2*n_chi+1) = 0;
    Bt((2*n_chi+2)*(3*i-2)+2*n_chi+1,(2*n_chi+2)*(3*i-2)+2*n_chi+1) = w2inv;
    % symmetry conditions of VI
    Bt((2*n_chi+2)*(3*i-2)+2,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-2)+2) = 0;
    Bt((2*n_chi+2)*(3*i-2)+2,(2*n_chi+2)*(3*i-2)+2) = w2inv;
    Bt((2*n_chi+2)*(3*i-2)+2*n_chi+2,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-2)+2*n_chi+2) = 0;
    Bt((2*n_chi+2)*(3*i-2)+2*n_chi+2,(2*n_chi+2)*(3*i-2)+2*n_chi+2) = w2inv;
    % symmetry conditions of YR
    Bt((2*n_chi+2)*(3*i-1)+1,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-1)+1) = 0;
    Bt((2*n_chi+2)*(3*i-1)+1,(2*n_chi+2)*(3*i-1)+1) = w2inv;
    Bt((2*n_chi+2)*(3*i-1)+2*n_chi+1,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-1)+2*n_chi+1) = 0;
    Bt((2*n_chi+2)*(3*i-1)+2*n_chi+1,(2*n_chi+2)*(3*i-1)+2*n_chi+1) = w2inv;
    % symmetry conditions of YI
    Bt((2*n_chi+2)*(3*i-1)+2,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-1)+2) = 0;
    Bt((2*n_chi+2)*(3*i-1)+2,(2*n_chi+2)*(3*i-1)+2) = w2inv;
    Bt((2*n_chi+2)*(3*i-1)+2*n_chi+2,:) = 0;
    Bt(:,(2*n_chi+2)*(3*i-1)+2*n_chi+2) = 0;
    Bt((2*n_chi+2)*(3*i-1)+2*n_chi+2,(2*n_chi+2)*(3*i-1)+2*n_chi+2) = w2inv;
end

subplot(2,4,3);spy(At);
subplot(2,4,7);spy(Bt);
disp([rank(At),rank(Bt)]);
%% regularity conditions at s = 0
% the procejure is same as symmetry conditions
for j = 1:n_chi+1
    % tilda A
    % XR
    At(2*(j-1)+1,:) = 0;
    At(:,2*(j-1)+1) = 0;
    At(2*(j-1)+1,2*(j-1)+1) = 1;
    % XI
    At(2*(j-1)+2,:) = 0;
    At(:,2*(j-1)+2) = 0;
    At(2*(j-1)+2,2*(j-1)+2) = 1;
    % tilda B
    % XR
    Bt(2*(j-1)+1,:) = 0;
    Bt(:,2*(j-1)+1) = 0;
    Bt(2*(j-1)+1,2*(j-1)+1) = w2inv;
    % XI
    Bt(2*(j-1)+2,:) = 0;
    Bt(:,2*(j-1)+2) = 0;
    Bt(2*(j-1)+2,2*(j-1)+2) = w2inv;
end

%% boundary conditions at plasma surface s = 1
% a wall straight on the plasma surface
%the procejure is same as symmetry conditions
for j = 1:n_chi+1
    % tilda A
    % XR
    At(3*n_s*(2*n_chi+2)+2*(j-1)+1,:) = 0;
    At(:,3*n_s*(2*n_chi+2)+2*(j-1)+1) = 0;
    At(3*n_s*(2*n_chi+2)+2*(j-1)+1,3*n_s*(2*n_chi+2)+2*(j-1)+1) = 1;
    % XI
    At(3*n_s*(2*n_chi+2)+2*(j-1)+2,:) = 0;
    At(:,3*n_s*(2*n_chi+2)+2*(j-1)+2) = 0;
    At(3*n_s*(2*n_chi+2)+2*(j-1)+2,3*n_s*(2*n_chi+2)+2*(j-1)+2) = 1;
    % tilda B
    % XR
    Bt(3*n_s*(2*n_chi+2)+2*(j-1)+1,:) = 0;
    Bt(:,3*n_s*(2*n_chi+2)+2*(j-1)+1) = 0;
    Bt(3*n_s*(2*n_chi+2)+2*(j-1)+1,3*n_s*(2*n_chi+2)+2*(j-1)+1) = w2inv;
    % XI
    Bt(3*n_s*(2*n_chi+2)+2*(j-1)+2,:) = 0;
    Bt(:,3*n_s*(2*n_chi+2)+2*(j-1)+2) = 0;
    Bt(3*n_s*(2*n_chi+2)+2*(j-1)+2,3*n_s*(2*n_chi+2)+2*(j-1)+2) = w2inv;
end
subplot(2,4,4);spy(At);
subplot(2,4,8);spy(Bt);
disp([rank(At),rank(Bt)]);
% At = sparse(At);
% Bt = sparse(Bt);
w2 = eig(Bt\At);
