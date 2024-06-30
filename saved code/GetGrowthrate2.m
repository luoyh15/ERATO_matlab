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
acoeff([1,end],:) = 1/2;
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
I1R = 1./pq(:).*DXDchiR-n*XI;
I1I = 1./pq(:).*DXDchiI+n*XR;
I2R = DXDsR+DVDchiR;
I2I = DXDsI+DVDchiI;
I3R = pH(:).*XR-pbetachi(:).*n.*pq(:).*XI...
    +n*pq(:).*VI+pbetachi(:).*DXDchiR+DXDsR;
I3I = pH(:).*XI+pbetachi(:).*n.*pq(:).*XR...
    -n*pq(:).*VR+pbetachi(:).*DXDchiI+DXDsI;
I4R = plnr2Ds(:).*XR+plnr2Dchi(:).*(VR+YR)...
    -n*pq(:).*YI+DXDsR+DVDchiR+DYDchiR;
I4I = plnr2Ds(:).*XI+plnr2Dchi(:).*(VI+YI)...
    +n*pq(:).*YR+DXDsI+DVDchiI+DYDchiI;
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
    sparse(i_Y+1,i_Y+3,-1,n_x,n_x)+...%YI
    sparse(i_Y+2*n_chi+1,i_Y+2*n_chi-1,-1,n_x,n_x);

% transformation using U
Uinv = inv(U);
At = Uinv'*A*Uinv;
Bt = Uinv'*B*Uinv;

% subplot(2,4,2);spy(At);
% subplot(2,4,6);spy(Bt);
% disp([rank(At),rank(Bt)]);
% force the new tilda variables to zero
w2inv = 1e-20; % which make the condition fullfilled when w2~=10^20

i_X = (2*n_chi+2)*(3*(1:n_s+1)-3)+1;
i_V = (2*n_chi+2)*(3*(1:n_s)-2)+1;
i_Y = (2*n_chi+2)*(3*(1:n_s)-1)+1;
i_x = [i_X,i_V,i_Y];
i_x = [i_x,i_x+1,i_x+2*n_chi,i_x+2*n_chi+1];
% symmetry conditions
%tilda A
At(i_x,:) = 0;
At(:,i_x) = 0;
At = At + sparse(i_x,i_x,1,n_x,n_x);
% tilda B
Bt(i_x,:) = 0;
Bt(:,i_x) = 0;
Bt = Bt + sparse(i_x,i_x,w2inv,n_x,n_x);

% the addition condition to make matrix not ill conditioned
i_x = [i_X+2,i_V+3,i_Y+3];
%tilda A
At(i_x,:) = 0;
At(:,i_x) = 0;
At = At + sparse(i_x,i_x,1,n_x,n_x);
% tilda B
Bt(i_x,:) = 0;
Bt(:,i_x) = 0;
Bt = Bt + sparse(i_x,i_x,w2inv,n_x,n_x);

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
