function gamma2 = GrowthRateOfModeN_save(quantities,n)

%% equilibrium state
gamma = 5/3;
Psi_s = quantities.psi_p;
%the coordinate
ps = quantities.ps;
pchi = quantities.pchi;
% psi
ppsi = ps.^2*Psi_s;
%r z
pr = quantities.pr;
pz = quantities.pz;
% the non-orthogonality betachi
pbetachi = quantities.pbetachi;
% lnr2Dchi
plnr2Dchi = quantities.plnr2Dchi;
% lnr2Ds
plnr2Ds = quantities.plnr2Ds;
% toroidal current density
pjphi = quantities.pjphi;
% pressure
pp = quantities.pp;
% psi derivative of p
ppDpsi = quantities.ppDpsi;
% flux function
pT = quantities.pT;
% s derivative of T
% pTDs = quantities.pTDs;
%safety factor
pq = quantities.pq;
%psi derivative of q
% pqDs = quantities.pqDs;
% psiGrad2
ppsiGrad2 = quantities.ppsiGrad2;
% H defined at the ERATO paper
pH = quantities.pH;
% lnpsiGradNormDpsi
plnpsiGradNormDpsi = quantities.plnpsiGradNormDpsi;
% lnrDpsi
plnrDpsi = quantities.plnrDpsi;
% mass density
prho = 1;
% mesh(pr,pz,pbetachi);
% pn = n./(1-alpha*(n*pq/n_chi).^2-beta*(n*pq/n_chi).^4);
% pn = n;

%% s and chi for the mesh
% n_s = 24;
% n_chi = 25;
% n_s = floor(n_s/2);
% n_chi = ceil(n_chi/2);
% ms = 0:1/n_s:1;
ms = quantities.ms;
mchi = quantities.mchi;
[n_chi,n_s] = size(ps);
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
mchistep = (mchi(2:end)-mchi(1:end-1))*ones(n_s,1)';
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
msstep = ones(n_chi,1)*(ms(2:end)-ms(1:end-1))';
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
a = 2*pq.^2.*ppsi.*pr.^4./(J.^3.*ppsiGrad2);
b = pT.^2.*pr.^2./(2*Psi_s*J);
c = ppsiGrad2.*pr.^2./(2*Psi_s*J);
d = pr.^4*gamma.*pp./(2*Psi_s*J);
e = 4*pr.^4.*ppsi./J.*(pjphi.^2./ppsiGrad2+...
    pjphi./pr.*plnpsiGradNormDpsi-ppDpsi.*plnrDpsi);
f = 2*prho.*ppsi.*pT.*pr.^2./(pq.*ppsiGrad2);
g = prho.*ppsiGrad2.*pq.*pr.^4./(2*pT*Psi_s);
h = prho.*pr.^4.*pT.*pq/(2*Psi_s);
%% the matrix form of Is
pn = n;
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
area = ones(n_chi,n_s).*mchistep.*msstep;
area([1,end],:) = 1/2*area([1,end],:);
A =   I1R'*(a(:)./ps(:).*area(:).*I1R) + I1I'*(a(:)./ps(:).*area(:).*I1I)...
    + I2R'*(b(:)./ps(:).*area(:).*I2R) + I2I'*(b(:)./ps(:).*area(:).*I2I)...
    + I3R'*(c(:)./ps(:).*area(:).*I3R) + I3I'*(c(:)./ps(:).*area(:).*I3I)...
    + I4R'*(d(:)./ps(:).*area(:).*I4R) + I4I'*(d(:)./ps(:).*area(:).*I4I)...
    - I5R'*(e(:)./ps(:).*area(:).*I5R) - I5I'*(e(:)./ps(:).*area(:).*I5I);
%% the matrix B
B =  XR'*(f(:)./ps(:).*area(:).*XR) + XI'*(f(:)./ps(:).*area(:).*XI) +...
    (VR+YR-pbetachi(:).*XR)'*(g(:)./ps(:).*area(:).*(VR+YR-pbetachi(:).*XR))+...
    (VI+YI-pbetachi(:).*XI)'*(g(:)./ps(:).*area(:).*(VI+YI-pbetachi(:).*XI))+...
    YR'*(h(:)./ps(:).*area(:).*YR) + YI'*(h(:)./ps(:).*area(:).*YI);
% subplot(2,4,1);spy(A);
% subplot(2,4,5);spy(B);
% disp([rank(A),rank(B)]);
% A = sparse(A);
% B = sparse(B);
%% symmetry conditions

i_X = (2*n_chi+2)*(3*(1:n_s+1)-3)+1;
i_V = (2*n_chi+2)*(3*(1:n_s)-2)+1;
i_Y = (2*n_chi+2)*(3*(1:n_s)-1)+1;
x_s = [i_X+1,i_V,i_Y];
x_as = [i_X,i_V+1,i_Y+1];
U = sparse(1:n_x,1:n_x,1,n_x,n_x)+...
    sparse(x_s,x_s+2,-1,n_x,n_x)+...
    sparse(x_s+2*n_chi,x_s+2*n_chi-2,-1,n_x,n_x)+...
    sparse(x_as,x_as+2,1,n_x,n_x)+...
    sparse(x_as+2*n_chi,x_as+2*n_chi-2,1,n_x,n_x);
% transformation using U
Uinv = inv(U);
At = Uinv'*A*Uinv;
Bt = Uinv'*B*Uinv;

% subplot(2,4,2);spy(At);
% subplot(2,4,6);spy(Bt);
% disp([rank(At),rank(Bt)]);
% force the new tilda variables to zero
w2inv = 1e-20; % which make the condition fullfilled when w2~=10^20

% symmetry conditions
i_X = (2*n_chi+2)*(3*(1:n_s+1)-3)+1;
i_V = (2*n_chi+2)*(3*(1:n_s)-2)+1;
i_Y = (2*n_chi+2)*(3*(1:n_s)-1)+1;
i_x = [i_X,i_V,i_Y];
i_x = [i_x,i_x+1,i_x+2*n_chi,i_x+2*n_chi+1];
%tilda A
At(i_x,:) = 0;
At(:,i_x) = 0;
At = At + sparse(i_x,i_x,1,n_x,n_x);
% tilda B
Bt(i_x,:) = 0;
Bt(:,i_x) = 0;
Bt = Bt + sparse(i_x,i_x,w2inv,n_x,n_x);

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
% eigVal = eig(full(At),full(Bt));
[~,p] = chol(At);
if p==0
    gamma2 = 0;
else   
        
    %     [eigFun,eigVal] = eig(full(At),full(Bt));
    %     w2 = diag(eigVal);
    %     [gamma2,Index] = min(w2);
    %     disp(gamma2);
    gamma2 = 0;
    w2guess = -1;
    while w2guess>-10^4    
        [eigFun,eigVal] = eigs(At-w2guess*Bt,Bt,1,'smallestabs');
        eigVal = eigVal+w2guess;
        if eigVal<0
            gamma2 = eigVal;
            break;
        end
        w2guess = w2guess*10;
    end
    % mu0 = -1.5;
    % [Fx,w2] = eigSolver(At,Bt,1e-8,mu0,rand(n_x,1));
    % w2 = w2+mu0;
    % disp(w2);
    Fx = eigFun;
    Fx = Uinv*Fx;
    Fx = Fx./(norm(Fx));
    %% plot the eigenFunction of the most ustable modes
    pXR = reshape(XR*Fx,[n_chi,n_s]);
    pXI = reshape(XI*Fx,[n_chi,n_s]);
    pVR = reshape(VR*Fx,[n_chi,n_s]);
    pVI = reshape(VI*Fx,[n_chi,n_s]);
    pYR = reshape(YR*Fx,[n_chi,n_s]);
    pYI = reshape(YI*Fx,[n_chi,n_s]);
    pX = complex(pXR,pXI);
    pV = complex(pVR,pVI);
    pY = complex(pYR,pYI);
    ppsiGradNorm = sqrt(ppsiGrad2);
    pksi1 = pT./(pq.*ppsiGradNorm).*pX;
    pksi2 = pr.*ppsiGradNorm./(2*ps*Psi_s).*(pV+pY-pbetachi.*pX);
    pksi3 = pr.*pT./(2*ps*Psi_s).*pY;
    % the compontant of ksi in r,z direction
    ppsiDr = quantities.ppsiDr;
    ppsiDz = quantities.ppsiDz;
    pksir = (pksi1.*ppsiDr-pksi2.*ppsiDz)./ppsiGradNorm;
    pksiz = (pksi1.*ppsiDz+pksi2.*ppsiDr)./ppsiGradNorm;
    %     quiver3(pr,zeros(size(pr)),pz,pksir,pksi3,pksiz);
    quiver(pr,pz,pksir,pksiz);
end
end



