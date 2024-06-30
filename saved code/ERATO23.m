function GrowthRate = ERATO23(E,a,q0,n,n_s)


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
% pn = n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct of matrixes A and B
% function GrowthRate = GetGrowthrate(n_s,n_chi)
% load('equilibrium.mat');
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

At = zeros(n_x,n_x);
Bt = zeros(n_x,n_x);
for j = 1:n_s
    for i = 1:n_chi
        % basis functions
        Pf = zeros(16,1);
        Pf(2) = 1/4; Pf(5) = 1/4; Pf(12) = 1/4; Pf(15) = 1/4;
        Pf(8) = 1/2; Pf(10) = 1/2; 
        dchi = mchi(i+1)-mchi(i);
        Pf(1) = -1/(2*dchi); Pf(4) = -Pf(1); Pf(11) = Pf(1); Pf(14) = Pf(4);
        Pf(7) = -1/dchi; Pf(9) = -Pf(7);
        ds = ms(j+1)-ms(j);
        Pf(3) = -1/(2*ds); Pf(6) = Pf(3); Pf(13) = -Pf(3); Pf(16) = Pf(13);
        % the coefficient a b c d e f g h
        Pc1 = zeros(16,1);
        Pc1(1) = a(i,j)/ps(i,j);Pc1(2) = Pc1(1);
        Pc1(3) = b(i,j)/ps(i,j);Pc1(4) = Pc1(3);
        Pc1(5) = c(i,j)/ps(i,j);Pc1(6) = Pc1(5);
        Pc1(7) = d(i,j)/ps(i,j);Pc1(8) = Pc1(7);
        Pc1(9) = -e(i,j)/ps(i,j);Pc1(10) = Pc1(9);
        Pc1(11) = f(i,j)/ps(i,j);Pc1(12) = Pc1(11);
        Pc1(13) = g(i,j)/ps(i,j);Pc1(14) = Pc1(13);
        Pc1(15) = h(i,j)/ps(i,j);Pc1(16) = Pc1(15);
        % quadratic terms 
        Pc = zeros(14,16);
        Pc(1,1) = 1/pq(i,j);Pc(4,1) = -pn(i,j);
        Pc(2,2) = 1/pq(i,j);Pc(3,2) = pn(i,j);
        
        Pc(5,3) = 1; Pc(7,3) = 1;
        Pc(6,4) = 1; Pc(8,4) = 1;
        
        Pc(1,5) = pbetachi(i,j); Pc(3,5) = pH(i,j);
        Pc(4,5) = -pn(i,j)*pq(i,j)*pbetachi(i,j); 
        Pc(5,5) = 1; Pc(10,5) = pn(i,j)*pq(i,j);
        Pc(2,6) = pbetachi(i,j); Pc(4,6) = pH(i,j);
        Pc(3,6) = pn(i,j)*pq(i,j)*pbetachi(i,j); 
        Pc(6,6) = 1; Pc(9,6) = -pn(i,j)*pq(i,j);
        
        Pc(3,7) = plnr2Ds(i,j); Pc(5,7) = 1;Pc(7,7) = 1;
        Pc(9,7) = plnr2Dchi(i,j);Pc(11,7) = 1;
        Pc(13,7) = plnr2Dchi(i,j);Pc(14,7) = -pn(i,j)*pq(i,j);
        Pc(4,8) = plnr2Ds(i,j); Pc(6,8) = 1;Pc(8,8) = 1;
        Pc(10,8) = plnr2Dchi(i,j);Pc(12,8) = 1;
        Pc(14,8) = plnr2Dchi(i,j);Pc(13,8) = pn(i,j)*pq(i,j);
        
        Pc(3,9) = 1; 
        Pc(4,10) = 1;
        
        Pc(3,11) = 1;
        Pc(4,12) = 1;
        Pc(3,13) = -pbetachi(i,j);Pc(9,13) = 1;Pc(13,13) = 1;
        Pc(4,14) = -pbetachi(i,j);Pc(10,14) = 1;Pc(14,14) = 1;
        Pc(13,15) = 1;
        Pc(14,16) = 1;
        
        Pv = zeros(16,16);
        Pv(1,:) = Pc(1,:)*Pf( 1)+Pc(3,:)*Pf( 2)+Pc(5,:)*Pf( 3);
        Pv(2,:) = Pc(2,:)*Pf( 1)+Pc(4,:)*Pf( 2)+Pc(6,:)*Pf( 3);
        Pv(3,:) = Pc(1,:)*Pf( 4)+Pc(3,:)*Pf( 5)+Pc(5,:)*Pf( 6);
        Pv(4,:) = Pc(2,:)*Pf( 4)+Pc(4,:)*Pf( 5)+Pc(6,:)*Pf( 6);
        Pv(5,:) = Pc(1,:)*Pf(11)+Pc(3,:)*Pf(12)+Pc(5,:)*Pf(13);
        Pv(6,:) = Pc(2,:)*Pf(11)+Pc(4,:)*Pf(12)+Pc(6,:)*Pf(13);
        Pv(7,:) = Pc(1,:)*Pf(14)+Pc(3,:)*Pf(15)+Pc(5,:)*Pf(16);
        Pv(8,:) = Pc(2,:)*Pf(14)+Pc(4,:)*Pf(15)+Pc(6,:)*Pf(16);
        
        Pv( 9,:) = Pc(7,:)*Pf(7)+Pc( 9,:)*Pf(8);
        Pv(10,:) = Pc(8,:)*Pf(7)+Pc(10,:)*Pf(8);
        Pv(11,:) = Pc(7,:)*Pf(9)+Pc( 9,:)*Pf(10);
        Pv(12,:) = Pc(8,:)*Pf(9)+Pc(10,:)*Pf(10);
        
        Pv(13,:) = Pc(11,:)*Pf(7)+Pc(13,:)*Pf(8);
        Pv(14,:) = Pc(12,:)*Pf(7)+Pc(14,:)*Pf(8);
        Pv(15,:) = Pc(11,:)*Pf(9)+Pc(13,:)*Pf(10);
        Pv(16,:) = Pc(12,:)*Pf(9)+Pc(14,:)*Pf(10);
        
        
        
        conA = Pv(:,1:10)*(Pc1(1:10).*Pv(:,1:10)');
        conB = Pv(:,11:end)*(Pc1(11:end).*Pv(:,11:end)');
        % symmetric conditions
        if i==1
            sign = -1;
            as_idx = [1,5,10,14];s_idx = [2,6,9,13];
            conA(:,as_idx+2) = conA(:,as_idx+2)+sign*conA(:,as_idx);
            conA(:, s_idx+2) = conA(:, s_idx+2)-sign*conA(:, s_idx);
            conA(as_idx+2,:) = conA(as_idx+2,:)+sign*conA(as_idx,:);
            conA( s_idx+2,:) = conA( s_idx+2,:)-sign*conA( s_idx,:);
            conA([as_idx,s_idx],:) = 0;
            conA(:,[as_idx,s_idx]) = 0;
            conA = conA+sparse([as_idx,s_idx],[as_idx,s_idx],1,16,16);
            
            conB(:,as_idx+2) = conB(:,as_idx+2)+sign*conB(:,as_idx);
            conB(:,s_idx+2) = conB(:,s_idx+2)-sign*conB(:,s_idx);
            conB(as_idx+2,:) = conB(as_idx+2,:)+sign*conB(as_idx,:);
            conB(s_idx+2,:) = conB(s_idx+2,:)-sign*conB(s_idx,:);
            conB([as_idx,s_idx],:) = 0;
            conB(:,[as_idx,s_idx]) = 0;
            conB = conB+sparse([as_idx,s_idx],[as_idx,s_idx],1,16,16);
            
            conA = 1/2*conA;
            conB = 1/2*conB;
        end
        if i==n_chi
            sign = -1;
            as_idx = [3,7,12,16];s_idx = [4,8,11,15];
            conA(:,as_idx-2) = conA(:,as_idx-2)+sign*conA(:,as_idx);
            conA(:, s_idx-2) = conA(:, s_idx-2)-sign*conA(:, s_idx);
            conA(as_idx-2,:) = conA(as_idx-2,:)+sign*conA(as_idx,:);
            conA( s_idx-2,:) = conA( s_idx-2,:)-sign*conA( s_idx,:);
            conA([as_idx,s_idx],:) = 0;
            conA(:,[as_idx,s_idx]) = 0;
            conA = conA+sparse([as_idx,s_idx],[as_idx,s_idx],1,16,16);
            
            conB(:,as_idx-2) = conB(:,as_idx-2)+sign*conB(:,as_idx);
            conB(:, s_idx-2) = conB(:, s_idx-2)-sign*conB(:, s_idx);
            conB(as_idx-2,:) = conB(as_idx-2,:)+sign*conB(as_idx,:);
            conB( s_idx-2,:) = conB( s_idx-2,:)-sign*conB( s_idx,:);
            conB([as_idx,s_idx],:) = 0;
            conB(:,[as_idx,s_idx]) = 0;
            conB = conB+sparse([as_idx,s_idx],[as_idx,s_idx],1,16,16);
            
            conA = 1/2*conA;
            conB = 1/2*conB;
        end
        xij = (2*n_chi+2)*(3*j-3)+2*(i-1);
        Map = xij+[1,2,3,4]'+[0,3,1,2]*(2*n_chi+2);
        Map = Map(:);
        At(Map,Map) = At(Map,Map)+ds*dchi*conA;
        Bt(Map,Map) = Bt(Map,Map)+ds*dchi*conB;
    end
end
%% symmetry conditions
% force the new tilda variables to zero
w2inv = 1; % which make the condition fullfilled when w2~=10^20

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
w2 = eig(full(At),full(Bt));
disp(min(w2));

GrowthRate =w2;


