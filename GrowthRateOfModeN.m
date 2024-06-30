function gamma2 = GrowthRateOfModeN(quantities,n)

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
%% construct the matrix A and B from each cell
dimcell = 16;
xR = sparse(1,[1,3,5,7],1/4,1,dimcell);
xI = sparse(1,[2,4,6,8],1/4,1,dimcell);
vR = sparse(1,[9,11],1/2,1,dimcell);
vI = sparse(1,[10,12],1/2,1,dimcell);
yR = sparse(1,[13,15],1/2,1,dimcell);
yI = sparse(1,[14,16],1/2,1,dimcell);
Atotal = sparse(n_x,n_x);
Btotal = sparse(n_x,n_x);

w2inv = 1e-20;
% transform matrix for symmetric condition
U1inv = sparse(1:dimcell,1:dimcell,1,dimcell,dimcell)-...
    sparse([1,5,10,14],[3,7,12,16],1,dimcell,dimcell)-...
    sparse([2,6,9,13],[4,8,11,15],-1,dimcell,dimcell);

U2inv = sparse(1:dimcell,1:dimcell,1,dimcell,dimcell)-...
    sparse([3,7,12,16],[1,5,10,14],1,dimcell,dimcell)-...
    sparse([4,8,11,15],[2,6,9,13],-1,dimcell,dimcell);

for i = 1:n_s
    for j = 1:n_chi
        msstep = ms(i+1)-ms(i);
        mchistep = mchi(j+1)-mchi(j);
        dxdscoef = [-1, -1, 1, 1]/(2*msstep);
        dxdsR = sparse(1,[1,3,5,7],dxdscoef,1,dimcell);
        dxdsI = sparse(1,[2,4,6,8],dxdscoef,1,dimcell);
        dxdccoef = [-1, 1, -1, 1]/(2*mchistep);
        dxdcR = sparse(1,[1,3,5,7],dxdccoef,1,dimcell);
        dxdcI = sparse(1,[2,4,6,8],dxdccoef,1,dimcell);
        dvdccoef = [-1,1]/mchistep;
        dvdcR = sparse(1,[9,11],dvdccoef,1,dimcell);
        dvdcI = sparse(1,[10,12],dvdccoef,1,dimcell);
        dydccoef = [-1,1]/mchistep;
        dydcR = sparse(1,[13,15],dydccoef,1,dimcell);
        dydcI = sparse(1,[14,16],dydccoef,1,dimcell);
        
        i1R = 1./pq(j,i).*dxdcR-n*xI;
        i1I = 1./pq(j,i).*dxdcI+n*xR;
        i2R = dxdsR+dvdcR;
        i2I = dxdsI+dvdcI;
        i3R = pH(j,i).*xR-pbetachi(j,i).*n.*pq(j,i).*xI...
            +n*pq(j,i).*vI+pbetachi(j,i).*dxdcR+dxdsR;
        i3I = pH(j,i).*xI+pbetachi(j,i).*n.*pq(j,i).*xR...
            -n*pq(j,i).*vR+pbetachi(j,i).*dxdcI+dxdsI;
        i4R = plnr2Ds(j,i).*xR+plnr2Dchi(j,i).*(vR+yR)...
            -n*pq(j,i).*yI+dxdsR+dvdcR+dydcR;
        i4I = plnr2Ds(j,i).*xI+plnr2Dchi(j,i).*(vI+yI)...
            +n*pq(j,i).*yR+dxdsI+dvdcI+dydcI;
        i5R = xR;
        i5I = xI;
        
        Acell =   i1R'*(a(j,i)./ps(j,i).*i1R) + i1I'*(a(j,i)./ps(j,i).*i1I)...
            + i2R'*(b(j,i)./ps(j,i).*i2R) + i2I'*(b(j,i)./ps(j,i).*i2I)...
            + i3R'*(c(j,i)./ps(j,i).*i3R) + i3I'*(c(j,i)./ps(j,i).*i3I)...
            + i4R'*(d(j,i)./ps(j,i).*i4R) + i4I'*(d(j,i)./ps(j,i).*i4I)...
            - i5R'*(e(j,i)./ps(j,i).*i5R) - i5I'*(e(j,i)./ps(j,i).*i5I);
        Acell = Acell*msstep*mchistep;
        Bcell =  xR'*(f(j,i)./ps(j,i).*xR) + xI'*(f(j,i)./ps(j,i).*xI) +...
            (vR+yR-pbetachi(j,i).*xR)'*(g(j,i)./ps(j,i).*(vR+yR-pbetachi(j,i).*xR))+...
            (vI+yI-pbetachi(j,i).*xI)'*(g(j,i)./ps(j,i).*(vI+yI-pbetachi(j,i).*xI))+...
            yR'*(h(j,i)./ps(j,i).*yR) + yI'*(h(j,i)./ps(j,i).*yI);
        Bcell = Bcell*msstep*mchistep;
        
        IndexMap = [(3*i-3)*(2*n_chi+2)+2*j-2+(1:4),(3*i)*(2*n_chi+2)+2*j-2+(1:4),...
            (3*i-2)*(2*n_chi+2)+2*j-2+(1:4),(3*i-1)*(2*n_chi+2)+2*j-2+(1:4)];
        % make matrixes symmetry
        Acell = (Acell+Acell')/2;
        Bcell = (Bcell+Bcell')/2;
%         % make small element to 0
%         Index_small = abs(Acell/max(abs(Acell(:))))<1e-15;
%         Acell(Index_small) = 0;
%         Acell = sparse(Acell);
%         Index_small = abs(Bcell/max(abs(Bcell(:))))<1e-15;
%         disp(min(abs(Bcell(Bcell~=0)/max(abs(Bcell(:))))));
%         Bcell(Index_small) = 0;
%         Bcell = sparse(Bcell);
        if(j==1)
            % symmetric condition
            Acell = Acell/2;
            Bcell = Bcell/2;
            Acell = U1inv'*Acell*U1inv;
            Bcell = U1inv'*Bcell*U1inv;
            Acell = SetDiag(Acell,[1,2,5,6,9,10,13,14],1);
            Bcell = SetDiag(Bcell,[1,2,5,6,9,10,13,14],w2inv);
        elseif(j==n_chi)
            Acell = Acell/2;
            Bcell = Bcell/2;
            Acell = U2inv'*Acell*U2inv;
            Bcell = U2inv'*Bcell*U2inv;
            Acell = SetDiag(Acell,[3,4,7,8,11,12,15,16],1);
            Bcell = SetDiag(Bcell,[3,4,7,8,11,12,15,16],w2inv);
        end
        if(i==1)
            Acell = SetDiag(Acell,[1,2,3,4],1);
            Bcell = SetDiag(Bcell,[1,2,3,4],w2inv);
        elseif(i==n_s)
            Acell = SetDiag(Acell,[5,6,7,8],1);
            Bcell = SetDiag(Bcell,[5,6,7,8],w2inv);
        end
%         % make small element to 0
%         Index_small = abs(Acell/max(abs(Acell(:))))<1e-15;
%         Acell(Index_small) = 0;
%         Acell = sparse(Acell);
%         Index_small = abs(Bcell/max(abs(Bcell(:))))<1e-15;
%         disp(min(abs(Bcell(Bcell~=0)/max(abs(Bcell(:))))));
%         Bcell(Index_small) = 0;
%         Bcell = sparse(Bcell);
        
        [IndexCol,IndexRow] = meshgrid(IndexMap,IndexMap);
        Atotal = Atotal + sparse(IndexRow(:),IndexCol(:),Acell(:),n_x,n_x);
        Btotal= Btotal + sparse(IndexRow(:),IndexCol(:),Bcell(:),n_x,n_x);
        
        % make small element to 0
%         Index_small = abs(Atotal)>0&&abs(Atotal/max(abs(Atotal(:))))<1e-15;
%         Atotal(Index_small) = 0;
%         Atotal = sparse(Atotal);
%         Index_small = abs(Btotal)>0&&abs(Btotal/max(abs(Btotal(:))))<1e-15;
%         disp(min(abs(Btotal(Btotal~=0)/max(abs(Btotal(:))))));
%         Btotal(Index_small) = 0;
%         Btotal = sparse(Btotal);
    end
end

[~,p] = chol(Atotal);
if p==0
    gamma2 = 0;
else
    gamma2 = 0;
    w2guess = -1;
    while w2guess>-10^3
%         [eigFun,eigVal] = eigs(Atotal,Btotal,1,w2guess);
        [eigFun,eigVal] = eigs(Atotal-w2guess*Btotal,Btotal,1,'smallestabs');
        eigVal = eigVal+w2guess;
        if eigVal<0
            gamma2 = eigVal;
            break;
        end
        w2guess = w2guess*10;
    end
    if gamma2<0
        plotEigenMode(eigFun);
    end
end


    function plotEigenMode(eigFun)
        Fx = eigFun;
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
        
        i_xij = i_xij + (2*n_chi+2);
        % VR
        V_coef = [1, 1]/2;
        VR = sparse(i_cell,i_xij  ,V_coef(1),n_cell,n_x)+...
            sparse(i_cell,i_xij+2,V_coef(2),n_cell,n_x);
        % VI
        VI = sparse(i_cell,i_xij+1,V_coef(1),n_cell,n_x)+...
            sparse(i_cell,i_xij+3,V_coef(2),n_cell,n_x);
        
        i_xij = i_xij + (2*n_chi+2);
        % YR
        Y_coef = [1, 1]/2;
        YR = sparse(i_cell,i_xij  ,Y_coef(1),n_cell,n_x)+...
            sparse(i_cell,i_xij+2,Y_coef(2),n_cell,n_x);
        % YI
        YI = sparse(i_cell,i_xij+1,Y_coef(1),n_cell,n_x)+...
            sparse(i_cell,i_xij+3,Y_coef(2),n_cell,n_x);
        %% retransform 
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
        
        Fx = U\Fx;
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


function A = SetDiag(A,index,value)
for i = 1:length(index)
    A(index(i),:) = 0;
    A(:,index(i)) = 0;
    A(index(i),index(i)) = value;
    A = sparse(A);
end
end



