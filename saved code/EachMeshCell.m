
dimcell = 16;
Mcell = zeros(dimcell,dimcell);
xR =  zeros(1,dimcell);
xI = zeros(1,dimcell);
yR =  zeros(1,dimcell);
yI = zeros(1,dimcell);
vR =  zeros(1,dimcell);
vI = zeros(1,dimcell);
dxdsR =  zeros(1,dimcell);
dxdsI = zeros(1,dimcell);
dxdcR =  zeros(1,dimcell);
dxdcI = zeros(1,dimcell);
dydcR =  zeros(1,dimcell);
dydcI = zeros(1,dimcell);
dvdcR =  zeros(1,dimcell);
dvdcI = zeros(1,dimcell);
xR([1,3,5,7]) = 1/4;
xI([2,4,6,8]) = 1/4;
yR([9,11]) = 1/2;
yI([10,12]) = 1/2;
vR([13,15]) = 1/2;
vI([14,16]) = 1/2;
Atotal = zeros(n_x,n_x);
Btotal = zeros(n_x,n_x);

% transform matrix for symmetric condition
U1 = diag(ones(1,dimcell));
U1(1,3) = 1;
U1(2,4) = -1;
U1(5,7) = 1;
U1(6,8) = -1;
U1(9,11) = -1;
U1(10,12) = 1;
U1(13,15) = -1;
U1(14,16) = 1;
U1inv = inv(U1);
U2 = diag(ones(1,dimcell));
U2(3,1) = 1;
U2(4,2) = -1;
U2(7,5) = 1;
U2(8,6) = -1;
U2(11,9) = -1;
U2(12,10) = 1;
U2(15,13) = -1;
U2(16,14) = 1;
U2inv = inv(U2);
for i = 1:n_s
    for j = 1:n_chi
        dxdscoef = [-1, -1, 1, 1]/(2*(ms(i+1)-ms(i)));
        dxdsR([1,3,5,7]) = dxdscoef;
        dxdsI([2,4,6,8]) = dxdscoef;
        dxdccoef = [-1, 1, -1, 1]/(2*(mchi(j+1)-mchi(j)));
        dxdcR([1,3,5,7]) = dxdccoef;
        dxdcI([2,4,6,8]) = dxdccoef;
        dydccoef = [-1,1]/(mchi(j+1)-mchi(j));
        dydcR([9,11]) = dydccoef;
        dydcI([10,12]) = dydccoef;
        dvdccoef = [-1,1]/(mchi(j+1)-mchi(j));
        dvdcR([13,15]) = dvdccoef;
        dvdcI([14,16]) = dvdccoef;
        
        i1R = 1./pq(j,i).*dxdcR-n*xI;
        i1I = 1./pq(j,i).*dxdcI+n*xR;
        i2R = dxdsR+dvdcR;
        i2I = dxdsI+dvdcI;
        i3R = pH(j,i).*xR-pbetachi(j,i).*n.*pq(j,i).*xI...
            -n*pq(j,i).*vI+pbetachi(j,i).*dxdcR+dxdsR;
        i3I = pH(j,i).*xI+pbetachi(j,i).*n.*pq(j,i).*xR...
            +n*pq(j,i).*vR+pbetachi(j,i).*dxdcI+dxdsI;
        i4R = plogr2Ds(j,i).*xR+plogr2Dchi(j,i).*(vR+pq(j,i).*yR)...
            -n*pq(j,i).^2.*yI+dxdsR+dvdcR+dydcR;
        i4I = plogr2Ds(j,i).*xI+plogr2Dchi(j,i).*(vI+pq(j,i).*yI)...
            +n*pq(j,i).^2.*yR+dxdsI+dvdcI+dydcI;
        i5R = xR;
        i5I = xI;
        
        Acell =   i1R'*(a(j,i)./ps(j,i).*i1R) + i1I'*(a(j,i)./ps(j,i).*i1I)...
            + i2R'*(b(j,i)./ps(j,i).*i2R) + i2I'*(b(j,i)./ps(j,i).*i2I)...
            + i3R'*(c(j,i)./ps(j,i).*i3R) + i3I'*(c(j,i)./ps(j,i).*i3I)...
            + i4R'*(d(j,i)./ps(j,i).*i4R) + i4I'*(d(j,i)./ps(j,i).*i4I)...
            - i5R'*(e(j,i)./ps(j,i).*i5R) - i5I'*(e(j,i)./ps(j,i).*i5I);
        
        Bcell =  xR'*(f(j,i)./ps(j,i).*xR) + xI'*(f(j,i)./ps(j,i).*xI) +...
            (vR+yR-pbetachi(j,i).*xR)'*(g(j,i)./ps(j,i).*(vR+yR-pbetachi(j,i).*xR))+...
            (vI+yI-pbetachi(j,i).*xI)'*(g(j,i)./ps(j,i).*(vI+yI-pbetachi(j,i).*xI))+...
            yR'*(h(j,i)./ps(j,i).*yR) + yI'*(h(j,i)./ps(j,i).*yI);
        IndexMap = [(3*i-3)*(2*n_chi+2)+2*j-2+(1:4),(3*i)*(2*n_chi+2)+2*j-2+(1:4),...
            (3*i-2)*(2*n_chi+2)+2*j-2+(1:4),(3*i-1)*(2*n_chi+2)+2*j-2+(1:4)];
        if(j==1)
            % symmetric condition
            Acell = Acell/2;
            Bcell = Bcell/2;
            Acell = U1inv'*Acell*U1inv;
            Bcell = U1inv'*Bcell*U1inv;
            Acell = SetDiag(Acell,[1,2,5,6,9,10,13,14],1);
            Bcell = SetDiag(Bcell,[1,2,5,6,9,10,13,14],1);            
        elseif(j==n_chi)
            Acell = Acell/2;
            Bcell = Bcell/2;
            Acell = U2inv'*Acell*U2inv;
            Bcell = U2inv'*Bcell*U2inv;
            Acell = SetDiag(Acell,[3,4,7,8,11,12,15,16],1);
            Bcell = SetDiag(Bcell,[3,4,7,8,11,12,15,16],1);
        end
%         if(i==1)
%             Acell = SetDiag(Acell,[1,2,3,4],1);
%             Bcell = SetDiag(Bcell,[1,2,3,4],1);
%         elseif(i==n_s)
%             Acell = SetDiag(Acell,[5,6,7,8],1);
%             Bcell = SetDiag(Bcell,[5,6,7,8],1);
%         end
        Atotal(IndexMap,IndexMap) = Atotal(IndexMap,IndexMap)+Acell;
        Btotal(IndexMap,IndexMap) = Btotal(IndexMap,IndexMap)+Bcell;
       
        
    end
end


function A = SetDiag(A,index,value)
for i = 1:length(index)
    A(index(i),:) = 0;
    A(:,index(i)) = 0;
    A(index(i),index(i)) = value;
end
end