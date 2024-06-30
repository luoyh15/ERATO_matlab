function I = IntegralPath(fpath, fint,r1,z1,r2,z2,zmax,eps)
% number of intervals
ninterval = 3;
% the precision of integral
ep = eps+1;
% discretization of r direction
rold = linspace(r1,r2,ninterval);
% the corresponding z
zold = zeros(size(rold));
zold(1) = z1; zold(end) = z2;
for i = 2:ninterval-1
    zold(i) = fzero(@(z) fpath(rold(i),z),[z1,zmax]);
end
dl = [0,sqrt((rold(2:end)-rold(1:end-1)).^2+(zold(2:end)-zold(1:end-1)).^2)];
T2 = sum(fint(rold,zold).*dl);
T1 = T2;
while (ep >= eps) 
    ninterval = 2*ninterval-1;     
    % discretization of r direction
    rmid = (rold(1:end-1)+rold(2:end))/2;
    % the corresponding z
    zmid = zeros(size(rmid));
    for i = 1:length(rmid)
        zmid(i) = fzero(@(z) fpath(rmid(i),z),[z1,zmax]);
    end
    rnew = zeros(1,ninterval) ;
    znew = zeros(size(rnew));
    rnew(1:2:end) = rold;
    rnew(2:2:end) = rmid;
    znew(1:2:end) = zold;
    znew(2:2:end) = zmid;
    dl = [0,sqrt((rnew(2:end)-rnew(1:end-1)).^2+(znew(2:end)-znew(1:end-1)).^2)];
    T2 = sum(fint(rnew,znew).*dl);
    ep = abs(T2 - T1);
    T1 = T2;
    rold = rnew;
    zold = znew;
    
end
I = T2;