q0=0.3;
n=2;
n_s = 14;
midx = 2*(1:2:n_s);
growthrate=zeros(size(midx));
for i = 1:length(midx)
    growthrate(i) = -min(ERATO18(2,1/3,q0,n,n_s,midx(i)));
end
plot(midx,growthrate);