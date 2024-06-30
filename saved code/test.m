n = 1:3;
q0 = 0.01:0.02:1.3;
n_s = 16;
% q0 = 0.1;
% n_s = (2:10)*10+4;
growthrate1 = zeros(length(q0),length(n));

for j=1:length(n)
    for i = 1:length(q0)
        growthrate1(i,j) = -min(ERATO21(2,1/3,q0(i),n(j),n_s));
    end
end
figure(4);hold on; 
for i=1:length(n)
    plot(q0,growthrate1(:,i)'.*q0.^2);
end
% plot(n_s,growthrate1.*q0.^2);